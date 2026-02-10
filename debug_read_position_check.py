#!/usr/bin/env python3
"""
Debug one read's indel coordinates through three paths:
1) raw numeric CIGAR tracking
2) scanner enum tracking (as_cigar + REF_CONSUMING_OPS)
3) actual scanner candidate output

This helps isolate where coordinate drift is introduced.
"""

import argparse
import yaml
import pysam
import pyfastx

from indel_scanner.configurator import PipelineConfig
from indel_scanner.scanner import ContigScanner, REF_CONSUMING_OPS
from indel_scanner.utils import as_cigar
from indel_scanner.STR_Classifier import STRClassifier
from indel_scanner.filters import IndelFilters
from indel_scanner.indel import FilterFlag


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--strdir", required=True)
    p.add_argument("--config", default="config.yaml")
    p.add_argument("--contig", required=True)
    p.add_argument("--read", required=True)
    p.add_argument("--output", default="/tmp/indel_debug_out")
    p.add_argument("--target-del-len", type=int, default=3)
    p.add_argument("--window-start", type=int, default=None)
    p.add_argument("--window-end", type=int, default=None)
    args = p.parse_args()

    print(
        "pysam op codes:",
        "CMATCH=", pysam.CMATCH,
        "CINS=", pysam.CINS,
        "CDEL=", pysam.CDEL,
        "CREF_SKIP=", pysam.CREF_SKIP,
        "CEQUAL=", pysam.CEQUAL,
        "CDIFF=", pysam.CDIFF,
    )

    bam = pysam.AlignmentFile(args.bam, "rb")
    read = None
    for r in bam.fetch(args.contig):
        if r.query_name == args.read and not r.is_secondary and not r.is_supplementary:
            read = r
            break
    if read is None:
        print("read not found")
        bam.close()
        return

    print("READ:", read.query_name)
    print("POS/CIGAR:", read.reference_start + 1, read.cigarstring, "MAPQ", read.mapping_quality)

    # A) Raw numeric tracking
    ref_num = read.reference_start
    for op, length in read.cigartuples or []:
        if op == pysam.CDEL and length == args.target_del_len:
            print("NUMERIC", f"{length}D", "at out_pos", ref_num + 1)
        if op in (pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF):
            ref_num += length

    # B) Scanner-style enum tracking
    ref_enum = read.reference_start
    none_ops = 0
    for op, length in read.cigartuples or []:
        cg = as_cigar(op)
        if cg is None:
            none_ops += 1
        if cg is not None and cg.name == "OP_D" and length == args.target_del_len:
            print("ENUM   ", f"{length}D", "at out_pos", ref_enum + 1)
        if cg in REF_CONSUMING_OPS:
            ref_enum += length
    print("as_cigar None ops:", none_ops)

    # C) Full scanner candidate output
    with open(args.config, "r", encoding="utf-8") as fh:
        cfg_data = yaml.safe_load(fh) or {}
    cfg_args = argparse.Namespace(
        bam=args.bam,
        fasta=args.fasta,
        output=args.output,
        strdir=args.strdir,
        config=args.config,
        contigs=args.contig,
        test_run=False,
        test_contigs_count=10,
    )
    cfg = PipelineConfig(cfg_args, cfg_data)
    scanner = ContigScanner(cfg)
    strc = STRClassifier(cfg, args.contig)
    fasta = pyfastx.Fasta(args.fasta)
    contig_seq = fasta[args.contig].seq

    cands = list(scanner._parse_cigar_for_candidates(read, strc, contig_seq))
    print("\nRAW_CANDIDATES")
    for c in cands:
        out_pos = c.ref_position if c.type.value == "ins" else c.ref_position + 1
        if args.window_start is not None and out_pos < args.window_start:
            continue
        if args.window_end is not None and out_pos > args.window_end:
            continue
        if c.type.value == "del" and c.length == args.target_del_len:
            print(
                "SCANNER",
                f"{c.length}D",
                "out_pos",
                out_pos,
                "ctx",
                c.sequencecontext_brackets(),
            )
    if args.window_start is None and args.window_end is None:
        print(f"(total candidates for read: {len(cands)})")

    # Build filter context the same way scanner does before IndelFilters.apply(...)
    location_counts = {}
    position_counts = {}
    for c in cands:
        key = (c.ref_position, c.type, c.length)
        location_counts[key] = location_counts.get(key, 0) + 1
        position_counts[c.ref_position] = position_counts.get(c.ref_position, 0) + 1

    target_positions = set()
    for p0 in position_counts.keys():
        for d in (-1, 0, 1):
            target_positions.add(p0 + d)

    coverage_by_pos = {}
    if target_positions:
        for rr in bam.fetch(contig=args.contig):
            if not scanner._valid_read(rr) or rr.mapping_quality < cfg.min_map_quality:
                continue
            for rp, refp in rr.get_aligned_pairs(matches_only=True):
                if rp is None or refp is None:
                    continue
                if refp in target_positions:
                    coverage_by_pos[refp] = coverage_by_pos.get(refp, 0) + 1

    filter_context = {
        "location_map": location_counts,
        "position_counts": position_counts,
        "coverage_by_pos": coverage_by_pos,
        "str_classifier": strc,
    }

    print("\nPOST_FILTER_STATUS")
    for c in cands:
        c.filter_mask = FilterFlag.NONE
        IndelFilters.apply(cfg.processor_filters, c, **filter_context)
        out_pos = c.ref_position if c.type.value == "ins" else c.ref_position + 1
        if args.window_start is not None and out_pos < args.window_start:
            continue
        if args.window_end is not None and out_pos > args.window_end:
            continue
        if c.filter_mask == FilterFlag.NONE:
            reason = "PASS"
        else:
            reasons = []
            for flag in FilterFlag:
                if flag is FilterFlag.NONE:
                    continue
                if c.filter_mask & flag:
                    reasons.append(flag.name)
            reason = ",".join(reasons) if reasons else str(int(c.filter_mask))
        print(
            c.type.value,
            out_pos,
            c.length,
            c.sequencecontext_brackets(),
            "MAPQ",
            c.map_quality,
            "FILTER",
            reason,
        )

    bam.close()


if __name__ == "__main__":
    main()
