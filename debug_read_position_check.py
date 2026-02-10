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
    for c in cands:
        if c.type.value == "del" and c.length == args.target_del_len:
            print(
                "SCANNER",
                f"{c.length}D",
                "out_pos",
                c.ref_position + 1,
                "ctx",
                c.sequencecontext_brackets(),
            )

    bam.close()


if __name__ == "__main__":
    main()
