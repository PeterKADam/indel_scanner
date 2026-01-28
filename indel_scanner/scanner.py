# indel_scanner/scanner.py
from pathlib import Path
from typing import Generator, Optional
import pysam
import pyfastx
import time
import logging

from .IO import write_passed_indels
from .indel import IndelRecord
from .STR_Classifier import STRClassifier
from .filters import IndelFilters
from .indel import INDEL_TYPE
from .utils import Cigar, as_cigar
from .configurator import PipelineConfig

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}


class ContigScanner:
    def __init__(self, config: PipelineConfig):
        self.config = config

    def _parse_cigar_for_candidates(
        self,
        read: pysam.AlignedSegment,
        str_classifier: STRClassifier,
        contig_seq: str,
    ) -> Generator[IndelRecord, None, None]:
        ref_pos_tracker = read.reference_start
        read_pos_tracker = 0
        if read.cigartuples is None:
            return
        qualities = read.query_qualities

        for op_int, length in read.cigartuples:
            op = as_cigar(op_int)

            if op == Cigar.OP_I:
                if length >= self.config.min_indel_size:
                    prefix_quality = None
                    indel_quality = None
                    suffix_quality = None
                    if qualities is not None:
                        flank_len = 5
                        prefix_start = max(0, read_pos_tracker - flank_len)
                        insertion_end = read_pos_tracker + length
                        suffix_end = min(len(qualities), insertion_end + flank_len)
                        prefix_quality = list(qualities[prefix_start:read_pos_tracker])
                        indel_quality = list(qualities[read_pos_tracker:insertion_end])
                        suffix_quality = list(qualities[insertion_end:suffix_end])
                    yield IndelRecord(
                        contig=read.reference_name,
                        ref_position=ref_pos_tracker,
                        type=INDEL_TYPE.INSERTION,
                        length=length,
                        prefix_context=read.query_sequence[
                            max(0, read_pos_tracker - 5) : read_pos_tracker
                        ],
                        suffix_context=read.query_sequence[
                            read_pos_tracker + length : read_pos_tracker + length + 5
                        ],
                        read_name=read.query_name,
                        in_STR=str_classifier.matches_rptrf_motif(
                            ref_pos_tracker, contig_seq
                        ),
                        map_quality=read.mapping_quality,
                        indel_content=read.query_sequence[
                            read_pos_tracker : read_pos_tracker + length
                        ],
                        prefix_quality=prefix_quality,
                        indel_quality=indel_quality,
                        suffix_quality=suffix_quality,
                    )
            elif op == Cigar.OP_D:
                if length >= self.config.min_indel_size:
                    prefix_quality = None
                    suffix_quality = None
                    if qualities is not None:
                        flank_len = 5
                        prefix_start = max(0, read_pos_tracker - flank_len)
                        suffix_end = min(len(qualities), read_pos_tracker + flank_len)
                        prefix_quality = list(qualities[prefix_start:read_pos_tracker])
                        suffix_quality = list(qualities[read_pos_tracker:suffix_end])
                    yield IndelRecord(
                        contig=read.reference_name,
                        ref_position=ref_pos_tracker,
                        type=INDEL_TYPE.DELETION,
                        length=length,
                        prefix_context=contig_seq[
                            max(0, ref_pos_tracker - 5) : ref_pos_tracker
                        ],
                        indel_content=contig_seq[
                            ref_pos_tracker : ref_pos_tracker + length
                        ],
                        suffix_context=contig_seq[
                            ref_pos_tracker + length : ref_pos_tracker + length + 5
                        ],
                        read_name=read.query_name,
                        in_STR=str_classifier.matches_rptrf_motif(
                            ref_pos_tracker, contig_seq
                        ),
                        map_quality=read.mapping_quality,
                        prefix_quality=prefix_quality,
                        suffix_quality=suffix_quality,
                    )

            if op in REF_CONSUMING_OPS:
                ref_pos_tracker += length
            if op in READ_CONSUMING_OPS:
                read_pos_tracker += length

    def _is_callable_at_position(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        str_classifier: STRClassifier,
        contig_seq: str,
    ) -> bool:
        if read_pos >= len(read.query_qualities) - 1:
            return False  # Cannot form a flank pair at the very end of a read.

        if (
            read.query_qualities[read_pos] < self.config.min_flank_quality
            or read.query_qualities[read_pos + 1] < self.config.min_flank_quality
        ):
            return False

        if IndelFilters.check_if_in_str(ref_pos, str_classifier):
            return False
        if IndelFilters.check_if_homopolymer_context(
            contig_seq, ref_pos, self.config.min_homopolymer_len
        ):
            return False

        return True

    def _is_snp_candidate(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        contig_seq: str,
    ) -> bool:
        if read_pos is None or ref_pos is None:
            return False
        if read.query_qualities is None:
            return False
        if read.query_qualities[read_pos] < self.config.min_base_quality:
            return False
        if read_pos >= len(read.query_sequence):
            return False
        if ref_pos >= len(contig_seq):
            return False
        read_base = read.query_sequence[read_pos]
        ref_base = contig_seq[ref_pos]
        if read_base == "N" or ref_base == "N":
            return False
        return read_base != ref_base

    def _classify_indel_bin(self, indel_length: int) -> Optional[str]:
        for bin_cfg in self.config.indel_bins:
            if bin_cfg["min"] <= indel_length <= bin_cfg["max"]:
                return bin_cfg["label"]
        return None

    def _bin_label(self, bin_cfg: dict, indel_type: INDEL_TYPE) -> str:
        prefix = "ins" if indel_type == INDEL_TYPE.INSERTION else "del"
        return f"{prefix}_{bin_cfg['label']}"

    def _get_processor_min_quality(self) -> int:
        for filter_cfg in self.config.processor_filters:
            if filter_cfg.get("name") == "low_minimum_indel_quality":
                params = filter_cfg.get("params", {})
                return params.get("min_quality", self.config.min_base_quality)
        return self.config.min_base_quality

    def _is_callable_for_insertion_length(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        contig_seq: str,
        str_classifier: STRClassifier,
        length: int,
    ) -> bool:
        if read_pos is None or ref_pos is None:
            return False
        qualities = read.query_qualities
        if qualities is None:
            return False
        if read_pos <= 0 or read_pos + length >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[read_pos] < self.config.min_flank_quality
        ):
            return False
        insertion_window = qualities[read_pos : read_pos + length]
        min_quality = self._get_processor_min_quality()
        if insertion_window and min(insertion_window) < min_quality:
            return False
        if IndelFilters.check_if_in_str(ref_pos, str_classifier):
            return False
        if IndelFilters.check_if_homopolymer_context(
            contig_seq, ref_pos, self.config.min_homopolymer_len
        ):
            return False
        return True

    def _is_callable_for_deletion_length(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        contig_seq: str,
        str_classifier: STRClassifier,
        length: int,
    ) -> bool:
        if read_pos is None or ref_pos is None:
            return False
        qualities = read.query_qualities
        if qualities is None:
            return False
        if read_pos <= 0 or read_pos >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[read_pos] < self.config.min_flank_quality
        ):
            return False
        if IndelFilters.check_if_in_str(ref_pos, str_classifier):
            return False
        if IndelFilters.check_if_homopolymer_context(
            contig_seq, ref_pos, self.config.min_homopolymer_len
        ):
            return False
        return True

    def scan_contig_streaming(
        self,
        contig_name: str,
        samfile: pysam.AlignmentFile,
        fasta: pyfastx.Fasta,
        in_memory: bool,
    ) -> tuple[Optional[list[IndelRecord]], int, dict]:
        start_time = time.time()
        str_classifier = STRClassifier(self.config, contig_name)
        contig_seq = fasta[contig_name].seq
        interrogated_bases_count = 0
        total_aligned_bases = 0
        sampling_bases_total = 0
        sampling_passable_by_type: dict[str, int] = {}
        sampling_active = (
            self.config.sampling_strategy == "largest_contig"
            and contig_name == self.config.sampling_contig
        )
        sampled_contig = sampling_active
        reads_processed = 0
        candidates: list[IndelRecord] = []
        location_counts = {}
        type_counts: dict[str, int] = {}

        for read in samfile.fetch(contig=contig_name):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.mapping_quality < self.config.min_map_quality
                or read.query_sequence is None
                or read.query_qualities is None
            ):
                continue
            reads_processed += 1

            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if ref_pos is None or read_pos is None:
                    continue
                total_aligned_bases += 1
                is_callable = self._is_callable_at_position(
                    read, read_pos, ref_pos, str_classifier, contig_seq
                )
                if is_callable:
                    interrogated_bases_count += 1

                if sampling_active and sampling_bases_total < self.config.sampling_bases:
                    sampling_bases_total += 1
                    if is_callable:
                        sampling_passable_by_type[self.config.snp_label] = (
                            sampling_passable_by_type.get(self.config.snp_label, 0) + 1
                        )
                    for bin_cfg in self.config.indel_bins:
                        length_for_bin = bin_cfg["max"]
                        if self._is_callable_for_insertion_length(
                            read,
                            read_pos,
                            ref_pos,
                            contig_seq,
                            str_classifier,
                            length_for_bin,
                        ):
                            ins_label = self._bin_label(bin_cfg, INDEL_TYPE.INSERTION)
                            sampling_passable_by_type[ins_label] = (
                                sampling_passable_by_type.get(ins_label, 0) + 1
                            )
                        if self._is_callable_for_deletion_length(
                            read,
                            read_pos,
                            ref_pos,
                            contig_seq,
                            str_classifier,
                            length_for_bin,
                        ):
                            del_label = self._bin_label(bin_cfg, INDEL_TYPE.DELETION)
                            sampling_passable_by_type[del_label] = (
                                sampling_passable_by_type.get(del_label, 0) + 1
                            )

                if is_callable and self._is_snp_candidate(
                    read, read_pos, ref_pos, contig_seq
                ):
                    type_counts[self.config.snp_label] = (
                        type_counts.get(self.config.snp_label, 0) + 1
                    )

            for candidate in self._parse_cigar_for_candidates(
                read, str_classifier, contig_seq
            ):
                candidates.append(candidate)
                key = (candidate.ref_position, candidate.type, candidate.length)
                location_counts[key] = location_counts.get(key, 0) + 1

        filter_context = {"location_map": location_counts, "str_classifier": str_classifier}
        passed_records: list[IndelRecord] = []
        for candidate in candidates:
            IndelFilters.apply(self.config.processor_filters, candidate, **filter_context)
            if not candidate.is_filtered():
                passed_records.append(candidate)
                indel_bin = self._classify_indel_bin(candidate.length)
                if indel_bin:
                    indel_type = candidate.type
                    if indel_type == INDEL_TYPE.INSERTION:
                        label = f"ins_{indel_bin}"
                    else:
                        label = f"del_{indel_bin}"
                    type_counts[label] = type_counts.get(label, 0) + 1

        if not in_memory:
            part_path = self.config.passed_parts_dir / f"{contig_name}.part.tsv"
            write_passed_indels(part_path, passed_records, self.config.write_buffer_size)

        end_time = time.time()
        stats = {
            "reads_processed": reads_processed,
            "candidates": len(candidates),
            "passed": len(passed_records),
            "total_aligned_bases": total_aligned_bases,
            "sampling_bases_total": sampling_bases_total,
            "sampling_passable_by_type": sampling_passable_by_type,
            "sampled_contig": sampled_contig,
            "type_counts": type_counts,
            "elapsed_s": end_time - start_time,
        }
        return (passed_records if in_memory else None, interrogated_bases_count, stats)


