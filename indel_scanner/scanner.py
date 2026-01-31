# indel_scanner/scanner.py
from pathlib import Path
from typing import Generator, Optional
from collections import defaultdict
import random
import pysam
import pyfastx
import time
import logging

from .IO import write_passed_indels
from .indel import IndelRecord
from .STR_Classifier import STRClassifier
from .filters import IndelFilters
from .homopolymer_classifier import HomopolymerClassifier
from .indel import INDEL_TYPE
from .utils import Cigar, as_cigar, flank_qualities_by_ref
from .configurator import PipelineConfig

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}


class ContigScanner:
    @staticmethod
    def _aligned_blocks(
        read: pysam.AlignedSegment,
    ) -> list[tuple[int, int, int]]:
        if read.cigartuples is None:
            return []
        blocks: list[tuple[int, int, int]] = []
        ref_pos = read.reference_start
        read_pos = 0
        for op_int, length in read.cigartuples:
            op = as_cigar(op_int)
            if op in (Cigar.OP_M, Cigar.OP_EQ, Cigar.OP_X):
                blocks.append((read_pos, ref_pos, length))
                read_pos += length
                ref_pos += length
            elif op in (Cigar.OP_I, Cigar.OP_S):
                read_pos += length
            elif op in (Cigar.OP_D, Cigar.OP_N):
                ref_pos += length
        return blocks
    @staticmethod
    def _repeat_unit_length(seq: str, min_units: int) -> Optional[int]:
        seq = seq.upper()
        if not seq or "N" in seq:
            return None
        seq_len = len(seq)
        if seq_len < min_units:
            return None
        max_motif_len = seq_len // min_units
        for motif_len in range(1, max_motif_len + 1):
            if seq_len % motif_len != 0:
                continue
            units = seq_len // motif_len
            if units < min_units:
                continue
            motif = seq[:motif_len]
            if motif * units == seq:
                return motif_len
        return None

    @staticmethod
    def _has_repeat_run(seq: str, min_units: int, max_motif_len: int) -> bool:
        seq = seq.upper()
        if not seq or "N" in seq:
            return False
        seq_len = len(seq)
        for motif_len in range(1, max_motif_len + 1):
            min_run_len = motif_len * min_units
            if seq_len < min_run_len:
                continue
            for start in range(0, seq_len - min_run_len + 1):
                motif = seq[start : start + motif_len]
                if "N" in motif:
                    continue
                run_len = 0
                idx = start
                while idx + motif_len <= seq_len and seq[idx : idx + motif_len] == motif:
                    run_len += motif_len
                    if run_len >= min_run_len:
                        return True
                    idx += motif_len
        return False
    def __init__(self, config: PipelineConfig):
        self.config = config

    @staticmethod
    def _valid_read(read: pysam.AlignedSegment) -> bool:
        if (
            read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or read.reference_name is None
            or read.query_name is None
            or read.query_sequence is None
            or read.query_qualities is None
        ):
            return False
        if read.cigartuples is None:
            return False
        for op_int, _ in read.cigartuples:
            op = as_cigar(op_int)
            if op is None or op in (Cigar.OP_S, Cigar.OP_N):
                return False
        return True

    def _parse_cigar_for_candidates(
        self,
        read: pysam.AlignedSegment,
        str_classifier: STRClassifier,
        contig_seq: str,
    ) -> Generator[IndelRecord, None, None]:
        ref_name = read.reference_name
        read_name = read.query_name
        seq = read.query_sequence
        qualities = read.query_qualities
        if ref_name is None or read_name is None or seq is None or qualities is None:
            return
        ref_pos_tracker = read.reference_start
        read_pos_tracker = 0
        if read.cigartuples is None:
            return

        for op_int, length in read.cigartuples:
            op = as_cigar(op_int)

            if op == Cigar.OP_I:
                if length >= self.config.min_indel_size:
                    indel_quality = None
                    flank_len = 5
                    insertion_end = read_pos_tracker + length
                    prefix_quality, suffix_quality = flank_qualities_by_ref(
                        read,
                        ref_pos_tracker,
                        length,
                        flank_len,
                        is_insertion=True,
                    )
                    if (
                        len(prefix_quality) < flank_len
                        or len(suffix_quality) < flank_len
                        or any(q is None for q in prefix_quality)
                        or any(q is None for q in suffix_quality)
                    ):
                        continue
                    indel_quality = list(qualities[read_pos_tracker:insertion_end])
                    indel_content = seq[read_pos_tracker : read_pos_tracker + length]
                    filter_cfg = getattr(self.config, "str_candidate_filter", {})
                    str_motif_match = str_classifier.matches_rptrf_motif_length(
                        ref_pos_tracker, contig_seq, length, indel_content
                    )
                    if filter_cfg.get("enabled", False) and not str_motif_match:
                        window_bp = int(filter_cfg.get("window_bp", 0))
                        if window_bp > 0 and str_classifier.is_within_str_window(
                            ref_pos_tracker, window_bp
                        ):
                            continue
                        check_positions = {ref_pos_tracker, max(0, ref_pos_tracker - 1)}
                        if any(
                            str_classifier.is_str_like(pos, contig_seq)
                            for pos in check_positions
                        ):
                            continue
                        motif_len = self._repeat_unit_length(
                            indel_content, int(filter_cfg.get("min_repeat_units", 3))
                        )
                        local_window = int(filter_cfg.get("local_window_bp", 30))
                        max_motif_len = int(filter_cfg.get("max_motif_len", 6))
                        win_start = max(0, ref_pos_tracker - local_window)
                        win_end = min(len(contig_seq), ref_pos_tracker + local_window)
                        local_seq = contig_seq[win_start:win_end]
                        if filter_cfg.get("micro_repeat_enabled", False):
                            micro_window = int(
                                filter_cfg.get("micro_repeat_window_bp", local_window)
                            )
                            micro_min_units = int(
                                filter_cfg.get("micro_repeat_min_units", 2)
                            )
                            micro_max_motif_len = int(
                                filter_cfg.get("micro_repeat_max_motif_len", max_motif_len)
                            )
                            micro_start = max(0, ref_pos_tracker - micro_window)
                            micro_end = min(len(contig_seq), ref_pos_tracker + micro_window)
                            micro_seq = contig_seq[micro_start:micro_end]
                            if self._has_repeat_run(
                                micro_seq, micro_min_units, micro_max_motif_len
                            ):
                                continue
                        has_repeat_in_indel = self._has_repeat_run(
                            indel_content,
                            int(filter_cfg.get("min_repeat_units", 3)),
                            max_motif_len,
                        )
                        if filter_cfg.get("aggressive", False):
                            aggressive_min_units = int(
                                filter_cfg.get("aggressive_min_repeat_units", 2)
                            )
                            aggressive_max_motif_len = int(
                                filter_cfg.get("aggressive_max_motif_len", max_motif_len)
                            )
                            if self._has_repeat_run(
                                local_seq, aggressive_min_units, aggressive_max_motif_len
                            ) or self._has_repeat_run(
                                indel_content,
                                aggressive_min_units,
                                aggressive_max_motif_len,
                            ):
                                continue
                        if (
                            motif_len is not None
                            or has_repeat_in_indel
                            or self._has_repeat_run(
                                local_seq,
                                int(filter_cfg.get("min_repeat_units", 3)),
                                max_motif_len,
                            )
                        ):
                            continue
                    motif_length = str_classifier.motif_length_at(ref_pos_tracker)
                    if motif_length == 1:
                        str_motif_match = False
                    yield IndelRecord(
                        contig=ref_name,
                        ref_position=ref_pos_tracker,
                        type=INDEL_TYPE.INSERTION,
                        length=length,
                        prefix_context=seq[
                            max(0, read_pos_tracker - 5) : read_pos_tracker
                        ],
                        suffix_context=seq[
                            read_pos_tracker + length : read_pos_tracker + length + 5
                        ],
                        read_name=read_name,
                        in_STR_region=self._span_has_str(
                            ref_pos_tracker, length, str_classifier
                        ),
                        in_STR=str_motif_match,
                        motif_length=motif_length,
                        map_quality=read.mapping_quality,
                        indel_content=indel_content,
                        prefix_quality=prefix_quality,
                        indel_quality=indel_quality,
                        suffix_quality=suffix_quality,
                    )
            elif op == Cigar.OP_D:
                if length >= self.config.min_indel_size:
                    flank_len = 5
                    prefix_quality, suffix_quality = flank_qualities_by_ref(
                        read,
                        ref_pos_tracker,
                        length,
                        flank_len,
                        is_insertion=False,
                    )
                    if (
                        len(prefix_quality) < flank_len
                        or len(suffix_quality) < flank_len
                        or any(q is None for q in prefix_quality)
                        or any(q is None for q in suffix_quality)
                    ):
                        continue
                    indel_content = contig_seq[
                        ref_pos_tracker : ref_pos_tracker + length
                    ]
                    filter_cfg = getattr(self.config, "str_candidate_filter", {})
                    str_motif_match = str_classifier.matches_rptrf_motif_length(
                        ref_pos_tracker, contig_seq, length, indel_content
                    )
                    if filter_cfg.get("enabled", False) and not str_motif_match:
                        window_bp = int(filter_cfg.get("window_bp", 0))
                        if window_bp > 0 and str_classifier.is_within_str_window(
                            ref_pos_tracker, window_bp
                        ):
                            continue
                        end_pos = ref_pos_tracker + max(0, length - 1)
                        if str_classifier.is_str_like(
                            ref_pos_tracker, contig_seq
                        ) or str_classifier.is_str_like(end_pos, contig_seq):
                            continue
                        motif_len = self._repeat_unit_length(
                            indel_content, int(filter_cfg.get("min_repeat_units", 3))
                        )
                        local_window = int(filter_cfg.get("local_window_bp", 30))
                        max_motif_len = int(filter_cfg.get("max_motif_len", 6))
                        win_start = max(0, ref_pos_tracker - local_window)
                        win_end = min(len(contig_seq), ref_pos_tracker + local_window)
                        local_seq = contig_seq[win_start:win_end]
                        if filter_cfg.get("micro_repeat_enabled", False):
                            micro_window = int(
                                filter_cfg.get("micro_repeat_window_bp", local_window)
                            )
                            micro_min_units = int(
                                filter_cfg.get("micro_repeat_min_units", 2)
                            )
                            micro_max_motif_len = int(
                                filter_cfg.get("micro_repeat_max_motif_len", max_motif_len)
                            )
                            micro_start = max(0, ref_pos_tracker - micro_window)
                            micro_end = min(len(contig_seq), ref_pos_tracker + micro_window)
                            micro_seq = contig_seq[micro_start:micro_end]
                            if self._has_repeat_run(
                                micro_seq, micro_min_units, micro_max_motif_len
                            ):
                                continue
                        has_repeat_in_indel = self._has_repeat_run(
                            indel_content,
                            int(filter_cfg.get("min_repeat_units", 3)),
                            max_motif_len,
                        )
                        if filter_cfg.get("aggressive", False):
                            aggressive_min_units = int(
                                filter_cfg.get("aggressive_min_repeat_units", 2)
                            )
                            aggressive_max_motif_len = int(
                                filter_cfg.get("aggressive_max_motif_len", max_motif_len)
                            )
                            if self._has_repeat_run(
                                local_seq, aggressive_min_units, aggressive_max_motif_len
                            ) or self._has_repeat_run(
                                indel_content,
                                aggressive_min_units,
                                aggressive_max_motif_len,
                            ):
                                continue
                        if (
                            motif_len is not None
                            or has_repeat_in_indel
                            or self._has_repeat_run(
                                local_seq,
                                int(filter_cfg.get("min_repeat_units", 3)),
                                max_motif_len,
                            )
                        ):
                            continue
                    motif_length = str_classifier.motif_length_at(ref_pos_tracker)
                    if motif_length == 1:
                        str_motif_match = False
                    yield IndelRecord(
                        contig=ref_name,
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
                        read_name=read_name,
                        in_STR_region=self._span_has_str(
                            ref_pos_tracker, length, str_classifier
                        ),
                        in_STR=str_motif_match,
                        motif_length=motif_length,
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
        qualities = read.query_qualities
        if qualities is None:
            return False
        if read_pos >= len(qualities) - 1:
            return False  # Cannot form a flank pair at the very end of a read.

        if (
            qualities[read_pos] < self.config.min_flank_quality
            or qualities[read_pos + 1] < self.config.min_flank_quality
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
        qualities = read.query_qualities
        if qualities is None:
            return False
        if qualities[read_pos] < self.config.min_base_quality:
            return False
        seq = read.query_sequence
        if seq is None:
            return False
        if read_pos >= len(seq):
            return False
        if ref_pos >= len(contig_seq):
            return False
        read_base = seq[read_pos]
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

    def _span_has_str(
        self, ref_pos: int, length: int, str_classifier: STRClassifier
    ) -> bool:
        if length <= 0:
            return False
        for offset in range(length):
            if IndelFilters.check_if_in_str(ref_pos + offset, str_classifier):
                return True
        return False

    def _span_has_homopolymer(self, contig_seq: str, ref_pos: int, length: int) -> bool:
        if length <= 0 or not contig_seq:
            return False
        min_len = self.config.min_homopolymer_len
        window_start = max(0, ref_pos - min_len)
        window_end = min(len(contig_seq), ref_pos + length + min_len)
        context_slice = contig_seq[window_start:window_end]
        interval_start = max(0, ref_pos - 1 - window_start)
        interval_end = min(len(context_slice) - 1, ref_pos + length - window_start)
        return HomopolymerClassifier.run_overlaps_interval(
            context_slice, (interval_start, interval_end), min_len=min_len
        )

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
        right_flank_idx = read_pos + length
        if read_pos <= 0 or right_flank_idx >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[right_flank_idx] < self.config.min_flank_quality
        ):
            return False
        insertion_window = qualities[read_pos : read_pos + length]
        min_quality = self._get_processor_min_quality()
        if insertion_window and min(insertion_window) < min_quality:
            return False
        if self._span_has_str(ref_pos, length, str_classifier):
            return False
        if self._span_has_homopolymer(contig_seq, ref_pos, length):
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
        right_flank_idx = read_pos + length
        if read_pos <= 0 or right_flank_idx >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[right_flank_idx] < self.config.min_flank_quality
        ):
            return False
        if self._span_has_str(ref_pos, length, str_classifier):
            return False
        if self._span_has_homopolymer(contig_seq, ref_pos, length):
            return False
        return True

    def _is_callable_for_insertion_length_str(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        contig_seq: str,
        str_classifier: STRClassifier,
        length: int,
    ) -> bool:
        motif_len = str_classifier.motif_length_at(ref_pos)
        if motif_len is None or motif_len <= 0:
            return False
        if length % motif_len != 0:
            return False
        if read_pos is None or ref_pos is None:
            return False
        qualities = read.query_qualities
        if qualities is None:
            return False
        right_flank_idx = read_pos + length
        if read_pos <= 0 or right_flank_idx >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[right_flank_idx] < self.config.min_flank_quality
        ):
            return False
        insertion_window = qualities[read_pos : read_pos + length]
        min_quality = self._get_processor_min_quality()
        if insertion_window and min(insertion_window) < min_quality:
            return False
        if not str_classifier.matches_rptrf_motif(ref_pos, contig_seq):
            return False
        if self._span_has_homopolymer(contig_seq, ref_pos, length):
            return False
        return True

    def _is_callable_for_deletion_length_str(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        contig_seq: str,
        str_classifier: STRClassifier,
        length: int,
    ) -> bool:
        motif_len = str_classifier.motif_length_at(ref_pos)
        if motif_len is None or motif_len <= 0:
            return False
        if length % motif_len != 0:
            return False
        if read_pos is None or ref_pos is None:
            return False
        qualities = read.query_qualities
        if qualities is None:
            return False
        right_flank_idx = read_pos + length
        if read_pos <= 0 or right_flank_idx >= len(qualities):
            return False
        if (
            qualities[read_pos - 1] < self.config.min_flank_quality
            or qualities[right_flank_idx] < self.config.min_flank_quality
        ):
            return False
        if not str_classifier.matches_rptrf_motif(ref_pos, contig_seq):
            return False
        if self._span_has_homopolymer(contig_seq, ref_pos, length):
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
        reads_processed = 0
        candidates: list[IndelRecord] = []
        location_counts = {}
        type_counts: dict[str, int] = {}

        for read in samfile.fetch(contig=contig_name):
            if not self._valid_read(read) or read.mapping_quality < self.config.min_map_quality:
                continue
            reads_processed += 1

            for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if ref_pos is None or read_pos is None:
                    continue
                is_callable = self._is_callable_at_position(
                    read, read_pos, ref_pos, str_classifier, contig_seq
                )
                if is_callable:
                    interrogated_bases_count += 1

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

        position_counts: dict[int, int] = defaultdict(int)
        for candidate in candidates:
            position_counts[candidate.ref_position] += 1

        target_positions = set()
        for pos in position_counts.keys():
            for offset in (-1, 0, 1):
                target_positions.add(pos + offset)

        coverage_by_pos: dict[int, int] = defaultdict(int)
        if target_positions:
            for read in samfile.fetch(contig=contig_name):
                if (
                    not self._valid_read(read)
                    or read.mapping_quality < self.config.min_map_quality
                ):
                    continue
                for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if ref_pos is None or read_pos is None:
                        continue
                    if ref_pos in target_positions:
                        coverage_by_pos[ref_pos] += 1

        filter_context = {
            "location_map": location_counts,
            "position_counts": position_counts,
            "coverage_by_pos": coverage_by_pos,
            "str_classifier": str_classifier,
        }
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
            "type_counts": type_counts,
            "elapsed_s": end_time - start_time,
        }
        return (passed_records if in_memory else None, interrogated_bases_count, stats)

    def scan_contig_callable_only(
        self,
        contig_name: str,
        samfile: pysam.AlignmentFile,
        fasta: pyfastx.Fasta,
    ) -> dict:
        str_classifier = STRClassifier(self.config, contig_name)
        contig_seq = fasta[contig_name].seq
        total_aligned_bases = 0
        sampling_bases_total = 0
        sampling_passable_by_type: dict[str, int] = {}
        sampling_totals_by_type: dict[str, int] = {}
        tract_counts_by_motif: dict[int, int] = {}
        callable_lengths = list(getattr(self.config, "callable_lengths", []))
        sampling_contigs = set(getattr(self.config, "sampling_contigs", []))
        sampling_targets = getattr(self.config, "sampling_contig_targets", {})
        contig_target = sampling_targets.get(contig_name, 0)
        sampling_active = contig_name in sampling_contigs and contig_target > 0
        contig_len = getattr(self.config, "sampling_contig_lengths", {}).get(
            contig_name, len(contig_seq)
        )
        sample_prob = 1.0 if contig_len <= 0 else min(1.0, contig_target / contig_len)
        rng = random.Random(self.config.sampling_random_seed + hash(contig_name) % 1000000)

        regions = str_classifier.iter_regions()
        region_motif_lens = [r[2] for r in regions]
        for motif_len in region_motif_lens:
            if motif_len:
                tract_counts_by_motif[motif_len] = (
                    tract_counts_by_motif.get(motif_len, 0) + 1
                )

        for read in samfile.fetch(contig=contig_name):
            if not self._valid_read(read) or read.mapping_quality < self.config.min_map_quality:
                continue
            blocks = self._aligned_blocks(read)
            aligned_len = sum(block[2] for block in blocks)
            total_aligned_bases += aligned_len
            if not sampling_active or sampling_bases_total >= contig_target:
                continue
            if aligned_len <= 0:
                continue
            expected = aligned_len * sample_prob
            sample_count = int(expected)
            if rng.random() < (expected - sample_count):
                sample_count += 1
            remaining = contig_target - sampling_bases_total
            if remaining <= 0:
                continue
            sample_count = min(sample_count, remaining, aligned_len)
            if sample_count <= 0:
                continue
            offsets = sorted(rng.sample(range(aligned_len), sample_count))
            block_idx = 0
            block_read_start, block_ref_start, block_len = blocks[block_idx]
            block_end = block_len
            for offset in offsets:
                while offset >= block_end and block_idx < len(blocks) - 1:
                    block_idx += 1
                    block_read_start, block_ref_start, block_len = blocks[block_idx]
                    block_end += block_len
                local_offset = offset - (block_end - block_len)
                read_pos = block_read_start + local_offset
                ref_pos = block_ref_start + local_offset
                sampling_bases_total += 1
                filter_cfg = getattr(self.config, "str_candidate_filter", {})
                non_str_repeat_block = False
                if filter_cfg.get("enabled", False):
                    local_window = int(filter_cfg.get("local_window_bp", 30))
                    max_motif_len = int(filter_cfg.get("max_motif_len", 6))
                    min_units = int(filter_cfg.get("min_repeat_units", 3))
                    win_start = max(0, ref_pos - local_window)
                    win_end = min(len(contig_seq), ref_pos + local_window)
                    local_seq = contig_seq[win_start:win_end]
                    non_str_repeat_block = self._has_repeat_run(
                        local_seq, min_units, max_motif_len
                    )
                    if filter_cfg.get("aggressive", False):
                        aggressive_min_units = int(
                            filter_cfg.get("aggressive_min_repeat_units", min_units)
                        )
                        aggressive_max_motif_len = int(
                            filter_cfg.get("aggressive_max_motif_len", max_motif_len)
                        )
                        non_str_repeat_block = non_str_repeat_block or self._has_repeat_run(
                            local_seq,
                            aggressive_min_units,
                            aggressive_max_motif_len,
                        )
                if self._is_callable_at_position(
                    read, read_pos, ref_pos, str_classifier, contig_seq
                ):
                    sampling_passable_by_type[self.config.snp_label] = (
                        sampling_passable_by_type.get(self.config.snp_label, 0) + 1
                    )
                for length_for_indel in callable_lengths:
                    span_in_str = self._span_has_str(
                        ref_pos, length_for_indel, str_classifier
                    )
                    if not span_in_str and non_str_repeat_block:
                        continue
                    motif_len = (
                        str_classifier.motif_length_at(ref_pos) if span_in_str else None
                    )
                    if span_in_str and motif_len and length_for_indel % motif_len == 0:
                        ins_label = f"ins_motif_{motif_len}bp"
                        del_label = f"del_motif_{motif_len}bp"
                        sampling_totals_by_type[ins_label] = (
                            sampling_totals_by_type.get(ins_label, 0) + 1
                        )
                        sampling_totals_by_type[del_label] = (
                            sampling_totals_by_type.get(del_label, 0) + 1
                        )
                        if self._is_callable_for_insertion_length_str(
                            read,
                            read_pos,
                            ref_pos,
                            contig_seq,
                            str_classifier,
                            length_for_indel,
                        ):
                            sampling_passable_by_type[ins_label] = (
                                sampling_passable_by_type.get(ins_label, 0) + 1
                            )
                        if self._is_callable_for_deletion_length_str(
                            read,
                            read_pos,
                            ref_pos,
                            contig_seq,
                            str_classifier,
                            length_for_indel,
                        ):
                            sampling_passable_by_type[del_label] = (
                                sampling_passable_by_type.get(del_label, 0) + 1
                            )
                        continue
                    if span_in_str:
                        continue
                    ins_label = f"ins_len_{length_for_indel}bp"
                    del_label = f"del_len_{length_for_indel}bp"
                    sampling_totals_by_type[ins_label] = (
                        sampling_totals_by_type.get(ins_label, 0) + 1
                    )
                    sampling_totals_by_type[del_label] = (
                        sampling_totals_by_type.get(del_label, 0) + 1
                    )
                    if self._is_callable_for_insertion_length(
                        read,
                        read_pos,
                        ref_pos,
                        contig_seq,
                        str_classifier,
                        length_for_indel,
                    ):
                        sampling_passable_by_type[ins_label] = (
                            sampling_passable_by_type.get(ins_label, 0) + 1
                        )
                    if self._is_callable_for_deletion_length(
                        read,
                        read_pos,
                        ref_pos,
                        contig_seq,
                        str_classifier,
                        length_for_indel,
                    ):
                        sampling_passable_by_type[del_label] = (
                            sampling_passable_by_type.get(del_label, 0) + 1
                        )
        return {
            "total_aligned_bases": total_aligned_bases,
            "sampling_bases_total": sampling_bases_total,
            "sampling_passable_by_type": sampling_passable_by_type,
            "sampling_totals_by_type": sampling_totals_by_type,
            "tract_counts_by_motif": tract_counts_by_motif,
            "sampled_contig": sampling_active,
        }


