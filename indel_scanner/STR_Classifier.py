import bisect
import math
from pathlib import Path
from typing import List, Tuple

from indel_scanner.configurator import PipelineConfig


class STRClassifier:
    """
    Classifies genomic positions as being inside or outside a
    Short Tandem Repeat (STR) region using binary search.

    This class is optimized for scenarios where the repeat regions are pre-sorted,
    but query positions arrive in an arbitrary (unsorted) order.
    """

    def __init__(self, config: PipelineConfig, contig: str):
        """
        Initializes the classifier by preparing lists for binary search.

        Args:
            repeat_regions: A list of (start, end) tuples, sorted by start position.
        """
        self.config = config
        self.contig = contig  # beautiful naming scheme right there m8
        self.imperfect_cfg = self.config.imperfect_str
        self.imperfect_enabled = bool(self.imperfect_cfg.get("enabled", False))
        self.imperfect_window_bp = int(self.imperfect_cfg.get("window_bp", 20))
        self.imperfect_motif_min = int(self.imperfect_cfg.get("motif_min", 2))
        self.imperfect_motif_max = int(self.imperfect_cfg.get("motif_max", 6))
        self.imperfect_max_mismatches = int(self.imperfect_cfg.get("max_mismatches", 2))
        self.imperfect_expand_bp = int(self.imperfect_cfg.get("expand_bp", 0))
        self.low_complexity_cfg = self.config.low_complexity
        self.low_complexity_enabled = bool(self.low_complexity_cfg.get("enabled", False))
        self.low_complexity_window_bp = int(self.low_complexity_cfg.get("window_bp", 32))
        self.low_complexity_entropy_threshold = float(
            self.low_complexity_cfg.get("entropy_threshold", 1.2)
        )
        self.repeat_run_cfg = self.config.repeat_run
        self.repeat_run_enabled = bool(self.repeat_run_cfg.get("enabled", False))
        self.repeat_run_min_bp = int(self.repeat_run_cfg.get("min_run_bp", 16))
        self.repeat_run_max_window_bp = int(self.repeat_run_cfg.get("max_window_bp", 80))
        self.repeat_run_max_mismatches = int(self.repeat_run_cfg.get("max_mismatches", 2))

        repeat_regions = self._read_repeat_regions(contig)

        if not repeat_regions:
            self.starts = ()
            self.ends = ()
        else:
            self.starts, self.ends = zip(*repeat_regions)

        self.num_regions = len(self.starts)

    def _read_repeat_regions(self, contig) -> List[Tuple[int, int]]:
        """
        Reads repeat regions from a file and returns a list of tuples
        (start, end, length, motif_size, motif_sequence), excluding homopolymer regions.

        Parameters:
        - file_path: Path to the file containing repeat regions.

        Returns:
        - List of tuples representing the start, end.
        """
        repeat_regions = []
        file_path: Path = self.config.str_directory / f"result-{contig}.txt"

        with open(file_path, "r") as file:
            lines = file.readlines()
            for line in lines:
                # Skip header and empty lines
                if (
                    line.startswith("*")
                    or line.strip() == ""
                    or line.startswith("Start")
                ):
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    start = int(parts[0])
                    end = int(parts[1])
                    if self.imperfect_expand_bp:
                        start = max(0, start - self.imperfect_expand_bp)
                        end = end + self.imperfect_expand_bp

                    # Extract motif info in the form of '4(ACCC)'
                    motif_info = parts[3]

                    # Split the motif info at the '(' and ')' to get the repeat number and sequence
                    _, motif_sequence = motif_info.split("(")
                    motif_sequence = motif_sequence.strip(")")

                    repeat_regions.append((start, end))

        return repeat_regions

    def is_in_str(self, position: int) -> bool:
        """
        Checks if a given genomic position falls within any STR region using binary search.

        Args:
            position: The genomic position to check.

        Returns:
            True if the position is within an STR, False otherwise.
        """
        if not self.num_regions:
            return False

        idx = bisect.bisect_right(self.starts, position)

        if idx == 0:
            # The position is before the start of all known regions.
            return False

        candidate_index = idx - 1

        candidate_start = self.starts[candidate_index]
        candidate_end = self.ends[candidate_index]

        return candidate_start <= position <= candidate_end

    def is_str_like(self, position: int, contig_seq: str) -> bool:
        if self.is_in_str(position):
            return True
        if not self.imperfect_enabled:
            if self._is_repeat_run(position, contig_seq):
                return True
            return self._is_low_complexity(position, contig_seq)
        if self._is_imperfect_repeat(position, contig_seq):
            return True
        if self._is_repeat_run(position, contig_seq):
            return True
        return self._is_low_complexity(position, contig_seq)

    def _is_imperfect_repeat(self, position: int, contig_seq: str) -> bool:
        if not contig_seq:
            return False
        if self.imperfect_window_bp <= 0:
            return False
        window = self._extract_window(contig_seq, position, self.imperfect_window_bp)
        return self._is_imperfect_repeat_window(window)

    def _is_low_complexity(self, position: int, contig_seq: str) -> bool:
        if not self.low_complexity_enabled:
            return False
        if not contig_seq:
            return False
        if self.low_complexity_window_bp <= 0:
            return False
        window = self._extract_window(contig_seq, position, self.low_complexity_window_bp)
        entropy = self._shannon_entropy(window)
        return entropy <= self.low_complexity_entropy_threshold

    def _is_repeat_run(self, position: int, contig_seq: str) -> bool:
        if not self.repeat_run_enabled:
            return False
        if not contig_seq:
            return False
        if self.repeat_run_max_window_bp <= 0:
            return False
        if self.repeat_run_min_bp <= 0:
            return False
        window = self._extract_window(contig_seq, position, self.repeat_run_max_window_bp)
        if len(window) < self.imperfect_motif_min * 2:
            return False
        center = len(window) // 2
        best_run = 0
        for motif_len in range(self.imperfect_motif_min, self.imperfect_motif_max + 1):
            if len(window) < motif_len * 2:
                continue
            for phase in range(motif_len):
                mismatches = 0
                run_left = 0
                idx = center - 1
                while idx >= 0:
                    expected = window[(idx - phase) % motif_len]
                    base = window[idx]
                    if base == "N" or base != expected:
                        mismatches += 1
                        if mismatches > self.repeat_run_max_mismatches:
                            break
                    run_left += 1
                    idx -= 1
                mismatches = 0
                run_right = 0
                idx = center
                while idx < len(window):
                    expected = window[(idx - phase) % motif_len]
                    base = window[idx]
                    if base == "N" or base != expected:
                        mismatches += 1
                        if mismatches > self.repeat_run_max_mismatches:
                            break
                    run_right += 1
                    idx += 1
                run_len = run_left + run_right
                if run_len > best_run:
                    best_run = run_len
                if run_len >= self.repeat_run_min_bp:
                    return True
        return False

    def _is_imperfect_repeat_window(self, window: str) -> bool:
        if len(window) < self.imperfect_motif_min * 2:
            return False
        for motif_len in range(self.imperfect_motif_min, self.imperfect_motif_max + 1):
            if len(window) < motif_len * 2:
                continue
            max_offset = min(motif_len, len(window) - motif_len + 1)
            for offset in range(max_offset):
                motif = window[offset : offset + motif_len]
                if "N" in motif:
                    continue
                mismatches = 0
                for idx, base in enumerate(window):
                    expected = motif[(idx - offset) % motif_len]
                    if base == "N" or base != expected:
                        mismatches += 1
                        if mismatches > self.imperfect_max_mismatches:
                            break
                if mismatches <= self.imperfect_max_mismatches:
                    return True
        return False

    @staticmethod
    def _shannon_entropy(seq: str) -> float:
        if not seq:
            return 0.0
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        total = 0
        for base in seq:
            if base in counts:
                counts[base] += 1
                total += 1
        if total == 0:
            return 0.0
        entropy = 0.0
        for count in counts.values():
            if count:
                p = count / total
                entropy -= p * math.log2(p)
        return entropy

    @staticmethod
    def _extract_window(contig_seq: str, position: int, window_bp: int) -> str:
        half_window = window_bp // 2
        start = max(0, position - half_window)
        end = min(len(contig_seq), start + window_bp)
        if end - start < window_bp:
            start = max(0, end - window_bp)
        return contig_seq[start:end]
