import bisect
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
            return False
        return self._is_imperfect_repeat(position, contig_seq)

    def _is_imperfect_repeat(self, position: int, contig_seq: str) -> bool:
        if not contig_seq:
            return False
        if self.imperfect_window_bp <= 0:
            return False
        window = self._extract_window(contig_seq, position, self.imperfect_window_bp)
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
    def _extract_window(contig_seq: str, position: int, window_bp: int) -> str:
        half_window = window_bp // 2
        start = max(0, position - half_window)
        end = min(len(contig_seq), start + window_bp)
        if end - start < window_bp:
            start = max(0, end - window_bp)
        return contig_seq[start:end]
