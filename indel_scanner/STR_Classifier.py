import bisect
from pathlib import Path
from typing import List, Tuple

from indel_scanner.indel_scanner.configurator import ScannerConfig


class STRClassifier:
    """
    Efficiently classifies genomic positions as being inside or outside a
    Short Tandem Repeat (STR) region using binary search.

    This class is optimized for scenarios where the repeat regions are pre-sorted,
    but query positions arrive in an arbitrary (unsorted) order.
    """

    def __init__(self, config: ScannerConfig, contig: str):
        """
        Initializes the classifier by preparing lists for binary search.

        Args:
            repeat_regions: A list of (start, end) tuples, sorted by start position.
        """
        self.config = config
        self.contig = contig  # beautiful naming scheme right here m8

        repeat_regions = self._read_repeat_regions(contig)

        if not repeat_regions:
            # Handle case with no repeat regions
            self.starts = []
            self.ends = []
        else:
            # Separate starts and ends for efficient lookup
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

                    # Extract motif info in the form of '4(ACCC)'
                    motif_info = parts[3]

                    # Split the motif info at the '(' and ')' to get the repeat number and sequence
                    _, motif_sequence = motif_info.split("(")
                    motif_sequence = motif_sequence.strip(")")
                    # The number of repeats (e.g., 4)

                    # Ensure the motif is not a homopolymer
                    if len(set(motif_sequence)) == 1:
                        continue

                    repeat_regions.append((start, end))

        return repeat_regions

    def is_in_str(self, position: int) -> bool:
        """
        Checks if a given genomic position falls within any STR region using binary search.
        This operation is very fast (O(log N)).

        Args:
            position: The genomic position to check.

        Returns:
            True if the position is within an STR, False otherwise.
        """
        if not self.num_regions:
            return False

        # Find the insertion point for 'position' in the sorted list of start coordinates.
        # bisect_right finds an insertion point which comes after (to the right of)
        # any existing entries of 'position'.
        # This tells us the index of the first region that *might* start after our position.
        idx = bisect.bisect_right(self.starts, position)

        # The candidate region that could contain our position must be the one
        # immediately to the left of our insertion point.
        # So we check the region at index `idx - 1`.

        if idx == 0:
            # The position is before the start of all known regions.
            return False

        candidate_index = idx - 1

        # Now we only need to perform one check:
        # Is the position within the bounds of this single candidate region?
        candidate_start = self.starts[candidate_index]
        candidate_end = self.ends[candidate_index]

        return candidate_start <= position <= candidate_end
