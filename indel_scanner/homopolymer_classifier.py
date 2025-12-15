import logging

from .indel import Indel

logger = logging.getLogger(__name__)

class HomopolymerClassifier():

    def __init__(self, indel_obj: Indel):
        self.indel_obj = indel_obj

    def is_homopolymer(self) -> bool:
        # Fast‑path: the indel itself is already a long enough run
        if len(set(self.indel_obj.indel_content)) == 1 and len(self.indel_obj.indel_content) >= 3:
            logger.debug(
                f"Filtering homopolymer (indel alone): {self.indel_obj.sequencecontext_brackets()}"
            )
            return True

        # Build the full context and compute the indel interval (inclusive)
        full_seq = self.indel_obj.sequencecontext()
        indel_start = len(self.indel_obj.prefix_context)  # first base of the indel
        indel_end = indel_start + len(self.indel_obj.indel_content) - 1  # last base of the indel

        if self.run_overlaps_interval(full_seq, (indel_start, indel_end)):
            logger.debug(
                f"Filtering homopolymer (spanning indel): {self.indel_obj.sequencecontext_brackets()}"
            )
            return True

        return False

    @staticmethod
    def run_overlaps_interval(
        seq: str,  # type: ignore
        interval: tuple[int, int],  # (start, end) inclusive, 0‑based
        min_len: int = 3,
    ) -> bool:
        """
        Return True if *seq* contains a stretch of the same base whose length is
        at least ``min_len`` **and** that stretch overlaps the interval
        ``(start, end)``.  The interval normally corresponds to the indel
        positions inside the concatenated string.
        """
        if not seq:
            return False

        start, end = interval  # inclusive indices of the indel
        cur_char = seq[0]
        cur_start = 0
        cur_len = 1

        for i in range(1, len(seq)):
            if seq[i] == cur_char:
                cur_len += 1
            else:
                # finish the previous run
                if cur_len >= min_len:
                    run_start, run_end = cur_start, cur_start + cur_len - 1
                    # does the run intersect the indel interval?
                    if not (run_end < start or run_start > end):
                        return True
                # start a new run
                cur_char = seq[i]
                cur_start = i
                cur_len = 1

        # check the very last run
        if cur_len >= min_len:
            run_start, run_end = cur_start, cur_start + cur_len - 1
            if not (run_end < start or run_start > end):
                return True

        return False

    def is_adjacent_to_homopolymer(self) -> bool:
        if len(self.indel_obj.prefix_context) >= 3:
            end_of_prefix = self.indel_obj.prefix_context[-3:]
            if len(set(end_of_prefix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (prefix): {end_of_prefix}"
                )
                return True
        if len(self.indel_obj.suffix_context) >= 3:
            start_of_suffix = self.indel_obj.suffix_context[:3]
            if len(set(start_of_suffix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (suffix): {start_of_suffix}"
                )
                return True
        return False

    def should_filter_indel(self) -> bool:
        """
        Return **True** if the indel must be removed because:
        • it participates in a homopolymer run (overlap), **or**
        • it sits directly next to a homopolymer of length ≥ min_len.

        """

        return self.is_adjacent_to_homopolymer() or self.is_homopolymer()