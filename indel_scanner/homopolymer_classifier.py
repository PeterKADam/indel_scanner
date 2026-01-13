# indel_scanner/homopolymer_classifier.py
import logging
from .indel import Indel

logger = logging.getLogger(__name__)


class HomopolymerClassifier:
    def __init__(self, indel_obj: Indel):
        self.indel_obj = indel_obj
        self.min_len = 3  # Default min_len for Indel object analysis

    @staticmethod
    def is_position_in_homopolymer_context(
        ref_seq: str, position: int, min_hp_length: int
    ) -> bool:
        """
        A fast, static method to check if a single genomic position is inside or
        immediately adjacent to a homopolymer run. This is optimized for the
        denominator calculation and does not require an Indel object.
        """
        # Define a search window that can capture any relevant homopolymer.
        # It needs to be large enough to contain a run of `min_hp_length` that
        # could be adjacent to our position.
        window_start = max(0, position - min_hp_length)
        window_end = min(len(ref_seq), position + 1 + min_hp_length)
        context_slice = ref_seq[window_start:window_end]

        # The interval to check is the position itself and its immediate neighbors,
        # translated into the coordinate system of our `context_slice`.
        relative_pos = position - window_start
        check_interval_start = max(0, relative_pos - 1)
        check_interval_end = relative_pos + 1

        # Reuse the core logic to see if any homopolymer run overlaps this small check interval.
        return HomopolymerClassifier.run_overlaps_interval(
            context_slice,
            (check_interval_start, check_interval_end),
            min_len=min_hp_length,
        )

    def is_homopolymer(self) -> bool:
        # Fast-path: the indel itself is already a long enough run
        if (
            len(set(self.indel_obj.indel_content)) == 1
            and len(self.indel_obj.indel_content) >= self.min_len
        ):
            logger.debug(
                f"Filtering homopolymer (indel alone): {self.indel_obj.sequencecontext_brackets()}"
            )
            return True

        full_seq = self.indel_obj.sequencecontext()
        indel_start = len(self.indel_obj.prefix_context)
        indel_end = indel_start + len(self.indel_obj.indel_content) - 1

        if HomopolymerClassifier.run_overlaps_interval(
            full_seq, (indel_start, indel_end), min_len=self.min_len
        ):
            logger.debug(
                f"Filtering homopolymer (spanning indel): {self.indel_obj.sequencecontext_brackets()}"
            )
            return True

        return False

    @staticmethod
    def run_overlaps_interval(
        seq: str, interval: tuple[int, int], min_len: int = 3
    ) -> bool:
        """
        Return True if *seq* contains a stretch of the same base whose length is
        at least `min_len` **and** that stretch overlaps the `interval`.
        """
        if not seq:
            return False

        start, end = interval  # inclusive indices of the interval to check
        cur_char = seq[0]
        cur_start = 0
        cur_len = 1

        for i in range(1, len(seq)):
            if seq[i] == cur_char:
                cur_len += 1
            else:
                if cur_len >= min_len:
                    run_start, run_end = cur_start, cur_start + cur_len - 1
                    if not (run_end < start or run_start > end):
                        return True

                cur_char = seq[i]
                cur_start = i
                cur_len = 1

        if cur_len >= min_len:
            run_start, run_end = cur_start, cur_start + cur_len - 1
            if not (run_end < start or run_start > end):
                return True

        return False

    def is_adjacent_to_homopolymer(self) -> bool:
        if len(self.indel_obj.prefix_context) >= self.min_len:
            end_of_prefix = self.indel_obj.prefix_context[-self.min_len :]
            if len(set(end_of_prefix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (prefix): {end_of_prefix}"
                )
                return True

        if len(self.indel_obj.suffix_context) >= self.min_len:
            start_of_suffix = self.indel_obj.suffix_context[: self.min_len]
            if len(set(start_of_suffix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (suffix): {start_of_suffix}"
                )
                return True

        return False

    def should_filter_indel(self) -> bool:
        return self.is_adjacent_to_homopolymer() or self.is_homopolymer()
