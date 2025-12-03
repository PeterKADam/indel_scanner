from dataclasses import dataclass
from enum import Enum, StrEnum
from abc import ABC, abstractmethod
import logging
from typing import List, Optional

logger = logging.getLogger(__name__)


def _format_quality_scores(scores: Optional[List[int]]) -> str:
    """Converts a list of integers into a comma-separated string."""
    if scores is None:
        return "NA"  # Use "NA" to signify missing data
    return ",".join(map(str, scores))


class INDEL_TYPE(StrEnum):
    INSERTION = "ins"
    DELETION = "del"


class TSV_HEADERS(Enum):
    SCANNER = [
        "contig",
        "ref_position",
        "type",
        "length",
        "[sequence]_context",
        "read_name",
        "in_STR",
    ]
    PROCESSOR = SCANNER + [
        "prefix_quality",
        "insertion_quality",
        "suffix_quality",
        "map_quality",
    ]


@dataclass
class Indel(ABC):
    contig: str
    ref_position: int
    length: int
    prefix_context: str
    indel_content: str
    suffix_context: str
    read_name: str
    type: INDEL_TYPE
    in_STR: bool
    map_quality: Optional[int] = None

    def tsv_sequence_field(self) -> str:
        return self.indel_content

    def to_scanner_tsv_row(self) -> List[str]:
        return [
            self.contig,
            str(self.ref_position),
            self.type.value,
            str(self.length),
            self.sequencecontext(),
            self.read_name,
            str(self.in_STR),
        ]

    def sequencecontext(self) -> str:
        return f"{self.prefix_context}[{self.tsv_sequence_field}]{self.suffix_context}"

    @abstractmethod
    def to_processor_tsv_row(self) -> List[str]:
        pass

    def _is_homopolymer(self) -> bool:
        # Fast‑path: the indel itself is already a long enough run
        if len(set(self.indel_content)) == 1 and len(self.indel_content) >= 3:
            logger.debug(
                f"Filtering homopolymer (indel alone): {self.prefix_context}[{self.indel_content}]{self.suffix_context}"
            )
            return True

        # Build the full context and compute the indel interval (inclusive)
        full_seq = f"{self.prefix_context}{self.indel_content}{self.indel_content}"
        indel_start = len(self.prefix_context)  # first base of the indel
        indel_end = indel_start + len(self.indel_content) - 1  # last base of the indel

        if self.run_overlaps_interval(full_seq, (indel_start, indel_end)):
            logger.debug(
                f"Filtering homopolymer (spanning indel): {self.prefix_context}[{self.indel_content}]{self.suffix_context}"
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

    def _is_adjacent_to_homopolymer(self) -> bool:
        if len(self.prefix_context) >= 3:
            end_of_prefix = self.prefix_context[-3:]
            if len(set(end_of_prefix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (prefix): {end_of_prefix}"
                )
                return True
        if len(self.suffix_context) >= 3:
            start_of_suffix = self.suffix_context[:3]
            if len(set(start_of_suffix)) == 1:
                logger.debug(
                    f"Filtering adjacent homopolymer (suffix): {start_of_suffix}"
                )
                return True
        return False

    def _should_filter_indel(self) -> bool:
        """
        Return **True** if the indel must be removed because:
        • it participates in a homopolymer run (overlap), **or**
        • it sits directly next to a homopolymer of length ≥ min_len.

        """

        return self._is_adjacent_to_homopolymer() or self._is_homopolymer()


@dataclass
class Insertion(Indel):
    indel_content: str  # type: ignore

    insertion_quality: Optional[List[int]] = None
    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    def to_processor_tsv_row(self) -> List[str]:
        base_row = self.to_scanner_tsv_row()
        base_row.extend(
            [
                _format_quality_scores(self.prefix_quality),
                _format_quality_scores(self.insertion_quality),
                _format_quality_scores(self.suffix_quality),
                str(self.map_quality),  # should never be None, #pray
            ]
        )
        return base_row


@dataclass
class Deletion(Indel):
    indel_content: str  # type: ignore

    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    def to_processor_tsv_row(self) -> List[str]:
        base_row = self.to_scanner_tsv_row()
        base_row.extend(
            [
                _format_quality_scores(self.prefix_quality),
                "NA",  # Placeholder for insertion_quality
                _format_quality_scores(self.suffix_quality),
                str(self.map_quality),
            ]
        )
        return base_row
