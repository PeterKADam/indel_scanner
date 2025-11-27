from dataclasses import dataclass
from enum import Enum, StrEnum
from abc import ABC, abstractmethod
from typing import List, Optional


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
    suffix_context: str
    read_name: str
    type: INDEL_TYPE
    in_STR: bool
    map_quality: Optional[int] = None

    @property
    @abstractmethod
    def tsv_sequence_field(self) -> str:
        pass

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


@dataclass
class Insertion(Indel):
    inserted_seq: str  # type: ignore

    insertion_quality: Optional[List[int]] = None
    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    @property
    def tsv_sequence_field(self) -> str:
        return self.inserted_seq

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
    reference_seq: str  # type: ignore

    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    @property
    def tsv_sequence_field(self) -> str:
        return self.reference_seq

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
