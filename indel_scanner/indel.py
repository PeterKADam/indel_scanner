from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, StrEnum
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


@dataclass #if you squint, this is still a dataclass
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

    @classmethod
    def create(
            cls,
            type: INDEL_TYPE,
            contig: str,
            ref_position: int,
            length: int,
            prefix_context: str,
            indel_content: str,
            suffix_context: str,
            read_name: str,
            in_STR: bool,
            map_quality: Optional[int],
            # Subclass-specific arguments
            insertion_quality: Optional[List[int]] = None,
            prefix_quality: Optional[List[int]] = None,
            suffix_quality: Optional[List[int]] = None,
    ) -> Indel:
        """
        Factory method to create an Insertion or Deletion object.
        """
        common_args = {
            "contig": contig,
            "ref_position": ref_position,
            "length": length,
            "prefix_context": prefix_context,
            "indel_content": indel_content,
            "suffix_context": suffix_context,
            "read_name": read_name,
            "in_STR": in_STR,
            "map_quality": map_quality,
        }

        if type == INDEL_TYPE.INSERTION:
            return Insertion(
                type=INDEL_TYPE.INSERTION,
                **common_args,
                insertion_quality=insertion_quality,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality
            )
        elif type == INDEL_TYPE.DELETION:
            return Deletion(
                type=INDEL_TYPE.DELETION,
                **common_args,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality
            )
        else:
            raise ValueError(f"Unknown INDEL_TYPE: {type}")

    def to_scanner_tsv_row(self) -> List[str]:
        return [
            self.contig,
            str(self.ref_position),
            self.type.value,
            str(self.length),
            self.sequencecontext_brackets(),
            self.read_name,
            str(self.in_STR),
        ]

    def sequencecontext_brackets(self) -> str:
        return f"{self.prefix_context}[{self.indel_content}]{self.suffix_context}"

    def sequencecontext(self) -> str:
        return f"{self.prefix_context}{self.indel_content}{self.suffix_context}"

    @abstractmethod
    def to_processor_tsv_row(self) -> List[str]:
        pass

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

