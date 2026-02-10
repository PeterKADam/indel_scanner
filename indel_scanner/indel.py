from __future__ import annotations
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum, StrEnum, IntFlag
from typing import List, Optional, Dict, Any

logger = logging.getLogger(__name__)

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
        "in_STR_region",
        "in_STR",
        "str_motif_length",
        "filter_reason",
    ]
    PROCESSOR = SCANNER + [
        "prefix_quality",
        "insertion_quality",
        "suffix_quality",
        "map_quality",
    ]


class FilterFlag(IntFlag):
    NONE = 0
    IS_IN_STR = 1 << 0
    IS_HOMOPOLYMER_CONTEXT = 1 << 1
    SIMILAR_IN_OTHER_READS = 1 << 2
    LOW_MINIMUM_INDEL_QUALITY = 1 << 3
    LOW_SINGLEBASE_FLANKING_QUALITY = 1 << 4


@dataclass
class IndelRecord:
    contig: str
    ref_position: int
    length: int
    prefix_context: str
    indel_content: str
    suffix_context: str
    read_name: str
    in_STR_region: bool
    in_STR: bool
    map_quality: int
    type: INDEL_TYPE
    filter_mask: FilterFlag = FilterFlag.NONE
    motif_length: Optional[int] = None

    def is_filtered(self) -> bool:
        return self.filter_mask != FilterFlag.NONE
    indel_quality: Optional[List[int]] = None
    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    def sequencecontext_brackets(self) -> str:
        """Returns the sequence context with brackets around the indel content."""
        return f"{self.prefix_context}[{self.indel_content}]{self.suffix_context}"

    def sequencecontext(self) -> str:
        """Returns the full sequence context without brackets."""
        return f"{self.prefix_context}{self.indel_content}{self.suffix_context}"


@dataclass
class Indel(IndelRecord, ABC):
    """
    An abstract factory class representing a generic Insertion or Deletion.
    """

    @classmethod
    def create(
        cls,
        type: INDEL_TYPE,
        *,
        contig: str,
        ref_position: int,
        length: int,
        prefix_context: str,
        indel_content: str,
        suffix_context: str,
        read_name: str,
        in_STR_region: bool,
        in_STR: bool,
        motif_length: Optional[int] = None,
        map_quality: Optional[int],
        # Subclass-specific arguments
        indel_quality: Optional[List[int]] = None,
        prefix_quality: Optional[List[int]] = None,
        suffix_quality: Optional[List[int]] = None,
    ) -> Indel:
        """Factory method to create an Insertion or Deletion object."""
        common_args = {
            "contig": contig,
            "ref_position": ref_position,
            "length": length,
            "prefix_context": prefix_context,
            "indel_content": indel_content,
            "suffix_context": suffix_context,
            "read_name": read_name,
            "in_STR_region": in_STR_region,
            "in_STR": in_STR,
            "motif_length": motif_length,
            "map_quality": map_quality,
        }
        if type == INDEL_TYPE.INSERTION:
            return Insertion(
                type=INDEL_TYPE.INSERTION,
                **common_args,
                indel_quality=indel_quality,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality,
            )
        elif type == INDEL_TYPE.DELETION:
            return Deletion(
                type=INDEL_TYPE.DELETION,
                **common_args,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality,
            )
        else:
            raise ValueError(f"Unknown INDEL_TYPE: {type}")

    def sequencecontext_brackets(self) -> str:
        """Returns the sequence context with brackets around the indel content."""
        return f"{self.prefix_context}[{self.indel_content}]{self.suffix_context}"

    def sequencecontext(self) -> str:
        """Returns the full sequence context without brackets."""
        return f"{self.prefix_context}{self.indel_content}{self.suffix_context}"

    @abstractmethod
    def to_processor_dict(self) -> Dict[str, Any]:
        """
        Converts the object's data into a dictionary for converting back to a Polars DataFrame.
        """
        raise NotImplementedError


@dataclass
class Insertion(Indel):
    type: INDEL_TYPE = INDEL_TYPE.INSERTION
    indel_quality: Optional[List[int]] = None
    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    def to_processor_dict(self) -> Dict[str, Any]:
        return {
            "contig": self.contig,
            "ref_position": self.ref_position,
            "type": self.type.value,
            "length": self.length,
            "sequence_context_brackets": self.sequencecontext_brackets(),
            "read_name": self.read_name,
            "in_STR_region": self.in_STR_region,
            "in_STR": self.in_STR,
            "str_motif_length": self.motif_length,
            "prefix_quality": self.prefix_quality,
            "insertion_quality": self.indel_quality,
            "suffix_quality": self.suffix_quality,
            "map_quality": self.map_quality,
            "filter_reason": "",
        }


@dataclass
class Deletion(Indel):
    type: INDEL_TYPE = INDEL_TYPE.DELETION
    prefix_quality: Optional[List[int]] = None
    suffix_quality: Optional[List[int]] = None

    def to_processor_dict(self) -> Dict[str, Any]:
        return {
            "contig": self.contig,
            "ref_position": self.ref_position,
            "type": self.type.value,
            "length": self.length,
            "sequence_context_brackets": self.sequencecontext_brackets(),
            "read_name": self.read_name,
            "in_STR_region": self.in_STR_region,
            "in_STR": self.in_STR,
            "str_motif_length": self.motif_length,
            "prefix_quality": self.prefix_quality,
            "insertion_quality": None,  # Deletions explicitly have no insertion quality
            "suffix_quality": self.suffix_quality,
            "map_quality": self.map_quality,
            "filter_reason": "",
        }


class IndelTsvFormatter:
    """
    Handles formatting of Indel objects into string-based TSV rows.
    """

    @staticmethod
    def _format_quality_scores(scores: Optional[List[int]]) -> str:
        """Converts a list of integers into a comma-separated string."""
        return ",".join(map(str, scores)) if scores is not None else "NA"

    @staticmethod
    def format_scanner_row(indel: Indel) -> List[str]:
        """Generates a TSV row for the 'scanner' output format."""
        return [
            indel.contig,
            str(indel.ref_position + 1),
            indel.type.value,
            str(indel.length),
            indel.sequencecontext_brackets(),
            indel.read_name,
            str(indel.in_STR_region),
            str(indel.in_STR),
            str(indel.motif_length) if indel.motif_length is not None else "NA",
            "",
        ]

    @staticmethod
    def format_processor_row(indel: Indel) -> List[str]:
        """
        Generates a TSV row for the 'processor' output format
        """
        base_row = IndelTsvFormatter.format_scanner_row(indel)

        if isinstance(indel, Insertion):
            base_row.extend(
                [
                    IndelTsvFormatter._format_quality_scores(indel.prefix_quality),
                    IndelTsvFormatter._format_quality_scores(indel.indel_quality),
                    IndelTsvFormatter._format_quality_scores(indel.suffix_quality),
                    str(indel.map_quality),
                ]
            )
        elif isinstance(indel, Deletion):
            base_row.extend(
                [
                    IndelTsvFormatter._format_quality_scores(indel.prefix_quality),
                    "NA",  # Placeholder for insertion_quality
                    IndelTsvFormatter._format_quality_scores(indel.suffix_quality),
                    str(indel.map_quality),
                ]
            )
        else:
            raise TypeError(f"Unsupported Indel type for formatting: {type(indel)}")

        return base_row
