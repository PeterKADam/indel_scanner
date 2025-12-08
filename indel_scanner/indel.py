from __future__ import annotations
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, StrEnum
from typing import List, Optional, Dict, Any

logger = logging.getLogger(__name__)

# --- Core Data Definitions ---

class INDEL_TYPE(StrEnum):
    INSERTION = "ins"
    DELETION = "del"

class TSV_HEADERS(Enum):
    SCANNER = [
        "contig", "ref_position", "type", "length",
        "[sequence]_context", "read_name", "in_STR",
    ]
    PROCESSOR = SCANNER + [
        "prefix_quality", "insertion_quality", "suffix_quality", "map_quality",
    ]

# --- Pure Data Classes for Indels ---

@dataclass
class Indel(ABC):
    """
    An abstract base class representing a generic Insertion or Deletion.
    This class is a pure data container. All formatting logic is externalized.
    """
    # Core attributes
    contig: str
    ref_position: int
    length: int
    prefix_context: str
    indel_content: str
    suffix_context: str
    read_name: str
    type: INDEL_TYPE
    in_STR: bool

    # Attributes populated from BAM file
    map_quality: Optional[int] = None

    @classmethod
    def create(
        cls,
        type: INDEL_TYPE, *,
        contig: str, ref_position: int, length: int,
        prefix_context: str, indel_content: str, suffix_context: str,
        read_name: str, in_STR: bool, map_quality: Optional[int],
        # Subclass-specific arguments
        insertion_quality: Optional[List[int]] = None,
        prefix_quality: Optional[List[int]] = None,
        suffix_quality: Optional[List[int]] = None,
    ) -> Indel:
        """Factory method to create an Insertion or Deletion object."""
        common_args = {
            "contig": contig, "ref_position": ref_position, "length": length,
            "prefix_context": prefix_context, "indel_content": indel_content,
            "suffix_context": suffix_context, "read_name": read_name,
            "in_STR": in_STR, "map_quality": map_quality,
        }
        if type == INDEL_TYPE.INSERTION:
            return Insertion(
                type=INDEL_TYPE.INSERTION, **common_args,
                insertion_quality=insertion_quality,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality
            )
        elif type == INDEL_TYPE.DELETION:
            return Deletion(
                type=INDEL_TYPE.DELETION, **common_args,
                prefix_quality=prefix_quality,
                suffix_quality=suffix_quality
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
        Converts the object's data into a dictionary.
        This is the primary mechanism for converting back to a Polars DataFrame.
        """
        raise NotImplementedError

@dataclass
class Insertion(Indel):
    type: INDEL_TYPE = INDEL_TYPE.INSERTION
    insertion_quality: Optional[List[int]] = None
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
            "in_STR": self.in_STR,
            "prefix_quality": self.prefix_quality,
            "insertion_quality": self.insertion_quality,
            "suffix_quality": self.suffix_quality,
            "map_quality": self.map_quality,
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
            "in_STR": self.in_STR,
            "prefix_quality": self.prefix_quality,
            "insertion_quality": None,  # Deletions explicitly have no insertion quality
            "suffix_quality": self.suffix_quality,
            "map_quality": self.map_quality,
        }

# --- Externalized Formatter Class ---

class IndelTsvFormatter:
    """
    Handles formatting of Indel objects into string-based TSV rows.
    Note: This class is useful for legacy code or debugging that needs to
    generate single rows, but the high-performance path uses the
    `.to_processor_dict()` method for bulk Polars conversion.
    """
    @staticmethod
    def _format_quality_scores(scores: Optional[List[int]]) -> str:
        """Converts a list of integers into a comma-separated string."""
        return ",".join(map(str, scores)) if scores is not None else "NA"

    def format_scanner_row(self, indel: Indel) -> List[str]:
        """Generates a TSV row for the 'scanner' output format."""
        return [
            indel.contig,
            str(indel.ref_position),
            indel.type.value,
            str(indel.length),
            indel.sequencecontext_brackets(),
            indel.read_name,
            str(indel.in_STR),
        ]

    def format_processor_row(self, indel: Indel) -> List[str]:
        """
        Generates a TSV row for the 'processor' output format,
        handling both Insertion and Deletion types.
        """
        base_row = self.format_scanner_row(indel)
        if isinstance(indel, Insertion):
            base_row.extend([
                self._format_quality_scores(indel.prefix_quality),
                self._format_quality_scores(indel.insertion_quality),
                self._format_quality_scores(indel.suffix_quality),
                str(indel.map_quality),
            ])
        elif isinstance(indel, Deletion):
            base_row.extend([
                self._format_quality_scores(indel.prefix_quality),
                "NA",  # Placeholder for insertion_quality
                self._format_quality_scores(indel.suffix_quality),
                str(indel.map_quality),
            ])
        else:
            raise TypeError(f"Unsupported Indel type for formatting: {type(indel)}")
        return base_row

