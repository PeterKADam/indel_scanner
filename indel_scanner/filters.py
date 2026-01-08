# indel_scanner/filters.py
import logging
from typing import Dict, Callable, Any, DefaultDict, List
from .homopolymer_classifier import HomopolymerClassifier
from .STR_Classifier import STRClassifier
from .indel import INDEL_TYPE, Indel

logger = logging.getLogger(__name__)


class IndelFilters:

    @staticmethod
    def check_if_in_str(position: int, str_classifier: STRClassifier) -> bool:
        """Checks if a genomic position is within a classified STR."""
        return str_classifier.is_in_str(position)

    @staticmethod
    def check_if_homopolymer_context(ref_seq: str, position: int) -> bool:
        """
        Checks if a genomic position is in or adjacent to a homopolymer context.
        This reuses the HomopolymerClassifier logic for a single source of truth.
        """
        return HomopolymerClassifier.is_position_in_homopolymer_context(ref_seq, position)

    @staticmethod
    def _is_in_str(indel: Indel, **kwargs) -> None:
        if indel.in_STR:
            indel.filter_reason.append("is_in_str")

    @staticmethod
    def _is_adjacent_to_homopolymer(indel: Indel, **kwargs) -> None:
        if HomopolymerClassifier(indel).is_adjacent_to_homopolymer():
            indel.filter_reason.append("is_adjacent_to_homopolymer")

    @staticmethod
    def _is_homopolymer(indel: Indel, **kwargs) -> None:
        if HomopolymerClassifier(indel).is_homopolymer():
            indel.filter_reason.append("is_homopolymer")

    @staticmethod
    def _poor_mapping_quality(indel: Indel, min_mapq: int = 30, **kwargs) -> None:
        if indel.map_quality <= min_mapq:
            indel.filter_reason.append("poor_mapping_quality")

    @staticmethod
    def _similar_indels_in_other_reads(
        indel: Indel, location_map: DefaultDict, **kwargs
    ) -> None:
        """Tags indels that have identical counterparts in other reads."""
        key = (indel.ref_position, indel.type, indel.length)

        if location_map.get(key, 0) > 1:
            indel.filter_reason.append("similar_indels_in_other_reads")

    @staticmethod
    def _low_minimum_indel_quality(
        indel: Indel, min_quality: int = 93, **kwargs
    ) -> None:
        if indel.type != INDEL_TYPE.INSERTION:
            return
        if min(indel.indel_quality) < min_quality:
            indel.filter_reason.append("low_minimum_indel_quality")

    @staticmethod
    def _low_minimum_flanking_quality(
        indel: Indel, min_flank_quality: int = 93, **kwargs
    ) -> None:
        if (
            min(indel.prefix_quality) < min_flank_quality
            or min(indel.suffix_quality) < min_flank_quality
        ):
            indel.filter_reason.append("low_minimum_flanking_quality")

    @staticmethod
    def _low_singlebase_flanking_quality(
        indel: Indel, min_flank_quality: int = 93, **kwargs
    ) -> None:
        if (int(indel.prefix_quality[-1]) < min_flank_quality) or (
            int(indel.suffix_quality[0]) < min_flank_quality
        ):
            indel.filter_reason.append("low_singlebase_flanking_quality")

    @staticmethod
    def apply(
        active_filters: List[Dict[str, Any]], indel: Indel, **context: Any
    ) -> None:
        for config in active_filters:
            filter_name = config.get("name")
            params = config.get("params", {})
            if filter_name in IndelFilters_all_filters:
                filter_func = IndelFilters_all_filters[filter_name]
                try:
                    filter_func(indel, **params, **context)
                except TypeError as e:
                    logger.error(
                        f"Error calling filter '{filter_name}'. Check parameters. Error: {e}"
                    )
            else:
                logger.warning(f"Filter '{filter_name}' not found. Skipping.")


IndelFilters_all_filters: Dict[str, Callable] = {
    "is_in_str": IndelFilters._is_in_str,
    "is_adjacent_to_homopolymer": IndelFilters._is_adjacent_to_homopolymer,
    "is_homopolymer": IndelFilters._is_homopolymer,
    "poor_mapping_quality": IndelFilters._poor_mapping_quality,
    "similar_indels_in_other_reads": IndelFilters._similar_indels_in_other_reads,
    "low_minimum_indel_quality": IndelFilters._low_minimum_indel_quality,
    "low_minimum_flanking_quality": IndelFilters._low_minimum_flanking_quality,
    "low_singlebase_flanking_quality": IndelFilters._low_singlebase_flanking_quality,
}
