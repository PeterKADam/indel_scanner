import logging
from functools import partial
from typing import List, Dict, Callable, Any, DefaultDict

# Assuming these are defined elsewhere
from .homopolymer_classifier import HomopolymerClassifier
from .indel import Indel

logger = logging.getLogger(__name__)


class IndelFilters:
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
            indel.filter_reason.append(f"poor_mapping_quality_<{min_mapq}")

    @staticmethod
    def _similar_indels_in_other_reads(
        indel: Indel, location_map: DefaultDict, **kwargs
    ) -> None:
        """Tags indels that have identical counterparts in other reads."""
        key = (indel.ref_position, indel.type, indel.length)

        if location_map.get(key, 0) > 1:
            indel.filter_reason.append("similar_indels_in_other_reads")

    @staticmethod
    def apply(
        active_filters: List[Dict[str, Any]], indel: Indel, **context: Any
    ) -> None:
        for config in active_filters:
            filter_name = config.get("name")
            params = config.get("params", {})
            if filter_name in IndelFilters._all_filters:
                filter_func = IndelFilters._all_filters[filter_name]
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
}
