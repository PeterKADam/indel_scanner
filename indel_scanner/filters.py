# indel_scanner/filters.py
import logging
from typing import Dict, Callable, Any, DefaultDict, List, Protocol, Optional
from .homopolymer_classifier import HomopolymerClassifier
from .STR_Classifier import STRClassifier
from .indel import INDEL_TYPE, FilterFlag


class FilterableIndel(Protocol):
    type: INDEL_TYPE
    ref_position: int
    length: int
    in_STR: bool
    filter_mask: FilterFlag
    indel_quality: Optional[List[int]]
    prefix_quality: Optional[List[int]]
    suffix_quality: Optional[List[int]]

logger = logging.getLogger(__name__)


class IndelFilters:
    # ========================================================================
    # SECTION 1: CORE LOGIC (STATIC METHODS FOR DENOMINATOR)
    # These methods are called from the scanner and operate on raw positions.
    # ========================================================================
    @staticmethod
    def check_if_in_str(position: int, str_classifier: STRClassifier) -> bool:
        return str_classifier.is_in_str(position)

    @staticmethod
    def check_if_homopolymer_context(
        ref_seq: str, position: int, min_hp_length: int
    ) -> bool:
        """
        Checks if a genomic position is in or adjacent to a homopolymer context.
        This uses the HomopolymerClassifier logic
        """
        return HomopolymerClassifier.is_position_in_homopolymer_context(
            ref_seq, position, min_hp_length
        )

    # ========================================================================
    # SECTION 2: INDEL-OBJECT FILTERS (WRAPPERS FOR PROCESSOR)
    # These methods operate on Indel objects for the 'process' stage.
    # ========================================================================
    @staticmethod
    def _is_in_str(indel: FilterableIndel, str_classifier: STRClassifier, **kwargs) -> None:
        if indel.in_STR:
            indel.filter_mask |= FilterFlag.IS_IN_STR

    @staticmethod
    def _is_homopolymer_or_adjacent(
        indel: FilterableIndel, min_homopolymer_len: int = 3, **kwargs
    ) -> None:
        classifier = HomopolymerClassifier(indel, min_len=min_homopolymer_len)
        if classifier.should_filter_indel():
            indel.filter_mask |= FilterFlag.IS_HOMOPOLYMER_CONTEXT

    @staticmethod
    def _similar_indels_in_other_reads(
        indel: FilterableIndel, location_map: DefaultDict, **kwargs
    ) -> None:
        key = (indel.ref_position, indel.type, indel.length)
        # We check for > 1 because the current indel is already in the map
        if location_map.get(key, 0) > 1:
            indel.filter_mask |= FilterFlag.SIMILAR_IN_OTHER_READS

    @staticmethod
    def _low_minimum_indel_quality(
        indel: FilterableIndel, min_quality: int = 93, **kwargs
    ) -> None:
        if indel.type != INDEL_TYPE.INSERTION or not indel.indel_quality:
            return
        if min(indel.indel_quality) < min_quality:
            indel.filter_mask |= FilterFlag.LOW_MINIMUM_INDEL_QUALITY

    @staticmethod
    def _low_singlebase_flanking_quality(
        indel: FilterableIndel, min_flank_quality: int = 93, **kwargs
    ) -> None:
        if not indel.prefix_quality or not indel.suffix_quality:
            return
        if (
            int(indel.prefix_quality[-1]) < min_flank_quality
            or int(indel.suffix_quality[0]) < min_flank_quality
        ):
            indel.filter_mask |= FilterFlag.LOW_SINGLEBASE_FLANKING_QUALITY

    # ========================================================================
    # SECTION 3: APPLY FUNCTION (UNCHANGED)
    # ========================================================================
    _all_filters: Dict[str, Callable] = {
        "is_in_str": _is_in_str,
        "is_homopolymer_context": _is_homopolymer_or_adjacent,
        "similar_indels_in_other_reads": _similar_indels_in_other_reads,
        "low_minimum_indel_quality": _low_minimum_indel_quality,
        "low_singlebase_flanking_quality": _low_singlebase_flanking_quality,
    }

    @staticmethod
    def apply(
        active_filters: List[Dict[str, Any]], indel: FilterableIndel, **context: Any
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
