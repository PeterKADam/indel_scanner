import logging
import time
import pysam
from enum import IntEnum


logger = logging.getLogger(__name__)


class Cigar(IntEnum):
    OP_I = pysam.CINS
    OP_D = pysam.CDEL
    OP_M = pysam.CMATCH
    OP_EQ = pysam.CEQUAL
    OP_X = pysam.CDIFF
    OP_N = pysam.CREF_SKIP
    OP_S = pysam.CSOFT_CLIP


def as_cigar(op: int) -> Cigar:
    try:
        return Cigar(op)
    except ValueError:
        logger.error(f"Unrecognized CIGAR operation: {op}")
        return None  # type: ignore


class Timer:
    def __init__(self, label: str, log_func=logger.info) -> None:
        self.label = label
        self.log_func = log_func
        self.start = 0.0

    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        duration = time.perf_counter() - self.start
        self.log_func(f"{self.label} completed in {duration:.2f}s")


def build_refpos_to_read_index(read: pysam.AlignedSegment) -> dict[int, int]:
    mapping: dict[int, int] = {}
    # Use aligned pairs to avoid insertion-adjacent index drift from
    # get_reference_positions(full_length=True) in some pysam versions.
    for read_idx, ref_pos in read.get_aligned_pairs(matches_only=True):
        if read_idx is None:
            continue
        if ref_pos is None:
            continue
        if ref_pos not in mapping:
            mapping[ref_pos] = read_idx
    return mapping


def flank_qualities_by_ref(
    read: pysam.AlignedSegment,
    ref_pos: int,
    length: int,
    flank_len: int,
    is_insertion: bool,
) -> tuple[list[int | None], list[int | None]]:
    qualities = read.query_qualities or []
    if flank_len <= 0:
        return [], []
    mapping = build_refpos_to_read_index(read)
    left_positions = range(ref_pos - flank_len, ref_pos)
    right_start = ref_pos if is_insertion else ref_pos + length
    right_positions = range(right_start, right_start + flank_len)

    prefix_quality: list[int | None] = []
    for pos in left_positions:
        read_idx = mapping.get(pos)
        if read_idx is None or not (0 <= read_idx < len(qualities)):
            prefix_quality.append(None)
            continue
        prefix_quality.append(int(qualities[read_idx]))

    suffix_quality: list[int | None] = []
    for pos in right_positions:
        read_idx = mapping.get(pos)
        if read_idx is None or not (0 <= read_idx < len(qualities)):
            suffix_quality.append(None)
            continue
        suffix_quality.append(int(qualities[read_idx]))

    return prefix_quality, suffix_quality
