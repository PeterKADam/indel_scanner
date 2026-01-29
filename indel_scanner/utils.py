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
