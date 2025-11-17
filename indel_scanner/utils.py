
import logging
from pathlib import Path
import pysam
from enum import Enum, IntEnum


logger = logging.getLogger(__name__)



class Cigar(IntEnum):
	OP_I = pysam.CINS
	OP_D = pysam.CDEL
	OP_M = pysam.CMATCH
	OP_EQ = pysam.CEQUAL
	OP_X = pysam.CDIFF
	OP_N = pysam.CREF_SKIP
	OP_S = pysam.CSOFT_CLIP
	




		