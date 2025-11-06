import csv
import os
import sys
import logging
import shutil
from pathlib import Path
import pysam
from enum import Enum
from typing import Any, Dict, Optional, Tuple, Self

import yaml

class Cigar(Enum):
	OP_I = pysam.CINS
	OP_D = pysam.CDEL
	OP_M = pysam.CMATCH
	OP_EQ = pysam.CEQUAL
	OP_X = pysam.CDIFF
	OP_N = pysam.CREF_SKIP
	OP_S = pysam.CSOFT_CLIP
	

logger = logging.getLogger(__name__)



def cleanup_temp_dir(temp_dir):
	logger.debug(f"Cleaning up temporary files in {temp_dir}")	
	shutil.rmtree(temp_dir)



		