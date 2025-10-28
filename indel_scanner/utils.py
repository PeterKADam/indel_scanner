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

def setup_output(args) -> tuple[Path, Path]:
	"""
	Sets up the output and temporary directories.
	If the temporary directory already exists from a previous run, it will be cleaned.
	"""
	output_dir = Path(args.output)
	temp_dir = output_dir / "temp_indel_parts"

	try:
		# 1. Ensure the main output directory exists.
		output_dir.mkdir(parents=True, exist_ok=True)
		logger.debug(f"Output directory set up at {output_dir}")

		# 2. Check if the temporary directory exists and clean it if necessary.
		if temp_dir.exists():
			logger.info(f"Found existing temporary directory at {temp_dir}. Cleaning it before use.")
			shutil.rmtree(temp_dir)

		# 3. Create a fresh, empty temporary directory.
		temp_dir.mkdir(parents=True, exist_ok=True)
		logger.debug(f"Temporary directory for partial results is ready at {temp_dir}")

	except Exception as e:
		logger.error(f"Error setting up output directories: {e}")
		sys.exit(1)
		
	return output_dir, temp_dir

def aggregate_partial_results(temp_dir, final_output_path):
	logger.info("\nAll contigs processed. Merging results...")
	logger.debug(f"Aggregating partial results from {temp_dir} into {final_output_path}")
	
	with open(final_output_path, 'w', newline='') as f_out:
		writer = csv.writer(f_out, delimiter='\t')
		
		# header
		writer.writerow(['Contig', 'Position', 'Type', 'Length',  'Sequence', 'Read_Name'])
		
		logger.debug(f"Writing ({len(os.listdir(temp_dir))}) partial results  to final output file...")
		for part_file in os.listdir(temp_dir):
			if part_file.endswith('.part.tsv'):
				part_path = os.path.join(temp_dir, part_file)
				with open(part_path, 'r') as f_in:
					reader = csv.reader(f_in, delimiter='\t')
					for row in reader:
						writer.writerow(row)

	logger.debug(f"Aggregation complete. Final output written to {final_output_path}")
	cleanup_temp_dir(temp_dir)

def cleanup_temp_dir(temp_dir):
	logger.debug(f"Cleaning up temporary files in {temp_dir}")	
	shutil.rmtree(temp_dir)

class BaseConfig:
	"""
	Base Configuration class, responsible for shared data, common properties,
	and the factory loading mechanism.
	"""
	raw_args: Any
	config_data: Dict[str, Any]
	mode: str
	contig: Optional[str]
	max_threads: int

	def __init__(self, raw_args: Any, config_data: Dict[str, Any], mode: str) -> None:
		"""Initializes common attributes using loaded config data."""
		self.raw_args = raw_args
		self.config_data = config_data
		self.mode = mode

		# Load universally common properties
		self.contig = config_data.get('contig')
		self.max_threads = config_data.get('max_threads', 1)

		# Post-initialization hook: Delegate control to the subclass for specific setup
		self._initialize_mode_specific_attributes()

	@classmethod
	def load_config(cls, args: Any, mode: str) -> Self:
		"""
		Generic method to load the appropriate configuration section from the file.
		This acts as the factory helper, performing file I/O before instantiation.
		"""
		config_path = getattr(args, 'config', None)
		if not config_path:
			logger.error("Configuration file path not provided.")
			sys.exit(1)

		try:
			with open(config_path, 'r') as f:
				full_config = yaml.safe_load(f)
			
			# Extract the section specific to the requested mode
			config_data = full_config.get(mode, {})
			if not config_data:
				logger.warning(f"Configuration section '{mode}' missing or empty in {config_path}.")

			# Instantiate the correct subclass (cls) using the loaded data
			return cls(args, config_data, mode)

		except Exception as e:
			logger.error(f"Error loading configuration file {config_path}: {e}")
			sys.exit(1)

	def _initialize_mode_specific_attributes(self):
		"""Placeholder for subclasses to define and initialize their unique fields."""
		# Base class implements nothing, relying on subclasses to override this
		pass


class ScannerConfig(BaseConfig):
	"""Configuration specific to the scanner mode."""
	bamfile:Path
	fastafile:Path
	min_length: int
	parallel_jobs: int
	output_dir: Path
	temp_dir: Path

	def _setup_output(self) -> Tuple[Path, Path]:
		"""Performs critical directory setup and cleanup for the scanner."""
		output_dir = Path(getattr(self.raw_args, 'output'))
		temp_dir = output_dir / "temp_indel_parts"
		
		try:
			output_dir.mkdir(parents=True, exist_ok=True)
			logger.debug(f"Output directory set up at {output_dir}")

			# CRITICAL SCANNER LOGIC: Clean up the temp directory before starting
			if temp_dir.exists():
				logger.info(f"Found existing temporary directory at {temp_dir}. Cleaning it before use.")
				shutil.rmtree(temp_dir)
			
			temp_dir.mkdir(parents=True, exist_ok=True)
			logger.debug(f"Temporary directory for partial results is ready at {temp_dir}")

		except Exception as e:
			logger.error(f"Error setting up output directories: {e}")
			sys.exit(1)
			
		return output_dir, temp_dir

	def _initialize_mode_specific_attributes(self):
		"""Sets scanner-specific fields and executes directory setup (cleanup is included)."""
		logger.info("Configuration Mode: Scanner")

		# Scanner-specific settings loaded from config_data
		self.min_length = self.config_data.get('min_length', 10)
		self.parallel_jobs = self.config_data.get('parallel_jobs', self.max_threads)
		self.bamfile = Path(self.raw_args.bam)
		self.fastafile = Path(self.raw_args.fasta)

		# Run Scanner-specific action: Directory setup and assign paths
		self.output_dir, self.temp_dir = self._setup_output()
		logger.info(f"Scanner output directories configured and cleaned.")



class ProcessorConfig(BaseConfig):
	"""Configuration specific to the processor mode."""
	output_format: str
	post_filter_threshold: float
	output_dir: Path
	temp_dir: Path

	def _initialize_mode_specific_attributes(self):
		"""Sets processor-specific fields and derives path information (no cleanup)."""
		logger.info("Configuration Mode: Processor")

		# Processor-specific settings loaded from config_data
		self.output_format = self.config_data.get('output_format', 'vcf')
		self.post_filter_threshold = self.config_data.get('post_filter_threshold', 0.95)
		
		self.bamfile = Path(self.raw_args.bam)
		self.fastafile = Path(self.raw_args.fasta)

		# Processor derives paths from raw arguments, assuming they exist
		self.output_dir = Path(getattr(self.raw_args, 'output'))
		self.temp_dir = self.output_dir / "temp_indel_parts"
		
		logger.info(f"Processor output configured (Format: {self.output_format}). No directory cleanup performed.")


# --- 4. Public Factory Interface ---

class Config:
	"""
	Public static factory interface to load specialized configuration classes.
	"""
	@staticmethod
	def ScannerConfig(args: Any) -> ScannerConfig:
		"""Factory method tailored for Scanner mode."""
		return ScannerConfig.load_config(args, 'scanner')

	@staticmethod
	def ProcessorConfig(args: Any) -> ProcessorConfig:
		"""Factory method tailored for Processor mode."""
		return ProcessorConfig.load_config(args, 'processor')

		