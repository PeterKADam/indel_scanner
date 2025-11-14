from pathlib import Path
import sys
import logging
import csv
from collections import defaultdict

import pyfastx
import pysam
from pysam import AlignedSegment
from typing import List, Optional, Tuple, Union, Dict

from .indel import INDEL_TYPE, TSV_HEADERS
from .indel import Insertion, Deletion
from .utils import Cigar
from indel_scanner.configurator import ProcessorConfig

logger = logging.getLogger(__name__)

class Processor:
	def __init__(self, config: ProcessorConfig):
		self.config = config
		self.reference = None  # Will be initialized in __enter__
		self.bam_handle = None  # Will be initialized in __enter__
		
		# Structure: {contig: {read_name: [list_of_indel_objects]}}
		self.indels_by_contig: Dict[str, Dict[str, List[Union[Insertion, Deletion]]]] = defaultdict(lambda: defaultdict(list))

		
		self.insertion_output: List[Insertion] = []
		self.deletion_output: List[Deletion] = []
		self.total_records_processed: int = 0

	def __enter__(self):
		logger.debug("Entering Processor context manager...")
		try:
			logger.debug("Opening BAM file for processing...")
			self.bam_handle = pysam.AlignmentFile(str(self.config.bamfile), "rb")
			logger.info(f"Successfully opened BAM file: {self.config.bamfile}")

			logger.debug("Opening reference FASTA for processing...")
			self.reference = pyfastx.Fasta(str(self.config.fastafile))
			logger.info(f"Successfully opened reference FASTA: {self.config.fastafile}")
			return self
		except Exception as e:
			logger.error(f"Error opening files: {e}")
			raise IOError(f"Failed to open input files.") from e

	def __exit__(self, exc_type, exc_value, traceback):
		if self.bam_handle:
			self.bam_handle.close()
			logger.debug(f"Closed BAM file handle for {self.config.bamfile}")

	def _load_and_filter_indels(self):
  
		logger.info(f"Loading and filtering indels from {self.config.input_file}...")
		
		try:
			with open(self.config.input_file, 'r', newline='') as infile:
				reader = csv.DictReader(infile, delimiter='\t')
				for i, row in enumerate(reader):
					
					try:
						contig = row['Contig']
						position = int(row['Position'])
						length = int(row['Length'])
						indel_type = INDEL_TYPE(row['Type'])
						context_string = row['Sequence']
						read_name = row['Read_Name']
						
					except (ValueError, KeyError, TypeError) as e:
						logger.warning(f"Skipping row #{i+2} due to malformed data: {row} | Error: {e}")
						continue

					prefix, indel_seq, suffix = self._parse_sequence_context(context_string)
					
					if self._is_homopolymer(prefix, indel_seq, suffix) or \
					   self._is_adjacent_to_homopolymer(prefix, indel_seq, suffix):
						continue  # Skip this record if it is a homopolymer

					indel_obj = None
					if indel_type == INDEL_TYPE.INSERTION:
						indel_obj = Insertion(
							contig=contig, ref_position=position, length=length,
							prefix_context=prefix, suffix_context=suffix, read_name=read_name,
							type=indel_type, inserted_seq=indel_seq
						)
					elif indel_type == INDEL_TYPE.DELETION:
						ref_seq = self.reference.fetch(contig, (position, position + length - 1)) # type: ignore
						indel_obj = Deletion(
							contig=contig, ref_position=position, length=length,
							prefix_context=prefix, suffix_context=suffix, read_name=read_name,
							type=indel_type, reference_seq=ref_seq
						)

					if indel_obj:
						self.indels_by_contig[contig][read_name].append(indel_obj)
						logger.info(f"added indel {contig}:{position}:{indel_type}")

		except FileNotFoundError:
			logger.error(f"Input file not found: {self.config.input_file}")
			raise
		
		num_indels = sum(len(indels) for reads in self.indels_by_contig.values() for indels in reads.values())
		logger.info(f"Loaded {num_indels} indels across {len(self.indels_by_contig)} contigs.")


	def _process_bam_file(self):
		if not self.bam_handle:
			raise IOError("BAM file handle is not open. Cannot process.")

		logger.info("Streaming BAM file and processing reads...")
		
		for contig, indels_on_this_contig in self.indels_by_contig.items():
			logger.info(f"Processing contig: {contig}...")
			
			for bam_read in self.bam_handle.fetch(contig):
				if bam_read.is_secondary or bam_read.is_supplementary:
					continue
				
				if bam_read.query_name in indels_on_this_contig:
					
					target_indels = indels_on_this_contig[bam_read.query_name]
					for indel_obj in target_indels:
						if indel_obj.type == INDEL_TYPE.INSERTION:
							self._populate_insertion_quality(bam_read, indel_obj) # type: ignore
							self.insertion_output.append(indel_obj) # type: ignore
						elif indel_obj.type == INDEL_TYPE.DELETION:
							self._populate_deletion_flanking_quality(bam_read, indel_obj) # type: ignore
							self.deletion_output.append(indel_obj) # type: ignore
						
						self.total_records_processed += 1
		
		logger.info(f"Processed {self.total_records_processed} total records.")


	def process(self):
		if not self.bam_handle or not self.reference:
			logger.error("Cannot run processing: BAM or FASTA file not accessible.")
			sys.exit(1)
		
		
		self._load_and_filter_indels()
		
		
		self._process_bam_file()
		
		logger.info(f"Total records processed and passed filters: {self.total_records_processed}")
		return self.insertion_output, self.deletion_output

	def _populate_insertion_quality(self, target_read: AlignedSegment, insertion: Insertion, flank_len: int = 5):
		""" Populates quality fields of an Insertion object using a pre-fetched read. """
		qualities = target_read.query_qualities
		if qualities is None:
			logger.debug(f"Read '{target_read.query_name}' has no quality scores ('*').")
			return

		ref_pos_tracker = target_read.reference_start
		read_pos_tracker = 0
		found_insertion = False
		
		REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
		READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}

		for op, length in target_read.cigartuples: # type: ignore
			if (op == Cigar.OP_I and ref_pos_tracker == insertion.ref_position and length == insertion.length):
				prefix_start = max(0, read_pos_tracker - flank_len)
				insertion_end = read_pos_tracker + length
				suffix_end = min(len(qualities), insertion_end + flank_len)

				insertion.prefix_quality = list(qualities[prefix_start:read_pos_tracker])
				insertion.insertion_quality = list(qualities[read_pos_tracker:insertion_end])
				insertion.suffix_quality = list(qualities[insertion_end:suffix_end])
				found_insertion = True
				break

			if op in REF_CONSUMING_OPS:
				ref_pos_tracker += length
			if op in READ_CONSUMING_OPS:
				read_pos_tracker += length
		
		if not found_insertion:
			logger.debug(f"CIGAR parsing failed to find a matching insertion for {insertion.read_name}")


	def _populate_deletion_flanking_quality(self, target_read: AlignedSegment, deletion: Deletion, flank_len: int = 5):
		""" Populates quality fields of a Deletion object using a pre-fetched read. """
		qualities = target_read.query_qualities
		if qualities is None:
			logger.debug(f"Read '{target_read.query_name}' has no quality scores ('*').")
			return

		ref_pos_tracker = target_read.reference_start
		read_pos_tracker = 0
		found_deletion = False

		REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
		READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}

		for op, length in target_read.cigartuples: # type: ignore
			
			if (op == Cigar.OP_D and ref_pos_tracker == deletion.ref_position and length == deletion.length):
			  
				prefix_start = max(0, read_pos_tracker - flank_len)
				suffix_end = min(len(qualities), read_pos_tracker + flank_len)
				
				deletion.prefix_quality = list(qualities[prefix_start:read_pos_tracker])
				deletion.suffix_quality = list(qualities[read_pos_tracker:suffix_end])
				found_deletion = True
				break

			if op in REF_CONSUMING_OPS:
				ref_pos_tracker += length
			if op in READ_CONSUMING_OPS:
				read_pos_tracker += length

		if not found_deletion:
			logger.debug(f"CIGAR parsing failed to find a matching deletion for {deletion.read_name}")


	
	def _parse_sequence_context(self, context_string: str) -> Tuple[str, str, str]:
	   
		try:
			start_bracket = context_string.index('[')
			end_bracket = context_string.index(']')
			if start_bracket < 0 or end_bracket < 0 or start_bracket >= end_bracket:
				raise ValueError("Invalid sequence context format.")
			prefix = context_string[:start_bracket]
			indel_seq = context_string[start_bracket + 1 : end_bracket]
			suffix = context_string[end_bracket + 1 :]
			return prefix, indel_seq, suffix
		except ValueError as e:
			logger.error(f"Failed to parse sequence context '{context_string}': {e}")
			return "", "", ""

	def _is_homopolymer(self, prefix: str, indel_seq: str, suffix: str) -> bool:
		
		if not indel_seq or len(set(indel_seq)) != 1:
			return False
		indel_base = indel_seq[0]
		run_length = len(indel_seq)
		for char in reversed(prefix):
			if char == indel_base: run_length += 1
			else: break
		for char in suffix:
			if char == indel_base: run_length += 1
			else: break
		if run_length >= 3:
			logger.debug(f"Filtering homopolymer: {prefix}[{indel_seq}]{suffix}")
			return True
		return False

	def _is_adjacent_to_homopolymer(self, prefix: str, indel_seq: str, suffix: str, min_len: int = 3) -> bool:
		
		if len(prefix) >= min_len:
			end_of_prefix = prefix[-min_len:]
			if len(set(end_of_prefix)) == 1:
				logger.debug(f"Filtering adjacent homopolymer (prefix): {end_of_prefix}")
				return True
		if len(suffix) >= min_len:
			start_of_suffix = suffix[:min_len]
			if len(set(start_of_suffix)) == 1:
				logger.debug(f"Filtering adjacent homopolymer (suffix): {start_of_suffix}")
				return True
		return False

	def _write_indels_to_file(self, output_path: Path, records: List[Union[Insertion, Deletion]]):
		if not records:
			logger.info(f"No records to write for {output_path.name}, skipping.")
			return
		logger.info(f"Writing {len(records)} enriched records to {output_path}...")
		try:
			sorted_records = sorted(records, key=lambda r: (r.contig, r.ref_position))
			with open(output_path, 'w', newline='') as outfile:
				writer = csv.writer(outfile, delimiter='\t')
				writer.writerow(TSV_HEADERS.PROCESSOR.value)
				for record in sorted_records:
					writer.writerow(record.to_processor_tsv_row())
			logger.info(f"Successfully wrote to {output_path}.")
		except IOError as e:
			logger.error(f"Failed to write output file at {output_path}: {e}")

	def write_output(self, insertions_path: Optional[Path] = None, deletions_path: Optional[Path] = None):
		if not any([insertions_path, deletions_path]):
			logger.warning("write_output called, but no output paths were provided.")
			return
		if insertions_path:
			self._write_indels_to_file(insertions_path, self.insertion_output) # type: ignore
		if deletions_path:
			self._write_indels_to_file(deletions_path, self.deletion_output) # type: ignore


def run_processor(config):
	logger.info("Starting indel processing task...")
	try:
		logger.info("Initializing Processor...")
		with Processor(config) as processor:
			logger.info("Processor initialized successfully.")
			processor.process()
			processor.write_output(config.insertions_path, config.deletions_path)
		logger.info("Indel processing finished successfully.")
	except IOError as e:
		logger.error(f"Indel processing failed: {e}")
		sys.exit(1)
	except Exception as e:
		logger.error(f"An unexpected error occurred during indel processing: {e}", exc_info=True)
		sys.exit(1)
