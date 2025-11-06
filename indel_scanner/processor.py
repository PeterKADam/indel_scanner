
import sys
import logging
import csv
import pyfastx
import pysam

from typing import List, Optional, Tuple
from .indel import INDEL_TYPE

from indel_scanner.configurator import ProcessorConfig


logger = logging.getLogger(__name__)

class Processor:
	def __init__(self, config:ProcessorConfig):
		self.config = config
		self.reference = pyfastx.Fasta(self.config.fastafile)
	

		 # Open BAM file (must be indexed)
		try:
			self.bam_handle = pysam.AlignmentFile(str(config.bamfile), "rb")
			logger.info(f"Successfully opened BAM file: {config.bamfile}")
		except Exception as e:
			logger.error(f"Error opening BAM file {config.bamfile}: {e}")
			self.bam_handle = None

		self.processed_records: List[List[str]] = []
		self.insertion_output: List[List[str]] = []
		self.deletion_output: List[List[str]] = []

	def __del__(self):
		"""Ensure the BAM file handle is closed when the object is destroyed."""
		if hasattr(self, 'bam_handle') and self.bam_handle:
			self.bam_handle.close()
			logger.debug(f"Closed BAM file handle for {self.config.bamfile}")
	
			
	def _parse_sequence_context(self, context_string: str) -> Tuple[str, str, str]:
		"""
		Parses the unified sequence context string: prefix[indel_sequence]suffix
		

		Args:
			context_string: prefix[indel_sequence]suffix.
			
		Returns:
			Tuple[prefix, indel_seq, suffix]

		"""
		try:
			# 1. Find the opening bracket
			start_bracket = context_string.index('[')
			# 2. Find the closing bracket
			end_bracket = context_string.index(']')
			
			# The structure must be correct, or it's an error
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
		"""
		Checks if an indel is part of a homopolymer run of 3 or more identical bases.

		This function is designed to filter sequencing artifacts where an indel of any
		length (e.g., 'A', 'CC') occurs within or extends a homopolymer sequence.

		Args:
			prefix: The sequence immediately before the indel.
			indel_seq: The sequence of the indel itself (e.g., 'A', 'GG').
			suffix: The sequence immediately after the indel.

		Returns:
			True if the indel is part of a homopolymer run of 3+, False otherwise.
		"""
		# 1. The indel sequence must exist and must itself be a homopolymer
		#    (e.g., 'A', 'GG', 'TTT'). It cannot be a mixed sequence like 'AG'.
		if not indel_seq or len(set(indel_seq)) != 1:
			return False

		# 2. Get the single base that makes up the homopolymer indel.
		indel_base = indel_seq[0]

		# 3. Calculate the total length of the continuous homopolymer run
		#    at the breakpoint by combining the indel and its flanking bases.

		# Start with the length of the indel itself.
		run_length = len(indel_seq)

		# Count matching bases by looking backwards from the end of the prefix.
		for char in reversed(prefix):
			if char == indel_base:
				run_length += 1
			else:
				# Stop counting as soon as the run is broken.
				break

		# Count matching bases by looking forwards from the start of the suffix.
		for char in suffix:
			if char == indel_base:
				run_length += 1
			else:
				# Stop counting as soon as the run is broken.
				break

		# 4. Final decision: Is the total contiguous run length 3 or more?
		if run_length >= 3:
			logger.debug(
				f"Filtering homopolymer: {prefix}[{indel_seq}]{suffix} "
				f"(run of '{indel_base}' with total length {run_length})"
			)
			return True

		return False

	
	def _is_adjacent_to_homopolymer(self, prefix: str, indel_seq: str, suffix: str, min_len: int = 3) -> bool:
		"""
		Checks if an indel is immediately adjacent to a homopolymer run of a minimum length.

		This filter is broader than a simple stutter filter. It removes any indel,
		regardless of its own sequence, if it touches a homopolymer run at the breakpoint.

		Args:
			prefix: The sequence immediately before the indel.
			indel_seq: The sequence of the indel itself.
			suffix: The sequence immediately after the indel.
			min_len: The minimum length to define a homopolymer run (e.g., 3 for 'AAA').

		Returns:
			True if the indel is adjacent to a homopolymer run of min_len or more.
		"""
		# Check the prefix: Does it end in a homopolymer run of at least min_len?
		if len(prefix) >= min_len:
			# Check if the last `min_len` characters of the prefix are all the same
			end_of_prefix = prefix[-min_len:]
			if len(set(end_of_prefix)) == 1:
				logger.debug(
					f"Filtering: Indel {prefix}[{indel_seq}]{suffix} is adjacent "
					f"to a homopolymer run in the prefix: '{end_of_prefix}'"
				)
				return True

		# Check the suffix: Does it start with a homopolymer run of at least min_len?
		if len(suffix) >= min_len:
			# Check if the first `min_len` characters of the suffix are all the same
			start_of_suffix = suffix[:min_len]
			if len(set(start_of_suffix)) == 1:
				logger.debug(
					f"Filtering: Indel {prefix}[{indel_seq}]{suffix} is adjacent "
					f"to a homopolymer run in the suffix: '{start_of_suffix}'"
				)
				return True

		return False

		
	def process(self):
		
		"""
		Main method to read the input file, filter, and separate outputs.
		"""

		if not self.bam_handle:
			logger.error("Cannot run processing: BAM file not accessible.")
			sys.exit(1)

		logger.info(f"Processing indel file: {self.config.input_file}")

		# Define expected output header (for clarity, though we write data only)
		OUTPUT_HEADER = ["contig", "position", "indel_type", "length", "prefix", "indel_seq", "suffix", "read_name"]

		with open(self.config.input_file, 'r', newline='') as infile:
			reader = csv.reader(infile, delimiter='\t')
			
			for i, row in enumerate(reader):
				if len(row) < 6:
					logger.warning(f"Skipping malformed row #{i+1}: {row}")
					continue
				
				contig, pos_str, indel_type, length_str, context_string, read_name = row
				
				# 1. Parse the sequence context
				prefix, indel_seq, suffix = self._parse_sequence_context(context_string)
				
				# 2. Filter homopolymers
				if self._is_homopolymer(prefix, indel_seq, suffix) or self._is_adjacent_to_homopolymer(prefix, indel_seq, suffix):
					continue # Skip this record if it is a homopolymer

				# 3. Store and separate outputs
				# Note: 'indel_seq' is the inserted sequence for INS, and the DELETED sequence for DEL.
				processed_row = [contig, pos_str, indel_type, length_str, prefix, indel_seq, suffix, read_name]
				self.processed_records.append(processed_row)

				if indel_type == INDEL_TYPE.INSERTION:
					self.insertion_output.append(processed_row)
				elif indel_type == INDEL_TYPE.DELETION:
					self.deletion_output.append(processed_row)

		logger.info(f"Total records processed: {len(self.processed_records)}")
		logger.info(f"Insertions found (post-filter): {len(self.insertion_output)}")
		logger.info(f"Deletions found (post-filter): {len(self.deletion_output)}")

		return self.insertion_output, self.deletion_output
	
	def query_insertion_quality(self, contig: str, position: int, read_name: str, inserted_sequence: str) -> Optional[List[int]]:
		"""
		Queries the original BAM file for the base quality scores of the inserted bases.
		
		Args:
			contig: Chromosome name.
			position: Genomic position of the insertion start (0-based).
			read_name: The query name of the read containing the insertion.
			inserted_sequence: The sequence of the insertion (e.g., 'C' or 'GAT').
			
		Returns:
			List of Phred quality scores corresponding to the inserted_sequence bases, 
			or None if the read/insertion cannot be located.
		"""
		if not self.bam_handle:
			logger.error("BAM handle is closed or unavailable for quality query.")
			return None
		
		try:
		
			# This is generally faster than iterating all reads in the file.
			# Position here is the 0-based coordinate immediately *before* the insertion site.
			for read in self.bam_handle.fetch(contig, position, position + 1):
				if read.query_name == read_name:
					# Found the target read. Now locate the insertion in the CIGAR.
					
					query_pos = 0  # Position along the read sequence
					ref_pos = read.reference_start  # Position along the reference
					
					target_qualities = []
					assert read is not None, logger.error("Read is None?") 
					assert read.query_sequence is not None, logger.error("Read has no query sequence")  
					assert read.cigartuples is not None,logger.error("No CIGAR information available")
					assert read.query_qualities is not None, logger.error("No quality scores available for read") 
					
					for operation, length in read.cigartuples:
						if operation == pysam.CINS: # Insertion relative to reference
							if ref_pos == position:
								# found the insertion at the correct reference coordinate.
								
								# Check if the bases in the read match the reported insertion
								read_insertion = read.query_sequence[query_pos:query_pos + length]
								if read_insertion == inserted_sequence:
									
									target_qualities = read.query_qualities[query_pos:query_pos + length]
									return target_qualities # type: ignore
								else:
									# Insertion sequence mismatch (shouldn't happen if TSV is accurate)
									logger.warning(f"CIGAR insertion mismatch for read {read_name} at {contig}:{position}.")
									return None
							
							# Inserted bases move the query position but not the reference position
							query_pos += length
							
						elif operation == pysam.CMATCH or operation == pysam.CDIFF: # M/=/X
							# Matched/Mismatched bases move both reference and query positions
							query_pos += length
							ref_pos += length
							
						elif operation == pysam.CDEL or operation == pysam.CREF_SKIP: # D/N
							# Deletions/Ref skips move reference position but not query position
							ref_pos += length
							
						elif operation == pysam.CSOFT_CLIP: # S
							# Soft clipping moves query position but not reference position in the aligned region
							query_pos += length

						

			logger.warning(f"Target read '{read_name}' not found near {contig}:{position} or insertion CIGAR incorrect.")
			return None

		except ValueError as e:
			logger.error(f"Pysam error during query: {e}. Check if BAM is indexed.")
			return None
		except Exception as e:
			logger.error(f"An unexpected error occurred during BAM query: {e}")
			return None


	def query_deletion_flanking_quality(self, contig: str, position: int, length: int = 5, deleted_length: int = 1) -> str:
		
		if not self.reference:
			return "Reference_Unavailable"
		
		try:
			position = int(position)
			
			# 1. Prefix (Upstream) region: [position - length, position)
			prefix_start = max(0, position - length)
			prefix_end = position
			
			# 2. Deleted region: [position, position + deleted_length)
			deleted_start = position
			deleted_end = position + deleted_length
			
			# 3. Suffix (Downstream) region: [deleted_end, deleted_end + length)
			suffix_start = deleted_end
			suffix_end = deleted_end + length
			
			contig_seq_obj = self.reference[contig]
			
			prefix_seq = contig_seq_obj.seq(prefix_start, prefix_end)
			deleted_seq = contig_seq_obj.seq(deleted_start, deleted_end)
			suffix_seq = contig_seq_obj.seq(suffix_start, suffix_end)
			
			full_context = f"{prefix_seq}[{deleted_seq}]{suffix_seq}"
			
			if 'N' in full_context.upper():
				logger.warning(f"Low quality reference bases ('N') found in context: {full_context}")
				
			return full_context.upper()
			
		except Exception as e:
			logger.error(f"Error querying FASTA for deletion context at {contig}:{position}: {e}")
			return "FASTA_Query_Failed"

	def write_output(self, base_output_path: str):
		
		
		OUTPUT_HEADER = ["contig", "position", "indel_type", "length", "prefix", "indel_seq", "suffix", "read_name"] #wrong

		ins_path = base_output_path.replace(".tsv", "_insertions.tsv")
		del_path = base_output_path.replace(".tsv", "_deletions.tsv")

		with open(ins_path, 'w', newline='') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(OUTPUT_HEADER)
			writer.writerows(self.insertion_output)
			logger.info(f"Wrote {len(self.insertion_output)} insertions to {ins_path}")

		with open(del_path, 'w', newline='') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(OUTPUT_HEADER)
			writer.writerows(self.deletion_output)
			logger.info(f"Wrote {len(self.deletion_output)} deletions to {del_path}")


# --- Helper for integration into main.py 
def run_processor(config):
	
	logger.info("Starting indel processing task...")

	# Initialize the processor
	processor = Processor(
		config=config
	)

	processor.process()

	processor.write_output(config.output_path)

	logger.info("Indel processing finished successfully.")