# indel_scanner/scanner_class.py
import csv
import os
from pathlib import Path
from typing import Generator, List
import pysam
import pyfastx
import time
import logging
from .indel import INDEL_TYPE, TSV_HEADERS, Insertion, Deletion
from .utils import Cigar, cleanup_temp_dir
from .configurator import ScannerConfig

logger = logging.getLogger(__name__)

class ContigScanner:
	"""
	A class to encapsulate the configuration and logic for scanning contigs.
	An instance of this class can be used as a worker in a multiprocessing pool.
	"""
	def __init__(self, config:ScannerConfig):
		"""
		Initializes the scanner with shared configuration that will be
		available to all calls within this worker process.
		"""
		self.config = config

	def scan_contig(self, contig_name:str) -> tuple[Path, str] | None:
	   
		try:
			with pysam.AlignmentFile(str(self.config.bamfile), "rb") as samfile:
				temp_output_path = self.config.temp_dir / f"{contig_name}.part.tsv"
				
				fasta = pyfastx.Fasta(str(self.config.fastafile)) # dosnt support context manager????
				logger.info(f"Scanning contig: {contig_name} @ {temp_output_path} \n with {(str(self.config.bamfile))} and {(str(self.config.fastafile))}")
				return self._process_reads(samfile, fasta, contig_name, temp_output_path)
			
		except Exception as e:
			logger.info(f"Error opening or processing files for contig {contig_name}: {e}")
			logger.error(f"Error opening or processing files for contig {contig_name}: {e}")
			return None

		
	def _parse_cigar(self, read: pysam.AlignedSegment, contig_seq: str) -> Generator[Insertion|Deletion]:
		"""
		Generator function to parse a CIGAR string and yield indel records.
		This isolates the core indel detection logic.
		"""
		ref_pos = read.reference_start
		query_pos = 0

		assert read is not None, logger.error("Read is None?")  # Read is None
		assert read.query_sequence is not None, logger.error("Read has no query sequence")  # Read has no query sequence
		assert read.cigartuples is not None,logger.error("No CIGAR information available")
		assert read.reference_name is not None,logger.error("No reference_name information available")
		assert read.query_name is not None,logger.error("No query_name information available")

		for op, length in read.cigartuples:

			logger.debug(f"Processing CIGAR operation {op} with length {length} at ref pos {ref_pos}, query pos {query_pos} in read {read.query_name}")
			if op == Cigar.OP_I:
				logger.debug(f"Found insertion of length {length} at ref pos {ref_pos}, query pos {query_pos} in read {read.query_name}")
				if length >= self.config.min_indel_size:
					logger.debug(f"insertion passed min size {self.config.min_indel_size}")
					prefix_context = read.query_sequence[max(0, query_pos - 5): query_pos]
					indel_seq = read.query_sequence[query_pos: query_pos + length]
					suffix_context = read.query_sequence[query_pos + length: query_pos + length + 5]
					

					yield Insertion(
									contig=read.reference_name,
									ref_position=ref_pos,
									type=INDEL_TYPE.INSERTION,
									length=length,
									prefix_context=prefix_context,
									suffix_context=suffix_context,
									read_name=read.query_name,
									inserted_seq=indel_seq
									)

				query_pos += length

			elif op == Cigar.OP_D:
				if length >= self.config.min_indel_size:

					# preloaded contig sequence
					ref_seq = contig_seq[ref_pos: ref_pos + length]
					prefix_context = contig_seq[max(0, ref_pos - 5): ref_pos]
					suffix_context = contig_seq[ref_pos+ length:ref_pos+length + 5]
					
					yield Deletion(
									contig=read.reference_name,
									ref_position=ref_pos,
									type=INDEL_TYPE.DELETION,
									length=length,
									prefix_context=prefix_context,
									reference_seq=ref_seq,
									suffix_context=suffix_context,
									read_name=read.query_name,
									)
				ref_pos += length
			elif op in [Cigar.OP_M, Cigar.OP_EQ, Cigar.OP_X]:
				ref_pos += length
				query_pos += length
			elif op == Cigar.OP_N:  # N is for introns, skips reference
				ref_pos += length
			elif op == Cigar.OP_S:  # S is soft-clip, consumes query
				query_pos += length	

	def _process_reads(self, samfile: pysam.AlignmentFile, fasta: pyfastx.Fasta, contig_name: str, temp_output_path: Path):
		"""
		Iterates through reads in a contig, uses a generator to find indels,
		and writes them to a file using a buffer.
		"""
		start_time = time.time()
		
		contig_seq = fasta[contig_name].seq
		
		logger.debug(f"Loaded contig sequence for {contig_name}, length: {len(contig_seq)}")

		BUFFER_SIZE = self.config.buffer_size
		results_buffer = []

		with open(temp_output_path, 'w') as f_out:
			writer = csv.writer(f_out, delimiter='\t')

			for read in samfile.fetch(contig=contig_name):
				if (read.is_unmapped or
					read.is_secondary or
					read.is_supplementary or
					read.query_sequence is None):
					#logger.info(f"Skipping read {read.query_name}: unmapped{read.is_unmapped}/secondary{read.is_secondary}/supplementary{read.is_supplementary}/no sequence{bool(read.query_sequence is None)}")
					continue

				for indel in self._parse_cigar(read, contig_seq):
					logger.debug(f"Detected indel: {indel}")
					results_buffer.append(indel.to_scanner_tsv_row())
				
				if len(results_buffer) >= BUFFER_SIZE:
					# Optimized writing without the csv module
					writer.writerows(results_buffer)
					results_buffer.clear()
			
			# Write any remaining results in the buffer
			if results_buffer:
				writer.writerows(results_buffer)

		end_time = time.time()
		return (temp_output_path, f"Finished {contig_name} in {end_time - start_time:.2f} s")


def run_scan(scanner:ContigScanner, contig_name):
	return scanner.scan_contig(contig_name)

def aggregate_partial_results(temp_dir, final_output_path):
	logger.info("\nAll contigs processed. Merging results...")
	logger.debug(f"Aggregating partial results from {temp_dir} into {final_output_path}")
	
	with open(final_output_path, 'w', newline='') as f_out:
		writer = csv.writer(f_out, delimiter='\t')
		
		# header
		writer.writerow(TSV_HEADERS.SCANNER.value)
		
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