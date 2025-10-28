from functools import partial
import os
import pysam
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn, TimeElapsedColumn
from multiprocessing import Pool, cpu_count
import logging

from .utils import Config, ScannerConfig, cleanup_temp_dir, aggregate_partial_results
from indel_scanner.scanner import run_scan
from .scanner import ContigScanner

logger = logging.getLogger(__name__)

def parallel_scan(args):

	scannerconfig: ScannerConfig = Config.ScannerConfig(args)
	scanner = ContigScanner(scannerconfig)

	logger.info(f"Using {scannerconfig.max_threads or cpu_count()} parallel processes.")
	
	logger.info("Starting indel scanning process...")

	with pysam.AlignmentFile(str(scannerconfig.bamfile), "rb") as samfile:
		contigs = samfile.header.references
		logger.info(f"found {len(contigs)} contigs.") 

	

	# Define the columns for the progress bar for a nice look
	progress_columns = [
		TextColumn("[progress.description]{task.description}"),
		BarColumn(),
		TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
		TextColumn("•"),
		TimeElapsedColumn(),
	]

	partial_files = []
	partialfunc = partial(run_scan, scanner)

	# Use the Progress context manager
	with Progress(*progress_columns, transient=False) as progress:
		# Add a task to the progress display
		task_id = progress.add_task("[green]Scanning contigs...", total=len(contigs))
		
		with Pool(processes=scannerconfig.max_threads or cpu_count()) as pool:
			# Use imap_unordered to get results as soon as they are ready
			results_iterator = pool.imap_unordered(partialfunc, contigs)
			
			# Iterate over the results from the workers
			for result_file, status_message in results_iterator: # pyright: ignore[reportGeneralTypeIssues]
				# If the worker sent back a status message, print it
				if status_message:
					# Use progress.print() instead of the standard print()
					progress.print(status_message)

				if result_file:
					partial_files.append(result_file)

				# Advance the progress bar by one step for each completed contig
				progress.update(task_id, advance=1)

	aggregate_partial_results(scanner.temp_dir, args.output)