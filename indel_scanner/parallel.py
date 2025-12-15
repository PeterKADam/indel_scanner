from functools import partial
import pysam
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from multiprocessing import Pool, cpu_count
import logging

from .scanner import ContigScanner, run_scan, aggregate_partial_results


logger = logging.getLogger(__name__)


def parallel_scan(scannerconfig):
    scanner = ContigScanner(scannerconfig)

    parallel_n = scannerconfig.num_processes or cpu_count()
    logger.info(f"Using {parallel_n} parallel processes.")

    logger.info("Starting indel scanning process...")

    with pysam.AlignmentFile(str(scannerconfig.bamfile), "rb") as samfile:
        contigs = samfile.header.references
        contigs = contigs[:10]  # --- IGNORE ---
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

        with Pool(processes=parallel_n) as pool:
            # Use imap_unordered to get results as soon as they are ready
            results_iterator = pool.imap_unordered(partialfunc, contigs)

            for result_file, status_message in results_iterator:  # pyright: ignore[reportGeneralTypeIssues]
                if status_message:
                    progress.print(status_message)

                if result_file:
                    partial_files.append(result_file)

                # Advance the progress bar by one step for each completed contig
                progress.update(task_id, advance=1)

    aggregate_partial_results(scannerconfig.temp_dir, scannerconfig.output_path)
