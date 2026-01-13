
from functools import partial
from pathlib import Path
import pysam
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from multiprocessing import Pool, cpu_count
import logging
from .configurator import ScannerConfig
from .scanner import ContigScanner, run_scan, aggregate_partial_results

logger = logging.getLogger(__name__)


def parallel_scan(scannerconfig: ScannerConfig) -> tuple[Path, int]:
    scanner = ContigScanner(scannerconfig)
    parallel_n = scannerconfig.num_processes or cpu_count()
    logger.info(f"Using {parallel_n} parallel processes.")
    logger.info("Starting indel scanning process...")

    with pysam.AlignmentFile(str(scannerconfig.bamfile), "rb") as samfile:
        contigs = list(samfile.header.references)
    logger.info(f"Found {len(contigs)} contigs.")

    progress_columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("•"),
        TimeElapsedColumn(),
    ]

    partial_files = []
    total_interrogated_bases = 0
    partialfunc = partial(run_scan, scanner)

    with Progress(*progress_columns, transient=False) as progress:
        task_id = progress.add_task("[green]Scanning contigs...", total=len(contigs))
        with Pool(processes=parallel_n) as pool:
            results_iterator = pool.imap_unordered(partialfunc, contigs)
            for result in results_iterator:
                if result:
                    result_file, status_message, contig_base_count = result
                    if status_message:
                        progress.print(status_message)
                    if result_file:
                        partial_files.append(result_file)
                    total_interrogated_bases += contig_base_count
                progress.update(task_id, advance=1)

    final_scanner_output = scannerconfig.output_path / "indel_scanner_results.tsv"
    aggregate_partial_results(scannerconfig.temp_dir, scannerconfig.output_path)

    return final_scanner_output, total_interrogated_bases
