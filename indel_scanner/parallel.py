
from pathlib import Path
import pysam
import pyfastx
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from multiprocessing import Pool, cpu_count, current_process, Manager
import logging
from queue import Empty
from .configurator import PipelineConfig
from .scanner import ContigScanner
from .IO import aggregate_tsv_parts, cleanup_temp_dir
from .indel import TSV_HEADERS
from .utils import Timer

logger = logging.getLogger(__name__)

_WORKER_BAM: pysam.AlignmentFile | None = None
_WORKER_FASTA: pyfastx.Fasta | None = None
_WORKER_CONFIG: PipelineConfig | None = None
_WORKER_STATUS_QUEUE = None


def _init_worker(
    bam_path: str,
    fasta_path: str,
    config: PipelineConfig,
    sampling_contig: str,
    status_queue,
):
    global _WORKER_BAM, _WORKER_FASTA, _WORKER_CONFIG, _WORKER_STATUS_QUEUE
    _WORKER_BAM = pysam.AlignmentFile(bam_path, "rb")
    _WORKER_FASTA = pyfastx.Fasta(fasta_path)
    _WORKER_CONFIG = config
    _WORKER_CONFIG.sampling_contig = sampling_contig
    _WORKER_STATUS_QUEUE = status_queue


def _process_contig_streaming(contig_name: str):
    if _WORKER_CONFIG is None or _WORKER_BAM is None or _WORKER_FASTA is None:
        raise RuntimeError("Worker not initialized with BAM/FASTA/config.")
    if _WORKER_STATUS_QUEUE is not None:
        _WORKER_STATUS_QUEUE.put(
            ("start", current_process().name, contig_name))
    scanner = ContigScanner(_WORKER_CONFIG)
    result = scanner.scan_contig_streaming(
        contig_name, _WORKER_BAM, _WORKER_FASTA, _WORKER_CONFIG.in_memory
    )
    if _WORKER_STATUS_QUEUE is not None:
        _WORKER_STATUS_QUEUE.put(("done", current_process().name, contig_name))
    return (contig_name, *result)


def parallel_pipeline(scannerconfig: PipelineConfig) -> tuple[Path, int, dict]:
    parallel_n = scannerconfig.num_processes or cpu_count()
    logger.info(f"Using {parallel_n} parallel processes.")
    logger.info("Starting streaming scan+process pipeline...")

    with pysam.AlignmentFile(str(scannerconfig.bamfile), "rb") as samfile:
        contigs = list(samfile.header.references)
        contig_lengths = list(samfile.header.lengths)
    contig_filter = getattr(scannerconfig, "contigs", None)
    if isinstance(contig_filter, (list, tuple, set)) and contig_filter:
        contig_lookup = dict(zip(contigs, contig_lengths))
        contigs = [c for c in contigs if c in contig_filter]
        contig_lengths = [contig_lookup[c] for c in contigs]
    logger.info(f"Found {len(contigs)} contigs.")
    if scannerconfig.sampling_strategy == "largest_contig":
        largest_idx = max(range(len(contigs)), key=lambda i: contig_lengths[i])
        sampling_contig = contigs[largest_idx]
    else:
        sampling_contig = contigs[0]

    progress_columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("•"),
        TimeElapsedColumn(),
        TextColumn("{task.fields[worker_status]}"),
    ]

    total_interrogated_bases = 0
    total_stats = {
        "reads": 0,
        "candidates": 0,
        "passed": 0,
        "total_aligned_bases": 0,
        "sampling_bases_total": 0,
        "sampling_passable_by_type": {},
        "type_counts": {},
    }
    in_memory_records = []

    with Manager() as manager:
        status_queue = manager.Queue()
        active_workers: dict[str, str] = {}

        def format_worker_status() -> str:
            if not active_workers:
                return "Workers: idle"
            entries = ", ".join(f"{name}={contig}" for name,
                                contig in sorted(active_workers.items()))
            return f"Workers: {entries}"

        def drain_status_queue():
            while True:
                try:
                    event, worker_name, contig_name = status_queue.get_nowait()
                except Empty:
                    break
                if event == "start":
                    active_workers[worker_name] = contig_name
                elif event == "done":
                    active_workers.pop(worker_name, None)

        with Progress(*progress_columns, transient=False) as progress:
            task_id = progress.add_task(
                "[green]Scanning contigs...",
                total=len(contigs),
                worker_status="Workers: idle",
            )
            with Pool(
                processes=parallel_n,
                initializer=_init_worker,
                initargs=(
                    str(scannerconfig.bamfile),
                    str(scannerconfig.fastafile),
                    scannerconfig,
                    sampling_contig,
                    status_queue,
                ),
            ) as pool:
                results_iterator = pool.imap_unordered(
                    _process_contig_streaming, contigs)
                with Timer("Parallel pipeline"):
                    for result in results_iterator:
                        drain_status_queue()
                        progress.update(
                            task_id, worker_status=format_worker_status())
                        if result:
                            contig_name, contig_records, contig_base_count, stats = result
                            progress.print(
                                f"[cyan]Completed contig[/cyan] {contig_name}")
                            total_interrogated_bases += contig_base_count
                            total_stats["reads"] += stats["reads_processed"]
                            total_stats["candidates"] += stats["candidates"]
                            total_stats["passed"] += stats["passed"]
                            total_stats["total_aligned_bases"] += stats["total_aligned_bases"]
                            if stats.get("sampled_contig"):
                                total_stats["sampling_bases_total"] += stats["sampling_bases_total"]
                                for k, v in stats["sampling_passable_by_type"].items():
                                    total_stats["sampling_passable_by_type"][k] = (
                                        total_stats["sampling_passable_by_type"].get(
                                            k, 0) + v
                                    )
                            for k, v in stats["type_counts"].items():
                                total_stats["type_counts"][k] = (
                                    total_stats["type_counts"].get(k, 0) + v
                                )

                            if scannerconfig.in_memory:
                                in_memory_records.extend(contig_records or [])
                                if (
                                    scannerconfig.max_in_memory_records
                                    and len(in_memory_records) > scannerconfig.max_in_memory_records
                                ):
                                    logger.warning(
                                        "In-memory limit exceeded; consider disabling in-memory mode."
                                    )
                            progress.update(task_id, advance=1)
                    drain_status_queue()
                    progress.update(
                        task_id, worker_status=format_worker_status())

    if scannerconfig.in_memory:
        from .IO import write_passed_indels

        write_passed_indels(
            scannerconfig.passed_indels_path,
            in_memory_records,
            scannerconfig.write_buffer_size,
        )
    else:
        aggregate_tsv_parts(
            scannerconfig.passed_parts_dir,
            scannerconfig.passed_indels_path,
            TSV_HEADERS.PROCESSOR.value,
        )
        cleanup_temp_dir(scannerconfig.passed_parts_dir)

    return scannerconfig.passed_indels_path, total_interrogated_bases, total_stats
