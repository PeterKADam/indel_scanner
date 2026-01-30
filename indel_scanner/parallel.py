
from pathlib import Path
import csv
import pysam
import pyfastx
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.live import Live
from rich.panel import Panel
from rich.table import Table
from rich.console import Console, Group
from multiprocessing import Pool, cpu_count, current_process, Queue
import logging
import time
from threading import Event, Lock, Thread
from queue import Empty
import sys
import select
import termios
import tty
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


def _process_contig_callable(contig_name: str):
    if _WORKER_CONFIG is None or _WORKER_BAM is None or _WORKER_FASTA is None:
        raise RuntimeError("Worker not initialized with BAM/FASTA/config.")
    scanner = ContigScanner(_WORKER_CONFIG)
    stats = scanner.scan_contig_callable_only(
        contig_name, _WORKER_BAM, _WORKER_FASTA
    )
    return contig_name, stats


def _collect_indel_lengths(passed_indels_path: Path) -> list[int]:
    if not passed_indels_path.exists():
        return []
    lengths = set()
    try:
        with open(passed_indels_path, "r") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                return []
            try:
                length_idx = header.index("length")
            except ValueError:
                return []
            for row in reader:
                if len(row) <= length_idx:
                    continue
                try:
                    lengths.add(int(row[length_idx]))
                except ValueError:
                    continue
    except Exception:
        return []
    return sorted(lengths)


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
    if contigs:
        contig_pairs = sorted(
            zip(contigs, contig_lengths), key=lambda x: x[1], reverse=True
        )
        contigs, contig_lengths = zip(*contig_pairs)
        contigs = list(contigs)
        contig_lengths = list(contig_lengths)
    logger.info(f"Found {len(contigs)} contigs.")
    if scannerconfig.sampling_strategy == "largest_contig":
        largest_idx = max(range(len(contigs)), key=lambda i: contig_lengths[i])
        sampling_contig = contigs[largest_idx]
        scannerconfig.sampling_contigs = [sampling_contig]
        scannerconfig.sampling_contig_lengths = {
            sampling_contig: contig_lengths[largest_idx]
        }
        scannerconfig.sampling_contig_targets = {
            sampling_contig: scannerconfig.sampling_target_aligned_bases
        }
    else:
        sampling_contig = contigs[0]
        if scannerconfig.sampling_strategy == "top_contigs_random":
            contig_pairs = sorted(
                zip(contigs, contig_lengths), key=lambda x: x[1], reverse=True
            )
            top_n = max(1, min(scannerconfig.sampling_top_n_contigs, len(contig_pairs)))
            sampling_contigs = [c for c, _ in contig_pairs[:top_n]]
            sampling_lengths = {c: l for c, l in contig_pairs[:top_n]}
            total_len = sum(sampling_lengths.values()) or 1
            targets = {}
            remaining = scannerconfig.sampling_target_aligned_bases
            for idx, contig in enumerate(sampling_contigs):
                if idx == len(sampling_contigs) - 1:
                    target = remaining
                else:
                    target = int(
                        round(
                            scannerconfig.sampling_target_aligned_bases
                            * (sampling_lengths[contig] / total_len)
                        )
                    )
                    remaining -= target
                targets[contig] = max(0, target)
            scannerconfig.sampling_contigs = sampling_contigs
            scannerconfig.sampling_contig_lengths = sampling_lengths
            scannerconfig.sampling_contig_targets = targets
        else:
            scannerconfig.sampling_contigs = [sampling_contig]
            scannerconfig.sampling_contig_lengths = {sampling_contig: contig_lengths[0]}
            scannerconfig.sampling_contig_targets = {
                sampling_contig: scannerconfig.sampling_target_aligned_bases
            }

    progress_columns = [
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("•"),
        TimeElapsedColumn(),
    ]

    total_interrogated_bases = 0
    total_stats = {
        "reads": 0,
        "candidates": 0,
        "passed": 0,
        "total_aligned_bases": 0,
        "sampling_bases_total": 0,
        "sampling_passable_by_type": {},
        "sampling_passable_by_motif": {},
        "tract_counts_by_motif": {},
        "type_counts": {},
    }
    in_memory_records = []

    status_queue = Queue()
    ui_queue = Queue()
    active_workers: dict[str, str] = {}
    worker_start_times: dict[str, float] = {}
    active_lock = Lock()

    def format_worker_table() -> Table:
        max_workers = 20
        console_width = Console().size.width
        min_cell_width = 24
        columns = max(1, min(4, console_width // (min_cell_width + 3)))
        display_workers = min(parallel_n, max_workers)
        cell_width = max(min_cell_width, (console_width // columns) - 3)
        contig_width = max(8, cell_width - 12)
        table = Table.grid(padding=(0, 1))
        for _ in range(columns):
            table.add_column(no_wrap=True, width=cell_width)
        with active_lock:
            entries = []
            for worker_id in range(1, display_workers + 1):
                worker_key = str(worker_id)
                contig = active_workers.get(worker_key, "idle")
                start_time = worker_start_times.get(worker_key)
                elapsed = int(time.monotonic() - start_time) if start_time else 0
                contig_label = contig[:contig_width]
                entries.append(
                    f"W{worker_id:<2} {contig_label:<{contig_width}} {elapsed:>5d}s"
                )
        for idx in range(0, len(entries), columns):
            row_cells = [cell.ljust(cell_width) for cell in entries[idx : idx + columns]]
            while len(row_cells) < columns:
                row_cells.append("")
            table.add_row(*row_cells)
        return table

    def drain_status_queue():
        while True:
            try:
                event, worker_name, contig_name = status_queue.get_nowait()
            except Empty:
                break
            with active_lock:
                worker_id = worker_name.split("-")[-1]
                if event == "start":
                    active_workers[worker_id] = contig_name
                    worker_start_times[worker_id] = time.monotonic()
                elif event == "done":
                    active_workers.pop(worker_id, None)
                    worker_start_times.pop(worker_id, None)

    progress = Progress(*progress_columns, transient=False, auto_refresh=False)
    stop_event = Event()
    quit_event = Event()

    def listen_for_quit():
        stdin = sys.__stdin__
        if stdin is None or not stdin.isatty():
            return
        fd = stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setcbreak(fd)
            while not stop_event.is_set():
                rlist, _, _ = select.select([stdin], [], [], 0.2)
                if rlist:
                    char = stdin.read(1)
                    if char.lower() == "q":
                        quit_event.set()
                        stop_event.set()
                        break
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)

    def status_updater():
        while not stop_event.is_set():
            total_advance = 0
            while True:
                try:
                    event, value = ui_queue.get_nowait()
                except Empty:
                    break
                if event == "advance":
                    total_advance += int(value)
            drain_status_queue()
            progress.update(task_id, advance=total_advance)
            worker_title = f"Workers ({min(parallel_n, 20)}/{parallel_n})"
            group = Group(
                Panel(progress, title="Progress", padding=(0, 1)),
                Panel.fit(format_worker_table(), title=worker_title, padding=(0, 1)),
            )
            live.update(group, refresh=True)
            stop_event.wait(1.0)

    worker_title = f"Workers ({min(parallel_n, 20)}/{parallel_n})"
    group = Group(
        Panel(progress, title="Progress", padding=(0, 1)),
        Panel.fit(format_worker_table(), title=worker_title, padding=(0, 1)),
    )
    with Live(group, refresh_per_second=1, transient=False) as live:
        task_id = progress.add_task(
            "[green]Scanning contigs...",
            total=len(contigs),
        )
        status_thread = Thread(target=status_updater, daemon=True)
        key_thread = Thread(target=listen_for_quit, daemon=True)
        status_thread.start()
        key_thread.start()
        try:
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
                        if quit_event.is_set():
                            pool.terminate()
                            break
                        if result:
                            contig_name, contig_records, contig_base_count, stats = result
                            total_interrogated_bases += contig_base_count
                            total_stats["reads"] += stats["reads_processed"]
                            total_stats["candidates"] += stats["candidates"]
                            total_stats["passed"] += stats["passed"]
                            for k, v in stats["type_counts"].items():
                                total_stats["type_counts"][k] = (
                                    total_stats["type_counts"].get(
                                        k, 0) + v
                                )

                            if scannerconfig.in_memory:
                                in_memory_records.extend(
                                    contig_records or [])
                                if (
                                    scannerconfig.max_in_memory_records
                                    and len(in_memory_records)
                                    > scannerconfig.max_in_memory_records
                                ):
                                    logger.warning(
                                        "In-memory limit exceeded; consider disabling in-memory mode."
                                    )
                            ui_queue.put(("advance", 1))
        finally:
            stop_event.set()
            status_thread.join()
            key_thread.join(timeout=0.5)
            drain_status_queue()
            worker_title = f"Workers ({min(parallel_n, 20)}/{parallel_n})"
            group = Group(
                Panel(progress, title="Progress", padding=(0, 1)),
                Panel.fit(format_worker_table(), title=worker_title, padding=(0, 1)),
            )
            live.update(group, refresh=True)

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

    observed_lengths = _collect_indel_lengths(scannerconfig.passed_indels_path)
    scannerconfig.callable_lengths = observed_lengths

    callable_workers = min(2, parallel_n)
    logger.info("Starting callable-bases pass with %s workers...", callable_workers)
    with Pool(
        processes=callable_workers,
        initializer=_init_worker,
        initargs=(
            str(scannerconfig.bamfile),
            str(scannerconfig.fastafile),
            scannerconfig,
            sampling_contig,
            None,
        ),
    ) as pool:
        for result in pool.imap_unordered(_process_contig_callable, contigs):
            if result:
                _, stats = result
                total_stats["total_aligned_bases"] += stats["total_aligned_bases"]
                if stats.get("sampled_contig"):
                    total_stats["sampling_bases_total"] += stats["sampling_bases_total"]
                    for k, v in stats["sampling_passable_by_type"].items():
                        total_stats["sampling_passable_by_type"][k] = (
                            total_stats["sampling_passable_by_type"].get(k, 0) + v
                        )
                    for k, v in stats.get("sampling_passable_by_motif", {}).items():
                        total_stats["sampling_passable_by_motif"][k] = (
                            total_stats["sampling_passable_by_motif"].get(k, 0) + v
                        )
                    for k, v in stats.get("tract_counts_by_motif", {}).items():
                        total_stats["tract_counts_by_motif"][k] = (
                            total_stats["tract_counts_by_motif"].get(k, 0) + v
                        )

    return scannerconfig.passed_indels_path, total_interrogated_bases, total_stats
