from __future__ import annotations

import argparse
import csv
import json
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Iterable

import pysam
import yaml

from indel_scanner.configurator import PipelineConfig
from indel_scanner.log import setup_logging
from indel_scanner.parallel import _init_worker, _process_contig_callable
from indel_scanner.reporting import (
    build_callable_bases_by_type,
    write_callable_bases_report,
)


logger = logging.getLogger(__name__)


def _select_sampling_contigs(contigs: list[str], lengths: list[int], config) -> str:
    if not contigs:
        raise ValueError("No contigs found in BAM.")

    if config.sampling_strategy == "largest_contig":
        largest_idx = max(range(len(contigs)), key=lambda i: lengths[i])
        sampling_contig = contigs[largest_idx]
        config.sampling_contigs = [sampling_contig]
        config.sampling_contig_lengths = {sampling_contig: lengths[largest_idx]}
        config.sampling_contig_targets = {
            sampling_contig: config.sampling_target_aligned_bases
        }
        return sampling_contig

    sampling_contig = contigs[0]
    if config.sampling_strategy == "top_contigs_random":
        contig_pairs = sorted(
            zip(contigs, lengths), key=lambda x: x[1], reverse=True
        )
        top_n = max(1, min(config.sampling_top_n_contigs, len(contig_pairs)))
        sampling_contigs = [c for c, _ in contig_pairs[:top_n]]
        sampling_lengths = {c: l for c, l in contig_pairs[:top_n]}
        total_len = sum(sampling_lengths.values()) or 1
        targets: dict[str, int] = {}
        remaining = config.sampling_target_aligned_bases
        for idx, contig in enumerate(sampling_contigs):
            if idx == len(sampling_contigs) - 1:
                target = remaining
            else:
                target = int(
                    round(
                        config.sampling_target_aligned_bases
                        * (sampling_lengths[contig] / total_len)
                    )
                )
                remaining -= target
            targets[contig] = max(0, target)
        config.sampling_contigs = sampling_contigs
        config.sampling_contig_lengths = sampling_lengths
        config.sampling_contig_targets = targets
    else:
        config.sampling_contigs = [sampling_contig]
        config.sampling_contig_lengths = {sampling_contig: lengths[0]}
        config.sampling_contig_targets = {
            sampling_contig: config.sampling_target_aligned_bases
        }
    return sampling_contig


def _load_config_data(config_path: Path) -> dict:
    if not config_path.exists():
        logger.warning(
            "Configuration file not found at %s. Using defaults.", config_path
        )
        return {}
    try:
        with open(config_path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    except Exception as exc:
        logger.error("Error loading configuration file %s: %s", config_path, exc)
        raise


def _build_pipeline_config(args: argparse.Namespace) -> PipelineConfig:
    config_data = _load_config_data(Path(args.config))
    return PipelineConfig(args, config_data)


def _iter_manifest_rows(manifest_path: Path) -> Iterable[dict]:
    with open(manifest_path, "r", newline="", encoding="utf-8") as handle:
        sample = handle.read(1024)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.DictReader(handle, dialect=dialect)
        for row in reader:
            yield row


def run_callable_only(
    config: PipelineConfig,
    workers: int | None,
    output_dir: Path | None,
) -> None:
    setup_logging(config.log_dir, enable_console=True)

    with pysam.AlignmentFile(str(config.bamfile), "rb") as samfile:
        contigs = list(samfile.header.references)
        contig_lengths = list(samfile.header.lengths)

    contig_filter = getattr(config, "contigs", None)
    if isinstance(contig_filter, (list, tuple, set)) and contig_filter:
        contig_lookup = dict(zip(contigs, contig_lengths))
        contigs = [c for c in contigs if c in contig_filter]
        contig_lengths = [contig_lookup[c] for c in contigs]

    if getattr(config, "test_run", False) and contigs:
        contig_pairs = list(zip(contigs, contig_lengths))
        target_count = min(
            len(contig_pairs), max(1, int(config.test_contigs_count))
        )
        contigs = [c for c, _ in contig_pairs[:target_count]]
        contig_lengths = [l for _, l in contig_pairs[:target_count]]

    config.callable_lengths = list(range(1, 11))
    sampling_contig = _select_sampling_contigs(contigs, contig_lengths, config)

    parallel_n = workers or config.num_processes or cpu_count()
    total_stats = {
        "total_aligned_bases": 0,
        "sampling_bases_total": 0,
        "sampling_passable_by_type": {},
        "sampling_totals_by_type": {},
        "tract_counts_by_motif": {},
    }

    logger.info(
        "Starting callable-only pass for 1-10bp with %s workers...",
        parallel_n,
    )

    with Pool(
        processes=parallel_n,
        initializer=_init_worker,
        initargs=(
            str(config.bamfile),
            str(config.fastafile),
            config,
            sampling_contig,
            None,
        ),
    ) as pool:
        for result in pool.imap_unordered(_process_contig_callable, contigs):
            if not result:
                continue
            _, stats = result
            total_stats["total_aligned_bases"] += stats.get(
                "total_aligned_bases", 0
            )
            if stats.get("sampled_contig"):
                total_stats["sampling_bases_total"] += stats.get(
                    "sampling_bases_total", 0
                )
                for k, v in stats.get("sampling_passable_by_type", {}).items():
                    total_stats["sampling_passable_by_type"][k] = (
                        total_stats["sampling_passable_by_type"].get(k, 0) + v
                    )
                for k, v in stats.get("sampling_totals_by_type", {}).items():
                    total_stats["sampling_totals_by_type"][k] = (
                        total_stats["sampling_totals_by_type"].get(k, 0) + v
                    )
            for k, v in stats.get("tract_counts_by_motif", {}).items():
                total_stats["tract_counts_by_motif"][k] = (
                    total_stats["tract_counts_by_motif"].get(k, 0) + v
                )

    callable_bases_by_type, str_region_counts_by_type = build_callable_bases_by_type(
        total_stats, config.snp_label
    )

    output_root = output_dir or config.output_path
    output_root.mkdir(parents=True, exist_ok=True)

    stats_path = output_root / "callable_stats_1_10.json"
    with open(stats_path, "w") as handle:
        json.dump(total_stats, handle, indent=2, sort_keys=True)

    callable_path = write_callable_bases_report(
        callable_bases_by_type=callable_bases_by_type,
        str_region_counts_by_type=str_region_counts_by_type,
        output_dir=output_root,
        report_filename="callable_bases_1_10.tsv",
    )

    logger.info("Callable-only stats written to %s", stats_path)
    logger.info("Callable-only bases report written to %s", callable_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Callable-only estimation for 1-10bp classes."
    )
    parser.add_argument(
        "-b", "--bam", required=False, help="Input BAM file (must be indexed)."
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=False,
        help="Reference FASTA file (must be indexed).",
    )
    parser.add_argument(
        "-o", "--output", required=False, help="Path to the output directory."
    )
    parser.add_argument(
        "-s",
        "--strdir",
        required=False,
        help="Directory containing STR results for the sample.",
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config.yaml",
        help="Path to the configuration file.",
    )
    parser.add_argument(
        "--contigs",
        default=None,
        help="Comma-separated list of contigs to process (default: all).",
    )
    parser.add_argument(
        "--test-run",
        action="store_true",
        help="Run on ~10 contigs around median length for quick testing.",
    )
    parser.add_argument(
        "--test-contigs-count",
        type=int,
        default=10,
        help="Number of contigs to include when --test-run is set.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes (defaults to config/CPU count).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override output directory (defaults to run output path).",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="CSV/TSV with columns: bam,fasta,output,strdir,config (optional).",
    )
    args = parser.parse_args()

    setup_logging(enable_console=True)

    if args.manifest:
        for row in _iter_manifest_rows(args.manifest):
            row_args = argparse.Namespace(
                bam=row.get("bam"),
                fasta=row.get("fasta"),
                output=row.get("output"),
                strdir=row.get("strdir"),
                config=row.get("config") or args.config,
                contigs=args.contigs,
                test_run=args.test_run,
                test_contigs_count=args.test_contigs_count,
            )
            config = _build_pipeline_config(row_args)
            run_callable_only(config, args.workers, args.output_dir)
        return

    missing = [k for k in ("bam", "fasta", "output", "strdir") if not getattr(args, k)]
    if missing:
        parser.error(f"Missing required args: {', '.join(missing)}")

    config = _build_pipeline_config(args)
    run_callable_only(config, args.workers, args.output_dir)


if __name__ == "__main__":
    main()
