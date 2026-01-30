import logging
from indel_scanner.configurator import Config
from indel_scanner.parallel import parallel_pipeline
from indel_scanner.reporting import (
    build_callable_bases_by_type,
    write_callable_bases_report,
    write_per_type_mutation_report,
)
from indel_scanner.log import setup_logging

logger = logging.getLogger(__name__)


def main():
    setup_logging(enable_console=False)
    config = Config.load()
    setup_logging(config.log_dir, enable_console=False)
    config.write_settings_file()
    logger.info(f"Starting pipeline with {config.num_processes} workers...")
    passed_indels_path, total_interrogated_bases, stats = parallel_pipeline(config)
    logger.info(
        "Pipeline complete. Total L_interrogated (denominator) = "
        f"{total_interrogated_bases}"
    )
    logger.info(
        "Pipeline stats: reads=%s candidates=%s passed=%s",
        stats["reads"],
        stats["candidates"],
        stats["passed"],
    )

    callable_bases_by_type = build_callable_bases_by_type(stats, config.snp_label)
    type_counts = dict(stats["type_counts"])
    type_counts.setdefault(config.snp_label, 0)

    logger.info(
        "Sampling total=%s; callable_by_type=%s",
        stats["sampling_bases_total"],
        callable_bases_by_type,
    )
    logger.info(
        "STR tracts by motif length=%s (RPTRF tract counts)",
        stats.get("tract_counts_by_motif", {}),
    )

    write_per_type_mutation_report(
        passed_indels_path=passed_indels_path,
        callable_bases_by_type=callable_bases_by_type,
        output_dir=config.output_path,
        report_filename=config.per_type_report_filename,
        indel_bins=config.indel_bins,
    )
    write_callable_bases_report(
        callable_bases_by_type=callable_bases_by_type,
        output_dir=config.output_path,
        report_filename=config.callable_bases_filename,
    )


if __name__ == "__main__":
    main()
