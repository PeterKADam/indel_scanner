import logging
from indel_scanner.configurator import Config
from indel_scanner.parallel import parallel_pipeline
from indel_scanner.reporting import (
    write_mutation_frequency_report,
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

    logger.info("Calculating mutation frequency...")
    write_mutation_frequency_report(
        passed_indels_path=passed_indels_path,
        total_interrogated_bases=total_interrogated_bases,
        output_dir=config.output_path,
        report_filename=config.report_filename,
    )

    sampling_total = stats["sampling_bases_total"]
    callable_bases_by_type: dict[str, float] = {}
    for mutation_type, passable_count in stats["sampling_passable_by_type"].items():
        if sampling_total > 0:
            callable_bases_by_type[mutation_type] = (
                passable_count / sampling_total
            ) * stats["total_aligned_bases"]
        else:
            callable_bases_by_type[mutation_type] = 0.0

    type_counts = dict(stats["type_counts"])
    type_counts.setdefault(config.snp_label, 0)
    if config.snp_label not in callable_bases_by_type:
        callable_bases_by_type[config.snp_label] = 0.0
    for bin_cfg in config.indel_bins:
        ins_label = f"ins_{bin_cfg['label']}"
        del_label = f"del_{bin_cfg['label']}"
        type_counts.setdefault(ins_label, 0)
        type_counts.setdefault(del_label, 0)
        callable_bases_by_type.setdefault(ins_label, 0.0)
        callable_bases_by_type.setdefault(del_label, 0.0)

    logger.info(
        "Sampling total=%s; callable_by_type=%s",
        sampling_total,
        callable_bases_by_type,
    )

    write_per_type_mutation_report(
        type_counts=type_counts,
        callable_bases_by_type=callable_bases_by_type,
        output_dir=config.output_path,
        report_filename=config.per_type_report_filename,
    )


if __name__ == "__main__":
    main()
