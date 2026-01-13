import logging
import csv
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
import pysam
import polars as pl
from indel_scanner.configurator import Config, ScannerConfig, ProcessorConfig
from indel_scanner.parallel import parallel_scan
from indel_scanner.processor import run_processor
from indel_scanner.log import setup_logging

logger = logging.getLogger(__name__)


def main():

    setup_logging()
    config = Config.load()

    if isinstance(config, (ScannerConfig, ProcessorConfig)):

        logger.info(f"Starting SCANNING stage with {config.args.workers} workers...")
        scanner_output_path, total_interrogated_bases = parallel_scan(config)
        logger.info(
            f"SCANNING stage complete. Total L_interrogated (denominator) = {total_interrogated_bases}"
        )


        if config.args.command == "process" or config.args.process:
            logger.info("Starting PROCESSING stage...")

            processor_config_data = config.yaml

            config.args.input = scanner_output_path

            processor_output_dir = config.output_path / "processed"
            config.args.output = processor_output_dir

            processor_config = ProcessorConfig(config.args, processor_config_data)
            run_processor(processor_config)

            logger.info("PROCESSING stage complete.")

            logger.info("Calculating mutation frequency...")
            passed_indels_path = processor_config.passed_indels_path
            n_indels = 0
            try:
                passed_df = pl.read_csv(passed_indels_path, separator="\t")
                n_indels = len(passed_df)
            except (FileNotFoundError, pl.exceptions.NoDataError):
                logger.warning(
                    f"{passed_indels_path} not found or is empty. Assuming 0 passed indels."
                )
                n_indels = 0

            mutation_frequency = 0.0
            mutation_frequency_str = "NA"
            if total_interrogated_bases > 0:
                mutation_frequency = n_indels / total_interrogated_bases
                mutation_frequency_str = f"{mutation_frequency:.10e}"
            else:
                logger.warning(
                    "Total interrogated bases is zero. Cannot calculate frequency."
                )

            logger.info(" FINAL RESULTS")
            logger.info(f"Final high-confidence indels (N_indels): {n_indels}")
            logger.info(
                f"Total interrogated bases (L_interrogated): {total_interrogated_bases}"
            )
            logger.info(f"De Novo Mutation Frequency: {mutation_frequency_str}")

            report_path = config.output_path / "final_mutation_frequency.tsv"
            try:
                with open(report_path, "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(
                        ["N_indels", "L_interrogated", "Mutation_Frequency"]
                    )
                    writer.writerow(
                        [n_indels, total_interrogated_bases, mutation_frequency_str]
                    )
                logger.info(f"Final results written to {report_path}")
            except Exception as e:
                logger.error(f"Failed to write final report file: {e}")

    else:
        logger.error(f"Unknown config: {config}")
        sys.exit(1)


if __name__ == "__main__":
    main()
