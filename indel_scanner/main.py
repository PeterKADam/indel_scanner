# indel_scanner/main.py
import argparse
import logging
import shutil
import subprocess
import csv
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
    """Main workflow orchestrator."""
    setup_logging()
    config = Config.load()

    if isinstance(config, (ScannerConfig, ProcessorConfig)):
        # --- 1. SCANNING STAGE ---
        # This stage is run for 'scan' and for 'scan --process'
        logger.info(f"Starting SCANNING stage with {config.args.workers} workers...")
        scanner_output_path, total_interrogated_bases = parallel_scan(config)
        logger.info(f"SCANNING stage complete. Total L_interrogated (denominator) = {total_interrogated_bases}")

        # --- 2. PROCESSING STAGE (Conditional) ---
        if config.args.command == "process" or config.args.process:
            logger.info("Starting PROCESSING stage...")

            # The processor needs a different config context, so we create it
            # ensuring input/output paths are correctly set.
            processor_config_data = config.yaml

            # The processor's input is the scanner's output
            config.args.input = scanner_output_path

            # Create a dedicated output directory for the processor
            processor_output_dir = config.output_path / "processed"
            config.args.output = processor_output_dir

            processor_config = ProcessorConfig(config.args, processor_config_data)
            run_processor(processor_config)

            logger.info("PROCESSING stage complete.")

            # --- 3. FINAL CALCULATION ---
            logger.info("Calculating final mutation frequency...")
            passed_indels_path = processor_config.passed_indels_path
            n_indels = 0
            try:
                # Count the number of rows in the final output file from the processor
                passed_df = pl.read_csv(passed_indels_path, separator='\t')
                n_indels = len(passed_df)
            except (FileNotFoundError, pl.exceptions.NoDataError):
                logger.warning(f"{passed_indels_path} not found or is empty. Assuming 0 passed indels.")
                n_indels = 0

            mutation_frequency = 0.0
            mutation_frequency_str = "NA"
            if total_interrogated_bases > 0:
                mutation_frequency = n_indels / total_interrogated_bases
                mutation_frequency_str = f"{mutation_frequency:.10e}"
            else:
                logger.warning("Total interrogated bases is zero. Cannot calculate frequency.")

            # --- 4. FINAL REPORTING ---
            logger.info(" FINAL RESULTS")
            logger.info(f"Final high-confidence indels (N_indels): {n_indels}")
            logger.info(f"Total interrogated bases (L_interrogated): {total_interrogated_bases}")
            logger.info(f"De Novo Mutation Frequency: {mutation_frequency_str}")

            # Write machine-readable output file
            report_path = config.output_path / "final_mutation_frequency.tsv"
            try:
                with open(report_path, "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(["N_indels", "L_interrogated", "Mutation_Frequency"])
                    writer.writerow([n_indels, total_interrogated_bases, mutation_frequency_str])
                logger.info(f"Final results written to {report_path}")
            except Exception as e:
                logger.error(f"Failed to write final report file: {e}")

    else:
        logger.error(f"Unknown config: {config}")
        sys.exit(1)

if __name__ == "__main__":
    main()
