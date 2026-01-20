import csv
import logging
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)


def write_mutation_frequency_report(
    passed_indels_path: Path,
    total_interrogated_bases: int,
    output_dir: Path,
    report_filename: str,
) -> dict:
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
        logger.warning("Total interrogated bases is zero. Cannot calculate frequency.")

    logger.info(" FINAL RESULTS")
    logger.info(f"Final high-confidence indels (N_indels): {n_indels}")
    logger.info(
        f"Total interrogated bases (L_interrogated): {total_interrogated_bases}"
    )
    logger.info(f"De Novo Mutation Frequency: {mutation_frequency_str}")

    report_path = output_dir / report_filename
    try:
        with open(report_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["N_indels", "L_interrogated", "Mutation_Frequency"])
            writer.writerow([n_indels, total_interrogated_bases, mutation_frequency_str])
        logger.info(f"Final results written to {report_path}")
    except Exception as e:
        logger.error(f"Failed to write final report file: {e}")

    return {
        "n_indels": n_indels,
        "total_interrogated_bases": total_interrogated_bases,
        "mutation_frequency": mutation_frequency,
        "mutation_frequency_str": mutation_frequency_str,
        "report_path": report_path,
    }


def write_per_type_mutation_report(
    type_counts: dict,
    callable_bases_by_type: dict,
    output_dir: Path,
    report_filename: str,
) -> Path:
    report_path = output_dir / report_filename
    try:
        with open(report_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["mutation_type", "count", "callable_bases", "frequency"])
            for mutation_type, count in sorted(type_counts.items()):
                callable_bases = callable_bases_by_type.get(mutation_type, 0.0)
                frequency = count / callable_bases if callable_bases > 0 else 0.0
                writer.writerow(
                    [mutation_type, count, f"{callable_bases:.0f}", f"{frequency:.10e}"]
                )
        logger.info(f"Per-type mutation report written to {report_path}")
    except Exception as e:
        logger.error(f"Failed to write per-type report: {e}")
    return report_path
