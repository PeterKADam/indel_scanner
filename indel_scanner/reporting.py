import csv
import logging
from pathlib import Path

import polars as pl
import pandas as pd

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


def build_callable_bases_by_type(
    stats: dict,
    snp_label: str,
) -> tuple[dict[str, float], dict[str, int]]:
    sampling_total = stats["sampling_bases_total"]
    sampling_totals_by_type = stats.get("sampling_totals_by_type")
    callable_bases_by_type: dict[str, float] = {}
    for mutation_type, passable_count in stats["sampling_passable_by_type"].items():
        sampling_total_for_type = (
            sampling_totals_by_type.get(mutation_type, sampling_total)
            if sampling_totals_by_type
            else sampling_total
        )
        if sampling_total_for_type > 0:
            callable_bases_by_type[mutation_type] = (
                passable_count / sampling_total_for_type
            ) * stats["total_aligned_bases"]
        else:
            callable_bases_by_type[mutation_type] = 0.0

    str_region_counts_by_type: dict[str, int] = {}
    tract_counts_by_motif = stats.get("tract_counts_by_motif", {})
    for motif_len, tract_count in tract_counts_by_motif.items():
        key_ins = f"ins_motif_{motif_len}bp"
        key_del = f"del_motif_{motif_len}bp"
        str_region_counts_by_type[key_ins] = int(tract_count)
        str_region_counts_by_type[key_del] = int(tract_count)

    if snp_label not in callable_bases_by_type:
        callable_bases_by_type[snp_label] = 0.0

    return callable_bases_by_type, str_region_counts_by_type


def write_callable_bases_report(
    callable_bases_by_type: dict[str, float],
    str_region_counts_by_type: dict[str, int],
    output_dir: Path,
    report_filename: str,
) -> Path:
    report_path = output_dir / report_filename
    try:
        with open(report_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["mutation_type", "callable_bases", "str_region_count"])
            all_types = set(callable_bases_by_type) | set(str_region_counts_by_type)
            for mutation_type in sorted(all_types):
                callable_bases = callable_bases_by_type.get(mutation_type, 0.0)
                str_regions = str_region_counts_by_type.get(mutation_type, 0)
                writer.writerow(
                    [mutation_type, f"{callable_bases:.0f}", str_regions]
                )
        logger.info(f"Callable bases report written to {report_path}")
    except Exception as e:
        logger.error(f"Failed to write callable bases report: {e}")
    return report_path


def write_per_type_mutation_report(
    passed_indels_path: Path,
    callable_bases_by_type: dict[str, float],
    str_region_counts_by_type: dict[str, int],
    output_dir: Path,
    report_filename: str,
    indel_bins: list[dict],
) -> Path:
    report_path = output_dir / report_filename
    passed_df = None
    if not passed_indels_path.exists():
        logger.warning(
            f"{passed_indels_path} not found or is empty. Assuming 0 passed indels."
        )
    else:
        try:
            passed_df = pd.read_csv(passed_indels_path, sep="\t")
        except Exception as e:
            logger.error(f"Failed to read {passed_indels_path}: {e}")

    def length_label(length: int) -> str:
        return f"len_{length}bp"

    try:
        with open(report_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(
                [
                    "mutation_type",
                    "str_class",
                    "count",
                    "unique_sites",
                    "callable_bases",
                    "str_region_count",
                    "frequency",
                ]
            )
            if passed_df is None or passed_df.empty:
                return report_path

            df = passed_df.copy()
            df["type_prefix"] = df["type"].map(lambda t: "ins_" if t == "ins" else "del_")
            df["length_label"] = df["length"].map(length_label)
            if "in_STR_region" in df.columns and "in_STR" in df.columns:
                df["str_motif"] = df["in_STR_region"] & df["in_STR"]
            else:
                df["str_motif"] = False
            if "str_motif_length" in df.columns:
                df["motif_length"] = pd.to_numeric(
                    df["str_motif_length"], errors="coerce"
                )
            else:
                df["motif_length"] = pd.NA
            df["mutation_type"] = df["type_prefix"] + df["length_label"]
            motif_mask = df["str_motif"] & df["motif_length"].notna()
            if motif_mask.any():
                df.loc[motif_mask, "mutation_type"] = (
                    df.loc[motif_mask, "type_prefix"]
                    + "motif_"
                    + df.loc[motif_mask, "motif_length"].astype(int).astype(str)
                    + "bp"
                )
            df["str_class"] = df["str_motif"].map(lambda v: "STR_motif" if v else "non_STR")

            grouped = df.groupby(["mutation_type", "str_class"]).size().reset_index()
            grouped = grouped.rename(columns={0: "count"})
            unique_df = df.drop_duplicates(
                subset=["contig", "ref_position", "type", "length", "str_class", "mutation_type"]
            )
            unique_counts = (
                unique_df.groupby(["mutation_type", "str_class"])
                .size()
                .to_dict()
            )

            for _, row in grouped.iterrows():
                mutation_type = str(row["mutation_type"])
                count = int(row["count"])
                str_class = str(row["str_class"])
                callable_bases = callable_bases_by_type.get(mutation_type, 0.0)
                str_region_count = str_region_counts_by_type.get(mutation_type, 0)
                denominator = callable_bases
                frequency = count / denominator if denominator > 0 else 0.0
                unique_sites = int(unique_counts.get((mutation_type, str_class), 0))
                writer.writerow(
                    [
                        mutation_type,
                        str_class,
                        count,
                        unique_sites,
                        f"{callable_bases:.0f}",
                        str_region_count,
                        f"{frequency:.10e}",
                    ]
                )
        logger.info(f"Per-type mutation report written to {report_path}")
    except Exception as e:
        logger.error(f"Failed to write per-type report: {e}")
    return report_path
