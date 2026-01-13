import csv
import logging
from pathlib import Path
import shutil
from typing import Callable, List, Any, Optional, Union, Generator
import polars as pl
from .indel import Insertion, Deletion, Indel

logger = logging.getLogger(__name__)


def write_records_with_polars(records: List[Union[Insertion, Deletion]], path: Path):
    """Converts records to a DataFrame and writes to a TSV file."""
    if not records:
        return

    data = [r.to_processor_dict() for r in records]
    df = pl.DataFrame(data)

    def format_list_col(col_name: str) -> pl.Expr:
        return (
            pl.col(col_name)
            .cast(pl.List(pl.String))
            .list.join(",")
            .fill_null("NA")
        )

    final_df = df.select(
        pl.col("contig"),
        pl.col("ref_position"),
        pl.col("type"),
        pl.col("length"),
        pl.col("sequence_context_brackets").alias("[sequence]_context"),
        pl.col("read_name"),
        pl.col("in_STR"),
        format_list_col("prefix_quality"),
        format_list_col(
            "insertion_quality"
        ),  # Will be null for deletions, handled by fill_null
        format_list_col("suffix_quality"),
        format_list_col("filter_reason").alias("filter_reason"),
        pl.col("map_quality"),
    )

    final_df.sort(["contig", "ref_position"]).write_csv(path, separator="\t")


def write_records_to_tsv(
    output_path: Path,
    records: Generator[Indel, None, None],
    row_converter: Callable[[Union[Insertion, Deletion]], List[Any]],
    header: Optional[List[str]] = None,
    buffer_size: int = 10000,
    sort_key: Optional[Callable[[Union[Insertion, Deletion]], Any]] = None,
) -> int:
    """
    Writes an iterable of records to a TSV file with buffering and optional sorting.

    Args:
            output_path: The path to the output file.
            records: An iterable (list or generator) of record objects.
            row_converter: A function that converts a record object to a list for the CSV writer.
            header: An optional list of strings for the header row.
            buffer_size: The number of rows to buffer in memory before writing.
            sort_key: An optional key function for sorting records before writing.
                              Note: Using this will load all records into memory.

    Returns:
            The total number of records written.
    """
    if sort_key:
        logger.debug(f"Sorting records for {output_path.name}...")
        records_to_write = sorted(records, key=sort_key)
    else:
        records_to_write = records

    record_count = 0
    try:
        with open(output_path, "w", newline="") as outfile:
            writer = csv.writer(outfile, delimiter="\t")
            if header:
                writer.writerow(header)

            buffer = []
            for record in records_to_write:
                buffer.append(row_converter(record))
                record_count += 1
                if len(buffer) >= buffer_size:
                    writer.writerows(buffer)
                    buffer.clear()

            # Write any remaining records in the buffer
            if buffer:
                writer.writerows(buffer)

        if record_count > 0:
            logger.info(f"Successfully wrote {record_count} records to {output_path}.")
        else:
            logger.info(
                f"No records to write for {output_path.name}, an empty file was created."
            )

        return record_count

    except IOError as e:
        logger.error(f"Failed to write output file at {output_path}: {e}")
        return 0


def cleanup_temp_dir(temp_dir):
    logger.debug(f"Cleaning up temporary files in {temp_dir}")
    shutil.rmtree(temp_dir)
