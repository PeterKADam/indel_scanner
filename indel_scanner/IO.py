import csv
import logging
from pathlib import Path
import shutil
from typing import Callable, List, Any, Optional, Union, Generator, Iterable
import polars as pl
from .indel import Insertion, Deletion, TSV_HEADERS, IndelRecord

logger = logging.getLogger(__name__)

# Coordinate contract:
# - Internal `ref_position` is 0-based and left-normalized.
# - TSV output `ref_position` is 1-based and left-normalized for IGV-friendly display.

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
        (pl.col("ref_position") + 1).alias("ref_position"),
        pl.col("type"),
        pl.col("length"),
        pl.col("sequence_context_brackets").alias("[sequence]_context"),
        pl.col("read_name"),
        pl.col("in_STR_region"),
        pl.col("in_STR"),
        pl.col("str_motif_length"),
        pl.col("filter_reason").alias("filter_reason"),
        format_list_col("prefix_quality"),
        format_list_col(
            "insertion_quality"
        ),  # Will be null for deletions, handled by fill_null
        format_list_col("suffix_quality"),
        pl.col("map_quality"),
    )

    final_df.sort(["contig", "ref_position"]).write_csv(path, separator="\t")


def write_records_to_tsv(
    output_path: Path,
    records: Iterable[Any],
    row_converter: Callable[[Any], List[Any]],
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


def candidate_to_scanner_row(record: IndelRecord) -> List[str]:
    return [
        record.contig,
        str(record.ref_position + 1),
        record.type.value,
        str(record.length),
        f"{record.prefix_context}[{record.indel_content}]{record.suffix_context}",
        record.read_name,
        str(record.in_STR_region),
        str(record.in_STR),
        str(record.motif_length) if record.motif_length is not None else "NA",
        "",
    ]


def passed_to_processor_row(record: IndelRecord) -> List[str]:
    def format_list(scores: Optional[List[int | None]]) -> str:
        if scores is None:
            return "NA"
        return ",".join("NA" if s is None else str(s) for s in scores)

    return [
        record.contig,
        str(record.ref_position + 1),
        record.type.value,
        str(record.length),
        f"{record.prefix_context}[{record.indel_content}]{record.suffix_context}",
        record.read_name,
        str(record.in_STR_region),
        str(record.in_STR),
        str(record.motif_length) if record.motif_length is not None else "NA",
        "",
        format_list(record.prefix_quality),
        format_list(record.indel_quality),
        format_list(record.suffix_quality),
        str(record.map_quality),
    ]


def write_candidate_records(
    output_path: Path, records: Iterable[IndelRecord], buffer_size: int
) -> int:
    return write_records_to_tsv(
        output_path=output_path,
        records=records,
        row_converter=candidate_to_scanner_row,
        header=TSV_HEADERS.SCANNER.value,
        buffer_size=buffer_size,
    )


def write_passed_indels(
    output_path: Path, records: Iterable[IndelRecord], buffer_size: int
) -> int:
    return write_records_to_tsv(
        output_path=output_path,
        records=records,
        row_converter=passed_to_processor_row,
        header=TSV_HEADERS.PROCESSOR.value,
        buffer_size=buffer_size,
    )


def aggregate_tsv_parts(
    parts_dir: Path, final_output_path: Path, header: List[str]
) -> None:
    logger.info(
        f"Aggregating partial results from {parts_dir} into {final_output_path}"
    )
    with open(final_output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(header)
        temp_dir_list = sorted(parts_dir.iterdir())
        logger.debug(
            f"Writing ({len(temp_dir_list)}) partial results to final output file..."
        )
        for part_path in temp_dir_list:
            if part_path.name.endswith(".part.tsv"):
                with open(part_path, "r") as f_in:
                    reader = csv.reader(f_in, delimiter="\t")
                    for row in reader:
                        if row == header:
                            continue
                        writer.writerow(row)


def iter_tsv_rows_in_batches(
    input_path: Path, batch_size: int
) -> Generator[List[List[str]], None, None]:
    with open(input_path, "r") as f_in:
        reader = csv.reader(f_in, delimiter="\t")
        next(reader, None)
        batch: List[List[str]] = []
        for row in reader:
            batch.append(row)
            if len(batch) >= batch_size:
                yield batch
                batch = []
        if batch:
            yield batch


def cleanup_temp_dir(temp_dir):
    logger.debug(f"Cleaning up temporary files in {temp_dir}")
    try:
        shutil.rmtree(temp_dir)
    except FileNotFoundError:
        pass
