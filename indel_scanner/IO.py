import csv
import logging
from pathlib import Path
import shutil
from typing import Iterable, Callable, List, Any, Optional, Union

from indel_scanner.utils import logger
from .indel import Insertion,Deletion

logger = logging.getLogger(__name__)

def _write_records_to_tsv(
	output_path: Path,
	records: Iterable[Union[Insertion, Deletion]],
	row_converter: Callable[[Union[Insertion, Deletion]], List[Any]],
	header: Optional[List[str]] = None,
	buffer_size: int = 10000,
	sort_key: Optional[Callable[[Union[Insertion, Deletion]], Any]] = None
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
		# Sorting requires consuming the entire iterable into a list.
		logger.debug(f"Sorting records for {output_path.name}...")
		records_to_write = sorted(list(records), key=sort_key)
	else:
		records_to_write = records

	record_count = 0
	try:
		with open(output_path, 'w', newline='') as outfile:
			writer = csv.writer(outfile, delimiter='\t')
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
			logger.info(f"No records to write for {output_path.name}, an empty file was created.")
		
		return record_count

	except IOError as e:
		logger.error(f"Failed to write output file at {output_path}: {e}")
		return 0


def cleanup_temp_dir(temp_dir):
	logger.debug(f"Cleaning up temporary files in {temp_dir}")
	shutil.rmtree(temp_dir)