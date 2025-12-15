import logging
import sys
from collections import defaultdict
from typing import List, Tuple, Dict

import polars as pl
import pyfastx
import pysam
from pysam import AlignedSegment

from .IO import write_records_with_polars
from .configurator import ProcessorConfig
from .filters import IndelFilters
from .indel import INDEL_TYPE, Indel
from .indel import Insertion, Deletion
from .utils import Cigar

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {
    Cigar.OP_M,
    Cigar.OP_D,
    Cigar.OP_N,
    Cigar.OP_EQ,
    Cigar.OP_X,
}

READ_CONSUMING_OPS = {
    Cigar.OP_M,
    Cigar.OP_I,
    Cigar.OP_S,
    Cigar.OP_EQ,
    Cigar.OP_X,
}


class Processor:
    def __init__(self, config: ProcessorConfig):
        self.config = config
        self.reference = None  # Will be initialized in __enter__
        self.bam_handle = None  # Will be initialized in __enter__

        # Structure: {contig: {read_name: [list_of_indel_objects]}}
        self.indels_by_contig: Dict[str, Dict[str, List[Indel]]] = defaultdict(
            lambda: defaultdict(list)
        )

        self.insertion_output: List[Insertion] = []
        self.deletion_output: List[Deletion] = []
        self.total_records_processed: int = 0

    def __enter__(self):
        logger.debug("Entering Processor context manager...")
        try:
            logger.debug("Opening BAM file for processing...")
            self.bam_handle = pysam.AlignmentFile(str(self.config.bamfile), "rb")
            logger.info(f"Successfully opened BAM file: {self.config.bamfile}")

            logger.debug("Opening reference FASTA for processing...")
            self.reference = pyfastx.Fasta(str(self.config.fastafile))
            logger.info(f"Successfully opened reference FASTA: {self.config.fastafile}")
            return self
        except Exception as e:
            logger.error(f"Error opening files: {e}")
            raise IOError("Failed to open input files.") from e

    def __exit__(self, exc_type, exc_value, traceback):
        if self.bam_handle:
            self.bam_handle.close()
            logger.debug(f"Closed BAM file handle for {self.config.bamfile}")

    def _load_indels(self):
        logger.info(f"Loading and filtering indels from {self.config.input_file}...")

        try:
            lazy_df = pl.scan_csv(
                self.config.input_file, separator="\t", has_header=True
            )
            df = lazy_df.collect()

            for row in df.iter_rows(named=True):
                contig = row["contig"]
                position = row["ref_position"]
                length = row["length"]
                indel_type = row["type"]
                read_name = row["read_name"]
                map_quality = row["map_quality"]

                prefix, indel_seq, suffix = self._parse_sequence_context(
                    row["[sequence]_context"]
                )

                indel_obj = Indel.create(
                    contig=contig,
                    ref_position=position,
                    length=length,
                    prefix_context=prefix,
                    suffix_context=suffix,
                    read_name=read_name,
                    type=indel_type,
                    indel_content=(
                        indel_seq
                        if indel_type == INDEL_TYPE.INSERTION
                        else self.reference.fetch(
                            contig,
                            (position, position + length - 1)
                            if indel_type == INDEL_TYPE.DELETION
                            else None,
                        )
                    ),
                    in_STR=row["in_STR"],
                    map_quality=map_quality,
                )

                if indel_obj:
                    self.indels_by_contig[contig][read_name].append(indel_obj)
                    logger.info(f"added indel {contig}:{position}:{indel_type}")

        except pl.exceptions.NoDataError:
            logger.warning(f"Input file is empty: {self.config.input_file}")

            return

        except Exception as e:
            logger.error(f"Failed to load indels with Polars: {e}")

            raise

        num_indels = sum(
            len(indels)
            for reads in self.indels_by_contig.values()
            for indels in reads.values()
        )
        logger.info(
            f"Loaded {num_indels} indels across {len(self.indels_by_contig)} contigs."
        )

    def _process_bam_file(self):
        if not self.bam_handle:
            raise IOError("BAM file handle is not open.")

        logger.info("Streaming BAM file and processing reads contig-by-contig...")

        filter_configuration = [
            {"name": "is_in_str"},
            {"name": "is_homopolymer"},
            {"name": "is_adjacent_to_homopolymer"},
            {
                "name": "poor_mapping_quality",
                "params": {"min_mapq": self.config.get("min_mapq", 30)},
            },
            {"name": "similar_indels_in_other_reads"},
        ]

        for contig, indels_on_this_contig in self.indels_by_contig.items():
            logger.info(f"Processing contig: {contig}...")

            # ========================================================================
            # Populate a TEMPORARY list of indels for only this contig.
            # ========================================================================
            indels_for_this_contig = []
            for bam_read in self.bam_handle.fetch(contig):
                if bam_read.is_secondary or bam_read.is_supplementary:
                    continue

                if bam_read.query_name in indels_on_this_contig:
                    for indel_obj in indels_on_this_contig[bam_read.query_name]:
                        # indel_obj.map_quality = bam_read.mapping_quality
                        if indel_obj.type == INDEL_TYPE.INSERTION:
                            self._populate_insertion_quality(bam_read, indel_obj)
                        elif indel_obj.type == INDEL_TYPE.DELETION:
                            self._populate_deletion_flanking_quality(
                                bam_read, indel_obj
                            )

                        indels_for_this_contig.append(indel_obj)

            if not indels_for_this_contig:
                logger.info(f"No indels found on contig {contig}. Skipping.")
                continue

            # ========================================================================
            # Build location map for this contig's indels.
            # ========================================================================
            indel_location_counts = defaultdict(int)
            for indel in indels_for_this_contig:
                key = (
                    indel.ref_position,
                    indel.type,
                    indel.length,
                )  # Contig is constant here
                indel_location_counts[key] += 1

            filter_context = {"location_map": indel_location_counts}

            # ========================================================================
            # Apply filters
            # ========================================================================
            for indel_obj in indels_for_this_contig:
                # Pass the contig-specific context to the filter function
                IndelFilters.apply(filter_configuration, indel_obj, **filter_context)

                # Append to the final, class-level output lists
                if indel_obj.type == INDEL_TYPE.INSERTION:
                    self.insertion_output.append(indel_obj)
                elif indel_obj.type == INDEL_TYPE.DELETION:
                    self.deletion_output.append(indel_obj)

        self.total_records_processed = len(self.insertion_output) + len(
            self.deletion_output
        )
        logger.info(
            f"Finished processing. Total records: {self.total_records_processed}."
        )

    def process(self):
        if not self.bam_handle or not self.reference:
            logger.error("Cannot run processing: BAM or FASTA file not accessible.")
            sys.exit(1)

        self._load_indels()

        self._process_bam_file()

        logger.info(
            f"Total records processed and passed filters: {self.total_records_processed}"
        )
        return self.insertion_output, self.deletion_output

    def _populate_insertion_quality(
        self, target_read: AlignedSegment, insertion: Insertion, flank_len: int = 5
    ):
        """Populates quality fields of an Insertion object using a pre-fetched read."""
        qualities = target_read.query_qualities
        if qualities is None:
            logger.debug(
                f"Read '{target_read.query_name}' has no quality scores ('*')."
            )
            return

        ref_pos_tracker = target_read.reference_start
        read_pos_tracker = 0
        found_insertion = False

        for op, length in target_read.cigartuples:  # type: ignore
            if (
                op == Cigar.OP_I
                and ref_pos_tracker == insertion.ref_position
                and length == insertion.length
            ):
                prefix_start = max(0, read_pos_tracker - flank_len)
                insertion_end = read_pos_tracker + length
                suffix_end = min(len(qualities), insertion_end + flank_len)

                insertion.prefix_quality = list(
                    qualities[prefix_start:read_pos_tracker]
                )
                insertion.insertion_quality = list(
                    qualities[read_pos_tracker:insertion_end]
                )
                insertion.suffix_quality = list(qualities[insertion_end:suffix_end])
                found_insertion = True
                break

            if op in REF_CONSUMING_OPS:
                ref_pos_tracker += length
            if op in READ_CONSUMING_OPS:
                read_pos_tracker += length

        if not found_insertion:
            logger.debug(
                f"CIGAR parsing failed to find a matching insertion for {insertion.read_name}"
            )

    def _populate_deletion_flanking_quality(
        self, target_read: AlignedSegment, deletion: Deletion, flank_len: int = 5
    ):
        """Populates quality fields of a Deletion object using a pre-fetched read."""
        qualities = target_read.query_qualities
        if qualities is None:
            logger.debug(
                f"Read '{target_read.query_name}' has no quality scores ('*')."
            )
            return

        ref_pos_tracker = target_read.reference_start
        read_pos_tracker = 0
        found_deletion = False

        for op, length in target_read.cigartuples:  # type: ignore
            if (
                op == Cigar.OP_D
                and ref_pos_tracker == deletion.ref_position
                and length == deletion.length
            ):
                prefix_start = max(0, read_pos_tracker - flank_len)
                suffix_end = min(len(qualities), read_pos_tracker + flank_len)

                deletion.prefix_quality = list(qualities[prefix_start:read_pos_tracker])
                deletion.suffix_quality = list(qualities[read_pos_tracker:suffix_end])
                found_deletion = True
                break

            if op in REF_CONSUMING_OPS:
                ref_pos_tracker += length
            if op in READ_CONSUMING_OPS:
                read_pos_tracker += length

        if not found_deletion:
            logger.debug(
                f"CIGAR parsing failed to find a matching deletion for {deletion.read_name}"
            )

    def _parse_sequence_context(self, context_string: str) -> Tuple[str, str, str]:
        try:
            start_bracket = context_string.index("[")
            end_bracket = context_string.index("]")
            if start_bracket < 0 or end_bracket < 0 or start_bracket >= end_bracket:
                raise ValueError("Invalid sequence context format.")
            prefix = context_string[:start_bracket]
            indel_seq = context_string[start_bracket + 1 : end_bracket]
            suffix = context_string[end_bracket + 1 :]
            return prefix, indel_seq, suffix
        except ValueError as e:
            logger.error(f"Failed to parse sequence context '{context_string}': {e}")
            return "", "", ""

    def write_output(self):
        if self.config.insertions_path and self.insertion_output:
            logger.info(
                f"Writing {len(self.insertion_output)} insertions to {self.config.insertions_path} using Polars..."
            )
            write_records_with_polars(
                self.insertion_output, self.config.insertions_path
            )

        if self.config.deletions_path and self.deletion_output:
            logger.info(
                f"Writing {len(self.deletion_output)} deletions to {self.config.deletions_path} using Polars..."
            )
            write_records_with_polars(self.deletion_output, self.config.deletions_path)


def run_processor(config):
    logger.info("Starting indel processing task...")
    try:
        logger.info("Initializing Processor...")
        with Processor(config) as processor:
            logger.info("Processor initialized successfully.")
            processor.process()
            processor.write_output()
        logger.info("Indel processing finished successfully.")
    except IOError as e:
        logger.error(f"Indel processing failed: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(
            f"An unexpected error occurred during indel processing: {e}", exc_info=True
        )
        sys.exit(1)
