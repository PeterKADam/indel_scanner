# indel_scanner/processor.py
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
from .indel import INDEL_TYPE, Indel, Insertion, Deletion
from .utils import Cigar

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}


class Processor:
    def __init__(self, config: ProcessorConfig):
        self.config = config
        self.reference = None
        self.bam_handle = None
        self.indels_by_contig: Dict[str, Dict[str, List[Indel]]] = defaultdict(lambda: defaultdict(list))
        self.passed_indels: List[Indel] = []
        self.total_records_processed = 0

    def __enter__(self):
        try:
            self.bam_handle = pysam.AlignmentFile(str(self.config.bamfile), "rb")
            self.reference = pyfastx.Fasta(str(self.config.fastafile))
            return self
        except Exception as e:
            logger.error(f"Error opening files: {e}")
            raise IOError("Failed to open input files.") from e

    def __exit__(self, exc_type, exc_value, traceback):
        if self.bam_handle:
            self.bam_handle.close()

    def _load_indels(self):
        logger.info(f"Loading and filtering indels from {self.config.input_file}...")
        try:
            lazy_df = pl.scan_csv(self.config.input_file, separator="\t", has_header=True)
            df = lazy_df.collect()
            for row in df.iter_rows(named=True):
                contig = row["contig"]
                position = row["ref_position"]
                length = row["length"]
                indel_type = row["type"]
                read_name = row["read_name"]
                map_quality = row["map_quality"]
                prefix, indel_seq, suffix = self._parse_sequence_context(row["[sequence]_context"])

                indel_content = ""
                if indel_type == INDEL_TYPE.INSERTION:
                    indel_content = indel_seq
                elif self.reference:
                    # For deletions, fetch content from reference
                    ref_end = position + length
                    indel_content = self.reference.fetch(contig, (position, ref_end))

                indel_obj = Indel.create(
                    contig=contig, ref_position=position, length=length,
                    prefix_context=prefix, suffix_context=suffix, read_name=read_name,
                    type=indel_type, indel_content=indel_content,
                    in_STR=row["in_STR"], map_quality=map_quality,
                )
                if indel_obj:
                    self.indels_by_contig[contig][read_name].append(indel_obj)

        except pl.exceptions.NoDataError:
            logger.warning(f"Input file is empty: {self.config.input_file}")
            return
        except Exception as e:
            logger.error(f"Failed to load indels with Polars: {e}")
            raise

        num_indels = sum(len(indels) for reads in self.indels_by_contig.values() for indels in reads.values())
        logger.info(f"Loaded {num_indels} indels across {len(self.indels_by_contig)} contigs.")

    def _process_bam_file(self):
        if not self.bam_handle:
            raise IOError("BAM file handle is not open.")

        logger.info("Streaming BAM file and applying filters...")
        filter_configuration = [
            {"name": "is_in_str"},
            {"name": "is_homopolymer"},
            {"name": "is_adjacent_to_homopolymer"},
            {"name": "similar_indels_in_other_reads"},
            {"name": "low_minimum_indel_quality", "params": {"min_quality": 93}},
            {"name": "low_singlebase_flanking_quality", "params": {"min_flank_quality": 93}},
        ]

        for contig, reads_on_contig in self.indels_by_contig.items():
            logger.info(f"Processing contig: {contig}...")

            # 1. First pass over reads in contig to populate quality and create a flat list
            indels_for_this_contig = []
            for bam_read in self.bam_handle.fetch(contig):
                if bam_read.is_secondary or bam_read.is_supplementary:
                    continue
                if bam_read.query_name in reads_on_contig:
                    for indel_obj in reads_on_contig[bam_read.query_name]:
                        if indel_obj.type == INDEL_TYPE.INSERTION:
                            self._populate_insertion_quality(bam_read, indel_obj)
                        elif indel_obj.type == INDEL_TYPE.DELETION:
                            self._populate_deletion_flanking_quality(bam_read, indel_obj)
                        indels_for_this_contig.append(indel_obj)

            if not indels_for_this_contig:
                continue

            # 2. Build location map for singleton filter
            indel_location_counts = defaultdict(int)
            for indel in indels_for_this_contig:
                key = (indel.ref_position, indel.type, indel.length)
                indel_location_counts[key] += 1
            filter_context = {"location_map": indel_location_counts}

            # 3. Apply all filters
            for indel_obj in indels_for_this_contig:
                IndelFilters.apply(filter_configuration, indel_obj, **filter_context)
                if not indel_obj.filter_reason: # No filters were triggered
                    self.passed_indels.append(indel_obj)

        self.total_records_processed = len(self.passed_indels)
        logger.info(f"Finished processing. Total records passed filters: {self.total_records_processed}")

    def process(self):
        if not self.bam_handle or not self.reference:
            logger.error("Cannot run processing: BAM or FASTA file not accessible.")
            sys.exit(1)
        self._load_indels()
        self._process_bam_file()
        return self.passed_indels

    def _populate_insertion_quality(self, target_read: AlignedSegment, insertion: Insertion, flank_len: int = 5):
        qualities = target_read.query_qualities
        if qualities is None: return
        ref_pos_tracker = target_read.reference_start
        read_pos_tracker = 0
        found = False
        for op, length in target_read.cigartuples: # type: ignore
            if op == Cigar.OP_I and ref_pos_tracker == insertion.ref_position and length == insertion.length:
                prefix_start = max(0, read_pos_tracker - flank_len)
                insertion_end = read_pos_tracker + length
                suffix_end = min(len(qualities), insertion_end + flank_len)
                insertion.prefix_quality = list(qualities[prefix_start:read_pos_tracker])
                insertion.indel_quality = list(qualities[read_pos_tracker:insertion_end])
                insertion.suffix_quality = list(qualities[insertion_end:suffix_end])
                found = True
                break
            if op in REF_CONSUMING_OPS: ref_pos_tracker += length
            if op in READ_CONSUMING_OPS: read_pos_tracker += length
        if not found: logger.debug(f"CIGAR parsing failed for insertion {insertion.read_name}")

    def _populate_deletion_flanking_quality(self, target_read: AlignedSegment, deletion: Deletion, flank_len: int = 5):
        qualities = target_read.query_qualities
        if qualities is None: return
        ref_pos_tracker = target_read.reference_start
        read_pos_tracker = 0
        found = False
        for op, length in target_read.cigartuples: # type: ignore
            if op == Cigar.OP_D and ref_pos_tracker == deletion.ref_position and length == deletion.length:
                prefix_start = max(0, read_pos_tracker - flank_len)
                suffix_end = min(len(qualities), read_pos_tracker + flank_len)
                deletion.prefix_quality = list(qualities[prefix_start:read_pos_tracker])
                deletion.suffix_quality = list(qualities[read_pos_tracker:suffix_end])
                found = True
                break
            if op in REF_CONSUMING_OPS: ref_pos_tracker += length
            if op in READ_CONSUMING_OPS: read_pos_tracker += length
        if not found: logger.debug(f"CIGAR parsing failed for deletion {deletion.read_name}")

    def _parse_sequence_context(self, context_string: str) -> Tuple[str, str, str]:
        try:
            start_bracket = context_string.index("[")
            end_bracket = context_string.index("]")
            prefix = context_string[:start_bracket]
            indel_seq = context_string[start_bracket + 1 : end_bracket]
            suffix = context_string[end_bracket + 1 :]
            return prefix, indel_seq, suffix
        except ValueError:
            return "", "", ""

    def write_output(self):
        if self.config.passed_indels_path and self.passed_indels:
            logger.info(f"Writing {len(self.passed_indels)} passed indels to {self.config.passed_indels_path}...")
            write_records_with_polars(self.passed_indels, self.config.passed_indels_path)

def run_processor(config: ProcessorConfig):
    logger.info("Starting indel processing task...")
    try:
        with Processor(config) as processor:
            processor.process()
            processor.write_output()
        logger.info("Indel processing finished successfully.")
    except Exception as e:
        logger.error(f"An unexpected error occurred during indel processing: {e}", exc_info=True)
        sys.exit(1)
