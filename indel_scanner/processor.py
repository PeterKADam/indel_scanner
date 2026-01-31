# indel_scanner/processor.py
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict
import pyfastx
import pysam
from pysam import AlignedSegment

from .IO import write_records_with_polars, iter_tsv_rows_in_batches
from .configurator import PipelineConfig
from .filters import IndelFilters
from .indel import INDEL_TYPE, Indel, Insertion, Deletion, TSV_HEADERS
from .utils import Cigar, flank_qualities_by_ref

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}


class Processor:
    def __init__(self, config: PipelineConfig, input_file: Path, output_path: Path):
        self.config = config
        self.input_file = input_file
        self.output_path = output_path
        self.reference = None
        self.bam_handle = None
        self.indels_by_contig: Dict[str, Dict[str, List[Indel]]] = defaultdict(lambda: defaultdict(list))
        self.passed_indels: List[Indel] = []
        self.total_records_processed = 0
        self.passed_indels_path = self.config.passed_indels_path
        self.filter_configuration = self.config.processor_filters

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
        logger.info(f"Loading and filtering indels from {self.input_file}...")
        try:
            for batch in iter_tsv_rows_in_batches(
                self.input_file, self.config.read_batch_size
            ):
                for row in batch:
                    if len(row) < len(TSV_HEADERS.SCANNER.value):
                        continue
                    contig = row[0]
                    position = int(row[1])
                    indel_type = row[2]
                    length = int(row[3])
                    sequence_context = row[4]
                    read_name = row[5]
                    in_str_region = row[6].lower() == "true"
                    in_str = row[7].lower() == "true"
                    motif_length = None
                    if len(row) > 8 and row[8] and row[8] != "NA":
                        try:
                            motif_length = int(row[8])
                        except ValueError:
                            motif_length = None
                    map_quality_idx = 13 if len(row) > 13 else 12
                    map_quality = (
                        int(row[map_quality_idx])
                        if len(row) > map_quality_idx and row[map_quality_idx]
                        else None
                    )
                    prefix, indel_seq, suffix = self._parse_sequence_context(
                        sequence_context
                    )

                    indel_content = ""
                    if indel_type == INDEL_TYPE.INSERTION:
                        indel_content = indel_seq
                    elif self.reference:
                        ref_end = position + length
                        indel_content = self.reference.fetch(contig, (position, ref_end))

                    indel_obj = Indel.create(
                        contig=contig,
                        ref_position=position,
                        length=length,
                        prefix_context=prefix,
                        suffix_context=suffix,
                        read_name=read_name,
                        type=indel_type,
                        indel_content=indel_content,
                        in_STR_region=in_str_region,
                        in_STR=in_str,
                        motif_length=motif_length,
                        map_quality=map_quality,
                    )
                    if indel_obj:
                        self.indels_by_contig[contig][read_name].append(indel_obj)

        except FileNotFoundError:
            logger.warning(f"Input file is empty or missing: {self.input_file}")
            return
        except Exception as e:
            logger.error(f"Failed to load indels from TSV: {e}")
            raise

        num_indels = sum(len(indels) for reads in self.indels_by_contig.values() for indels in reads.values())
        logger.info(f"Loaded {num_indels} indels across {len(self.indels_by_contig)} contigs.")

    def _process_bam_file(self):
        if not self.bam_handle:
            raise IOError("BAM file handle is not open.")

        logger.info("Streaming BAM file and applying filters...")
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
                IndelFilters.apply(self.filter_configuration, indel_obj, **filter_context)
                if not indel_obj.is_filtered(): # No filters were triggered
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
                insertion_end = read_pos_tracker + length
                insertion.prefix_quality, insertion.suffix_quality = flank_qualities_by_ref(
                    target_read,
                    insertion.ref_position,
                    insertion.length,
                    flank_len,
                    is_insertion=True,
                )
                insertion.indel_quality = list(qualities[read_pos_tracker:insertion_end])
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
                deletion.prefix_quality, deletion.suffix_quality = flank_qualities_by_ref(
                    target_read,
                    deletion.ref_position,
                    deletion.length,
                    flank_len,
                    is_insertion=False,
                )
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
        if self.passed_indels_path and self.passed_indels:
            logger.info(f"Writing {len(self.passed_indels)} passed indels to {self.passed_indels_path}...")
            write_records_with_polars(self.passed_indels, self.passed_indels_path)

def run_processor(config: PipelineConfig, input_file: Path, output_path: Path):
    logger.info("Starting indel processing task...")
    if not input_file.is_file():
        logger.error(f"Input file not found at: {input_file}")
        sys.exit(1)
    try:
        with Processor(config, input_file, output_path) as processor:
            processor.process()
            processor.write_output()
        logger.info("Indel processing finished successfully.")
    except Exception as e:
        logger.error(f"An unexpected error occurred during indel processing: {e}", exc_info=True)
        sys.exit(1)
