# indel_scanner/scanner.py
import csv
import os
from pathlib import Path
from typing import Generator, Optional
import pysam
import pyfastx
import time
import logging

from .indel import Indel, IndelTsvFormatter
from .IO import write_records_to_tsv
from .STR_Classifier import STRClassifier
from .filters import IndelFilters
from .indel import INDEL_TYPE, TSV_HEADERS
from .utils import Cigar, as_cigar
from .configurator import ScannerConfig

logger = logging.getLogger(__name__)

REF_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_D, Cigar.OP_N, Cigar.OP_EQ, Cigar.OP_X}
READ_CONSUMING_OPS = {Cigar.OP_M, Cigar.OP_I, Cigar.OP_S, Cigar.OP_EQ, Cigar.OP_X}


class ContigScanner:
    def __init__(self, config: ScannerConfig):
        self.config = config

    def scan_contig(self, contig_name: str) -> Optional[tuple[Path, str, int]]:
        try:
            with pysam.AlignmentFile(str(self.config.bamfile), "rb") as samfile:
                temp_output_path = self.config.temp_dir / f"{contig_name}.part.tsv"
                fasta = pyfastx.Fasta(str(self.config.fastafile))
                return self._process_reads_single_pass(
                    samfile, contig_name, temp_output_path, fasta
                )
        except Exception as e:
            logger.error(
                f"Error opening or processing files for contig {contig_name}: {e}"
            )
            return None

    def _parse_cigar_for_candidates(
        self,
        read: pysam.AlignedSegment,
        str_classifier: STRClassifier,
        contig_seq: str,
    ) -> Generator[Indel, None, None]:
        ref_pos_tracker = read.reference_start
        read_pos_tracker = 0
        if read.cigartuples is None:
            return

        for op_int, length in read.cigartuples:
            op = as_cigar(op_int)

            if op == Cigar.OP_I:
                if length >= self.config.min_indel_size:
                    yield Indel.create(
                        contig=read.reference_name,
                        ref_position=ref_pos_tracker,
                        type=INDEL_TYPE.INSERTION,
                        length=length,
                        prefix_context=read.query_sequence[
                            max(0, read_pos_tracker - 5) : read_pos_tracker
                        ],
                        suffix_context=read.query_sequence[
                            read_pos_tracker + length : read_pos_tracker + length + 5
                        ],
                        read_name=read.query_name,
                        in_STR=str_classifier.is_in_str(ref_pos_tracker),
                        map_quality=read.mapping_quality,
                        indel_content=read.query_sequence[
                            read_pos_tracker : read_pos_tracker + length
                        ],
                    )
            elif op == Cigar.OP_D:
                if length >= self.config.min_indel_size:
                    yield Indel.create(
                        contig=read.reference_name,
                        ref_position=ref_pos_tracker,
                        type=INDEL_TYPE.DELETION,
                        length=length,
                        prefix_context=contig_seq[
                            max(0, ref_pos_tracker - 5) : ref_pos_tracker
                        ],
                        indel_content=contig_seq[
                            ref_pos_tracker : ref_pos_tracker + length
                        ],
                        suffix_context=contig_seq[
                            ref_pos_tracker + length : ref_pos_tracker + length + 5
                        ],
                        read_name=read.query_name,
                        in_STR=str_classifier.is_in_str(ref_pos_tracker),
                        map_quality=read.mapping_quality,
                    )

            if op in REF_CONSUMING_OPS:
                ref_pos_tracker += length
            if op in READ_CONSUMING_OPS:
                read_pos_tracker += length

    def _is_callable_at_position(
        self,
        read: pysam.AlignedSegment,
        read_pos: int,
        ref_pos: int,
        str_classifier: STRClassifier,
        contig_seq: str,
    ) -> bool:
        if read_pos >= len(read.query_qualities) - 1:
            return False  # Cannot form a flank pair at the very end of a read.

        if (
            read.query_qualities[read_pos] < self.config.min_flank_quality
            or read.query_qualities[read_pos + 1] < self.config.min_flank_quality
        ):
            return False

        if IndelFilters.check_if_in_str(ref_pos, str_classifier):
            return False
        if IndelFilters.check_if_homopolymer_context(
            contig_seq, ref_pos, self.config.min_homopolymer_len
        ):
            return False

        return True

    def _process_reads_single_pass(
        self,
        samfile: pysam.AlignmentFile,
        contig_name: str,
        temp_output_path: Path,
        fasta: pyfastx.Fasta,
    ):
        start_time = time.time()
        str_classifier = STRClassifier(self.config, contig_name)
        contig_seq = fasta[contig_name].seq
        interrogated_bases_count = 0

        def indel_candidate_generator():
            nonlocal interrogated_bases_count
            for read in samfile.fetch(contig=contig_name):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.mapping_quality < self.config.min_map_quality
                    or read.query_sequence is None
                    or read.query_qualities is None
                ):
                    continue

                # --- Denominator Calculation ---
                for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if ref_pos is None:
                        continue
                    if self._is_callable_at_position(
                        read, read_pos, ref_pos, str_classifier, contig_seq
                    ):
                        interrogated_bases_count += 1

                # --- Numerator Candidate Generation ---
                yield from self._parse_cigar_for_candidates(
                    read, str_classifier, contig_seq
                )

        write_records_to_tsv(
            output_path=temp_output_path,
            records=indel_candidate_generator(),
            row_converter=lambda indel: IndelTsvFormatter.format_scanner_row(indel),
            buffer_size=self.config.buffer_size,
        )

        end_time = time.time()
        status_message = (
            f"Finished {contig_name} in {end_time - start_time:.2f} s. "
            f"Found {interrogated_bases_count} interrogated bases."
        )
        return (temp_output_path, status_message, interrogated_bases_count)


def run_scan(scanner: ContigScanner, contig_name):
    return scanner.scan_contig(contig_name)


def aggregate_partial_results(temp_dir, final_output_path):
    logger.info("\nAll contigs processed. Merging results...")
    final_output_path = Path(final_output_path) / "indel_scanner_results.tsv"
    logger.info(f"Aggregating partial results from {temp_dir} into {final_output_path}")

    with open(final_output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(TSV_HEADERS.SCANNER.value)
        temp_dir_list = os.listdir(temp_dir)
        logger.debug(
            f"Writing ({len(temp_dir_list)}) partial results to final output file..."
        )
        for part_file in temp_dir_list:
            if part_file.endswith(".part.tsv"):
                part_path = os.path.join(temp_dir, part_file)
                with open(part_path, "r") as f_in:
                    reader = csv.reader(f_in, delimiter="\t")
                    for row in reader:
                        writer.writerow(row)
    logger.debug(f"Aggregation complete. Final output written to {final_output_path}")
