import csv
import os
from pathlib import Path
from typing import Generator, Optional
import pysam
import pyfastx
import time
import logging

from .indel import Indel, IndelTsvFormatter
from .IO import write_records_to_tsv, cleanup_temp_dir
from .STR_Classifier import STRClassifier
from .indel import INDEL_TYPE, TSV_HEADERS
from .utils import Cigar, as_cigar
from .configurator import ScannerConfig

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

class ContigScanner:
    """
    A class to encapsulate the configuration and logic for scanning contigs.
    An instance of this class can be used as a worker in a multiprocessing pool.
    """

    def __init__(self, config: ScannerConfig):
        """
        Initializes the scanner with shared configuration that will be
        available to all calls within this worker process.
        """
        self.config = config

    def scan_contig(self, contig_name: str) -> Optional[tuple[Path, str]]:
        try:
            with pysam.AlignmentFile(str(self.config.bamfile), "rb") as samfile:
                temp_output_path = self.config.temp_dir / f"{contig_name}.part.tsv"

                fasta = pyfastx.Fasta(
                    str(self.config.fastafile)
                )  # dosnt support context manager????
                logger.info(
                    f"Scanning contig: {contig_name} @ {temp_output_path} \n with {(str(self.config.bamfile))} and {(str(self.config.fastafile))}"
                )
                return self._process_reads(
                    samfile, fasta, contig_name, temp_output_path
                )

        except Exception as e:
            logger.error(
                f"Error opening or processing files for contig {contig_name}: {e}"
            )
            return None

    def _parse_cigar(
        self,
        read: pysam.AlignedSegment,
        fasta: pyfastx.Fasta,
        contig_name: str,
        str_classifier: STRClassifier,
        contig_seq,
    ) -> Generator[Indel, None, None]:
        ref_pos_tracker = read.reference_start
        read_pos_tracker = 0

        if read is None or read.query_sequence is None or read.cigartuples is None:
            logger.warning(f"Skipping read with missing data: {read.query_name}")
            return

        for op_int, length in read.cigartuples:
            op = as_cigar(op_int)
            logger.debug(
                f"Processing CIGAR operation {op} with length {length} at ref pos {ref_pos_tracker}, query pos {read_pos_tracker} in read {read.query_name}"
            )

            if (op not in READ_CONSUMING_OPS) and (op not in REF_CONSUMING_OPS):
                logger.warning(f"op {op} was not caught")

            elif op == Cigar.OP_I:
                logger.debug(
                    f"Found insertion of length {length} at ref pos {ref_pos_tracker}, query pos {read_pos_tracker} in read {read.query_name}"
                )
                if length >= self.config.min_indel_size:
                    logger.debug(
                        f"insertion passed min size {self.config.min_indel_size}"
                    )

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

    def _generate_indels_from_contig(
        self,
        samfile: pysam.AlignmentFile,
        fasta: pyfastx.Fasta,
        contig_name: str,
        str_classifier: STRClassifier,
        contig_seq,
    ) -> Generator[Indel, None, None]:
        for read in samfile.fetch(contig=contig_name):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.query_sequence is None
            ):
                continue

            yield from self._parse_cigar(
                read, fasta, contig_name, str_classifier, contig_seq
            )

    def _process_reads(
        self,
        samfile: pysam.AlignmentFile,
        fasta: pyfastx.Fasta,
        contig_name: str,
        temp_output_path: Path,
    ):
        start_time = time.time()

        str_classifier = STRClassifier(self.config, contig_name)
        contig_seq = fasta[contig_name].seq

        indel_generator = self._generate_indels_from_contig(
            samfile, fasta, contig_name, str_classifier, contig_seq
        )

        write_records_to_tsv(
            output_path=temp_output_path,
            records=indel_generator,
            row_converter=lambda indel: IndelTsvFormatter.format_scanner_row(indel),
            buffer_size=self.config.buffer_size,
            # No header or sorting for this temporary file
        )

        end_time = time.time()
        return (
            temp_output_path,
            f"Finished {contig_name} in {end_time - start_time:.2f} s",
        )


def run_scan(scanner: ContigScanner, contig_name):
    return scanner.scan_contig(contig_name)


def aggregate_partial_results(temp_dir, final_output_path):
    logger.info("\nAll contigs processed. Merging results...")
    final_output_path = Path(final_output_path) / "indel_scanner_results.tsv"
    logger.info(f"Aggregating partial results from {temp_dir} into {final_output_path}")
    with open(final_output_path, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")

        # header
        writer.writerow(TSV_HEADERS.SCANNER.value)

        temp_dir_list = os.listdir(temp_dir)
        logger.debug(
            f"Writing ({len(temp_dir_list)}) partial results  to final output file..."
        )
        for part_file in temp_dir_list:
            if part_file.endswith(".part.tsv"):
                part_path = os.path.join(temp_dir, part_file)
                with open(part_path, "r") as f_in:
                    reader = csv.reader(f_in, delimiter="\t")
                    for row in reader:
                        if len(row) != 9:
                            logger.info(f"Skipping malformed row in {part_file}: {row}")
                            continue
                        writer.writerow(row)

    logger.debug(f"Aggregation complete. Final output written to {final_output_path}")
    # cleanup_temp_dir(temp_dir)
