# tests/test_scanner.py
import pytest
from unittest.mock import Mock, MagicMock, patch, mock_open, call
from pathlib import Path

# --- IMPORTANT: Adjust these imports to match your project structure ---
from indel_scanner.scanner import ContigScanner, aggregate_partial_results
from indel_scanner.indel import Insertion, Deletion
from indel_scanner.indel import TSV_HEADERS


# Fixtures from conftest.py are automatically available

class TestParseCigar:
    """ Unit tests for the core CIGAR parsing logic in _parse_cigar. """

    @pytest.fixture
    def mock_dependencies(self, mocker):
        """Mocks dependencies needed by _parse_cigar."""
        mock_fasta_contig = Mock()
        mock_fasta_contig.seq = "N" * 500
        mock_fasta = {"chr1": mock_fasta_contig}

        mock_str_classifier = Mock()
        mocker.patch('indel_scanner.scanner.STRClassifier', return_value=mock_str_classifier)
        return mock_fasta, mock_str_classifier

    def test_parses_simple_insertion(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_fasta, mock_str_classifier = mock_dependencies
        mock_str_classifier.is_in_str.return_value = False

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 50), (1, 7), (0, 50)],  # 50M, 7I, 50M
            query_sequence="A" * 50 + "GATTACA" + "C" * 50,
        )

        results = list(scanner._parse_cigar(read, mock_fasta, "chr1"))

        assert len(results) == 1
        indel = results[0]
        assert isinstance(indel, Insertion)
        assert indel.ref_position == 150
        assert indel.length == 7
        assert indel.indel_content == "GATTACA"
        mock_str_classifier.is_in_str.assert_called_once_with(150)

    def test_parses_simple_deletion(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_fasta, mock_str_classifier = mock_dependencies
        mock_str_classifier.is_in_str.return_value = True
        mock_fasta["chr1"].seq = "A" * 150 + "GATTACA" + "C" * 150

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 50), (2, 7), (0, 50)],  # 50M, 7D, 50M
            query_sequence="A" * 50 + "C" * 50,
        )

        results = list(scanner._parse_cigar(read, mock_fasta, "chr1"))

        assert len(results) == 1
        indel = results[0]
        assert isinstance(indel, Deletion)
        assert indel.ref_position == 150
        assert indel.length == 7
        assert indel.indel_content == "GATTACA"
        mock_str_classifier.is_in_str.assert_called_once_with(150)

    def test_ignores_indels_below_min_size(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_fasta, _ = mock_dependencies
        read = read_factory(
            cigartuples=[(0, 50), (1, 2), (0, 20), (2, 1), (0, 30)],
            query_sequence="A" * 102
        )
        results = list(scanner._parse_cigar(read, mock_fasta, "chr1"))
        assert len(results) == 0


class TestContigScannerOrchestration:
    """Tests the high-level orchestration methods of ContigScanner."""

    def test_scan_contig_happy_path(self, mocker, mock_scanner_config):
        m_pysam = mocker.patch('indel_scanner.scanner.pysam.AlignmentFile')
        m_pyfastx = mocker.patch('indel_scanner.scanner.pyfastx.Fasta')
        m_process_reads = mocker.patch('indel_scanner.scanner.ContigScanner._process_reads')

        mock_samfile = MagicMock()
        m_pysam.return_value.__enter__.return_value = mock_samfile
        mock_fasta_instance = Mock()
        m_pyfastx.return_value = mock_fasta_instance

        scanner = ContigScanner(mock_scanner_config)
        contig_name = "chr1"
        expected_temp_path = mock_scanner_config.temp_dir / f"{contig_name}.part.tsv"

        scanner.scan_contig(contig_name)

        m_pysam.assert_called_once_with(str(mock_scanner_config.bamfile), "rb")
        m_pyfastx.assert_called_once_with(str(mock_scanner_config.fastafile))
        m_process_reads.assert_called_once_with(
            mock_samfile,
            mock_fasta_instance,
            contig_name,
            expected_temp_path
        )

    def test_scan_contig_handles_exception(self, mocker, mock_scanner_config):
        mocker.patch('indel_scanner.scanner.pysam.AlignmentFile', side_effect=IOError("File not found"))
        m_logger = mocker.patch('indel_scanner.scanner.logger')
        scanner = ContigScanner(mock_scanner_config)

        result = scanner.scan_contig("chr1")

        assert result is None
        m_logger.error.assert_called_once()


class TestAggregation:
    """Tests the final aggregation and cleanup step."""

    def test_aggregate_partial_results(self, mocker, tmp_path):
        """
        GIVEN a temporary directory with several partial.tsv files
        WHEN aggregate_partial_results is called
        THEN it should merge their contents into a final output file and clean up.
        """
        # Arrange
        m_listdir = mocker.patch('indel_scanner.scanner.os.listdir',
                                 return_value=['p1.part.tsv', 'p2.part.tsv', 'other.log'])
        m_open = mocker.patch('builtins.open', new_callable=mock_open)
        m_csv_writer = mocker.patch('indel_scanner.scanner.csv.writer')
        m_cleanup = mocker.patch('indel_scanner.scanner.cleanup_temp_dir')

        temp_dir = tmp_path / "temp"
        temp_dir.mkdir()
        final_output_path = tmp_path / "final.tsv"

        # The keys here are strings, which is correct because the side_effect uses str(path)
        file_content_map = {
            str(temp_dir / 'p1.part.tsv'): "row1\tdata1\nrow2\tdata2",
            str(temp_dir / 'p2.part.tsv'): "row3\tdata3",
        }

        # This side effect dynamically returns the correct mock file handle based on path
        m_open.side_effect = lambda path, *args, **kwargs: mock_open(
            read_data=file_content_map.get(str(path), "")
        ).return_value

        mock_writer_instance = Mock()
        m_csv_writer.return_value = mock_writer_instance

        # Act
        aggregate_partial_results(temp_dir, final_output_path)

        # Assert
        m_listdir.assert_called_once_with(temp_dir)

        # Check that the final output file and the partial files were opened
        m_open.assert_any_call(final_output_path, "w", newline="")

        # --- FIX: Convert Path objects to strings to match the actual call ---
        m_open.assert_any_call(str(temp_dir / 'p1.part.tsv'), "r")
        m_open.assert_any_call(str(temp_dir / 'p2.part.tsv'), "r")

        # Check that the CSV writer wrote the header and all data rows
        mock_writer_instance.writerow.assert_any_call(TSV_HEADERS.SCANNER.value)
        mock_writer_instance.writerow.assert_any_call(['row1', 'data1'])
        mock_writer_instance.writerow.assert_any_call(['row2', 'data2'])
        mock_writer_instance.writerow.assert_any_call(['row3', 'data3'])
        assert mock_writer_instance.writerow.call_count == 4

        m_cleanup.assert_called_once_with(temp_dir)

