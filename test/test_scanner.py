# tests/test_scanner.py
import pytest
from unittest.mock import Mock

# --- IMPORTANT: Adjust these imports to match your project structure ---
from indel_scanner.scanner import ContigScanner
from indel_scanner.indel import INDEL_TYPE


# Fixtures from conftest.py are automatically available

class TestParseCigar:
    """ Unit tests for the core CIGAR parsing logic in _parse_cigar. """

    @pytest.fixture
    def mock_dependencies(self):
        """Mocks dependencies needed by _parse_cigar."""
        mock_str_classifier = Mock()
        return mock_str_classifier

    def test_parses_simple_insertion(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = mock_dependencies
        mock_str_classifier.is_in_str.return_value = False

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 50), (1, 7), (0, 50)],  # 50M, 7I, 50M
            query_sequence="A" * 50 + "GATTACA" + "C" * 50,
        )

        contig_seq = "N" * 500
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )

        assert len(results) == 1
        indel = results[0]
        assert indel.type == INDEL_TYPE.INSERTION
        assert indel.ref_position == 150
        assert indel.length == 7
        assert indel.indel_content == "GATTACA"
        mock_str_classifier.is_in_str.assert_called_once_with(150)

    def test_parses_simple_deletion(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = mock_dependencies
        mock_str_classifier.is_in_str.return_value = True

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 50), (2, 7), (0, 50)],  # 50M, 7D, 50M
            query_sequence="A" * 50 + "C" * 50,
        )

        contig_seq = "A" * 150 + "GATTACA" + "C" * 150
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )

        assert len(results) == 1
        indel = results[0]
        assert indel.type == INDEL_TYPE.DELETION
        assert indel.ref_position == 150
        assert indel.length == 7
        assert indel.indel_content == "GATTACA"
        mock_str_classifier.is_in_str.assert_called_once_with(150)

    def test_ignores_indels_below_min_size(self, mock_scanner_config, read_factory, mock_dependencies):
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = mock_dependencies
        read = read_factory(
            cigartuples=[(0, 50), (1, 2), (0, 20), (2, 1), (0, 30)],
            query_sequence="A" * 102
        )
        contig_seq = "N" * 200
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert len(results) == 0


class TestIndelBinning:
    def test_classify_indel_bins(self, mock_scanner_config):
        scanner = ContigScanner(mock_scanner_config)
        assert scanner._classify_indel_bin(1) == "indel_1bp"
        assert scanner._classify_indel_bin(2) == "indel_2_3bp"
        assert scanner._classify_indel_bin(3) == "indel_2_3bp"
        assert scanner._classify_indel_bin(4) == "indel_4_10bp"
        assert scanner._classify_indel_bin(10) == "indel_4_10bp"
        assert scanner._classify_indel_bin(11) is None

