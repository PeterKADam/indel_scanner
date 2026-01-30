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


class TestRepeatUnitDetection:
    def test_valid_read_rejects_soft_clips(self, mock_scanner_config, read_factory):
        scanner = ContigScanner(mock_scanner_config)
        read = read_factory(
            reference_start=100,
            cigartuples=[(4, 5), (0, 10)],  # 5S, 10M
            query_sequence="A" * 15,
        )
        read.query_qualities = [93] * len(read.query_sequence)
        assert not scanner._valid_read(read)

    def test_aligned_blocks_from_cigar(self, mock_scanner_config, read_factory):
        scanner = ContigScanner(mock_scanner_config)
        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 5), (1, 2), (0, 3), (2, 4), (0, 2)],
            query_sequence="A" * 12,
        )
        blocks = scanner._aligned_blocks(read)
        assert blocks == [(0, 100, 5), (7, 105, 3), (10, 112, 2)]

    def test_repeat_unit_length_detects_exact_repeats(self, mock_scanner_config):
        scanner = ContigScanner(mock_scanner_config)
        assert scanner._repeat_unit_length("ATATAT", 3) == 2
        assert scanner._repeat_unit_length("TATATA", 3) == 2
        assert scanner._repeat_unit_length("ACACACAC", 3) == 2

    def test_repeat_unit_length_rejects_interrupted(self, mock_scanner_config):
        scanner = ContigScanner(mock_scanner_config)
        assert scanner._repeat_unit_length("ATATGCAT", 3) is None
        assert scanner._repeat_unit_length("ATATATC", 3) is None

    def test_has_repeat_run_detects_local_repeat(self, mock_scanner_config):
        scanner = ContigScanner(mock_scanner_config)
        seq = "ATATCATATATATCATATATATC"
        assert scanner._has_repeat_run(seq, 3, 6) is True
        assert scanner._has_repeat_run("ACGTTGCAAC", 3, 6) is False


class TestStrCandidateFilter:
    def test_skips_str_like_insertion_near_rptrf(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = True
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (1, 6), (0, 10)],  # 10M, 6I, 10M
            query_sequence="A" * 10 + "ATATAT" + "C" * 10,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "N" * 200
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert results == []

    def test_keeps_non_str_like_insertion(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = False
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (1, 6), (0, 10)],  # 10M, 6I, 10M
            query_sequence="A" * 10 + "GATTAC" + "C" * 10,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "N" * 200
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert len(results) == 1

    def test_filters_insertion_with_repeat_units(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
            "local_window_bp": 30,
            "max_motif_len": 6,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = False
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        repeat_insert = "TACAA" * 3
        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (1, len(repeat_insert)), (0, 10)],
            query_sequence="A" * 10 + repeat_insert + "C" * 10,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "ACGT" * 40
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert results == []

    def test_skips_str_like_deletion_near_rptrf(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = True
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (2, 6), (0, 10)],  # 10M, 6D, 10M
            query_sequence="A" * 20,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "A" * 50 + "ATATCATATATATCATATATATC" + "A" * 50
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert results == []

    def test_keeps_deletion_when_not_near_rptrf(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
            "local_window_bp": 30,
            "max_motif_len": 6,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = False
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (2, 6), (0, 10)],  # 10M, 6D, 10M
            query_sequence="A" * 20,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "ACGTTGCAAGCTTGACCGTACGATCGATGCACTG"
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert len(results) == 1

    def test_filters_short_indel_if_local_repeat_run(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
            "local_window_bp": 30,
            "max_motif_len": 6,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = True
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=100,
            cigartuples=[(0, 10), (2, 2), (0, 10)],  # 10M, 2D, 10M
            query_sequence="A" * 20,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "A" * 50 + "ATATCATATATATCATATATATC" + "A" * 50
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert results == []

    def test_filters_specific_repeat_deletion_near_rptrf(self, mock_scanner_config, read_factory):
        mock_scanner_config.min_indel_size = 1
        mock_scanner_config.str_candidate_filter = {
            "enabled": True,
            "min_repeat_units": 3,
            "window_bp": 100,
            "local_window_bp": 30,
            "max_motif_len": 6,
        }
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier = Mock()
        mock_str_classifier.is_within_str_window.return_value = True
        mock_str_classifier.is_str_like.return_value = False
        mock_str_classifier.is_in_str.return_value = False
        mock_str_classifier.matches_rptrf_motif_length.return_value = False
        mock_str_classifier.motif_length_at.return_value = None

        read = read_factory(
            reference_start=95,
            cigartuples=[(0, 10), (2, 2), (0, 10)],  # 10M, 2D, 10M
            query_sequence="A" * 20,
        )
        read.query_qualities = [93] * len(read.query_sequence)

        contig_seq = "A" * 100 + "ATATCATATATATCATATATATC" + "A" * 100
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert results == []

        mock_scanner_config.str_candidate_filter = {"enabled": False}
        scanner = ContigScanner(mock_scanner_config)
        mock_str_classifier.is_within_str_window.return_value = False
        mock_str_classifier.is_str_like.return_value = False
        results = list(
            scanner._parse_cigar_for_candidates(read, mock_str_classifier, contig_seq)
        )
        assert len(results) == 1
        assert results[0].sequencecontext_brackets() == "ATATC[AT]ATATA"

