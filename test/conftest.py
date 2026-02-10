# tests/conftest.py
import pytest
from unittest.mock import Mock, MagicMock
from pathlib import Path

# Mock the application's config object
@pytest.fixture
def mock_scanner_config():
    """Provides a mock PipelineConfig object for tests."""
    config = Mock()
    config.bamfile = Path("/fake/test.bam")
    config.fastafile = Path("/fake/test.fasta")
    config.temp_dir = Path("/fake/temp")
    config.output_path = Path("/fake/output.tsv")
    config.min_indel_size = 3  # Set a specific value for predictable tests
    config.write_buffer_size = 1000
    config.min_flank_quality = 93
    config.min_base_quality = 20
    config.min_map_quality = 0
    config.indel_bins = [
        {"label": "indel_1bp", "min": 1, "max": 1},
        {"label": "indel_2_3bp", "min": 2, "max": 3},
        {"label": "indel_4_10bp", "min": 4, "max": 10},
    ]
    config.str_candidate_filter = {
        "enabled": False,
        "min_repeat_units": 3,
        "window_bp": 100,
        "local_window_bp": 30,
        "max_motif_len": 6,
    }
    config.sampling_strategy = "largest_contig"
    config.sampling_bases = 1000000
    config.breakpoint_coherence_window_bp = 60
    config.snp_label = "snp"
    return config

@pytest.fixture
def read_factory():
    """
    A factory fixture to create mock pysam.AlignedSegment (read) objects.
    This is the most critical test double for testing the CIGAR parsing logic.
    """
    def _create_read(
        query_name="read1",
        reference_name="chr1",
        reference_start=100,
        cigartuples=None,
        query_sequence=None,
    ):
        def _aligned_pairs(
            matches_only=False,
        ):
            pairs = []
            read_pos = 0
            ref_pos = reference_start
            for op, length in (cigartuples or []):
                # M, =, X
                if op in (0, 7, 8):
                    for i in range(length):
                        pairs.append((read_pos + i, ref_pos + i))
                    read_pos += length
                    ref_pos += length
                # I
                elif op == 1:
                    if not matches_only:
                        for i in range(length):
                            pairs.append((read_pos + i, None))
                    read_pos += length
                # D, N
                elif op in (2, 3):
                    if not matches_only:
                        for i in range(length):
                            pairs.append((None, ref_pos + i))
                    ref_pos += length
                # S
                elif op == 4:
                    read_pos += length
            return pairs

        read = MagicMock()
        read.query_name = query_name
        read.reference_name = reference_name
        read.reference_start = reference_start
        read.cigartuples = cigartuples or []
        read.query_sequence = query_sequence
        read.query_qualities = (
            [93] * len(query_sequence) if query_sequence is not None else None
        )
        read.get_aligned_pairs.side_effect = _aligned_pairs
        # Ensure the mock doesn't complain about missing attributes
        read.is_unmapped = False
        read.is_secondary = False
        read.is_supplementary = False
        return read

    return _create_read

