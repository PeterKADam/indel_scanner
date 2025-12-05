# tests/conftest.py
import pytest
from unittest.mock import Mock, MagicMock
from pathlib import Path

# Mock the application's config object
@pytest.fixture
def mock_scanner_config():
    """Provides a mock ScannerConfig object for tests."""
    config = Mock()
    config.bamfile = Path("/fake/test.bam")
    config.fastafile = Path("/fake/test.fasta")
    config.temp_dir = Path("/fake/temp")
    config.output_path = Path("/fake/output.tsv")
    config.min_indel_size = 3  # Set a specific value for predictable tests
    config.buffer_size = 1000
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
        read = MagicMock()
        read.query_name = query_name
        read.reference_name = reference_name
        read.reference_start = reference_start
        read.cigartuples = cigartuples or []
        read.query_sequence = query_sequence
        # Ensure the mock doesn't complain about missing attributes
        read.is_unmapped = False
        read.is_secondary = False
        read.is_supplementary = False
        return read

    return _create_read

