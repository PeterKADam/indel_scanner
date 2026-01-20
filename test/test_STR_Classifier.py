# tests/test_str_classifier.py
from pathlib import Path
from unittest.mock import MagicMock, mock_open
import pytest
from indel_scanner.configurator import PipelineConfig
from indel_scanner.STR_Classifier import STRClassifier


# A fixture to provide a mock PipelineConfig, avoiding repetition.
@pytest.fixture
def mock_config(tmp_path):
    """Provides a mock PipelineConfig pointing to a temporary directory."""
    config = MagicMock(spec=PipelineConfig)
    config.str_directory = tmp_path  # Use pytest's tmp_path for a realistic Path object
    return config


class TestSTRClassifier:
    def test_initialization_and_parsing_with_real_format(self, mocker, mock_config):
        """
        GIVEN a file content string that matches the user's real file format
        WHEN STRClassifier is initialized
        THEN it should correctly parse the file, skipping headers and filtering out homopolymers.
        """
        # Arrange: This data is synthetic but  mimics the real file format.
        realistic_file_content = """
        **********************************TRs Found by RPTRF**********************************
Start End Len Motif Size( Sequence )
77 85 9 2(TA) TATATATAT
98 115 18 6(ATTTAT) ATTTATATTTATATTTAT
200 220 20 20(G) GGGGGGGGGGGGGGGGGGGG
128 147 20 6(ATTTAT) ATTTATATTTATATTTATAT
"""

        # Mock the built-in 'open' to return our fake content
        mocker.patch("builtins.open", mock_open(read_data=realistic_file_content))

        # Act: Initialize the classifier
        classifier = STRClassifier(config=mock_config, contig="h1tg000001l")

        # Assert: Check that the internal lists were populated correctly.
        # The line starting with '*' and 'Start' should be skipped.
        # The empty line should be skipped.
        # The homopolymer '20(G)' at position 200 should be skipped.
        assert classifier.starts == (77, 98, 128)
        assert classifier.ends == (85, 115, 147)
        assert classifier.num_regions == 3

    def test_initialization_with_no_valid_regions(self, mocker, mock_config):
        """
        GIVEN a file containing only headers and homopolymers
        WHEN STRClassifier is initialized
        THEN it should handle the empty case gracefully.
        """
        fake_file_content = """
Start   End     Length  Motif
200     250     50      5(A)
300     350     50      4(G)
"""
        mocker.patch("builtins.open", mock_open(read_data=fake_file_content))
        classifier = STRClassifier(config=mock_config, contig="chr1")

        assert classifier.starts == ()
        assert classifier.ends == ()
        assert classifier.num_regions == 0

    def test_is_in_str_with_no_regions(self, mocker, mock_config):
        """
        GIVEN a classifier initialized with no STR regions
        WHEN is_in_str is called
        THEN it should always return False.
        """
        mocker.patch("builtins.open", mock_open(read_data=""))  # Empty file
        classifier = STRClassifier(config=mock_config, contig="chr1")

        assert not classifier.is_in_str(100)
        assert not classifier.is_in_str(500)

    # Use parametrize to efficiently test many boundary conditions
    @pytest.mark.parametrize("position, expected_result", [
        # Before all regions
        (99, False),
        # Exactly on the start of the first region
        (100, True),
        # Inside the first region
        (150, True),
        # Exactly on the end of the first region
        (200, True),
        # Between the two regions
        (250, False),
        # Exactly on the start of the second region
        (300, True),
        # Inside the second region
        (350, True),
        # Exactly on the end of the second region
        (400, True),
        # After all regions
        (401, False),
    ])
    def test_is_in_str_logic(self, mocker, mock_config, position, expected_result):
        """
        GIVEN a classifier with known regions
        WHEN is_in_str is called with various positions
        THEN it should return the correct boolean value.
        """
        # For this test, it's easier to mock the method that produces the data
        # rather than the file read itself. This isolates the `is_in_str` logic.
        mocker.patch.object(STRClassifier, '_read_repeat_regions', return_value=[
            (100, 200),
            (300, 400)
        ])

        classifier = STRClassifier(config=mock_config, contig="chr1")

        # We need to manually set these since we bypassed part of __init__
        classifier.starts, classifier.ends = zip(*[(100, 200), (300, 400)])
        classifier.num_regions = 2

        assert classifier.is_in_str(position) == expected_result

    def test_file_not_found(self, mocker, mock_config):
        """
        GIVEN a contig that does not have a corresponding STR file
        WHEN STRClassifier is initialized
        THEN it should raise a FileNotFoundError.
        """
        # The mock_open doesn't raise FileNotFoundError by default.
        # So we patch 'open' to raise the error instead.
        mocker.patch("builtins.open").side_effect = FileNotFoundError

        with pytest.raises(FileNotFoundError):
            STRClassifier(config=mock_config, contig="non_existent_contig")

