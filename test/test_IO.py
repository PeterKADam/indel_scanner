# tests/test_io.py
import pytest
from unittest.mock import Mock, mock_open
from pathlib import Path
# Updated to import from your specified module path.
from indel_scanner.IO import _write_records_to_tsv, cleanup_temp_dir

@pytest.fixture
def mock_records():
    """
    Provides a list of mock record objects for testing.
    """
    record1 = Mock(name="record1")
    record1.contig = "chr1"
    record1.pos = 200
    record2 = Mock(name="record2")
    record2.contig = "chr1"
    record2.pos = 100
    record3 = Mock(name="record3")
    record3.contig = "chr2"
    record3.pos = 50
    return [record1, record2, record3]


def mock_row_converter(record):
    """A test double for the real row_converter function."""
    return [record.contig, record.pos]


class TestIOWriter:
    """Unit tests for the _write_records_to_tsv function using mocking."""

    def test_write_records_streaming_happy_path(self, mocker, mock_records):
        """
        GIVEN: A list of records without a sort key.
        WHEN: _write_records_to_tsv is called.
        THEN: It should open the correct file and write the header and all records in original order.
        """
        # Arrange
        m_open = mocker.patch("indel_scanner.IO.open", mock_open())
        m_logger = mocker.patch("indel_scanner.IO.logger")

        mock_writer_instance = Mock()
        m_csv_writer = mocker.patch("indel_scanner.IO.csv.writer")
        m_csv_writer.return_value = mock_writer_instance

        output_path = Path("/fake/output.tsv")
        header = ["Contig", "Position"]

        # Act
        count = _write_records_to_tsv(
            output_path=output_path,
            records=mock_records,
            row_converter=mock_row_converter,
            header=header,
            sort_key=None
        )

        # Assert
        assert count == 3
        m_open.assert_called_once_with(output_path, "w", newline="")

        mock_writer_instance.writerow.assert_called_once_with(header)

        # Correctly capture the argument passed to writerows
        mock_writer_instance.writerows.assert_called_once()
        written_rows = mock_writer_instance.writerows.call_args[0][0]

        expected_rows = [
            ['chr1', 200],
            ['chr1', 100],
            ['chr2', 50],
        ]
        assert written_rows == expected_rows
        m_logger.info.assert_called_with(f"Successfully wrote 3 records to {output_path}.")

    def test_write_records_with_sorting(self, mocker, mock_records):
        """
        GIVEN: A list of records and a sort key.
        WHEN: _write_records_to_tsv is called.
        THEN: It should log the sorting action and write all records in the correct sorted order.
        """
        # Arrange
        mocker.patch("indel_scanner.IO.open", mock_open())
        m_logger = mocker.patch("indel_scanner.IO.logger")

        mock_writer_instance = Mock()
        m_csv_writer = mocker.patch("indel_scanner.IO.csv.writer")
        m_csv_writer.return_value = mock_writer_instance

        output_path = Path("/fake/output_sorted.tsv")
        sort_func = lambda r: r.pos

        # Act
        count = _write_records_to_tsv(
            output_path=output_path,
            records=mock_records,
            row_converter=mock_row_converter,
            sort_key=sort_func
        )

        # Assert
        assert count == 3
        m_logger.debug.assert_called_with(f"Sorting records for {output_path.name}...")

        written_rows = mock_writer_instance.writerows.call_args[0][0]

        expected_rows = [
            ['chr2', 50],  # First
            ['chr1', 100],  # Second
            ['chr1', 200],  # Third
        ]
        assert written_rows == expected_rows

    def test_write_no_records(self, mocker):
        """
        GIVEN: An empty list of records.
        WHEN: _write_records_to_tsv is called.
        THEN: It should create an empty file with only a header and log the appropriate message.
        """
        # Arrange
        mocker.patch("indel_scanner.IO.open", mock_open())
        m_logger = mocker.patch("indel_scanner.IO.logger")

        mock_writer_instance = Mock()
        m_csv_writer = mocker.patch("indel_scanner.IO.csv.writer")
        m_csv_writer.return_value = mock_writer_instance

        output_path = Path("/fake/empty.tsv")
        header = ["Contig", "Position"]

        # Act
        count = _write_records_to_tsv(
            output_path=output_path,
            records=[],  # Empty list
            row_converter=mock_row_converter,
            header=header
        )

        # Assert
        assert count == 0

        mock_writer_instance.writerow.assert_called_once_with(header)
        mock_writer_instance.writerows.assert_not_called()

        m_logger.info.assert_called_with(
            f"No records to write for {output_path.name}, an empty file was created."
        )

    def test_io_error_handling(self, mocker, mock_records):
        """
        GIVEN: The open() call will raise an IOError.
        WHEN: _write_records_to_tsv is called.
        THEN: It should log the error and return 0.
        """
        # Arrange
        mocker.patch("indel_scanner.IO.open", side_effect=IOError("Permission denied"))
        m_logger = mocker.patch("indel_scanner.IO.logger")
        output_path = Path("/restricted/output.tsv")

        # Act
        count = _write_records_to_tsv(
            output_path=output_path,
            records=mock_records,
            row_converter=mock_row_converter
        )

        # Assert
        assert count == 0
        m_logger.error.assert_called_once_with(
            f"Failed to write output file at {output_path}: Permission denied"
        )


class TestCleanup:
    """Unit tests for the cleanup_temp_dir function."""

    def test_cleanup_temp_dir_calls_rmtree(self, mocker):
        """
        GIVEN: A directory path.
        WHEN: cleanup_temp_dir is called.
        THEN: It should call shutil.rmtree with the correct path.
        """
        # Arrange
        mock_rmtree = mocker.patch("indel_scanner.IO.shutil.rmtree")
        mock_logger = mocker.patch("indel_scanner.IO.logger")
        temp_dir = Path("/fake/temp/dir")

        # Act
        cleanup_temp_dir(temp_dir)

        # Assert
        mock_rmtree.assert_called_once_with(temp_dir)
        mock_logger.debug.assert_called_once_with(
            f"Cleaning up temporary files in {temp_dir}"
        )
