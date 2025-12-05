
import pytest
from unittest.mock import Mock, mock_open, call
from pathlib import Path

# Adjust this import path to match your project structure.
# This assumes your functions are in a file at 'your_app/io_utils.py'.
from indel_scanner.IO import _write_records_to_tsv, cleanup_temp_dir


@pytest.fixture
def mock_records():
    """
    Provides a list of mock record objects for testing.
    This avoids a dependency on the real Insertion/Deletion classes,
    keeping the test focused only on the I/O function's logic.
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
        # Arrange: Mock the dependencies (open and logger).
        m_open = mocker.patch("app.io_utils.open", mock_open())
        m_logger = mocker.patch("app.io_utils.logger")

        output_path = Path("/fake/output.tsv")
        header = ["Contig", "Position"]

        # Act: Call the function under test.
        count = _write_records_to_tsv(
            output_path=output_path,
            records=mock_records,
            row_converter=mock_row_converter,
            header=header,
            sort_key=None  # Explicitly test the non-sorting path
        )

        # Assert: Verify the results.
        assert count == 3
        m_open.assert_called_once_with(output_path, "w", newline="")

        handle = m_open()
        # Check that the header was written
        handle.writerow.assert_called_once_with(header)

        # Check that the data rows were passed to writerows
        # Note: We check the call to writerows, which is what csv.writer does internally.
        written_rows = handle.writerows.call_args[0][0]
        expected_rows = [
            ['chr1', 200],
            ['chr1', 100],
            ['chr2', 50],
        ]
        assert written_rows == expected_rows

        # Verify the success log message
        m_logger.info.assert_called_with(f"Successfully wrote 3 records to {output_path}.")

    def test_write_records_with_sorting(self, mocker, mock_records):
        """
        GIVEN: A list of records and a sort key.
        WHEN: _write_records_to_tsv is called.
        THEN: It should log the sorting action and write all records in the correct sorted order.
        """
        # Arrange
        m_open = mocker.patch("app.io_utils.open", mock_open())
        m_logger = mocker.patch("app.io_utils.logger")

        output_path = Path("/fake/output_sorted.tsv")
        sort_func = lambda r: r.pos  # Sort by position

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

        # Get the arguments passed to the CSV writer's writerows method
        written_rows = m_open().writerows.call_args[0][0]

        # Check that the data rows were written in the *sorted* order (pos 50, 100, 200)
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
        m_open = mocker.patch("app.io_utils.open", mock_open())
        m_logger = mocker.patch("app.io_utils.logger")

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

        handle = m_open()
        handle.writerow.assert_called_once_with(header)
        handle.writerows.assert_not_called()  # No data rows should be written

        # Check the specific log message for no records
        m_logger.info.assert_called_with(
            f"No records to write for {output_path.name}, an empty file was created."
        )

    def test_io_error_handling(self, mocker, mock_records):
        """
        GIVEN: The open() call will raise an IOError.
        WHEN: _write_records_to_tsv is called.
        THEN: It should log the error and return 0.
        """
        # Arrange: Configure the mock_open to raise an error when called
        mocker.patch("app.io_utils.open", side_effect=IOError("Permission denied"))
        m_logger = mocker.patch("app.io_utils.logger")

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
        mock_rmtree = mocker.patch("app.io_utils.shutil.rmtree")
        mock_logger = mocker.patch("app.io_utils.logger")
        temp_dir = Path("/fake/temp/dir")

        # Act
        cleanup_temp_dir(temp_dir)

        # Assert
        mock_rmtree.assert_called_once_with(temp_dir)
        mock_logger.debug.assert_called_once_with(
            f"Cleaning up temporary files in {temp_dir}"
        )

