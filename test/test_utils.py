import pytest
from unittest.mock import Mock

from indel_scanner.utils import (
    as_cigar,
    Cigar,
    build_refpos_to_read_index,
    flank_qualities_by_ref,
)

### positive
def test_as_cigar_with_valid_op():
    """
    Tests if as_cigar correctly converts a valid integer to the Cigar enum.
    """
    # Arrange: Define the input and expected output
    op_integer = 0  # In pysam, CMATECH is 0
    expected_cigar_member = Cigar.OP_M

    # Act: Call the function we are testing
    result = as_cigar(op_integer)

    # Assert: Check if the result is what we expect
    assert result == expected_cigar_member


@pytest.mark.parametrize("op_int, expected_cigar", [
    (0, Cigar.OP_M),      # CMATCH
    (1, Cigar.OP_I),      # CINS
    (2, Cigar.OP_D),      # CDEL
    (4, Cigar.OP_S),      # CSOFT_CLIP
    (7, Cigar.OP_EQ),     # CEQUAL
    (8, Cigar.OP_X),      # CDIFF
])
def test_as_cigar_with_multiple_valid_ops(op_int, expected_cigar):
    assert as_cigar(op_int) == expected_cigar

### negative
def test_as_cigar_with_invalid_op():
    # Arrange
    invalid_op_integer = 99

    # Act
    result = as_cigar(invalid_op_integer)

    # Assert
    assert result is None

###logging?
def test_as_cigar_logs_error_on_invalid_op(caplog):
    """
    Tests if an error is logged when an invalid CIGAR operation is provided.
    `caplog` is a special pytest "fixture" that captures logging output.
    """
    # Arrange
    invalid_op_integer = 99

    # Act
    as_cigar(invalid_op_integer)

    # Assert
    assert "Unrecognized CIGAR operation: 99" in caplog.text
    # You can also be more specific about the log level
    assert len(caplog.records) == 1
    assert caplog.records[0].levelname == "ERROR"


def test_build_refpos_to_read_index_uses_aligned_pairs_mapping():
    read = Mock()
    read.get_aligned_pairs.return_value = [
        (0, 100),
        (1, 101),
        (2, 102),
        (3, None),  # insertion base
        (4, 103),
        (5, 104),
    ]

    mapping = build_refpos_to_read_index(read)

    assert mapping == {100: 0, 101: 1, 102: 2, 103: 4, 104: 5}
    read.get_aligned_pairs.assert_called_once_with(matches_only=True)


def test_flank_qualities_by_ref_handles_insertion_suffix_from_reference_anchor():
    read = Mock()
    read.query_qualities = [10, 11, 12, 13, 14, 15]
    # Simulate 3M1I2M around ref 100..104.
    read.get_aligned_pairs.return_value = [
        (0, 100),
        (1, 101),
        (2, 102),
        (3, None),
        (4, 103),
        (5, 104),
    ]

    prefix, suffix = flank_qualities_by_ref(
        read=read,
        ref_pos=103,
        length=1,
        flank_len=2,
        is_insertion=True,
    )

    # Prefix looks at ref 101,102; suffix looks at 103,104.
    assert prefix == [11, 12]
    assert suffix == [14, 15]

