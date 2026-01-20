import pytest

from indel_scanner.homopolymer_classifier import HomopolymerClassifier
from indel_scanner.indel import (
    IndelTsvFormatter,
    Insertion,
    Deletion,
    INDEL_TYPE,
    FilterFlag,
)


@pytest.mark.parametrize(
    "scores, expected",
    [
        ([10, 20, 30], "10,20,30"),  # Test with a typical list of integers
        ([], ""),  # Test with an empty list
        (None, "NA"),  # Test the special 'None' case
    ],
)
def test_format_quality_scores(scores, expected):
    """Tests the _format_quality_scores helper function with various inputs."""
    assert IndelTsvFormatter._format_quality_scores(scores) == expected


@pytest.mark.parametrize(
    "seq, interval, min_len, expected",
    [
        # Cases that should return True (overlap exists)
        ("AAACCCGGG", (0, 2), 3, True),  # Run is the interval: [AAA]CCCGGG
        ("AAACCCGGG", (1, 3), 3, True),  # Interval overlaps run: A[AAC]CCGGG
        ("AAACCCGGG", (2, 4), 3, True),  # Interval overlaps two runs: AA[ACC]CGGG
        ("AAACCCGGG", (3, 5), 3, True),  # Interval is the run: AAA[CCC]GGG
        ("AAACCCGGG", (0, 8), 3, True),  # Interval contains runs: [AAACCCGGG]
        ("TTT", (0, 0), 3, True),  # Overlap at the start: [T]TT
        ("GGGGG", (1, 3), 3, True),  # Interval is inside run: G[GGG]G
        ("AAACCCGGG", (4, 6), 3, True),  # Extension on prefix: AAAC[CCG]GG
        (
            "AAACCCGGG",
            (0, 1),
            3,
            True,
        ),  # indel starts spanning homopolymer into suffix: [AA]ACCCGGG
        # Cases that should return False (no overlap or no valid run)
        ("TATCCGAGT", (4, 6), 3, False),  # no homopolymers: TAT[CCG]AGT
        ("AAACCCGGG", (0, 1), 4, False),  # Run 'AAA' is not long enough for min_len=4
        ("AAGGCC", (0, 5), 3, False),  # No runs of min_len: [AAGGCC]
        ("TTTAA", (3, 4), 3, False),  # Run is before the interval: TTT[AA]
        ("", (0, 1), 3, False),  # Empty sequence
        # Case where run length is less than min_len
        (
            "AAACCCGGG",
            (0, 2),
            4,
            False,
        ),  # Run exists, but is shorter than min_len=4: [AAA]CCCGGG
    ],
)
def test_run_overlaps_interval(seq, interval, min_len, expected):
    """Tests the static method Indel.run_overlaps_interval with various scenarios."""
    assert (
        HomopolymerClassifier.run_overlaps_interval(seq, interval, min_len) == expected
    )


# ==============================================================================
# Part 3: Test Objects with Fixtures
# ==============================================================================


@pytest.fixture
def insertion_instance() -> Insertion:
    """Provides a standard Insertion object for use in multiple tests."""
    return Insertion(
        contig="chr1",
        ref_position=100,
        length=3,
        prefix_context="GATT",
        indel_content="ACA",
        suffix_context="TACC",
        read_name="read_1",
        type=INDEL_TYPE.INSERTION,
        in_STR=False,
        filter_mask=FilterFlag.NONE,
        map_quality=60,
        indel_quality=[30, 31, 32],
        prefix_quality=[25, 26, 27, 28],
        suffix_quality=[35, 36, 37, 38],
    )


@pytest.fixture
def deletion_instance() -> Deletion:
    """Provides a standard Deletion object for use in multiple tests."""
    return Deletion(
        contig="chr2",
        ref_position=200,
        length=2,
        prefix_context="AAT",
        indel_content="GG",
        suffix_context="CCA",
        read_name="read_2",
        type=INDEL_TYPE.DELETION,
        in_STR=True,
        filter_mask=FilterFlag.NONE,
        map_quality=50,
        prefix_quality=[20, 21, 22],
        suffix_quality=[28, 29, 30],
    )


# ==============================================================================
# Part 4: Testing the Instance Methods
# ==============================================================================


def test_sequence_context(insertion_instance):
    """Tests the sequencecontext method."""
    # Arrange: The fixture provides the object
    # Act: Call the method
    result = insertion_instance.sequencecontext_brackets()
    # Assert: Check the output
    assert result == "GATT[ACA]TACC"


def test_to_scanner_tsv_row(insertion_instance, deletion_instance):
    """Tests the to_scanner_tsv_row method for both Indel types."""
    # Test Insertion
    ins_row = IndelTsvFormatter.format_scanner_row(insertion_instance)
    expected_ins = ["chr1", "100", "ins", "3", "GATT[ACA]TACC", "read_1", "False", ""]
    assert ins_row == expected_ins

    # Test Deletion
    del_row = IndelTsvFormatter.format_scanner_row(deletion_instance)
    expected_del = ["chr2", "200", "del", "2", "AAT[GG]CCA", "read_2", "True", ""]
    assert del_row == expected_del


def test_insertion_to_processor_tsv_row(insertion_instance):
    """Tests the processor-specific TSV row for an Insertion."""
    proc_row = IndelTsvFormatter.format_processor_row(insertion_instance)
    expected_proc = [
        "chr1",
        "100",
        "ins",
        "3",
        "GATT[ACA]TACC",
        "read_1",
        "False",
        "",  # Scanner part
        "25,26,27,28",
        "30,31,32",
        "35,36,37,38",
        "60",  # Processor part
    ]
    assert proc_row == expected_proc


def test_deletion_to_processor_tsv_row(deletion_instance):
    """Tests the processor-specific TSV row for a Deletion."""
    proc_row = IndelTsvFormatter.format_processor_row(deletion_instance)
    expected_proc = [
        "chr2",
        "200",
        "del",
        "2",
        "AAT[GG]CCA",
        "read_2",
        "True",
        "",  # Scanner part
        "20,21,22",
        "NA",
        "28,29,30",
        "50",  # Processor part
    ]
    assert proc_row == expected_proc


# ==============================================================================
# Part 5: Testing the Filtering Logic
# ==============================================================================


@pytest.mark.parametrize(
    "prefix, indel, suffix, expected_filter",
    [
        # Should be filtered
        ("GATT", "AAA", "TACC", True),  # Indel itself is a homopolymer
        ("GAAA", "TT", "TACC", True),  # Indel extends a homopolymer run
        ("GATT", "CC", "CCGT", True),  # Indel is adjacent to suffix homopolymer
        ("GTTT", "AC", "TACC", True),  # Indel is adjacent to prefix homopolymer
        (
            "GATT",
            "AC",
            "CCT",
            True,
        ),  # Indel starts a spanning homopolymer into the suffix
        ("GAA", "AT", "TACC", True),  # Indel ends a spanning homopolymer
        # Should NOT be filtered
        (
            "GATT",
            "AC",
            "TACC",
            False,
        ),  # No homopolymers, TODO: filter more stricly when indel is in homopolymer
        ("GATT", "AA", "TACC", False),  # Indel is too short to be a homopolymer
        ("GTT", "AC", "TACC", False),  # Adjacent prefix is too short
    ],
)
def test_should_filter_indel(prefix, indel, suffix, expected_filter):
    """Tests the _should_filter_indel logic with various contexts."""
    # We can create a minimal object for this test
    test_indel = Insertion(
        contig="chr1",
        ref_position=1,
        length=len(indel),
        prefix_context=prefix,
        indel_content=indel,
        suffix_context=suffix,
        read_name="r1",
        type=INDEL_TYPE.INSERTION,
        in_STR=False,
        filter_mask=FilterFlag.NONE,
        map_quality=50,
    )
    assert HomopolymerClassifier(test_indel).should_filter_indel() == expected_filter
