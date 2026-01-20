
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from indel_scanner import main
from indel_scanner.configurator import PipelineConfig


@pytest.fixture
def mock_main_dependencies(mocker):
    """A pytest fixture to mock all dependencies of the main() function."""
    mocks = {
        "setup_logging": mocker.patch("indel_scanner.main.setup_logging"),
        "Config_load": mocker.patch("indel_scanner.main.Config.load"),
        "parallel_pipeline": mocker.patch("indel_scanner.main.parallel_pipeline"),
        "report": mocker.patch("indel_scanner.main.write_mutation_frequency_report"),
        "per_type_report": mocker.patch("indel_scanner.main.write_per_type_mutation_report"),
        "logger": mocker.patch("indel_scanner.main.logger"),
    }
    return mocks


def test_main_runs_full_pipeline(mock_main_dependencies):
    """
    GIVEN a valid pipeline config
    WHEN main() is called
    THEN scanning, processing, and reporting should run in sequence.
    """
    mock_config = MagicMock(spec=PipelineConfig)
    mock_config.num_processes = 4
    mock_config.passed_indels_path = Path("out/processed/final_passed_indels.tsv")
    mock_config.output_path = Path("out")
    mock_config.report_filename = "final_mutation_frequency.tsv"
    mock_config.per_type_report_filename = "per_type.tsv"
    mock_config.snp_label = "snp"
    mock_config.indel_bins = [
        {"label": "indel_1bp", "min": 1, "max": 1},
        {"label": "indel_2_3bp", "min": 2, "max": 3},
        {"label": "indel_4_10bp", "min": 4, "max": 10},
    ]

    mock_main_dependencies["Config_load"].return_value = mock_config
    mock_main_dependencies["parallel_pipeline"].return_value = (
        Path("out/processed/final_passed_indels.tsv"),
        42,
        {
            "reads": 1,
            "candidates": 2,
            "passed": 1,
            "total_aligned_bases": 100,
            "sampling_bases_total": 10,
            "sampling_passable_by_type": {"snp": 2, "ins_indel_1bp": 1},
            "type_counts": {"snp": 1},
        },
    )

    main.main()

    mock_main_dependencies["parallel_pipeline"].assert_called_once_with(mock_config)
    mock_main_dependencies["report"].assert_called_once_with(
        passed_indels_path=Path("out/processed/final_passed_indels.tsv"),
        total_interrogated_bases=42,
        output_dir=mock_config.output_path,
        report_filename=mock_config.report_filename,
    )
    mock_main_dependencies["per_type_report"].assert_called_once_with(
        type_counts={
            "snp": 1,
            "ins_indel_1bp": 0,
            "del_indel_1bp": 0,
            "ins_indel_2_3bp": 0,
            "del_indel_2_3bp": 0,
            "ins_indel_4_10bp": 0,
            "del_indel_4_10bp": 0,
        },
        callable_bases_by_type={
            "snp": 20.0,
            "ins_indel_1bp": 10.0,
            "del_indel_1bp": 0.0,
            "ins_indel_2_3bp": 0.0,
            "del_indel_2_3bp": 0.0,
            "ins_indel_4_10bp": 0.0,
            "del_indel_4_10bp": 0.0,
        },
        output_dir=mock_config.output_path,
        report_filename=mock_config.per_type_report_filename,
    )

