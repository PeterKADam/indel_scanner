from unittest.mock import MagicMock
from datetime import datetime

import pytest
import yaml

from indel_scanner.configurator import Config, PipelineConfig


class TestConfigFactory:
    """Tests the main Config.load() factory method."""

    def test_load_returns_pipeline_config(self, mocker):
        """
        GIVEN a valid YAML file
        WHEN Config.load() is called
        THEN it should instantiate and return a PipelineConfig object.
        """
        mock_args = MagicMock(
            config="path/to/config.yaml",
            bam="in.bam",
            fasta="in.fasta",
            output="out_dir",
            strdir="str_dir",
        )
        mocker.patch.object(Config, "_parse_args", return_value=mock_args)

        mock_yaml_data = {"scanner": {"num_processes": 8}}
        mocker.patch("pathlib.Path.exists", return_value=True)
        mocker.patch("builtins.open", mocker.mock_open(read_data=yaml.dump(mock_yaml_data)))

        cfg = Config.load()

        assert isinstance(cfg, PipelineConfig)
        assert cfg.num_processes == 8

    def test_load_handles_missing_config_file_gracefully(self, mocker):
        """
        GIVEN the config file path does not exist
        WHEN Config.load() is called
        THEN it should proceed with an empty config dict and not raise an error.
        """
        mock_args = MagicMock(
            config="nonexistent.yaml",
            bam="in.bam",
            fasta="in.fasta",
            output="out",
            strdir="str",
        )
        mocker.patch.object(Config, "_parse_args", return_value=mock_args)

        mocker.patch("pathlib.Path.exists", return_value=False)
        mock_logger_warning = mocker.patch("indel_scanner.configurator.logger.warning")

        cfg = Config.load()

        mock_logger_warning.assert_called_with(
            "Configuration file not found at nonexistent.yaml. Using command-line args and defaults."
        )
        assert isinstance(cfg, PipelineConfig)


class TestPipelineConfig:
    """Tests the specific logic of the PipelineConfig class."""

    def test_setup_output_creates_directories_correctly(self, tmp_path):
        """
        GIVEN a path to a non-existent output directory
        WHEN PipelineConfig is initialized
        THEN it should create the output and temp directories.
        """
        output_dir = tmp_path / "scanner_output"
        args = MagicMock(
            output=str(output_dir),
            bam="in.bam",
            fasta="in.fasta",
            strdir="path/to/str",
        )
        config_data = {}

        cfg = PipelineConfig(args, config_data)

        assert output_dir.exists() and output_dir.is_dir()
        assert cfg.output_path.exists() and cfg.output_path.is_dir()
        assert cfg.temp_dir.exists() and cfg.temp_dir.is_dir()

    def test_setup_output_cleans_existing_temp_directory(self, tmp_path, mocker):
        """
        GIVEN an existing temp directory with a file inside it
        WHEN PipelineConfig is initialized
        THEN it should delete the old temp directory and create a new, empty one.
        """
        output_dir = tmp_path / "scanner_output"
        mock_datetime = mocker.patch("indel_scanner.configurator.datetime")
        mock_datetime.now.return_value = datetime(2026, 1, 20, 0, 0, 0)
        run_dir = output_dir / "2026-01-20_000000_in"
        temp_dir = run_dir / "temp_indel_parts"
        temp_dir.mkdir(parents=True, exist_ok=True)
        dummy_file = temp_dir / "stale_data.txt"
        dummy_file.write_text("old data")
        assert dummy_file.exists()

        args = MagicMock(
            output=str(output_dir),
            bam="in.bam",
            fasta="in.fasta",
            strdir="path/to/str",
        )
        config_data = {}

        cfg = PipelineConfig(args, config_data)

        assert output_dir.exists()
        assert cfg.temp_dir.exists()
        assert not dummy_file.exists()

