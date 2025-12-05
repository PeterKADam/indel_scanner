import pytest
import yaml
from unittest.mock import MagicMock, call
from pathlib import Path

# --- FIX: Ensure this import path matches your file structure ---
# Since your file is configurator.py, we import from there.
from indel_scanner.configurator import Config, ScannerConfig, ProcessorConfig


class TestConfigFactory:
    """Tests the main Config.load() factory method."""

    def test_load_returns_scanner_config_for_scan_command(self, mocker):
        """
        GIVEN the command-line command is 'scan'
        WHEN Config.load() is called
        THEN it should instantiate and return a ScannerConfig object.
        """
        # --- Arrange ---
        mock_args = MagicMock(
            command="scan",
            config="path/to/config.yaml",
            bam="in.bam",
            fasta="in.fasta",
            output="out_dir",
            strdir="str_dir"
        )
        mocker.patch.object(Config, '_parse_args', return_value=mock_args)

        mock_yaml_data = {"scan": {"num_processes": 8, "min_indel_size": 2, "strdir": "yaml/str/dir"}}
        mocker.patch('pathlib.Path.exists', return_value=True)
        mocker.patch('builtins.open', mocker.mock_open(read_data=yaml.dump(mock_yaml_data)))

        # --- FIX: Correct the patch path to use the 'configurator' module name ---
        mock_scanner_init = mocker.patch('indel_scanner.configurator.ScannerConfig.__init__', return_value=None)

        # --- Act ---
        Config.load()

        # --- Assert ---
        mock_scanner_init.assert_called_once_with(mock_args,
                                                  {"num_processes": 8, "min_indel_size": 2, "strdir": "yaml/str/dir"})

    def test_load_returns_processor_config_for_process_command(self, mocker):
        """
        GIVEN the command-line command is 'process'
        WHEN Config.load() is called
        THEN it should instantiate and return a ProcessorConfig object.
        """
        # --- Arrange ---
        mock_args = MagicMock(
            command="process",
            config="path/to/config.yaml",
            input="in.tsv",
            bam="in.bam",
            fasta="in.fasta",
            output="out_dir"
        )
        mocker.patch.object(Config, '_parse_args', return_value=mock_args)

        mock_yaml_data = {"process": {"output_format": "json"}}
        mocker.patch('pathlib.Path.exists', return_value=True)
        mocker.patch('builtins.open', mocker.mock_open(read_data=yaml.dump(mock_yaml_data)))

        # --- FIX: Correct the patch path to use the 'configurator' module name ---
        mock_processor_init = mocker.patch('indel_scanner.configurator.ProcessorConfig.__init__', return_value=None)

        # --- Act ---
        Config.load()

        # --- Assert ---
        mock_processor_init.assert_called_once_with(mock_args, {"output_format": "json"})

    def test_load_handles_missing_config_file_gracefully(self, mocker):
        """
        GIVEN the config file path does not exist
        WHEN Config.load() is called
        THEN it should proceed with an empty config dict and not raise an error.
        """
        # --- Arrange ---
        mock_args = MagicMock(command="scan", config="nonexistent.yaml", bam="in.bam", fasta="in.fasta", output="out",
                              strdir="str")
        mocker.patch.object(Config, '_parse_args', return_value=mock_args)

        mocker.patch('pathlib.Path.exists', return_value=False)

        # --- FIX: Correct the patch path to use the 'configurator' module name ---
        mock_logger_warning = mocker.patch('indel_scanner.configurator.logger.warning')
        mock_scanner_init = mocker.patch('indel_scanner.configurator.ScannerConfig.__init__', return_value=None)

        # --- Act ---
        Config.load()

        # --- Assert ---
        mock_logger_warning.assert_called_with(
            "Configuration file not found at nonexistent.yaml. Using command-line args and defaults."
        )
        mock_scanner_init.assert_called_once_with(mock_args, {})


class TestScannerConfig:
    """Tests the specific logic of the ScannerConfig class, especially file I/O."""

    def test_setup_output_creates_directories_correctly(self, tmp_path):
        """
        GIVEN a path to a non-existent output directory
        WHEN ScannerConfig is initialized
        THEN it should create both the main output directory and the temp directory.
        """
        # --- Arrange ---
        output_dir = tmp_path / "scanner_output"

        args = MagicMock(command="scan", output=str(output_dir), bam="in.bam", fasta="in.fasta")

        # --- FIX: Add the missing 'strdir' key to the config data ---
        config_data = {"min_indel_size": 5, "strdir": "path/to/str"}

        # --- Act ---
        ScannerConfig(args, config_data)

        # --- Assert ---
        temp_dir = output_dir / "temp_indel_parts"
        assert output_dir.exists() and output_dir.is_dir()
        assert temp_dir.exists() and temp_dir.is_dir()

    def test_setup_output_cleans_existing_temp_directory(self, tmp_path):
        """
        GIVEN an existing temp directory with a file inside it
        WHEN ScannerConfig is initialized
        THEN it should delete the old temp directory and create a new, empty one.
        """
        # --- Arrange ---
        output_dir = tmp_path / "scanner_output"
        temp_dir = output_dir / "temp_indel_parts"

        temp_dir.mkdir(parents=True, exist_ok=True)
        dummy_file = temp_dir / "stale_data.txt"
        dummy_file.write_text("old data")
        assert dummy_file.exists()

        args = MagicMock(command="scan", output=str(output_dir), bam="in.bam", fasta="in.fasta")

        # --- FIX: Add the missing 'strdir' key to the config data ---
        config_data = {"strdir": "path/to/str"}

        # --- Act ---
        ScannerConfig(args, config_data)

        # --- Assert ---
        assert output_dir.exists()
        assert temp_dir.exists()
        assert not dummy_file.exists()


class TestProcessorConfig:
    """Tests the specific logic of the ProcessorConfig class."""

    def test_setup_output_fails_if_input_file_is_missing(self, tmp_path, mocker):
        """
        GIVEN a path to a non-existent input file
        WHEN ProcessorConfig is initialized
        THEN it should log an error and call sys.exit(1).
        """
        # --- Arrange ---
        input_file = tmp_path / "nonexistent_input.tsv"
        output_dir = tmp_path / "processor_output"
        assert not input_file.exists()

        mock_sys_exit = mocker.patch('sys.exit')

        # --- FIX: Correct the patch path to use the 'configurator' module name ---
        mock_logger_error = mocker.patch('indel_scanner.configurator.logger.error')

        args = MagicMock(
            command="process",
            input=str(input_file),
            output=str(output_dir),
            bam="in.bam",
            fasta="in.fasta"
        )
        config_data = {}

        # --- Act ---
        ProcessorConfig(args, config_data)

        # --- Assert ---
        mock_logger_error.assert_has_calls([
            call(f"Input file not found at: {input_file}"),
            call("The 'process' command requires a valid input file generated by the 'scan' command.")
        ])
        mock_sys_exit.assert_called_once_with(1)

    def test_setup_output_succeeds_if_input_file_exists(self, tmp_path):
        """
        GIVEN a valid, existing input file
        WHEN ProcessorConfig is initialized
        THEN it should create the output directory and not exit.
        """
        # --- Arrange ---
        input_file = tmp_path / "real_input.tsv"
        input_file.write_text("header\ndata")
        output_dir = tmp_path / "processor_output"
        assert not output_dir.exists()

        args = MagicMock(
            command="process",
            input=str(input_file),
            output=str(output_dir),
            bam="in.bam",
            fasta="in.fasta"
        )
        config_data = {}

        # --- Act ---
        ProcessorConfig(args, config_data)

        # --- Assert ---
        assert output_dir.exists() and output_dir.is_dir()

