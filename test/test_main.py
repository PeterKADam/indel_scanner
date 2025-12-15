
import pytest
from unittest.mock import MagicMock
from indel_scanner import main
from indel_scanner.configurator import ScannerConfig, ProcessorConfig

@pytest.fixture
def mock_main_dependencies(mocker):
    """A pytest fixture to mock all dependencies of the main() function."""
    mocks = {
        'setup_logging': mocker.patch('indel_scanner.main.setup_logging'),
        'Config_load': mocker.patch('indel_scanner.main.Config.load'),
        'run_scan': mocker.patch('indel_scanner.main.run_scan'),
        'process': mocker.patch('indel_scanner.main.process'),
        'sys_exit': mocker.patch('sys.exit'),
        'logger': mocker.patch('indel_scanner.main.logger')
    }
    return mocks

def test_main_with_scanner_config_and_process_flag(mock_main_dependencies):
    """
    GIVEN a ScannerConfig with args.process = True
    WHEN main() is called
    THEN run_scan() and process() should be called and the program should not exit with an error.
    """
    mock_args = MagicMock()
    mock_args.process = True
    mock_config = MagicMock(spec=ScannerConfig)
    mock_config.args = mock_args
    mock_main_dependencies['Config_load'].return_value = mock_config

    main.main()

    mock_main_dependencies['run_scan'].assert_called_once_with(mock_config)
    mock_main_dependencies['process'].assert_called_once_with(mock_config)
    mock_main_dependencies['sys_exit'].assert_not_called()

def test_main_with_scanner_config_no_process_flag(mock_main_dependencies):
    """
    GIVEN a ScannerConfig with args.process = False
    WHEN main() is called
    THEN only run_scan() should be called and the program should not exit with an error.
    """
    mock_args = MagicMock()
    mock_args.process = False
    mock_config = MagicMock(spec=ScannerConfig)
    mock_config.args = mock_args
    mock_main_dependencies['Config_load'].return_value = mock_config

    main.main()

    mock_main_dependencies['run_scan'].assert_called_once_with(mock_config)
    mock_main_dependencies['process'].assert_not_called()
    mock_main_dependencies['sys_exit'].assert_not_called()

def test_main_with_processor_config(mock_main_dependencies):
    """
    GIVEN a ProcessorConfig
    WHEN main() is called
    THEN only process() should be called and the program should not exit with an error.
    """
    mock_config = MagicMock(spec=ProcessorConfig)
    mock_main_dependencies['Config_load'].return_value = mock_config

    main.main()

    mock_main_dependencies['run_scan'].assert_not_called()
    mock_main_dependencies['process'].assert_called_once_with(mock_config)

    # According to the new logic, this path is also a success case.
    mock_main_dependencies['sys_exit'].assert_not_called()

