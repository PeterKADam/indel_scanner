import pytest
from unittest.mock import MagicMock, call, ANY

# --- IMPORTANT: Adjust the import path to match your project structure ---
from indel_scanner.parallel import parallel_pipeline


class TestParallelPipeline:
    """
    Unit tests for the orchestration logic in parallel.py.
    """

    @pytest.fixture
    def mock_scanner_config(self):
        """
        Provides a mock configuration object for the scanner.
        This fixture is reused by all tests in this class.
        """
        config = MagicMock()
        config.bamfile = "path/to/mock.bam"
        config.temp_dir = "mock/temp/dir"
        config.passed_parts_dir = "mock/passed/parts"
        config.passed_indels_path = "mock/output/final.tsv"
        config.num_processes = 4  # Explicitly set a number of processes
        config.in_memory = False
        config.max_in_memory_records = 0
        config.fastafile = "path/to/mock.fasta"
        config.sampling_strategy = "largest_contig"
        return config

    def test_parallel_pipeline_happy_path(self, mocker, mock_scanner_config):
        """
        GIVEN a valid scanner configuration
        WHEN parallel_pipeline is called
        THEN it should correctly orchestrate the entire workflow:
             1. Initialize ContigScanner
             2. Open the BAM file to get contigs
             3. Set up and run a multiprocessing Pool
             4. Call the final aggregation function with the correct paths
        """
        # --- Arrange: Mock all external dependencies ---
        # Mock the components from your own application
        m_aggregate = mocker.patch('indel_scanner.parallel.aggregate_tsv_parts')
        m_logger = mocker.patch('indel_scanner.parallel.logger')

        # Mock external libraries (pysam, multiprocessing, rich)
        m_pysam = mocker.patch('indel_scanner.parallel.pysam.AlignmentFile')
        m_pool = mocker.patch('indel_scanner.parallel.Pool')
        m_progress = mocker.patch('indel_scanner.parallel.Progress')

        # --- Configure the behavior of the mocks ---
        # 1. Configure pysam mock to return a fake SAM file with two contigs
        mock_samfile = MagicMock()
        mock_samfile.header.references = ['chr1', 'chr2']
        mock_samfile.header.lengths = [1000, 2000]
        m_pysam.return_value.__enter__.return_value = mock_samfile

        # 2. Configure the multiprocessing Pool mock
        mock_pool_instance = MagicMock()
        # Simulate the pool returning results for each contig
        mock_pool_instance.imap_unordered.return_value = [
            (
                None,
                10,
                {
                    "reads_processed": 1,
                    "candidates": 2,
                    "passed": 1,
                    "total_aligned_bases": 50,
                    "sampling_bases_total": 5,
                    "sampling_passable_by_type": {"snp": 2, "ins_indel_1bp": 1},
                    "sampled_contig": True,
                    "type_counts": {"snp": 1},
                },
            ),
            (
                None,
                12,
                {
                    "reads_processed": 2,
                    "candidates": 3,
                    "passed": 2,
                    "total_aligned_bases": 75,
                    "sampling_bases_total": 0,
                    "sampling_passable_by_type": {},
                    "sampled_contig": False,
                    "type_counts": {"del_indel_1bp": 2},
                },
            ),
        ]
        m_pool.return_value.__enter__.return_value = mock_pool_instance

        # 3. Configure the Progress bar mock to avoid side effects
        mock_progress_instance = MagicMock()
        m_progress.return_value.__enter__.return_value = mock_progress_instance

        # --- Act: Call the function we are testing ---
        parallel_pipeline(mock_scanner_config)

        # --- Assert: Verify that the orchestration happened as expected ---
        # Verify initial setup
        m_logger.info.assert_has_calls([
            call(f"Using {mock_scanner_config.num_processes} parallel processes."),
            call("Starting streaming scan+process pipeline..."),
            call("Found 2 contigs.")
        ])
        m_pysam.assert_called_once_with(str(mock_scanner_config.bamfile), "rb")

        # Verify multiprocessing setup and execution
        m_pool.assert_called_once()
        mock_pool_instance.imap_unordered.assert_called_once()

        # Verify progress bar was used
        mock_progress_instance.add_task.assert_called_once_with("[green]Scanning contigs...", total=2)
        assert mock_progress_instance.update.call_count == 2

        # Verify final aggregation step
        m_aggregate.assert_called_once_with(
            mock_scanner_config.passed_parts_dir,
            mock_scanner_config.passed_indels_path,
            ANY,
        )

    def test_parallel_pipeline_uses_cpu_count_when_processes_not_specified(self, mocker, mock_scanner_config):
        """
        GIVEN a scanner configuration where num_processes is None
        WHEN parallel_pipeline is called
        THEN it should use cpu_count() to determine the number of processes for the Pool.
        """
        # --- Arrange ---
        # Modify the config for this specific test case
        mock_scanner_config.num_processes = None

        # Mock dependencies.
        mocker.patch('indel_scanner.parallel.aggregate_tsv_parts')
        mocker.patch('indel_scanner.parallel.logger')
        m_pysam = mocker.patch('indel_scanner.parallel.pysam.AlignmentFile')

        # --- FIX: Patch 'Pool' where it is looked up, not where it is defined ---
        m_pool = mocker.patch('indel_scanner.parallel.Pool')

        mocker.patch('indel_scanner.parallel.Progress')
        m_cpu_count = mocker.patch('indel_scanner.parallel.cpu_count', return_value=8)

        # Basic setup for pysam to allow the function to run
        mock_samfile = MagicMock()
        mock_samfile.header.references = ['chr1']
        mock_samfile.header.lengths = [1000]
        m_pysam.return_value.__enter__.return_value = mock_samfile

        # Basic setup for Pool context manager
        m_pool.return_value.__enter__.return_value = MagicMock()

        # --- Act ---
        parallel_pipeline(mock_scanner_config)

        # --- Assert ---
        # The key assertion: verify the Pool was created with the value from cpu_count
        m_cpu_count.assert_called_once()
        m_pool.assert_called_once()
