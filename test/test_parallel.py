import pytest
from queue import Empty
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
        config.sampling_target_aligned_bases = 10000000
        config.sampling_top_n_contigs = 10
        config.sampling_random_seed = 1
        config.callable_lengths = [1, 2, 3]
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
        mocker.patch('indel_scanner.parallel._collect_indel_lengths', return_value=[1])
        m_logger = mocker.patch('indel_scanner.parallel.logger')
        queue_mock = mocker.patch('indel_scanner.parallel.Queue')
        queue_instances = []
        for _ in range(4):
            instance = MagicMock()
            instance.get_nowait.side_effect = Empty
            queue_instances.append(instance)
        queue_mock.side_effect = queue_instances

        # Mock external libraries (pysam, multiprocessing, rich)
        m_pysam = mocker.patch('indel_scanner.parallel.pysam.AlignmentFile')
        m_pool = mocker.patch('indel_scanner.parallel.Pool')
        m_progress = mocker.patch('indel_scanner.parallel.Progress')
        mocker.patch('indel_scanner.parallel.Live')
        mocker.patch('indel_scanner.parallel.Thread')
        event_mock = mocker.patch('indel_scanner.parallel.Event')
        event_instance = MagicMock()
        event_instance.is_set.return_value = False
        event_instance.wait.return_value = None
        event_mock.return_value = event_instance

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
                "chr1",
                [],
                10,
                {
                    "reads_processed": 1,
                    "candidates": 2,
                    "passed": 1,
                    "type_counts": {"snp": 1},
                },
            ),
            (
                "chr2",
                [],
                12,
                {
                    "reads_processed": 2,
                    "candidates": 3,
                    "passed": 2,
                    "type_counts": {"del_indel_1bp": 2},
                },
            ),
        ]
        pre_callable_pool = MagicMock()
        pre_async = MagicMock()
        pre_async.get.return_value = [
            (
                "chr1",
                {
                    "total_aligned_bases": 50,
                    "sampling_bases_total": 5,
                    "sampling_passable_by_type": {"snp": 2, "ins_len_1bp": 1},
                    "tract_counts_by_motif": {},
                    "sampled_contig": True,
                },
            )
        ]
        pre_callable_pool.map_async.return_value = pre_async
        mock_pool_instance.__enter__.return_value = mock_pool_instance
        m_pool.side_effect = [pre_callable_pool, mock_pool_instance]

        # 3. Configure the Progress bar mock to avoid side effects
        mock_progress_instance = MagicMock()
        m_progress.return_value = mock_progress_instance

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
        assert m_pool.call_count == 2
        mock_pool_instance.imap_unordered.assert_called_once()

        # Verify progress bar was used
        mock_progress_instance.add_task.assert_called_once_with(
            "[green]Scanning contigs...", total=2
        )

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
        mocker.patch('indel_scanner.parallel._collect_indel_lengths', return_value=[1])
        mocker.patch('indel_scanner.parallel.logger')
        queue_mock = mocker.patch('indel_scanner.parallel.Queue')
        queue_instances = []
        for _ in range(4):
            instance = MagicMock()
            instance.get_nowait.side_effect = Empty
            queue_instances.append(instance)
        queue_mock.side_effect = queue_instances
        m_pysam = mocker.patch('indel_scanner.parallel.pysam.AlignmentFile')

        # --- FIX: Patch 'Pool' where it is looked up, not where it is defined ---
        m_pool = mocker.patch('indel_scanner.parallel.Pool')

        mocker.patch('indel_scanner.parallel.Progress')
        mocker.patch('indel_scanner.parallel.Live')
        mocker.patch('indel_scanner.parallel.Thread')
        event_mock = mocker.patch('indel_scanner.parallel.Event')
        event_instance = MagicMock()
        event_instance.is_set.return_value = False
        event_instance.wait.return_value = None
        event_mock.return_value = event_instance
        m_cpu_count = mocker.patch('indel_scanner.parallel.cpu_count', return_value=8)

        # Basic setup for pysam to allow the function to run
        mock_samfile = MagicMock()
        mock_samfile.header.references = ['chr1']
        mock_samfile.header.lengths = [1000]
        m_pysam.return_value.__enter__.return_value = mock_samfile

        # Basic setup for Pool context manager
        mock_pool_instance = MagicMock()
        pre_callable_pool = MagicMock()
        pre_async = MagicMock()
        pre_async.get.return_value = []
        pre_callable_pool.map_async.return_value = pre_async
        mock_pool_instance.__enter__.return_value = MagicMock()
        m_pool.side_effect = [pre_callable_pool, mock_pool_instance]

        # --- Act ---
        parallel_pipeline(mock_scanner_config)

        # --- Assert ---
        # The key assertion: verify the Pool was created with the value from cpu_count
        m_cpu_count.assert_called_once()
        assert m_pool.call_count == 2
