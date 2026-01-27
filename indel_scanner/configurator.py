import argparse
import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

import yaml

logger = logging.getLogger(__name__)

DEFAULT_SCANNER: Dict[str, Any] = {
    "num_processes": 4,
    "min_indel_size": 1,
    "preload_contigs": True,
    "min_map_quality": 30,
    "min_base_quality": 20,
    "min_flank_quality": 93,
    "min_homopolymer_len": 3,
    "temp_dir_name": "temp_indel_parts",
    "output_filename": "indel_scanner_results.tsv",
}

DEFAULT_PROCESSOR: Dict[str, Any] = {
    "output_subdir": "processed",
    "passed_indels_filename": "final_passed_indels.tsv",
    "filters": [
        {"name": "is_in_str"},
        {"name": "is_homopolymer_context", "params": {"min_homopolymer_len": 3}},
        {"name": "similar_indels_in_other_reads"},
        {"name": "low_minimum_indel_quality", "params": {"min_quality": 93}},
        {
            "name": "low_singlebase_flanking_quality",
            "params": {"min_flank_quality": 93},
        },
    ],
}

DEFAULT_REPORTING: Dict[str, Any] = {
    "mutation_frequency_filename": "final_mutation_frequency.tsv",
    "per_type_frequency_filename": "per_type_mutation_frequency.tsv",
}

DEFAULT_IO: Dict[str, Any] = {
    "write_buffer_size": 1000,
    "read_batch_size": 10000,
}

DEFAULT_PIPELINE: Dict[str, Any] = {
    "in_memory": False,
    "max_in_memory_records": 500000,
    "passed_parts_dir_name": "temp_passed_indel_parts",
    "sampling_bases": 1000000,
    "sampling_strategy": "largest_contig",
    "imperfect_str": {
        "enabled": False,
        "expand_bp": 0,
        "window_bp": 20,
        "motif_min": 2,
        "motif_max": 6,
        "max_mismatches": 2,
    },
    "low_complexity": {
        "enabled": False,
        "window_bp": 32,
        "entropy_threshold": 1.2,
    },
    "repeat_run": {
        "enabled": False,
        "min_run_bp": 16,
        "max_window_bp": 80,
        "max_mismatches": 2,
    },
    "indel_bins": [
        {"label": "indel_1bp", "min": 1, "max": 1},
        {"label": "indel_2_3bp", "min": 2, "max": 3},
        {"label": "indel_4_10bp", "min": 4, "max": 10},
    ],
    "snp_label": "snp",
}


class PipelineConfig:
    def __init__(self, args: argparse.Namespace, config_data: dict) -> None:
        self.args = args
        self.yaml = config_data

        self.bamfile: Path = Path(self.args.bam)
        self.fastafile: Path = Path(self.args.fasta)
        self.sample_name = self.bamfile.stem
        self.base_output_path: Path = Path(self.args.output)
        run_stamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        self.run_id = run_stamp
        self.output_path: Path = self.base_output_path / self.sample_name / run_stamp
        self.log_dir: Path = self.output_path / "logs"
        self.str_directory: Path = Path(self.args.strdir)

        self.scanner = {**DEFAULT_SCANNER, **(config_data.get("scanner") or {})}
        self.processor = {**DEFAULT_PROCESSOR, **(config_data.get("processor") or {})}
        self.reporting = {**DEFAULT_REPORTING, **(config_data.get("reporting") or {})}
        self.io = {**DEFAULT_IO, **(config_data.get("io") or {})}
        self.pipeline = {**DEFAULT_PIPELINE, **(config_data.get("pipeline") or {})}
        self.pipeline["imperfect_str"] = {
            **DEFAULT_PIPELINE["imperfect_str"],
            **(self.pipeline.get("imperfect_str") or {}),
        }
        self.pipeline["low_complexity"] = {
            **DEFAULT_PIPELINE["low_complexity"],
            **(self.pipeline.get("low_complexity") or {}),
        }
        self.pipeline["repeat_run"] = {
            **DEFAULT_PIPELINE["repeat_run"],
            **(self.pipeline.get("repeat_run") or {}),
        }

        self.num_processes: int = self.scanner["num_processes"]
        self.min_indel_size: int = self.scanner["min_indel_size"]
        self.preload: bool = self.scanner["preload_contigs"]
        scanner_write_buffer = (config_data.get("scanner") or {}).get("write_buffer_size")
        if (config_data.get("io") or {}).get("write_buffer_size") is not None:
            self.write_buffer_size = self.io["write_buffer_size"]
        elif scanner_write_buffer is not None:
            self.write_buffer_size = scanner_write_buffer
        else:
            self.write_buffer_size = self.io["write_buffer_size"]

        self.read_batch_size: int = self.io["read_batch_size"]
        self.min_map_quality: int = self.scanner["min_map_quality"]
        self.min_base_quality: int = self.scanner["min_base_quality"]
        self.min_flank_quality: int = self.scanner["min_flank_quality"]
        self.min_homopolymer_len: int = self.scanner["min_homopolymer_len"]

        self.temp_dir = self.output_path / self.scanner["temp_dir_name"]
        self.scan_output_path = self.output_path / self.scanner["output_filename"]

        self.processor_output_dir = self.output_path / self.processor["output_subdir"]
        self.passed_indels_path = (
            self.processor_output_dir / self.processor["passed_indels_filename"]
        )
        self.processor_filters = self.processor["filters"]

        self.report_filename = self.reporting["mutation_frequency_filename"]
        self.per_type_report_filename = self.reporting["per_type_frequency_filename"]
        self.in_memory: bool = self.pipeline["in_memory"]
        self.max_in_memory_records: int = self.pipeline["max_in_memory_records"]
        self.passed_parts_dir = self.output_path / self.pipeline["passed_parts_dir_name"]
        self.sampling_bases: int = self.pipeline["sampling_bases"]
        self.sampling_strategy: str = self.pipeline["sampling_strategy"]
        self.indel_bins = self.pipeline["indel_bins"]
        self.snp_label: str = self.pipeline["snp_label"]
        self.imperfect_str = self.pipeline["imperfect_str"]
        self.low_complexity = self.pipeline["low_complexity"]
        self.repeat_run = self.pipeline["repeat_run"]
        self.sampling_contig: str = ""
        self.contigs = self._parse_contigs(self.args.contigs, self.scanner.get("contigs"))

        self._setup_output()

    def _parse_contigs(self, cli_value, yaml_value):
        raw = cli_value if cli_value is not None else yaml_value
        if raw is None:
            return None
        if isinstance(raw, list):
            return [str(c).strip() for c in raw if str(c).strip()]
        if isinstance(raw, str):
            parts = [p.strip() for p in raw.split(",")]
            return [p for p in parts if p]
        return None

    def _setup_output(self) -> None:
        try:
            self.output_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Output directory set up at {self.output_path}")
            if self.temp_dir.exists():
                logger.info(
                    f"Found existing temporary directory at {self.temp_dir}. Cleaning it before use."
                )
                shutil.rmtree(self.temp_dir, ignore_errors=True)
            self.temp_dir.mkdir(parents=True, exist_ok=True)
            self.processor_output_dir.mkdir(parents=True, exist_ok=True)
            self.log_dir.mkdir(parents=True, exist_ok=True)
            if self.passed_parts_dir.exists():
                shutil.rmtree(self.passed_parts_dir, ignore_errors=True)
            self.passed_parts_dir.mkdir(parents=True, exist_ok=True)
            logger.debug(
                f"Temporary directory for partial results is ready at {self.temp_dir}"
            )
        except Exception as e:
            logger.error(f"Error setting up output directories: {e}")
            sys.exit(1)

    def write_settings_file(self) -> None:
        settings_path = self.output_path / "settings_used.yaml"
        settings_payload = {
            "args": self._serialize_for_yaml(vars(self.args)),
            "input": {
                "bam": str(self.bamfile),
                "fasta": str(self.fastafile),
                "str_directory": str(self.str_directory),
            },
            "output": {
                "base_output_path": str(self.base_output_path),
                "sample_name": self.sample_name,
                "run_id": self.run_id,
                "output_path": str(self.output_path),
                "log_dir": str(self.log_dir),
            },
            "config": self._serialize_for_yaml(self.yaml),
            "resolved": self._serialize_for_yaml(
                {
                    "scanner": self.scanner,
                    "processor": self.processor,
                    "reporting": self.reporting,
                    "io": self.io,
                    "pipeline": self.pipeline,
                }
            ),
        }
        try:
            with open(settings_path, "w", encoding="utf-8") as handle:
                yaml.safe_dump(settings_payload, handle, sort_keys=False)
            logger.info("Wrote settings file to %s", settings_path)
        except Exception as exc:
            logger.error("Failed to write settings file at %s: %s", settings_path, exc)

    @staticmethod
    def _serialize_for_yaml(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {key: PipelineConfig._serialize_for_yaml(val) for key, val in value.items()}
        if isinstance(value, list):
            return [PipelineConfig._serialize_for_yaml(item) for item in value]
        if isinstance(value, tuple):
            return [PipelineConfig._serialize_for_yaml(item) for item in value]
        return value


class Config:
    @staticmethod
    def _parse_args() -> argparse.Namespace:
        parser = argparse.ArgumentParser(
            description="A toolkit for scanning and processing indels from BAM files."
        )
        parser.add_argument(
            "-b", "--bam", required=True, help="Input BAM file (must be indexed)."
        )
        parser.add_argument(
            "-f",
            "--fasta",
            required=True,
            help="Reference FASTA file (must be indexed).",
        )
        parser.add_argument(
            "-o", "--output", required=True, help="Path to the output directory."
        )
        parser.add_argument(
            "-s",
            "--strdir",
            required=True,
            help="Directory containing STR results for the sample.",
        )
        parser.add_argument(
            "-c",
            "--config",
            default="config.yaml",
            help="Path to the configuration file.",
        )
        parser.add_argument(
            "--contigs",
            default=None,
            help="Comma-separated list of contigs to process (default: all).",
        )
        return parser.parse_args()

    @staticmethod
    def load() -> PipelineConfig:
        args = Config._parse_args()
        config_path = Path(args.config)
        full_config = {}
        if config_path.exists():
            try:
                with open(config_path, "r") as f:
                    full_config = yaml.safe_load(f) or {}
            except Exception as e:
                logger.error(f"Error loading configuration file {config_path}: {e}")
                sys.exit(1)
        else:
            logger.warning(
                f"Configuration file not found at {config_path}. Using command-line args and defaults."
            )

        return PipelineConfig(args, full_config)
