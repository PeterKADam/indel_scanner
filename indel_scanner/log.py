import datetime
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def setup_logging(log_dir: Path | None = None) -> None:
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    for handler in list(root_logger.handlers):
        root_logger.removeHandler(handler)

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
    root_logger.addHandler(console)

    if log_dir is not None:
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%dT%H.%M.%S")
        info_log = logging.FileHandler(
            log_dir / f"{timestamp}-indel_scanner.log",
            mode="w",
            encoding="utf-8",
        )
        info_log.setLevel(logging.INFO)
        info_log.setFormatter(formatter)

        debug_log = logging.FileHandler(
            log_dir / "indel_scanner.log",
            mode="w",
            encoding="utf-8",
        )
        debug_log.setLevel(logging.DEBUG)
        debug_log.setFormatter(formatter)

        root_logger.addHandler(info_log)
        root_logger.addHandler(debug_log)

    logger.setLevel(logging.DEBUG)
