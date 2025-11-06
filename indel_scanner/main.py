# indel_scanner/indel_scanner/main.py

import sys
import logging
import time

from .parallel import parallel_scan
from .log import setup_logging
from .processor import run_processor
from .configurator import Config, ScannerConfig, ProcessorConfig
logger = logging.getLogger(__name__)


	

def run_scan(config):

	start_time = time.time()
	parallel_scan(config)
	end_time = time.time()

	logger.info(f"SCANNER: Total time taken: {end_time - start_time:.2f} seconds.")
	
def process(config):

	start_time = time.time()
	run_processor(config)
	end_time = time.time()

	logger.info(f"PROCESSOR: Total time taken: {end_time - start_time:.2f} seconds.")

def main():
	setup_logging()

	config = Config.load()

	# --- Dispatch to the correct function based on the command ---
	if isinstance(config,ScannerConfig):
		run_scan(config)
		if config.args.process:
			process(config)

	elif isinstance(config,ProcessorConfig):
		process(config)

	sys.exit(1)

if __name__ == '__main__':
	main()
