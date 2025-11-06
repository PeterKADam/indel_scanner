import logging
logger = logging.getLogger(__name__)
def setup_logging():
	logging.basicConfig(level=logging.INFO,
						format='%(asctime)s - %(levelname)s - %(message)s',
						filename='indel_scanner.log',
						filemode='w',
						encoding='utf-8')
	
	logger.setLevel(logging.DEBUG)
	logfile = logging.FileHandler('indel_scanner.log',mode='w',encoding='utf-8',)
	logfile.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
	logfile.setLevel(logging.DEBUG)

	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	console.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))

	logger.addHandler(logfile)
	logger.addHandler(console)