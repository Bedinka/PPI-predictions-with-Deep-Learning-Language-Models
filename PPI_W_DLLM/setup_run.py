import os
import sys
import datetime
import logging

class LoggerWriter:
    def __init__(self, level):
        self.level = level

    def write(self, message):
        if message.strip():
            self.level(message)

    def flush(self):
        pass

def setup_run_directory():
    # Create a timestamped directory for this run
    run_directory = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    if not os.path.exists(run_directory):
        os.makedirs(run_directory)
    
    # Set up logging to both a file and the console
    log_file = os.path.join(run_directory, f'{run_directory}.log')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ])
    
    # Redirect stdout and stderr to logging
    sys.stdout = LoggerWriter(logging.info)
    sys.stderr = LoggerWriter(logging.error)

    # Log the creation of the run directory
    logging.info("Run directory created at %s", run_directory)

    # Change the current working directory to the run directory
    os.chdir(run_directory)
    
    return run_directory
