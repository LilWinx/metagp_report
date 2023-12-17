import os
import sys
import logging
import pandas as pd

def check_files(file):
    """
    Check input files if they exist and have contents
    """

    if os.path.isfile(file) is True and os.stat(file).st_size != 0:
        truemsg = file + " exists and not empty, continuing..."
        logging.info(truemsg)
    else:
        msg = (
            file
            + " either file is does not exist or is empty, please check files. Exiting."
        )
        logging.critical(msg)
        sys.exit(1)


def check_folders(folder):
    """
    Check the output folder if it exists, if not make new directory.
    """
    if os.path.exists(folder) is True:
        truemsg = folder + " output folder exists"
        logging.info(truemsg)
    else:
        os.makedirs(folder)
        msg = folder + " does not exist, making output folder"
        logging.info(msg)

def check_input_folder(folder):
    """
    Check the input folder contains the correct number and types of files.
    """
    input_files = []
    for file in os.listdir(folder):
        if not file.startswith("."):
            input_files.append(file)
    count_input = len(input_files)
    if count_input > 10:
        logging.critical(f"There more than the expected number of files in {folder}, please double check")
        sys.exit(1)
    for file in input_files:
        if file.endswith((".html", ".txt", ".png", ".csv", ".cov_stats")):
            check_files(os.path.join(folder, file))
        else:
            logging.critical(f"Invalid file type detected please check or remove")
            sys.exit(1)
    return input_files

def load_csv(file):
    input_file = pd.read_csv(file)
    return input_file

def check_na_files(folder):
    keywords = ['CzDna', 'CzRna']
    file_extension = 'zscore.csv'
    found_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if any(keyword.lower() in file.lower() for keyword in keywords) and file.lower().endswith(file_extension.lower()):
                found_files.append(os.path.join(root, file))

    if len(found_files) == 2:
        logging.info(f"Found both files")
    else:
        logging.error(f"Either you have extra files or are missing files")
    for file in found_files:
        logging.info(f"{file} is one of the files")
    return found_files

    
        