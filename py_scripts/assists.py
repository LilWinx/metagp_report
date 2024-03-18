import os
import sys
import logging
import pandas as pd
import subprocess

def run_cmd(command):
    """
    Run commands with error outputs.
    """
    logging.info("Running command: %s", command)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    result.stdout = result.stdout.decode()
    result.stderr = result.stderr.decode()
    if result.returncode != 0:
        logging.critical("Failed to run command: %s", result.args)
        logging.critical("stdout: %s", result.stdout)
        logging.critical("stderr: %s", result.stderr)
        sys.exit(1)
    return result

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
    if count_input > 11:
        logging.critical(f"There more than the expected number of files in {folder}, please double check")
        sys.exit(1)
    for file in input_files:
        if file.endswith((".html", ".txt", ".png", ".csv", ".cov_stats", ".tsv")):
            check_files(os.path.join(folder, file))
        else:
            logging.critical(f"Invalid file type detected please check or remove")
            sys.exit(1)
    return input_files

def load_csv(file):
    if file.endswith(".csv"):
        input_file = pd.read_csv(file)
    elif file.endswith(".tsv"):
        input_file = pd.read_csv(file, sep="\t")
    else:
        logging.error("Invalid file type, how did you even get to this point???")
    return input_file

def check_na_files(folder, type):
    keywords = ['CzDna', 'CzRna']
    if type == "zscore":
        file_extension = 'zscore.csv'
    elif type == "krona":
        file_extension = "krona.tsv"
    else:
        logging.info("No type specified.")
    found_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if any(keyword.lower() in file.lower() for keyword in keywords) and file.lower().endswith(file_extension.lower()):
                found_files.append(os.path.join(root, file))

    if len(found_files) == 2:
        logging.info(f"Found both {type} files")
    else:
        logging.error(f"Either you have extra files or are missing files")
    for file in found_files:
        logging.info(f"{file} is one of the files")
    return found_files

    
        