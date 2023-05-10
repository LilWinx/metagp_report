import os
import sys
import logging

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
    input_files = os.listdir(folder)
    count_input = len(input_files)
    if count_input > 3:
        logging.critical(f"There more than the expected number of files in {folder}, please double check")
        sys.exit(1)
    for file in input_files:
        if file.endswith(".html") or file.endswith(".txt"):
            check_files(os.path.join(folder, file))
        else:
            logging.critical(f"Invalid file type detected please check or remove")
            sys.exit(1)
    return input_files
    
        