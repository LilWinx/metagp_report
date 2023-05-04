import os
import sys
import argparse
import logging
import warnings
from py_scripts import arguments
from py_scripts import assists

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)
warnings.simplefilter(action="ignore", category=FutureWarning)
formatter = logging.Formatter(
    "Meta-GP_Report:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p"
)


def main():
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    assists.check_folders(args.outdir)

    # set outdir defaults - if no outdir is set, it will default to either the fasta or R1 location
    if args.outdir is None:
        default = os.getcwd()
        outdir = default
    else:
        outdir = args.outdir

    # force creation of new folder if not already exists
    folder_exists = os.path.exists(outdir)
    if not folder_exists:
        os.makedirs(outdir)

    if args.manual != True:
        assists.check_files(args.text)
        input_data = args.text
    else:
        pn = input(f"Patient Name:")
        mrn = input(f"MRN:")
        accession = input(f"Accession:")
        doc = input(f"Date of Collection:")
        seqdate = input(f"Sequence Date:")
        runname = input(f"Run Name:")
        repdate = input(f"Date of Report:")

if __name__ == "__main__":
    main()