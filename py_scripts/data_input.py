import logging
import sys
from automate_html import Inputs


def auto_txt_input(pt_txt):
    with open(pt_txt, 'r') as pt_details:
        lines = pt_details.read().splitlines()
    line_count = len(lines)
    if line_count < 8:
        logging.error(f"Text file does not contain correct number of data fields, there should be 8")
        sys.exit(1)
    inputs = Inputs(*lines)
    return inputs

def manual_input():
    pn = input(f"Patient Name:"),
    mrn = input(f"MRN:"),
    accession = input(f"Accession:"),
    doc = input(f"Date of Collection:"),
    seqdate = input(f"Sequence Date:"),
    runname = input(f"Run Name:"),
    repdate = input(f"Date of Report:"),
    wgsid = input(f"WGS ID:")
    inputs = Inputs(pn, mrn, accession, doc, seqdate, runname, repdate, wgsid)
    return inputs


