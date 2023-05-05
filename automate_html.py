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
    
    template_html = os.path.join(os.path.dirname(__file__), "index.html")
    assists.check_files(template_html)

    inframe_krona = args.krona
    assists.check_files(inframe_krona)

    # set outdir defaults - if no outdir is set, it will default to krona location
    if args.outdir is None:
        outdir = os.path.dirname(args.krona)
    else:
        outdir = args.outdir

    # force creation of new folder if not already exists
    folder_exists = os.path.exists(outdir)
    if not folder_exists:
        os.makedirs(outdir)
    output_file = os.path.join(outdir, "report.html")
    assists.check_folders(outdir)

    if args.manual != True:
        assists.check_files(args.txt)
        pt_details = open(args.txt, 'r')
        lines = pt_details.readlines()
        line_count = len(lines)
        if line_count < 8:
            logging.error(f"Text file does not contain correct number of data fields, there should be 8")
        pn = lines[0].strip("\n")
        mrn = lines[1].strip("\n")
        accession = lines[2].strip("\n")
        doc = lines[3].strip("\n")
        seqdate = lines[4].strip("\n")
        runname = lines[5].strip("\n")
        repdate = lines[6].strip("\n")
        wgsid = lines[7].strip("\n")
        pt_details.close()
    else:
        pn = input(f"Patient Name:")
        mrn = input(f"MRN:")
        accession = input(f"Accession:")
        doc = input(f"Date of Collection:")
        seqdate = input(f"Sequence Date:")
        runname = input(f"Run Name:")
        repdate = input(f"Date of Report:")
        wgsid = input(f"WGS ID:")
    
    
    replace_dict = {
        "py_pn_ph": pn,
        "py_mrn_ph": mrn,
        "py_acc_ph": accession,
        "py_doc_ph": doc,
        "py_seqdate_ph": seqdate,
        "py_rundate_ph": runname,
        "py_repdate_ph": repdate,
        "py_wgsid_ph": wgsid,
        "py_krona_ph": inframe_krona
    }

    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        template = template.replace(key, value)

    with open(output_file, "w") as outfile:
        outfile.write(template)

if __name__ == "__main__":
    main()