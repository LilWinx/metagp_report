import os
import sys
import argparse
import logging
import warnings
from py_scripts import arguments
from py_scripts import assists
from py_scripts import base64_encode
from py_scripts import data_input

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)
warnings.simplefilter(action="ignore", category=FutureWarning)
logging.basicConfig(level=logging.INFO, format="Meta-GP Report:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p")

def main():
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    logger = logging.getLogger()
    template_html = os.path.join(os.path.dirname(__file__), "index.html")
    assists.check_files(template_html)
    

    iframe_krona = args.krona
    assists.check_files(iframe_krona)
    base64_krona = base64_encode.html_base64_encode(iframe_krona)

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
        [pn, mrn, accession, doc, seqdate, runname, repdate, wgsid] = data_input.auto_txt_input()
    else:
        [pn, mrn, accession, doc, seqdate, runname, repdate, wgsid] = data_input.manual_input()
    
    replace_dict = {
        "py_pn_ph": pn,
        "py_mrn_ph": mrn,
        "py_acc_ph": accession,
        "py_doc_ph": doc,
        "py_seqdate_ph": seqdate,
        "py_rundate_ph": runname,
        "py_repdate_ph": repdate,
        "py_wgsid_ph": wgsid,
        "py_krona_ph": base64_krona
    }

    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        template = template.replace(key, value)

    with open(output_file, "w") as outfile:
        outfile.write(template)

if __name__ == "__main__":
    main()