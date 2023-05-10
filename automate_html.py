import os
import sys
import argparse
import logging
import warnings
from py_scripts import arguments
from py_scripts import assists
from py_scripts import base64_encode
from py_scripts import data_input
from py_scripts import pathogen_db_search

__version__ = "0.0.1"
logging.getLogger().setLevel(logging.INFO)
warnings.simplefilter(action="ignore", category=FutureWarning)
logging.basicConfig(level=logging.INFO, format="Meta-GP Report:%(levelname)s:%(asctime)s: %(message)s", datefmt="%y/%m/%d %I:%M:%S %p")

class Inputs:
    def __init__(self, pn, mrn, accession, doc, seqdate, runname, repdate, wgsid):
        self.pn = pn
        self.mrn = mrn
        self.accession = accession
        self.doc = doc
        self.seqdate = seqdate
        self.runname = runname
        self.repdate = repdate
        self.wgsid = wgsid

def main():
    parser = arguments.create_parser()  # pylint: disable=E1101
    args = parser.parse_args()
    logger = logging.getLogger()
    template_html = os.path.join(os.path.dirname(__file__), "index.html")
    assists.check_files(template_html)
    
    # set outdir defaults - if no outdir is set, it will default to krona location
    if args.output is None:
        outdir = os.path.dirname(args.input)
    else:
        outdir = args.output

    # force creation of new folder if not already exists
    folder_exists = os.path.exists(outdir)
    if not folder_exists:
        os.makedirs(outdir)
    output_file = os.path.join(outdir, "report.html")
    assists.check_folders(outdir)

    iframe_krona = args.krona
    assists.check_files(iframe_krona)
    base64_krona = base64_encode.html_base64_encode(iframe_krona)

    if args.manual != True and args.input is None:
        assists.check_files(args.txt)
        inputs = data_input.auto_txt_input(args.txt)
    elif args.input is not None:
        assists.check_folders(args.input)
        file_list = assists.check_input_folder(args.input)
        for file in file_list:
            if file.endswith("pt.txt"):
                patient_txt = os.path.join(args.input, file)
            elif file.endswith("sp.txt"):
                species_details = open(os.path.join(args.input, file), 'r')
                species_list = species_details.readlines()
            elif file.endswith("krona.html"):
                iframe_krona = os.path.join(args.input, file)
                base64_krona = base64_encode.html_base64_encode(iframe_krona)
            elif file.endswith(".png"):
                coverage_png = os.path.join(args.input, file)
                base64_cov_png = base64_encode.html_base64_encode(coverage_png)
        inputs = data_input.auto_txt_input(patient_txt)
        db_accordion = pathogen_db_search.pathogen_search(species_list)

    else:
        inputs = data_input.manual_input()
    
    img_logo = os.path.join(os.path.dirname(__file__), "assets/nswhp-logo.png")
    base64_logo_png = base64_encode.html_base64_encode(img_logo)
    
    replace_dict = {
        "py_logo_ph": base64_logo_png,
        "py_pn_ph": inputs.pn,
        "py_mrn_ph": inputs.mrn,
        "py_acc_ph": inputs.accession,
        "py_doc_ph": inputs.doc,
        "py_seqdate_ph": inputs.seqdate,
        "py_rundate_ph": inputs.runname,
        "py_repdate_ph": inputs.repdate,
        "py_wgsid_ph": inputs.wgsid,
        "py_krona_ph": base64_krona,
        "py_coverageimg_ph": base64_cov_png,
    }

    replace_dict.update(db_accordion)
    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        logging.info(f"Replacing {key} to {value} in the HTML")
        template = template.replace(key, value)

    with open(output_file, "w") as outfile:
        outfile.write(template)

if __name__ == "__main__":
    main()