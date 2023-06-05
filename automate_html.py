import os
import re
import sys
import argparse
import logging
import warnings
from py_scripts import arguments
from py_scripts import assists
from py_scripts import base64_encode
from py_scripts import pathogen_db_search
from py_scripts import result_interpret

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

    assists.check_folders(args.input)
    file_list = assists.check_input_folder(args.input)
    for file in file_list:
        if file.startswith("TestProjectWinkieFon") and file.endswith(".csv"):
            patient_data = result_interpret.clinican_results(os.path.join(args.input, file), args.wgsid)
        elif file.endswith("sp.txt"):
            species_list = result_interpret.species_list(os.path.join(args.input, file))
        elif file.endswith("krona.html"):
            iframe_krona = os.path.join(args.input, file)
            base64_krona = base64_encode.html_base64_encode(iframe_krona)
        elif file.endswith("_coverageplot.png"):
            coverage_png = os.path.join(args.input, file)
            base64_cov_png = base64_encode.html_base64_encode(coverage_png)
            match = re.search(".*_([A-Z]+_[A-Z0-9]+.*[0-9]*)_.*", file) # thx jake "_".join(file.split("_")[1:3])
            used_reference = match.group(1)
        elif file.endswith("_hbar.png"):
            hbar_png = os.path.join(args.input, file)
            base64_hbar_png = base64_encode.html_base64_encode(hbar_png)
    db_accordion = pathogen_db_search.pathogen_search(species_list)

    img_logo = os.path.join(os.path.dirname(__file__), "assets/nswhp-logo.png")
    mgimg_logo = os.path.join(os.path.dirname(__file__), "assets/metagp-logo.png")
    base64_logo_png = base64_encode.html_base64_encode(img_logo)
    base64_metagp = base64_encode.html_base64_encode(mgimg_logo)

    replace_dict = {
        "py_logo_ph": base64_logo_png,
        "py_metagplogo_ph": base64_metagp,
        "py_finalpathref_ph": used_reference,
        "py_krona_ph": base64_krona,
        "py_coverageimg_ph": base64_cov_png,
        "py_hbar_ph": base64_hbar_png
    }
    dictionary_list = [patient_data, db_accordion]
    for dictionary in dictionary_list:
        replace_dict.update(dictionary)
    print(replace_dict)
    with open(template_html, "r") as html_template:
        template = html_template.read()

    for key, value in replace_dict.items():
        logging.info(f"Replacing {key} to {value} in the HTML")
        template = template.replace(key, str(value))

    with open(output_file, "w") as outfile:
        outfile.write(template)

if __name__ == "__main__":
    main()