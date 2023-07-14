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
from py_scripts import tpmAbundances_top10

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
        if file.startswith("MetaGP") and file.endswith(".csv"):
            patient_data = result_interpret.clinican_results(os.path.join(args.input, file), args.wgsid)
        elif file.endswith("kraken2.html"):
            iframe_krona = os.path.join(args.input, file)
            base64_krona = base64_encode.html_base64_encode(iframe_krona)
        elif file.endswith("_hbar.png"):
            hbar_png = os.path.join(args.input, file)
            base64_hbar_png = base64_encode.html_base64_encode(hbar_png)
        elif file.endswith("tpmAbundances.txt"):
            tpm_file = os.path.join(args.input, file)
            db_accordion = tpmAbundances_top10.read_in_tpm(tpm_file)
        if file.endswith("_coverageplot.png"):
            coverage_png = os.path.join(args.input, file)
            base64_cov_png = base64_encode.html_base64_encode(coverage_png)
            match = re.search(".*_([A-Z]+_[A-Z0-9]+.*[0-9]*)_.*", file) # thx jake "_".join(file.split("_")[1:3])
            used_reference = match.group(1)
            div_base64_cov_png = f'''
                <div class="coverage_plot">
                    <p style="padding-left: 10px; font-size: 12px">
                        <b>Figure 1:</b> Whole genome coverage map of sequencing reads to {used_reference}.
                    </p>
                    <div class="centre_img">
                        <img src="data:image/png;base64, {base64_cov_png}">
                        <br>
                    </div>
                </div>
            '''
        else:
            div_base64_cov_png = ""

    # base64 encode all set images
    img_names = ["nswhp-logo.png", 
                 "metagp-logo.png", 
                 "bacteria.png", 
                 "virus.png", 
                 "fungi.png", 
                 "parasite.png"]
    base64_images = []

    for img in img_names:
        img_path = os.path.join(os.path.dirname(__file__), "assets", img)
        base64_image = base64_encode.html_base64_encode(img_path)
        variable_name = "base64_" + os.path.splitext(img)[0]
        locals()[variable_name] = base64_image
        base64_images.append(variable_name)

    replace_dict = {
        "py_logo_ph": locals()[base64_images[0]],
        "py_metagplogo_ph": locals()[base64_images[1]],
        "py_bacteria_icon_ph": locals()[base64_images[2]],
        "py_virus_icon_ph": locals()[base64_images[3]],
        "py_fungi_icon_ph": locals()[base64_images[4]],
        "py_parasite_icon_ph": locals()[base64_images[5]],
        "py_krona_ph": base64_krona,
        "py_coverageimg_ph": div_base64_cov_png,
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