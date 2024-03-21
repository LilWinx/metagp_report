import os
import re
import sys
import argparse
import logging
import warnings
import gen_pdf
from py_scripts import arguments
from py_scripts import assists
from py_scripts import base64_encode
from py_scripts import pathogen_db_search
from py_scripts import result_interpret
#from py_scripts import tpmAbundances_top10
from py_scripts import zscore_top10
from py_scripts import merge_natype


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

    na_type = None
    if args.natype == "both":
        zscore_found_list = assists.check_na_files(args.input, "zscore")
        krona_found_list = assists.check_na_files(args.input, "krona")
        tpm_file = merge_natype.merge_dna_and_rna(zscore_found_list)
        krona_file = merge_natype.merge_dna_and_rna(krona_found_list)
        new_krona = args.input + "/merged_krona_input.tsv"
        krona_file.to_csv(new_krona, sep = '\t', index=False)
        krona_cmd = f"ktImportTaxonomy -m 4 -i {new_krona} -o {args.input}/merged-rpm-zscore.html"
        assists.run_cmd(krona_cmd)

    elif args.natype == "dna":
        na_type = "DNA"
    elif args.natype == "rna":
        na_type = "RNA"
    else:
        na_type = None

    coverage_png = None
    div_base64_cov_png = ""
    used_reference = ""
    fig_no = 1
    img_type = None
    for file in file_list:
        if file.startswith("MetaGP") and file.endswith(".csv"):
            patient_data = result_interpret.clinican_results(os.path.join(args.input, file), args.wgsid)
        elif file.endswith("_coverageplot.png"):
            coverage_png = os.path.join(args.input, file)
            covstats = os.path.join(args.input, "coverage_ref.cov_stats")
            with open(covstats, 'r') as cov_stats:
                lines = cov_stats.readlines()
                second_row = lines[1].rstrip('\n')  # Get the second row and remove newline character
                used_reference = second_row.split('\t')[0]  # Extract first instance before first tab
            base64_cov_png = base64_encode.html_base64_encode(coverage_png)
            div_base64_cov_png = f'''
                <div class="coverage_plot">
                    <p style="padding-left: 10px; font-size: 12px">
                        <b>Figure {fig_no}:</b> Whole genome coverage map of sequencing reads to {used_reference}.
                    </p>
                    <div class="centre_img">
                        <img src="data:image/png;base64, {base64_cov_png}" class="covplot">
                        <br>
                    </div>
                </div>
            '''
            fig_no += 1
        elif file.endswith("_hbar.png"):
            hbar_png = os.path.join(args.input, file)
            base64_hbar_png = base64_encode.html_base64_encode(hbar_png)
            div_hbar_cov_png = f'''
            <div>
                <p style="padding-left: 10px; font-size: 12px">
                    <b>Figure {fig_no}:</b> Relative abundance of the Top 10 species following filtering.
                </p>
                <br>
                <div class="centre_img">
                    <img src="data:image/png;base64, {base64_hbar_png}" class="covplot">
                </div>
                <br>      
            </div>
            '''
            fig_no += 1
            img_type = "hbar"
        elif file.endswith("_quaddonut.png"):
            hbar_png = os.path.join(args.input, file)
            base64_qd_png = base64_encode.html_base64_encode(hbar_png)
            div_hbar_cov_png = f'''
            <div>
                <p style="padding-left: 10px; font-size: 12px">
                    <b>Figure {fig_no}:</b> Relative abundance of the Top 10 species following filtering from both DNA and RNA mNGS.
                </p>
                <br>
                <div class="centre_img">
                    <img src="data:image/png;base64, {base64_qd_png}" class="donut">
                </div>
                <br>      
            </div>
            '''
            fig_no += 1
            img_type = "donut"
        elif file.endswith("zscore.csv"):
            if args.natype == "dna" or args.natype == "rna" and na_type != None:
                tpm_file = os.path.join(args.input, file)
                tpm_df = assists.load_csv(tpm_file)
                db_accordion = zscore_top10.read_in_tpm(tpm_df, na_type)
            elif args.natype == "both":
                db_accordion = zscore_top10.read_in_tpm(tpm_file, na_type)
            else:
                tpm_file = os.path.join(args.input, file)
                tpm_df = assists.load_csv(tpm_file)
                db_accordion = zscore_top10.read_in_tpm(tpm_df, na_type)
        elif file.endswith("rpm-zscore.html"):
            iframe_krona = os.path.join(args.input, file)
            if args.natype == "both":
                iframe_krona = os.path.join(args.input, "merged-rpm-zscore.html")
            base64_krona = base64_encode.html_base64_encode(iframe_krona)
            div_krona_html = f'''
            <div data-html2canvas-ignore="true" class="krona">
                <p style="padding-left: 10px;"><b>Figure {fig_no}:</b> Taxonomic Classification of contigs and unassembled reads post-pipeline.</p>
                <iframe src="data:text/html;base64,{base64_krona}">
                </iframe>
            </div>
            '''

    # base64 encode all set images
    img_names = ["nswhp-logo.png", 
                 "metagp-logo.png", 
                 "bacteria.png", 
                 "virus.png", 
                 "fungi.png", 
                 "parasite.png"]
    base64_images = []
    images = []

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
        "py_krona_ph": div_krona_html,
        "py_coverageimg_ph": div_base64_cov_png,
        "py_hbarordonuts_ph": div_hbar_cov_png
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

    # replace base64 back to original image for pdf gen
    for img in img_names:
        img_path = os.path.join(os.path.dirname(__file__), "assets", img)
        variable_name = "img_" + os.path.splitext(img)[0]
        locals()[variable_name] = img_path
        images.append(variable_name)

    update_dict = {
        "py_logo_ph": locals()[images[0]],
        "py_metagplogo_ph": locals()[images[1]],
        "py_bacteria_icon_ph": locals()[images[2]],
        "py_virus_icon_ph": locals()[images[3]],
        "py_fungi_icon_ph": locals()[images[4]],
        "py_parasite_icon_ph": locals()[images[5]],
        "py_coverageimg_ph": coverage_png,
        "py_hbar_ph": hbar_png
    }
    replace_dict.update(update_dict)
    gen_pdf.pdf_template(outdir, replace_dict, used_reference, img_type)


if __name__ == "__main__":
    main()