import logging
import sys
from py_scripts import arguments
parser = arguments.create_parser()  # pylint: disable=E1101
args = parser.parse_args()

pt_data_list = []
input_prompts = [
    f"Patient Name:",
    f"MRN:",
    f"Accession:",
    f"Date of Collection:",
    f"Sequence Date:",
    f"Run Name:",
    f"Date of Report:",
    f"WGS ID:",
]
def auto_txt_input():
    pt_details = open(args.txt, 'r')
    lines = pt_details.readlines()
    for line in lines:
        crop_line = line.strip("\n")
        pt_data_list.append(crop_line)
    line_count = len(lines)
    if line_count < 8:
        logging.error(f"Text file does not contain correct number of data fields, there should be 8")
        sys.exit(1)
    pt_details.close()
    return pt_data_list

def manual_input():
    for prompt in input_prompts:
        line = input(prompt)
        pt_data_list.append(line)
    return pt_data_list
