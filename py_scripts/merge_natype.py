import os
import logging
import pandas as pd
from py_scripts import assists




"""
Script specifically to control the merge of zscore csv from both RNA and DNA
"""


def merge_dna_and_rna(found_list):
    for file_path in found_list:
        file_name = os.path.basename(file_path)
        if file_name.startswith("CzDna"):
            dna_file = file_path
        elif file_name.startswith("CzRna"):
            rna_file = file_path
    logging.info(f"DNA file: {dna_file}")
    logging.info(f"RNA file: {rna_file}")
    dna_df = assists.load_csv(dna_file)
    rna_df = assists.load_csv(rna_file)

    dna_df['na_type'] = "DNA"
    rna_df['na_type'] = "RNA"
    
    merged_df = pd.concat([dna_df, rna_df])
    
    if "krona" in dna_file:
        merged_df = merged_df.groupby(['#queryID', '#taxID']).agg({'#score': 'mean', 'rpm_sample': 'mean'}).reset_index()
    else:
        merged_df.reset_index(drop=True)
    return merged_df

