import os
import logging
import pandas as pd
from py_scripts import assists




"""
Script specifically to control the merge of zscore csv from both RNA and DNA
"""


def merge_dna_and_rna(found_list):
    dna_file = next((file for file in found_list if 'CzDna' in file), None)
    logging.info(f"DNA file: {dna_file}")
    rna_file = next((file for file in found_list if 'CzRna' in file), None)
    logging.info(f"RNA file: {rna_file}")
    dna_df = assists.load_csv(dna_file)
    rna_df = assists.load_csv(rna_file)

    dna_df['na_type'] = "DNA"
    rna_df['na_type'] = "RNA"

    merged_df = pd.concat([dna_df, rna_df])
    merged_df.reset_index(drop=True)
    return merged_df
