import sys
import os
import pandas as pd
from py_scripts import pathogen_db_search

# summarise all the data from the pipeline
# this script will extract the top 10 of each domain/kingdom
# specifically pulling Bacteria, Viruses, Fungi and Kingdom

def read_in_tpm(tpm_file):
    accordion_dict = {}
    colnames = ["TPM", 
                "Kingdom", 
                "Phylum", 
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
    ]
    
    fungi_phyla = [
        "Ascomycota",
        "Basidiomycota",
        "Chytridiomycota",
        "Microsporidia",
        "Glomeromycota",
        "Zygomycota"
    ]

    parasite_phyla = [
        "Apicomplexa",
        "Ciliophora",
        "Bacillariophyta",
        "Cercozoa",
        "Euglenozoa",
        "Heterolobosea",
        "Parabasalia",
        "Fornicata",
        "Evosea",
        "Streptophyta",
    ]

    input_file = pd.read_csv(tpm_file, sep="\t", names=colnames)
    tpm_input = input_file[input_file["Phylum"].str.contains("Chordata")==False]
    tpm_file = tpm_input[:50].reset_index(drop=True)
    d_bacteria = tpm_file[
        tpm_file['Kingdom'].str.contains("Bacteria")][:10].reset_index(drop=True)
    d_viruses = tpm_file[tpm_file['Kingdom'].str.contains("Viruses")][:10].reset_index(drop=True)
    k_fungi = tpm_file[tpm_file['Phylum'].str.contains('|'.join(fungi_phyla))][:10].reset_index(drop=True)
    k_parasite = tpm_file[tpm_file['Phylum'].str.contains('|'.join(parasite_phyla))][:10].reset_index(drop=True)
    for tpm_df in [d_bacteria, d_viruses, k_fungi, k_parasite]:
        tpm_df.sort_values(by="TPM", ascending = False)
    bacteria_list = d_bacteria['Species'].to_list()
    viruses_list = d_viruses['Species'].to_list()
    fungi_list = k_fungi['Species'].to_list()
    parasite_list = k_parasite['Species'].to_list()
    
    kd_list = [
        (bacteria_list, "Bacteria"), 
        (viruses_list, "Viruses"), 
        (fungi_list, "Fungi"), 
        (parasite_list, "Parasites")
    ]

    for species_list, dk_status in kd_list:
        accordion_stuff = pathogen_db_search.pathogen_search(species_list, dk_status)
        accordion_dict.update(accordion_stuff)
    return accordion_dict