import sys
import os
import pandas as pd
from py_scripts import pathogen_db_search

# summarise all the data from the pipeline
# this script will extract the top 10 of each domain/kingdom
# specifically pulling Bacteria, Viruses, Fungi and Kingdom

def read_in_tpm(tpm_file):
    accordion_dict = {}

    fungi_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/fungi.txt"), sep = "\t", header = 0)
    parasites_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/parasites.txt"), sep = "\t" , header = 0)
    input_file = pd.read_csv(tpm_file)
    tpm_input = input_file[input_file["phylum"].str.contains("Chordata")==False]
    tpm_input = tpm_input[tpm_input['ranked_score'] > 1]
    d_bacteria = tpm_input[
        tpm_input['superkingdom'].str.contains("Bacteria")][:10].reset_index(drop=True)
    d_viruses = tpm_input[tpm_input['superkingdom'].str.contains("Viruses")][:10].reset_index(drop=True)
    k_fungi = tpm_input[tpm_input['species'].isin(fungi_db['#Organism/Name'])][:10].reset_index(drop=True)
    k_parasite = tpm_input[tpm_input['species'].isin(parasites_db['#Organism/Name'])][:10].reset_index(drop=True)
    for tpm_df in [d_bacteria, d_viruses, k_fungi, k_parasite]:
        tpm_df.sort_values(by="tpm", ascending = False)
    bacteria_list = d_bacteria['species'].to_list()
    viruses_list = d_viruses['species'].to_list()
    fungi_list = k_fungi['species'].to_list()
    parasite_list = k_parasite['species'].to_list()
    
    kd_list = [
        (bacteria_list, "Bacteria"), 
        (viruses_list, "Viruses"), 
        (fungi_list, "Fungi"), 
        (parasite_list, "Parasites")
    ]
    
    for species_list, dk_status in kd_list:
        accordion_stuff = pathogen_db_search.pathogen_search(species_list, dk_status)
        accordion_dict.update(accordion_stuff)
        score_stuff = get_ranked_score(tpm_input, species_list, dk_status)
        accordion_dict.update(score_stuff)
    return accordion_dict

def get_ranked_score(og_df, species_list, dk_status):
    score_dict = {}
    for i, species in enumerate(species_list, start = 1):
        if species != '-':
            match_row = og_df.loc[og_df['species'] == species]
            ranked_score = round(match_row['ranked_score'].values[0], 2)
        else: 
            ranked_score = '-'
        if dk_status == "Bacteria":
            zscore_ph = f"py_species{i}zscore_ph"
        if dk_status == "Viruses":
            zscore_ph = f"py_vspecies{i}zscore_ph"
        if dk_status == "Fungi":
            zscore_ph = f"py_fspecies{i}zscore_ph"
        if dk_status == "Parasites":
            zscore_ph = f"py_pspecies{i}zscore_ph"
        score_dict[zscore_ph] = ranked_score
    return score_dict