import sys
import os
import pandas as pd
from py_scripts import pathogen_db_search

# summarise all the data from the pipeline
# this script will extract the top 10 of each domain/kingdom
# specifically pulling Bacteria, Viruses, Fungi and Kingdom

def read_in_tpm(input_file, na_type):
    accordion_dict = {}
    fungi_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/fungi.txt"), sep = "\t", header = 0)
    parasites_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/parasites.txt"), sep = "\t" , header = 0)
    tpm_input = input_file[input_file["phylum"].str.contains("Chordata")==False]
    tpm_input = tpm_input[tpm_input['zscore'] > 1]
    tpm_input = tpm_input.sort_values(by=['zscore', 'rpm_sample',], ascending=[False, False]).dropna(subset='species').drop_duplicates(subset='species')
    top10only = top10only = tpm_input[tpm_input['zscore'] > 50].nlargest(10, 'zscore')['species'].reset_index(drop=True)
    if na_type == "DNA" or na_type == "RNA":
        top10only['na_type'] == na_type
    else:
        pass

    # Filter and select top 10 for Bacteria with zscore > 50
    d_bacteria = tpm_input[tpm_input['superkingdom'].str.contains("Bacteria") & (tpm_input['zscore'] > 50)].nlargest(10, 'zscore').reset_index(drop=True)

    # Filter and select top 10 for Viruses with zscore > 50
    d_viruses = tpm_input[tpm_input['superkingdom'].str.contains("Viruses") & (tpm_input['zscore'] > 50)].nlargest(10, 'zscore').reset_index(drop=True)

    # Filter and select top 10 for Fungi with zscore > 50
    k_fungi = tpm_input[tpm_input['species'].isin(fungi_db['#Organism/Name']) & (tpm_input['zscore'] > 50)].nlargest(10, 'zscore').reset_index(drop=True)

    # Filter and select top 10 for Parasites with zscore > 50
    k_parasite = tpm_input[tpm_input['species'].isin(parasites_db['#Organism/Name']) & (tpm_input['zscore'] > 50)].nlargest(10, 'zscore').reset_index(drop=True)

    for tpm_df in [d_bacteria, d_viruses, k_fungi, k_parasite]:
        tpm_df.sort_values(by="rpm_sample", ascending = False)
    bacteria_list = d_bacteria['species'].to_list()
    viruses_list = d_viruses['species'].to_list()
    fungi_list = k_fungi['species'].to_list()
    parasite_list = k_parasite['species'].to_list()
    top10only_list = top10only.to_list()
    
    kd_list = [
        (bacteria_list, "Bacteria"), 
        (viruses_list, "Viruses"), 
        (fungi_list, "Fungi"), 
        (parasite_list, "Parasites"),
        (top10only_list, "Top10")
    ]
    
    for species_list, dk_status in kd_list:
        accordion_stuff = pathogen_db_search.pathogen_search(species_list, dk_status)
        accordion_dict.update(accordion_stuff)
        score_stuff = get_ranked_score(tpm_input, species_list, dk_status)
        accordion_dict.update(score_stuff)
    return accordion_dict

def get_ranked_score(og_df, species_list, dk_status):
    score_dict = {}
    na_type = "-"
    for i, species in enumerate(species_list, start = 1):
        if species != '-':
            match_row = og_df.loc[og_df['species'] == species]
            ranked_score = round(match_row['zscore'].values[0], 2)
            na_type = match_row['na_type'].values[0]
        else: 
            ranked_score = '-'
            na_type = '-'
        if dk_status == "Bacteria":
            zscore_ph = f"py_species{i}zscore_ph"
        if dk_status == "Viruses":
            zscore_ph = f"py_vspecies{i}zscore_ph"
        if dk_status == "Fungi":
            zscore_ph = f"py_fspecies{i}zscore_ph"
        if dk_status == "Parasites":
            zscore_ph = f"py_pspecies{i}zscore_ph"
        if dk_status == "Top10":
            zscore_ph = f"py_tspecies{i}zscore_ph"
        score_dict[zscore_ph] = f"Z-score: {ranked_score} | Result NA Type: {na_type}"
    return score_dict