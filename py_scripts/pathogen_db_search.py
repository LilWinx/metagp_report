import os
import logging
import pandas as pd
import numpy as np

def pathogen_search(species_list, dk_status):
    """
    Fill in information for the accordion, using the input file.
    """
    species_dict = {}
    status_dict = {}    
    zscore_dict = {}
    pathogen_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/pathogen_list.csv"), header = 0)
    contam_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/known_contaminants.csv"), header = 0)
    logging.info("Reading Pathogen_List database to fill the species list with metadata")
    if len(species_list) < 10:
        logging.info(f"There was not enough species in {dk_status} detected, filling with dashes.")
        species_list += ['-'] * (10 - len(species_list))
    
    for i, species in enumerate(species_list, start = 1):
        match_species = pathogen_db[(pathogen_db['Species'].str.lower() == species.lower()) | (pathogen_db['AltNames'].str.lower() == species.lower())]
        contam_species = contam_db[(contam_db['Species'] == species.lower())]
        if not contam_species.empty:
            
            species = f"{species}"
            if contam_species['Reason'].isna().all() is True:
                contam_status = contam_status.iloc[0]['Statement']
            else:
                contam_status = ". " + species + " " + contam_species.iloc[0]['Statement'] + " " + contam_species.iloc[0]['Reason']
        else:
            contam_status = ""

        if not match_species.empty:
            species = f"{species} \u2757"
            status = match_species.iloc[0]['Status']
            additions = []
            if not pd.isnull(match_species.iloc[0]['Disease_type']):
                additions.append(" causing a " + str(match_species.iloc[0]['Disease_type']) + " infection")
            if not pd.isnull(match_species.iloc[0]['Disease_Name']):
                additions.append(" called " + str(match_species.iloc[0]['Disease_Name']))
            if not pd.isnull(match_species.iloc[0]['Commensal']):
                commensal_locations = match_species.iloc[0]['Commensal'].split('; ')
                if len(commensal_locations) > 1:
                    commensal_locations[-2] += " and " + commensal_locations[-1]
                    commensal_locations = commensal_locations[:-1]
                additions.append(". " + match_species.iloc[0]['Species'] + " is a commensal microbe found in the " + ', '.join(commensal_locations))
            if not pd.isnull(match_species.iloc[0]['NNDSS_Notifiable']):
                additions.append(", it is also a notifiable disease in NSW")
            if additions:
                status += "".join(additions)
        elif species == "-":
            status = ""
        else:
            status = " is not a known pathogen"
        status = status + contam_status

        if dk_status == "Bacteria":
            species_ph = f"py_species{i}_ph"
            status_ph = f"py_species{i}status_ph"
            zscore_ph = f"py_species{i}zscore_ph"
        if dk_status == "Viruses":
            species_ph = f"py_vspecies{i}_ph"
            status_ph = f"py_vspecies{i}status_ph"
            zscore_ph = f"py_vspecies{i}zscore_ph"
        if dk_status == "Fungi":
            species_ph = f"py_fspecies{i}_ph"
            status_ph = f"py_fspecies{i}status_ph"
            zscore_ph = f"py_fspecies{i}zscore_ph"
        if dk_status == "Parasites":
            species_ph = f"py_pspecies{i}_ph"
            status_ph = f"py_pspecies{i}status_ph"
            zscore_ph = f"py_pspecies{i}zscore_ph"
        if dk_status == "Top10":
            species_ph = f"py_tspecies{i}_ph"
            status_ph = f"py_tspecies{i}status_ph"
            zscore_ph = f"py_tspecies{i}zscore_ph"
        species_dict[species_ph] = species
        status_dict[status_ph] = species.replace(" \u2757", " ") + status
    
    list_of_dicts = [status_dict, zscore_dict]
    for dictionary in list_of_dicts:
        species_dict.update(dictionary)
    return species_dict
    


