import os
import logging
import pandas as pd
import numpy as np

species_dict = {}
status_dict = {}
def pathogen_search(species_list):
    """
    Fill in information for the accordion, using the input file.
    """
    species_list = list(map(lambda x:x.strip("\n"), species_list))
    pathogen_db = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/pathogen_list.csv"), header = 0)
    for i, species in enumerate(species_list, start = 1):
        match_species = pathogen_db[(pathogen_db['Species'] == species) | (pathogen_db['AltNames'] == species)]
        if not match_species.empty:
            status = match_species.iloc[0]['Status']
            additions = []
            if not pd.isnull(match_species.iloc[0]['Disease_type']):
                additions.append(" causing a " + str(match_species.iloc[0]['Disease_type']) + " infection")
            if not pd.isnull(match_species.iloc[0]['Disease_Name']):
                additions.append(" called " + str(match_species.iloc[0]['Disease_Name']))
            if not pd.isnull(match_species.iloc[0]['Commensal']):
                additions.append(". " + match_species.iloc[0]['Species'] + " is a commensal microbe")
            if not pd.isnull(match_species.iloc[0]['NNDSS_Notifiable']):
                additions.append(", it is also a notifiable disease in NSW")
            if additions:
                status += "".join(additions)
        else:
            status = "not a known pathogen"
        species_ph = f"py_species{i}_ph"
        status_ph = f"py_species{i}status_ph"
        species_dict[species_ph] = species
        status_dict[status_ph] = status
    species_dict.update(status_dict)
    return species_dict
    


