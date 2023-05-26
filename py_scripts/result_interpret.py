import os
import logging
import pandas as pd

def species_list(file):
    """
    Fill in information for the accordion, using the input file.
    """
    species_file = pd.read_csv(file, sep="\t", header=0)
    species_list = species_file['species'].to_list()
    return species_list

def clinican_results(file, mrn):
    """
    Get data like suspected pathogen, reason for clinican request and clinical information
    This data comes from the .csv file from RedCap.
    """
    ph_dictionary = {
        "mrn": "py_mrn_ph",
        "pn": "py_pn_ph",
        "seq_accession": "py_acc_ph",
        "dobage": "py_dob_ph",
        "date_req": "py_reqdate_ph",
        "date_collection": "py_doc_ph",
        "seq_date": "py_seqdate_ph",
        "seq_run_name": "py_runname_ph",
        "rep_gen_date": "py_repdate_ph",
        "wgs_id": "py_wgsid_ph",
        "comments": "py_clincomm_ph",
        "sus_pathogen": "py_suspathogen_ph",
        "sample_type": "py_sampletype_ph"
    }
    sample_type_dict = {
        1.0: "Whole Blood",
        2.0: "Serum",
        3.0: "Plasma",
        4.0: "CSF",
        5.0: "Nasopharyngeal Aspirate",
        6.0: "Nasopharyngeal Swab",
        11.0: "Sputum", 
        12.0: "Oral swab",
        13.0: "Wound Swab",
        7.0: "Saliva",
        8.0: "Urine",
        9.0: "Stool",
        10.0: "Rectal Swab",
    }

    pt_data = pd.read_csv(file, sep=",", header=0)
    filt_pt = pt_data[pt_data['mrn'] == int(mrn)]
    filt_pt.reset_index(drop=True, inplace=True)
    filt_pt['pn'] = filt_pt['first_name'] + ' ' + filt_pt['last_name']
    filt_pt['dobage'] = filt_pt['dob'] + ' (' + str(filt_pt['age'][0]) + " y/o)"
    filt_pt['sample_type'] = filt_pt['sample_type'].map(sample_type_dict)
    null_mask = filt_pt['sample_type'].isnull()
    filt_pt.loc[null_mask, 'sample_type'] = filt_pt.loc[null_mask, 'sample_other']
    clean_pt_data = filt_pt.drop(columns=[col for col in filt_pt.columns if col not in ph_dictionary.keys()])
    row = clean_pt_data.iloc[0]
    pt_dict = row.to_dict()
    pt_dict_map = {ph_dictionary.get(key, key): value for key, value in pt_dict.items()}
    print(pt_dict_map)
    return pt_dict_map


    
