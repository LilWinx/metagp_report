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

def clinican_results(file, wgsid):
    """
    Get data like suspected pathogen, reason for clinican request and clinical information
    This data comes from the .csv file from RedCap.
    """
    ph_dictionary = {
        "mrn": "py_mrn_ph",
        "pn": "py_pn_ph",
        "seq_accession": "py_acc_ph",
        "dobage": "py_dob_ph",
        "req_date": "py_reqdate_ph",
        "date_onset": "py_date_onset_ph",
        "date_collection": "py_doc_ph",
        "seq_date": "py_seqdate_ph",
        "seq_run_name": "py_runname_ph",
        "rep_gen_date": "py_repdate_ph",
        "wgs_id": "py_wgsid_ph",
        "comments": "py_clincomm_ph",
        "sus_pathogen": "py_suspathogen_ph",
        "fin_sample_type": "py_sampletype_ph",
        "na_type": "py_natype_ph",
        "ext_kit": "py_extkit_ph",
        "library_prep": "py_libraryprep_ph",
        "seq_platform": "py_seqplatform_ph",
        "final_result": "py_finalpathogen_ph",
        "binf_interpretation": "py_riinterpret_ph",
        "binf_reason": "py_rireason_ph",
    }

    na_type_dict = {
        0: "DNA",
        1: "RNA",
        2: "DNA & RNA",
    }

    extraction_kit_dict = {
        1.0: "Qiagen DNeasy Ultraclean Kit",
        2.0: "ZymoBIOMICS DNA/RNA Miniprep Kits",
        3.0: "Roche Diagnostics High Pure PCR Template Preparation Kit",
        4.0: "Roche Diagnostics MagNA Pure 96 DNA and Viral NA Small/Large Volume Kit",
        5.0: "bioMérieux EasyMag",
    }

    library_prep_dict = {
        1.0: "Illumina NexteraXT Library Preparation Kit",
        2.0: "Illumina Nextera DNA Library Preparation Kit",
        3.0: "Twist Comprehensive Panel",
        4.0: "ONT Rapid Barcoding Kit",
        5.0: "RAPIDPrep Protocol",
        6.0: "NEBNext Ultra II RNA Library Prep Kit for Illumina",
    }

    seq_platorm_dict = {
        1.0: "Illumina iSeq",
        2.0: "Illumina MiniSeq",
        3.0: "Illumina NextSeq 500",
        4.0: "Illumina NextSeq 2000",
        5.0: "ONT MinION",
        6.0: "ONT GridION",
    }

    sample_data = pd.read_csv(file, sep=",", header=0, encoding='latin-1')
    pt_data = sample_data.iloc[:,:13] # split dataframe to only have first 13 columns
    pt_data.dropna(how='all', inplace=True) # drop the empty rows 
    pt_data.reset_index(drop=True, inplace=True) # reset the index
    result_data = sample_data.iloc[:, 13:] # 2nd split to get other half
    filt_pt = result_data[result_data['wgs_id'] == str(wgsid)] # only have the row with the correct wgsid
    filt_pt.reset_index(drop=True, inplace=True) # reset the index
    combined_pt = pd.concat([pt_data, filt_pt], axis = 1) # recombine the two dataframes
    combined_pt['pn'] = combined_pt['first_name'] + ' ' + combined_pt['last_name'] # combine first and last name for patient name
    combined_pt['dobage'] = combined_pt['dob'] + ' (' + str(combined_pt['age'][0]) + " y/o)" # combine dob & age
    #MAP TIME
    map_dict = {'na_type': na_type_dict, 
                'ext_kit': extraction_kit_dict, 
                'library_prep': library_prep_dict,
                'seq_platform': seq_platorm_dict,
                }
    
    for key, value in map_dict.items():
        combined_pt[key] = combined_pt[key].map(value)
    
    # return to regularly scheduled programming
    clean_pt_data = combined_pt.drop(columns=[col for col in combined_pt.columns if col not in ph_dictionary.keys()]) # drop unnecessary columns not needed in report
    clean_pt_data['mrn'] = clean_pt_data["mrn"].astype('int64')
    clean_pt_data = clean_pt_data.fillna("-").replace(0, "-")
    row = clean_pt_data.iloc[0]
    pt_dict = row.to_dict()
    pt_dict_map = {ph_dictionary.get(key, key): value for key, value in pt_dict.items()}
    return pt_dict_map


    
