import re
import os
import logging
from Bio import Entrez
import pandas as pd
from urllib.request import HTTPError
import time

Entrez.email = "winkie.fong@sydney.edu.au"

viruses="/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db/viruses.txt"
fungi = os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/fungi.txt")
parasites = os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/parasites.txt")
bacteria="/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db/prokaryotes-refonly.txt"

def open_virus_file(file):
    db_file = pd.read_csv(file, sep= "\t", header=0)
    db_file = db_file[db_file["Segmemts"].str.contains("NC_")==True]
    db_file["Ref_List"] = db_file["Segmemts"].apply(extract_NC)
    db_file["#Organism/Name"].drop_duplicates(inplace= True)
    for index, row in db_file.iterrows():
        db_file.at[index, 'FTP Path'] = get_ftp_path(row["Ref_List"])
    return db_file[["#Organism/Name", 'FTP Path']]

def open_eukaryotes_file(file):
    db_file = pd.read_csv(file, sep= "\t", header=0)
    db_file["#Organism/Name"].drop_duplicates(inplace= True)
    for index, row in db_file.iterrows():
        db_file.at[index, 'FTP Path'] = get_ftp_path_eukaryotes(row["Assembly Accession"])
    return db_file[["#Organism/Name", 'FTP Path']]

def extract_NC(string):
    matches = re.findall(r'NC_\d+\.\d+', string)
    if matches:
        return matches
    else:
        return None

def get_ftp_path(terms):
    ftp_paths = []
    for term in terms:
        try:
            search_handle = Entrez.esearch(db="assembly", term=term)
            search_record = Entrez.read(search_handle)
            search_handle.close()

            assembly_id = search_record["IdList"][0]
            assembly_handle = Entrez.esummary(db="assembly", id=assembly_id)
            assembly_summary = Entrez.read(assembly_handle)
            assembly_handle.close()

            ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
            ftp_paths.append(ftp_path)
        except HTTPError:
            time.sleep(2)
            continue
    return ftp_paths[0]

def get_ftp_path_eukaryotes(term):
    try:
        search_handle = Entrez.esearch(db="assembly", term=term)
        logging.info(f"Searching for {term}")
        search_record = Entrez.read(search_handle)
        search_handle.close()

        assembly_id = search_record["IdList"][0]
        logging.info(f"Found {term}, getting Assembly ID {assembly_id}")
        assembly_handle = Entrez.esummary(db="assembly", id=assembly_id)
        assembly_summary = Entrez.read(assembly_handle)
        assembly_handle.close()

        ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
        logging.info(f"Found {assembly_id}, ftp path for this sample is {ftp_path}, adding to dataframe")
    except HTTPError:
        time.sleep(2)
        logging.info("Too fast, slowing down")
    return ftp_path

def combine_ref_lists():
    virus = open_virus_file(viruses)
    logging.info(f"Opening {virus}")
    fungus = open_eukaryotes_file(fungi)
    logging.info(f"Opening {fungi}")
    parasite = open_eukaryotes_file(parasites)
    logging.info(f"Opening {parasites}")
    bacterium_db = pd.read_csv(bacteria, sep= "\t", header=0)
    bacterium = bacterium_db[["#Organism/Name", "FTP Path"]]
    logging.info(f"Opening {bacteria}")
    combined = pd.concat([virus, fungus, parasite, bacterium])
    combined.to_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/ref_list.txt"), sep="\t")
    logging.info(f"Almost done! joining to csv")
