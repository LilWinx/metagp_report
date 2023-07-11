import re
import os
from Bio import Entrez
import pandas as pd
from urllib.request import HTTPError
import time

Entrez.email = "winkie.fong@sydney.edu.au"

viruses = os.path.join(os.path.dirname(os.path.dirname(__file__)), "database/virus_nconly.txt")
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
    result = db_file[["#Organism/Name", 'FTP Path']]
    result.to_csv("/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db/virus_ref_list.txt", sep="\t")
    return result

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
            print(f"searching for {term}")

            assembly_id = search_record["IdList"][0]
            assembly_handle = Entrez.esummary(db="assembly", id=assembly_id)
            assembly_summary = Entrez.read(assembly_handle)
            assembly_handle.close()
            print(f"searching for {assembly_id}")

            ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
            ftp_paths.append(ftp_path)
            print(f"found ftp path {ftp_path}")
        except HTTPError:
            time.sleep(2)
            print(f"too fast!")
            continue
    return ftp_paths[0]

def get_ftp_path_eukaryotes(term):
    try:
        search_handle = Entrez.esearch(db="assembly", term=term)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        assembly_id = search_record["IdList"][0]
        assembly_handle = Entrez.esummary(db="assembly", id=assembly_id)
        assembly_summary = Entrez.read(assembly_handle)
        assembly_handle.close()

        ftp_path = assembly_summary["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"]
    except HTTPError:
        time.sleep(2)
    return ftp_path

open_virus_file(viruses)