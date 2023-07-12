import re
import sys
import subprocess
from Bio import Entrez
import pandas as pd
from urllib.request import HTTPError
import time

Entrez.email = "winkie.fong@sydney.edu.au"
Entrez.api_key = "dfa1ec2df7ac97add34b8fdec4d128a58a09"

query = sys.argv[1]
outdir = sys.argv[2]
#db = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db/full-ncbi-list.txt"
db = "/project/MetaGP/ncbi_db/full-ncbi-list.txt"

def search_query(query):
    ftp_path = None
    ncbi_db = pd.read_csv(db, sep = "\t", header=0)
    match = ncbi_db[ncbi_db["#Organism/Name"].str.lower() == query.lower()]
    if len(match) > 0:
        bacteria_match = match[match["Kingdom/Phylum"] == "Bacteria"]
        if len(bacteria_match) > 0:
            ftp_path = bacteria_match["FTP Path"].iloc[0]
        virus_match = match[match["Kingdom/Phylum"] == "Viruses"]
        if len(virus_match) > 0:
            nc_list = match["NC_SearchTerm"].apply(extract_NC)   
            ftp_path = get_ftp_path(nc_list)
        eukaryote_match = match[match["Kingdom/Phylum"] == "Eukaryotes"]
        if len(eukaryote_match) > 0:
            ftp_path = get_ftp_path(eukaryote_match["Assembly Accession"])
    print(ftp_path)
    download_ref(ftp_path, outdir)

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

def download_ref(ftp_path):
    rsync_path = ftp_path.replace("ftp:", "rsync:")
    retries = 10
    for attempt in range(1, retries + 1):
        try:
            subprocess.call(["rsync", "--copy-links", "--recursive", "--times", "--verbose", rsync_path, outdir])
            print(f"Downloading {rsync_path} from NCBI")
        except subprocess.CalledProcessError as e:
            print(f"Download attempt {attempt} failed: {e}")
            print(f"Retrying...")
            time.sleep(1)
    print('Maximum number of retries reached. Download failed.')

search_query(query)