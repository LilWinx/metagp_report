import sys
import re
import os
import pandas as pd
from Bio import Entrez
import ssl
import certifi

ssl._create_default_https_context = ssl._create_unverified_context


#file = sys.argv[1] 
# needs --use-mpa-style for kraken2
wd = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/bor/kraken"
#input="CSF5703A.kraken2.nonhuman.nonrRNA.std.report"
input = "3160270484.kraken.mpareport"
file = wd + "/" + input
outfile = wd + "/species.txt"
db_path = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db"


def db_read_in(file):
    db = pd.read_csv(file, sep="\t", header=0)
    return db

def read_count(sorted_kraken, domain):
    filt_kraken = sorted_kraken[sorted_kraken['tree'].str.contains(domain)]
    read_count = filt_kraken[filt_kraken['tree'] == domain]['reads'].values[0]
    filt_kraken['taxon_RA'] = (filt_kraken['reads']/read_count)*100
    return filt_kraken

def total_read_count(sorted_kraken):
    k_domain_list = ["d__Bacteria", "d__Viruses", "d__Eukaryota", "d__Archaea"]
    k_domain_readcount = []
    for domain in k_domain_list:
        k_domain = sorted_kraken[sorted_kraken['tree'] == domain]['reads'].values[0]
        k_domain_readcount.append(k_domain)
    return k_domain_readcount

def genome_size(species_df):
    bacteria_db = db_read_in(os.path.join(db_path, "prokaryotes-refonly.txt"))
    bacteria_size = bacteria_db[['#Organism/Name', 'TaxID', 'Size (Mb)']]
    virus_db = db_read_in(os.path.join(db_path, "viruses.txt"))    
    bacteria_size['Size (bp)'] = bacteria_size['Size (Mb)']*100000
    virus_size = virus_db[['#Organism/Name', 'TaxID', 'Size (Kb)']].drop_duplicates(subset='#Organism/Name')
    virus_size['Size (bp)'] = virus_size['Size (Kb)']*1000
    db_merge = pd.concat([bacteria_size, virus_size], ignore_index=True)
    species_merge = pd.merge(species_df, db_merge, left_on='species', right_on='#Organism/Name')
    return species_merge
    
def read_tsv(file):
    reorder_list = [
        'species',
        'TaxID',
        'Size (bp)',
        'reads',
        'taxon_RA',
        'total_RA',
        'GPM'
    ]
    kraken_input = pd.read_csv(file, sep="\t", header=None)
    kraken_input.columns = ['tree', "reads"]
    kraken_sort = kraken_input.sort_values(by="reads", ascending = False)
    total_reads = sum(total_read_count(kraken_sort))

    d_bacteria = read_count(kraken_sort, "d__Bacteria") # extract only bacteria
    d_viruses = read_count(kraken_sort, "d__Viruses") # extract only viruses
    #d_eukaryota = read_count(kraken_sort, "d__Eukaryota") # TO DO extract only fungal
    # TO DO extract only parasitic
    d_list = [d_bacteria, d_viruses]
    merged_df = pd.concat(d_list, ignore_index=True)
    filt_again = merged_df[merged_df['tree'].str.contains('s__')]
    str_regex = r"\|s__([^|]+)"
    filt_again['species'] = filt_again["tree"].str.extract(str_regex, expand=False)
    species_merge = genome_size(filt_again)
    species_merge['GPM'] = total_reads/(species_merge['reads'] * species_merge['Size (bp)'])
    species_merge['total_RA'] = (species_merge['reads']/total_reads)*100
    re_sorted_kraken = species_merge.sort_values(by="total_RA", ascending = False)
    re_sorted_kraken.drop(columns=['tree', 'Size (Mb)', 'Size (Kb)', '#Organism/Name'], inplace=True)
    top50 = re_sorted_kraken[:50].reset_index(drop=True)
    kraken_top10 = top50.sort_values(by="taxon_RA", ascending = False)[:10].reset_index(drop=True)
    kraken_top10 = kraken_top10[reorder_list]
    kraken_top10[['TaxID', 'Size (bp)']] = kraken_top10[['TaxID', 'Size (bp)']].astype('int64')
    kraken_top10.to_csv(outfile, sep="\t", index=False)

    
    #kraken_top50 = re_sorted_kraken[:1].reset_index(drop=True)
    #kraken_top50 = re_sorted_kraken[:50].reset_index(drop=True)

    #kraken_top50.to_csv(wd+"/merged.txt", sep="\t", index=False)
    #kraken_top10 = re_sorted_kraken[:10].reset_index()
    
    
    #
    return re_sorted_kraken

read_tsv(file)                            
