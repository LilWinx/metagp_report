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
wd = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/test/input"
#input="CSF5703A.kraken2.nonhuman.nonrRNA.std.report"
input = "SRR14351804.kraken.mpareport"
file = wd + "/" + input
outfile = wd + "/species.txt"
db_path = "/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/db"


def db_read_in(file):
    db = pd.read_csv(file, sep="\t", header=0)
    db['#Organism/Name'] = db['#Organism/Name'].apply(remove_xtra_str)
    return db

def remove_xtra_str(value):
    words = value.split()
    if len(words) > 2 and 'subsp.' not in value:
        value = ' '.join(words[:2])
    return value

def read_count(sorted_kraken, taxon):
    filt_kraken = sorted_kraken[sorted_kraken['tree'].str.contains(taxon)]
    if taxon == 'k__Fungi':
        read_count = filt_kraken[filt_kraken['tree'] == "d__Eukaryota|k__Fungi"]['reads'].values[0]
    else:
        read_count = filt_kraken[filt_kraken['tree'] == taxon]['reads'].values[0]
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
    bacteria_size['Size (bp)'] = bacteria_size['Size (Mb)']*100000
    virus_db = db_read_in(os.path.join(db_path, "viruses.txt"))    
    virus_size = virus_db[['#Organism/Name', 'TaxID', 'Size (Kb)']].drop_duplicates(subset='#Organism/Name')
    virus_size['Size (bp)'] = virus_size['Size (Kb)']*1000
    fungi_protozoa_db = db_read_in(os.path.join(db_path, "eukaryotes.txt"))
    fp_size = fungi_protozoa_db[['#Organism/Name', 'TaxID', 'Size (Mb)']].drop_duplicates(subset='#Organism/Name')
    fp_size['Size (bp)'] = bacteria_size['Size (Mb)']*100000
    db_merge = pd.concat([bacteria_size, virus_size, fp_size], ignore_index=True)
    species_merge = pd.merge(species_df, db_merge, left_on='species', right_on='#Organism/Name', how='left')
    return species_merge

def gen_species_list(d_domain, total_reads, d_name):
    reorder_list = [
        'species',
        'TaxID',
        'Size (bp)',
        'reads',
        'taxon_RA',
        #'total_RA',
        'GPM'
    ]
    domain_filt = d_domain[d_domain['tree'].str.contains('s__')]
    str_regex = r"\|s__([^|]+)"
    domain_filt['species'] = domain_filt["tree"].str.extract(str_regex, expand=False)
    species_merge = genome_size(domain_filt)
    species_merge['GPM'] = total_reads/(species_merge['reads'] * species_merge['Size (bp)'])
    species_merge.drop_duplicates(subset='species', inplace=True)
    kraken_top10 = species_merge.sort_values(by="taxon_RA", ascending = False)[:10].reset_index(drop=True)
    kraken_top10 = kraken_top10[reorder_list]
    kraken_top10[['TaxID', 'Size (bp)']] = kraken_top10[['TaxID', 'Size (bp)']].fillna(0).astype(int)
    kraken_top10.to_csv(wd + "/" + d_name + ".txt", sep="\t", index=False)


def read_tsv(file):
    kraken_input = pd.read_csv(file, sep="\t", header=None)
    kraken_input.columns = ['tree', "reads"]
    kraken_sort = kraken_input.sort_values(by="reads", ascending = False)
    total_reads = sum(total_read_count(kraken_sort))

    d_bacteria = read_count(kraken_sort, "d__Bacteria") # extract only bacteria
    d_viruses = read_count(kraken_sort, "d__Viruses") # extract only viruses
    k_fungi = read_count(kraken_sort, "k__Fungi") # extract only fungi
    #k_protozoa = read_count(kraken_sort, "k__Protozoa") # extract only protozoa

    domain_list = [
        (d_bacteria, "d_bacteria"),
        (d_viruses, "d_viruses"),
        (k_fungi, "k_fungi"),
        #(k_protozoa, "k_protozoa"),
    ]
    for domain, dname in domain_list:
        gen_species_list(domain, total_reads, dname)


read_tsv(file)                            
