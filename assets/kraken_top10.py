import sys
import re
import pandas as pd

#file = sys.argv[1] 
# needs --use-mpa-style for kraken2
wd="/Users/wfon4473/Documents/Bioinformatics/all_testdirs/meta-gp_reports_tests/kraken"
input="CSF5703A.kraken2.nonhuman.nonrRNA.std.report"
file=wd+"/"+input
outfile=wd+"/species.txt"


def read_count(sorted_kraken, domain):
    filt_kraken = sorted_kraken[sorted_kraken['tree'].str.contains(domain)]
    read_count = filt_kraken[filt_kraken['tree'] == domain]['reads'].values[0]
    filt_kraken['perc'] = filt_kraken['reads']/read_count
    return filt_kraken

def read_tsv(file):
    kraken_input = pd.read_csv(file, sep="\t", header=None)
    kraken_input.columns = ['tree', "reads"]
    kraken_sort = kraken_input.sort_values(by="reads", ascending = False)
    d_bacteria = read_count(kraken_sort, "d__Bacteria") # extract only bacteria
    d_viruses = read_count(kraken_sort, "d__Viruses") # extract only virus
    #d_eukaryota = read_count(kraken_sort, "d__Eukaryota") # TO DO extract only fungal
    # TO DO extract only parasitic
    d_list = [d_bacteria, d_viruses]
    merged_df = pd.concat(d_list, ignore_index=True)
    filt_again = merged_df[merged_df['tree'].str.contains('s__')]
    re_sorted_kraken = filt_again.sort_values(by="perc", ascending = False)
    kraken_top10 = re_sorted_kraken[:10].reset_index()
    str_regex = r"\|s__([^|]+)"
    kraken_top10['species'] = kraken_top10["tree"].str.extract(str_regex, expand=False)
    kraken_top10['species'].to_csv(outfile, sep="\t", header=False, index=False)
    return kraken_top10

read_tsv(file)                            
