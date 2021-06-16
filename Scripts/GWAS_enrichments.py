#!/usr/bin/env python
# coding: utf-8
from collections import defaultdict, Counter
import json
import numpy as np
import pandas as pd
import os
from scipy import stats
from time import time
import timeit
import argparse
import sys

# # parsing experimental factor ontology for child traits

# method to verify if a given EFO trait/term is obsolete
def is_obsolete(lines,index):
    replaced_ids = []
    last_index = len(lines)
    while(True):
        if lines[index].startswith("[Term]") or index==last_index or lines[index].startswith("[Typedef]"):
            return False
        if lines[index].startswith("is_obsolete"):
            return True
        index = index + 1

def get_is_a(line):
    #is_a: GO:0000732 ! strand displacement
    # is_a: EFO:0003892 ! pulmonary function measurement
    # splits = line.split(":")
    parent = line.split(" ! ")[1]
    return parent.strip()

# this function returns list of alternate ID's for a specific term
def get_alt_id(lines,index):
    alt_ids = []
    while(True):
        if lines[index].startswith("alt_id"):
            splits = lines[index].strip().split(":")
            if len(splits)>2:
                term_id = splits[1] + ":" + splits[2]
            else:
                term_id = splits[1]
            alt_ids.append(term_id.strip())
            index = index+1
        else:
            break
    return alt_ids

# this function parses all the parent terms of a given term
def getParents(lines,index):
    parents = []
    last_index = len(lines)
    while(True): #looping till next term starts
        if index==last_index or lines[index].startswith("[Term]") or lines[index].startswith("[Typedef]"):
            break
        else:
            if lines[index].startswith("is_a"):
                parent = get_is_a(lines[index])
                if parent is not None:
                    parents.append(parent)
                index = index + 1
            else:
                index = index + 1
    return parents[::-1] # returning parents in the descending order of their distance from root

def get_name(line):
    #id: EFO:0000001
    splits = line.split(":")
    term_name = splits[1]+":"+splits[2]
    return term_name.strip()

# this function performs the fisher's exact test for overlap between two sets of genes
def fisher_test(geneset1, geneset2):
    overlap = set.intersection(set(geneset1), set(geneset2))
    if len(overlap) ==0:
        return set(), 1.0
    a = len(overlap)
    b = len(geneset1) - len(overlap)
    c = len(geneset2) - len(overlap)
    d = 20000 - (len(geneset1) + len(geneset2) - len(overlap))
    result = stats.fisher_exact([[a,b],[c,d]])
    fet = '{:0.3e}'.format(result[1])
    return overlap, float(fet)

# computes adjusted p-values
def calc_adjp(pvals):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    pvals = np.asfarray(pvals)
    descend_p = pvals.argsort()[::-1]
    orig = descend_p.argsort()
    steps = float(len(pvals)) / np.arange(len(pvals), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * pvals[descend_p]))
    return q[orig]

def SplitData(df,target_columns,separator):
    new_rows = []
    df.apply(splitListToRows,axis=1,args=(new_rows,target_columns,separator))
    new_df = pd.DataFrame(new_rows)
    return new_df

#method to split columns having multiple values to multiple rows
def splitListToRows(row,row_accumulator,target_columns,separator):
    split_rows = []
    for target_column in target_columns:
        split_rows.append(row[target_column].split(separator))
    #separate for multiple columns
    for i in range(len(split_rows[0])):#looping over each split value
        new_row = row.to_dict()
        for j in range(len(split_rows)):#looping over each column
            #print(split_rows[j][i].strip())
            new_row[target_columns[j]] = split_rows[j][i].strip()
        row_accumulator.append(new_row)

def enrichment_handler(cid):
    print("*" * 5, cid, "*" * 5)
    cluster_genes = node_stats[node_stats['cluster'] == cid]['gene'].to_list()
    enc_results_cluster = trait_enrichments(cluster_genes, cid)
    return enc_results_cluster

# Enrichments of trait-associated genes in a given DEG
def trait_enrichments(genes, geneset_name):
    enc_results = []

    for trait in all_traits:
        trait_list = [trait]
        if trait in child_dict:
            trait_list = trait_list + child_dict[trait]
        trait_genes = set(GWAS_df[GWAS_df['MAPPED_TRAIT'].map(lambda trait:
                                                                  trait in trait_list)]['MAPPED_GENE'].to_list())
        if len(trait_genes) == 0:
            continue
        overlap, fet = fisher_test(genes, trait_genes)
        if len(overlap) == 0:
            continue
        logp = -np.log10(fet)
        enc_results.append([trait, len(child_dict[trait]), len(trait_genes), geneset_name,
                                len(genes), len(overlap), ",".join(list(overlap)), fet, logp])
    enc_results = pd.DataFrame(enc_results, columns=["Trait", "# Child Traits", "# Mapped genes",
                                                         "DEG", "DEG_Size", "Overlap_Size", "Overlap",
                                                         "P-value(FET)", "logP"])
    enc_results['q-val'] = calc_adjp(enc_results['P-value(FET)'].to_list())
    enc_results['logq'] = -np.log10(enc_results['q-val'])
    return enc_results

# '/path/to/git/repository'
parser = argparse.ArgumentParser(description="command line arguments")
parser.add_argument('--obo_file', type=str,
                    help="Path to the EFO ontology OBO file (https://www.ebi.ac.uk/efo/efo.obo)",
                    required=True)
parser.add_argument('--assoc_file', type=str,
                    help="Path to the file containing GWAS Catalog associations (https://www.ebi.ac.uk/gwas/docs/file-downloads)",
                    required=True)
parser.add_argument('--cluster_file', type=str,
                    help="Path to tab-delimited file (two columns) containing genes and their module memberships",
                    required=True)
parser.add_argument('--remove_intergenic', action='store_true',
                    help="Path to tab-delimited file (two columns) containing genes and their module memberships",
                    required=False)
parser.add_argument('--outpath', type=str,
                    help="Output path where the result file needs to be stored",
                    required=False)
args = parser.parse_args()

start_time = timeit.default_timer()

if not os.path.exists(args.obo_file):
    print("***** The EFO OBO file is not available. Please provide the correct path. You can download from the below URL:")
    print("https://www.ebi.ac.uk/efo/efo.obo")
    sys.exit(1)
efo_lines = open(args.obo_file, "r", encoding="utf-8").readlines()
print("Number of lines the OBO file:", len(efo_lines))

# parsing EFO OBO file
print("*"*5, "Parsing EFO OBO file", "*"*5)
efo = {}
parent_dict = defaultdict(list)
child_dict = defaultdict(list)
i=0
while i<len(efo_lines):
    if efo_lines[i].startswith("[Term]") and not is_obsolete(efo_lines,i+1) and 'EFO' in efo_lines[i+1]:
        term_id = get_name(efo_lines[i+1])
        term_name = efo_lines[i+2].split(":")[1].strip()
        efo[term_id] = term_name
        parents = getParents(efo_lines,i+1)

        #updating the child dict of the parents
        for parent in parents:
            child_dict[parent].append(term_name)
            parent_dict[term_name].append(parent)
    i=i+1
print(len(efo), len(parent_dict), len(child_dict))

# reading GWAS Catalog data
if not os.path.exists(args.assoc_file):
    print("***** The GWAS Catalog association file does not exist. Please provide the correct path. You can download it from the below URL:")
    print("https://www.ebi.ac.uk/gwas/docs/file-downloads")
    sys.exit(1)

GWAS_df = pd.read_csv(args.assoc_file,sep = "\t")
print(GWAS_df.shape)

# removing unwanted columns
GWAS_df.drop(['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK', 'STUDY',
                 'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'GENOTYPING TECHNOLOGY', 'STUDY ACCESSION',
                 'PLATFORM [SNPS PASSING QC]'], axis = 1, inplace=True)
GWAS_df.drop_duplicates(inplace=True)
print(GWAS_df.shape)

GWAS_df['MAPPED_TRAIT'] = GWAS_df['MAPPED_TRAIT'].astype(str)
GWAS_df['MAPPED_GENE'] = GWAS_df['MAPPED_GENE'].astype(str)

GWAS_df = SplitData(GWAS_df,["MAPPED_GENE"],", ") #splitting multiple mapped genes
GWAS_df = SplitData(GWAS_df,['MAPPED_TRAIT'],', ')#splitting multiple mapped traits

# removing intergenic associations
if args.remove_intergenic:
    print("***** Removing intergenic associations")
    GWAS_df = GWAS_df[GWAS_df['INTERGENIC']==0]
    print(GWAS_df.shape)

all_traits = GWAS_df['MAPPED_TRAIT'].unique()
print(len(all_traits))

all_child_traits = []
for trait in all_traits:
    all_child_traits.extend(child_dict[trait])
all_child_traits = set(all_child_traits)
print(len(all_child_traits))

# # GWAS Catalog trait enrichments in SARS-CoV-2 gene modules
print("*"*5, "Computing trait enrichments in SARS-CoV-2 modules", "*"*5)
try:
    import multiprocessing as mp
except ImportError as e:
    print("Multiprocessing package doesn't exist. Please install it using the below command:")
    print(" pip install multiprocessing")
    sys.exit(1)

print("Number of processors:",mp.cpu_count())

node_stats = pd.read_csv(args.cluster_file, sep = "\t")
print(node_stats.head(2))

# removing unassigned genes
node_stats = node_stats[-node_stats['cluster'].isna()]
# candidate clusters
selected_clusters = list({cid:size for cid,size in dict(Counter(node_stats['cluster'])).items() if size > 4}.keys())
print("Number of candidate clusters:", len(selected_clusters))

#initiating multiprocessing pool object
pool = mp.Pool(mp.cpu_count())
results = pool.map_async(enrichment_handler, [cid for cid in selected_clusters]).get()
pool.close()

enc_results = pd.concat(results, axis = 0)
print(enc_results.shape)
end_time = timeit.default_timer()
print("Time elapsed (mins):", round((end_time - start_time)/60, 4))

outfile = os.path.join(args.outpath, "module_GWAS_enrichments.txt")
enc_results.to_csv(outfile, sep="\t")