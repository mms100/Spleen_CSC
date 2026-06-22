#!/usr/bin/env python
# coding: utf-8

# In[1]:


#load packages
import os
import scanpy as sc
import decoupler as dc
import argparse

# Only needed for processing
import numpy as np
import pandas as pd
import scipy
# Needed for some plotting
import matplotlib.pyplot as plt
import seaborn.objects as so
from scipy.sparse import csc_matrix, csr_matrix
from pathlib import Path



# In[2]:


# Create the ArgumentParser object
parser = argparse.ArgumentParser(description='Process some CSV files.')

# Define the expected command-line arguments
parser.add_argument('--input_dir_1', type=str, required=True, help='Path to the input directory')
parser.add_argument('--outdir_3', type=str, required=True, help='Path to the output directory')
parser.add_argument('--cond1', type=str, required=True, help='Path to the output directory')
parser.add_argument('--cond2', type=str, required=True, help='Path to the output directory')
parser.add_argument('--species', type=str, choices=['human', 'mouse'], required=True,
                    help='Species to use: "human" or "mouse"')
# Parse the arguments
args = parser.parse_args()

# Assign the parsed arguments to variables
input_dir = args.input_dir_1
output_dir = args.outdir_3
cond1 = args.cond1
cond2 = args.cond2
species = args.species.lower()

# In[4]:


# Initialize the directory path, CSV files, and names
directory_path = input_dir
csv_files = os.listdir(directory_path)  # Assumes all files in the directory are CSV files
#names_list = ['_'.join(file.split('.')[0].split('_')[1:]) for file in csv_files]


prefix = f"{args.cond1}vs{args.cond2}_"
names_list = [file.split('.')[0].replace(prefix, '') for file in csv_files]

print(names_list)
# Use a list comprehension to read CSV files into DataFrames
filtered_list = [pd.read_csv(os.path.join(directory_path, file), encoding='latin-1') [["GeneID", "t_stat"]] for file in csv_files]

# Create a dictionary to associate names with DataFrames 
named_dataframes = dict(zip(names_list, filtered_list))
named_dataframes


# In[7]:


#organised_dataframe = [named_dataframes[name].rename(columns={ "GeneID": "genesymbols", "t_stat": name}).set_index("genesymbols").T for name in names_list]
organised_dataframe = [
    named_dataframes[name]
    .replace([np.inf, -np.inf], np.nan)  # Replace Inf/-Inf with NaN
    .dropna(subset=["t_stat"])  # Drop rows where t_stat is NaN
    .rename(columns={"GeneID": "genesymbols", "t_stat": name})  # Rename columns
    .set_index("genesymbols")
    .T
    for name in names_list]


# In[5]:
# Import databases based on species
if species == 'human':
    print("Loading MSigDB human resource...")
    msigdb = pd.read_csv('/path/to/databases/MSigDB_human_database.csv')
    msigdb = msigdb[['genesymbol', 'collection', 'geneset']]
else:
    print("Loading mouse MSigDB database...")
    msigdb = pd.read_csv('/path/to/databases/Msg_mouse_database.csv')
    msigdb = msigdb[['genesymbol', 'collection', 'geneset']]

# In[6]:


msigdb["collection"].unique()


# In[7]:


# Assume msigdb is your existing DataFrame

collections_to_process = ['go_biological_process', 'reactome_pathways', "kegg_pathways", 'hallmark']


# Initialize an empty dictionary to store the results
msigdb_results = {}

for collection in collections_to_process:
    # Filter the msigdb DataFrame based on the current collection
    msigdb_collection = msigdb[msigdb['collection'] == collection]
    # Remove duplicated entries
    msigdb_collection = msigdb_collection[~msigdb_collection.duplicated(['geneset', 'genesymbol'])]
    #rename
    msigdb_collection.loc[:, 'geneset'] = ['_'.join(name.split('_')[1:]) for name in msigdb_collection['geneset']]

    geneset_size = msigdb_collection.groupby("geneset").size()
    gsea_genesets = geneset_size.index[(geneset_size > 15) & (geneset_size < 500)]

    msigdb_collection = msigdb_collection[msigdb_collection["geneset"].isin(gsea_genesets)]

    # Add a column of weights =1
    msigdb_collection["weight"] = 1

    msigdb_collection = msigdb_collection.rename(columns={"genesymbol": "target", "geneset": "source"})

    # Store the result in the dictionary
    msigdb_results[f"{collection}"] = msigdb_collection



# In[8]:


#apply decouplR
# Initialize dictionaries for all collections specified above
pathway_acts = {}
pathway_pvals = {}

for collection in collections_to_process:
    # Initialize dictionaries for each collection
    acts_dict = {}
    pvals_dict = {}

    for i, name in enumerate(names_list):
        acts_dict[name], pvals_dict[name] = dc.mt.ulm(data=organised_dataframe[i], net=msigdb_results[collection])

        
    pathway_acts[collection] = acts_dict
    pathway_pvals[collection] = pvals_dict


# In[9]:


#correct pvals
#for name in names_list:
    #pathway_pvals['reactome_pathways'][name].loc[name ,:] = scipy.stats.false_discovery_control(pathway_pvals['reactome_pathways'][name].loc[name ,:])



# In[10]:


#additaional step to make the dictionary keys fits the output directory names
key_mapping = {
    'hallmark': 'Hallmark',
    'go_biological_process': 'GO',
    'reactome_pathways': 'Reactome',
    'kegg_pathways': 'KEGG'
}

# Update the dictionary keys based on the mapping
for old_key, new_key in key_mapping.items():
    if old_key in pathway_pvals and old_key in pathway_acts:
        pathway_pvals[new_key] = pathway_pvals[old_key]
        pathway_acts[new_key] = pathway_acts[old_key]
        del pathway_pvals[old_key]
        del pathway_acts[old_key]

print(pathway_pvals.keys())
print(pathway_acts.keys())


# In[11]:


collections_to_process = ['GO', 'Reactome', "KEGG", 'Hallmark']



for collection in collections_to_process:
    # Create the output directory if it doesn't exist
    output_directory = f'{output_dir}/{collection}/'
    os.makedirs(output_directory, exist_ok=True)

    # Concatenate pathway_acts along columns and save to 'acts.csv'
    acts_df = pd.concat([df.T.fillna(0) for df in pathway_acts[collection].values()], axis=1)
    acts_df.to_csv(os.path.join(output_directory, 'acts.csv'), index=True)

    # Concatenate pathway_pvals along columns and save to 'pvals.csv'
    pvals_df = pd.concat([df.T.fillna(0) for df in pathway_pvals[collection].values()], axis=1)
    pvals_df.to_csv(os.path.join(output_directory, 'pvals.csv'), index=True)


print("end_of_script")
