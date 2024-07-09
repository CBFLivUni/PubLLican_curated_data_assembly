import os
import pandas as pd


# from the set of curated data, get a list of all unique PMIDs that contributed to the curation
# this list will be used to run Noah's pipeline using different LLMs to scorre their performance
# and to test how fine-tuning impacts this (how similar are the GO terms annotations of the LLMs)

PMIDs_data_dir = os.path.join("..", "output")
merged_gaf_file = os.path.join(PMIDs_data_dir, "merged_gaf.tsv")

# read in merged_gaf.tsv column with pubmed IDs as a list
df = pd.read_csv(merged_gaf_file, sep='\t', comment='!', header=None)

# Assuming the PubMed IDs are in the 6th column (index 5)
pubmed_ids = df[5].astype(str)  # Convert to string to handle potential float values

# Remove 'PMID:' prefix and strip whitespace
pubmed_ids = pubmed_ids.str.replace("PMID:", "").str.strip()

# Get unique PMIDs
unique_pmids = pubmed_ids.unique()

# Convert to list
unique_pmids_list = unique_pmids.tolist()

# Print the list of unique PMIDs
print(unique_pmids_list)
