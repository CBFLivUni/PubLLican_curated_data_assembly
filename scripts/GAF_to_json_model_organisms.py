import os
import pandas as pd
import json
import zipfile
import re
import tempfile
import gzip
import shutil

# get data paths
data_folder_path = os.path.join("..", "data", "GO_consortium_model_species")
output_json_file_path = os.path.join("..", "output", "model_species_merged_gaf.json")
output_tsv_file_path = os.path.join("..", "output", "model_species_merged_gaf.tsv")

# store GAF dataframes in a list
gaf_dataframes = []

# find all .gz files and read as data frame, add it to the ist of dfs
for root, dirs, files in os.walk(data_folder_path):
    for file in files:
        print(file)
        if file.endswith(".gaf.gz"):
            zip_path = os.path.join(root, file)
            with gzip.open(zip_path, 'rt') as gaf_file:
                try:
                    df = pd.read_csv(gaf_file, sep='\t', low_memory=False, comment='!', header=None,
                                     compression='infer')

                    gaf_dataframes.append(df)
                except Exception as e:
                    print(f"Failed to read GAF file {file}: {e}")



# merge all dataframes into one
if gaf_dataframes:
    merged_df = pd.concat(gaf_dataframes, ignore_index=True)

    # Save the DF as a TSV file
    merged_df.to_csv(output_tsv_file_path, sep='\t', index=False)
    print(f"TSV file written to {output_tsv_file_path}")


    # dictionary to hold the final JSON structure
    json_data = {}

    # Iterate over each row in the dataframe and get the relevant data for:
    for index, row in merged_df.iterrows():
        # paper
        paper_id = str(row[5])# PUBMED ID is in the 6th column

        species_taxon_id = str(row[12]) # 13th column, remove taxon

        # gene
        gene_shorthand = row[2]  # gene shorthand notation is in col 3
        gene_veupath_id = row[1]  # is in veupath db is in column 2
        gene_full_name = row[9]  # column 10 is the long gene name

        # GO terms
        go_term_id = row[4]  # GO term is in the 4th column
        evidence_code = row[6]  # curators code for what sort of evidence there was for this annotation

        # Create structure if paper_id not in json_data
        if paper_id not in json_data:
            json_data[paper_id] = {"PMID": paper_id, "species": []}

        # Find or create the species entry
        species_entry = next((species for species in json_data[paper_id]["species"] if species["taxon_id"] == species_taxon_id),
                             None)
        if not species_entry:
            species_entry = {"taxon_id": species_taxon_id, "genes": []}
            json_data[paper_id]["species"].append(species_entry)

        # Find or create the gene entry
        gene_entry = next((gene for gene in species_entry["genes"] if gene["shorthand_name"] == gene_shorthand), None)
        if not gene_entry:
            gene_entry = {
                "VEuPAthDB_ID": gene_veupath_id,
                "shorthand_name": gene_shorthand,
                "full_name": gene_full_name,
                "GO_terms": []
            }
            species_entry["genes"].append(gene_entry)

        # Add the GO term to the gene entry if not already present
        go_term_entry = {"GO_ID": go_term_id, "evidence_code": evidence_code}
        if go_term_entry not in gene_entry["GO_terms"]:
            gene_entry["GO_terms"].append(go_term_entry)

    # Convert the dictionary to a list
    json_output = list(json_data.values())

    # Write the dictionary to a JSON file
    with open(output_json_file_path, 'w', encoding="UTF-8") as json_file:
        json.dump(json_output, json_file, indent=4)
    print(f"Merged .gaf files written to {output_json_file_path}")
else:
    print("No .gaf files found to merge.")


#
# print(gaf_dataframes)
# # Check the structure of one of the DataFrames for debugging purposes
# if gaf_dataframes:
#     print("Column names and first few rows of one DataFrame for debugging:")
#     print(gaf_dataframes[0].head())
#     print(gaf_dataframes[0].columns)
