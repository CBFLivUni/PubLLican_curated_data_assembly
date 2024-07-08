import os
import pandas as pd
import json
import zipfile

# get data paths
data_folder_path = os.path.join("..", "data")
output_file_path = os.path.join("..", "output", "merged_gaf.json")

# store GAF dataframes in a list
gaf_dataframes = []

# find all .zip files ending with _curated_GOterms.zip
for root, dirs, files in os.walk(data_folder_path):
    for file in files:
        if file.endswith("_curated_GOterms.zip"):
            zip_path = os.path.join(root, file)
            # Open each zip file
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                # for each file witin
                for zip_info in zip_ref.infolist():
                    if zip_info.filename.endswith(".gaf"):
                        #TODO:
                        # # extract species from file name and add to column

                        # extract the GAF file to a temporary location
                        with zip_ref.open(zip_info) as gaf_file:
                            try:
                                # rad the GAF file and append it to the list of DFs
                                df = pd.read_csv(gaf_file, sep='\t', comment='!', header=None)
                                gaf_dataframes.append(df)
                            except Exception as e:
                                print(f"Failed to read {zip_info.filename} from {zip_path}: {e}")


# merge all dataframes into one if we
if gaf_dataframes:
    merged_df = pd.concat(gaf_dataframes, ignore_index=True)

    # dictionary to hold the final JSON structure
    json_data = {}

    # Iterate over each row in the dataframe and get the relevant data for:
    for index, row in merged_df.iterrows():
        # paper
        paper_id = row[4].replace("PMID:", "").strip()  # PUBMED ID is in the 5th column; remove PMID: from the start

        # species
        species_name = row[12] #13th column is the new column we created with species name extracted from file name
        species_taxon_id = row[9].replace("taxon:", "").strip()  # 10th column, remove taxon: from start; keep to make sure we havent got species name wrong

        # gene
        gene = row[2]  # gene shorthand notation is in row 3
        gene_id_veupath = row[1] # is in veupath db is in column 2
        gene_name = row[7] # column 8 is the long gene name

        # GO terms
        GOterm_id = row[3]  # GO term is in the 4th column
        evidence_code = row[5] # curators code for what sort of evidence there was for this annotation

        # Create structure if paper_id not in json_data
        #TODO:
        # update final json format
        if paper_id not in json_data:
            json_data[paper_id] = {"PMID": paper_id, "species": []}

        # Find or create the species entry
        species_entry = next((species for species in json_data[paper_id]["species"] if species["name"] == species_name), None)
        if not species_entry:
            species_entry = {"name": species_name, "genes": []}
            json_data[paper_id]["species"].append(species_entry)

        # Find or create the gene entry
        gene_entry = next((gene for gene in species_entry["genes"] if gene["name"] == gene_name), None)
        if not gene_entry:
            gene_entry = {"name": gene_name, "GO_terms": []}
            species_entry["genes"].append(gene_entry)

        # Add the GO term to the gene entry if not already present
        if go_term not in gene_entry["GO_terms"]:
            gene_entry["GO_terms"].append(go_term)

    # Convert the dictionary to a list
    json_output = list(json_data.values())

    # Write the dictionary to a JSON file
    with open(output_file_path, 'w') as json_file:
        json.dump(json_output, json_file, indent=4)
    print(f"Merged .gaf files written to {output_file_path}")
else:
    print("No .gaf files found to merge.")