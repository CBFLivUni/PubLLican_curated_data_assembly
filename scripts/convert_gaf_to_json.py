import os
import pandas as pd
import json
import zipfile
import re
import tempfile
import gzip
import shutil

# get data paths
data_folder_path = os.path.join("..", "data", "PFalciparum_3d7_test")
output_json_file_path = os.path.join("..", "output", "PFalciparum_3d7_gaf.json")
output_tsv_file_path = os.path.join("..", "output", "PFalciparum_3d7_gaf.tsv")

# store GAF dataframes in a list
gaf_dataframes = []

# function to convert the gaf files to df and add a column for species and strain
def process_gaf_file(gaf_file, species, strain, filename):
    try:
        df = pd.read_csv(gaf_file, sep='\t', comment='!', header=None, compression='infer')
        # Add species and strain columns
        df['Species'] = species
        df['Strain'] = strain
        gaf_dataframes.append(df)
    except Exception as e:
        print(f"Failed to read GAF file {filename}: {e}")


# function to get the species and strain from the file name
def extract_species_and_strain(filename_base):
    match = re.search(r'\d{2}_([^_]+?)_', filename_base)

    if match:
        species_strain = match.group(1)
        capital_indices = [i for i, c in enumerate(species_strain) if c.isupper()]
        if len(capital_indices) >= 2:
            species = f"{species_strain[0]}.{species_strain[1:capital_indices[1]].lower()}"
            strain = species_strain[capital_indices[1]:]
        else:
            species = species_strain
            strain = ""
    else:
        species = ""
        strain = ""
    return species, strain


def process_zip_file(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for zip_info in zip_ref.infolist():
            filename_base = os.path.basename(zip_info.filename)
            species, strain = extract_species_and_strain(filename_base)

            if zip_info.filename.endswith(".zip"):
                # Handle nested zip files
                with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                    tmp_file.write(zip_ref.read(zip_info.filename))
                    tmp_file_path = tmp_file.name
                process_zip_file(tmp_file_path)
                os.remove(tmp_file_path)
            elif zip_info.filename.endswith(".gaf"):
                # Directly process .gaf files
                with zip_ref.open(zip_info) as gaf_file:
                    process_gaf_file(gaf_file, species, strain, zip_info.filename)
            elif zip_info.filename.endswith(".gz"):
                # Temporarily extract .gz files
                with zip_ref.open(zip_info) as extracted_file:
                    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                        shutil.copyfileobj(extracted_file, tmp_file)
                        tmp_file_path = tmp_file.name

                    # Process the extracted .gz file
                    with gzip.open(tmp_file_path, 'rt') as gaf_file:
                        process_gaf_file(gaf_file, species, strain, zip_info.filename)

                    os.remove(tmp_file_path)


# find all .zip files ending with _curated_GOterms.zip
for root, dirs, files in os.walk(data_folder_path):
    for file in files:
        if file.endswith("_curated_GOterms.zip"):
            zip_path = os.path.join(root, file)
            process_zip_file(zip_path)

print(gaf_dataframes)
# Check the structure of one of the DataFrames for debugging purposes
if gaf_dataframes:
    print("Column names and first few rows of one DataFrame for debugging:")
    print(gaf_dataframes[0].head())
    print(gaf_dataframes[0].columns)
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
        paper_id = str(row[5]).replace("PMID:", "").strip()  # PUBMED ID is in the 6th column; remove PMID: from the start

        # species
        species_name = row['Species']  # 18h column is the new column we created with species name extracted from file name
        species_taxon_id = str(row[12]).replace("taxon:", "").strip()  # 13th column, remove taxon: from start; keep to make sure we havent got species name wrong
        strain = row['Strain']
        # gene
        gene_shorthand = row[2]  # gene shorthand notation is in row 3
        gene_veupath_id = row[1]  # is in veupath db is in column 2
        gene_full_name = row[9]  # column 10 is the long gene name

        # GO terms
        go_term_id = row[4]  # GO term is in the 4th column
        evidence_code = row[6]  # curators code for what sort of evidence there was for this annotation

        # Create structure if paper_id not in json_data
        if paper_id not in json_data:
            json_data[paper_id] = {"PMID": paper_id, "species": []}

        # Find or create the species entry
        species_entry = next((species for species in json_data[paper_id]["species"] if species["name"] == species_name),
                             None)
        if not species_entry:
            species_entry = {"name": species_name, "strain": strain, "taxon_id": species_taxon_id, "genes": []}
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
