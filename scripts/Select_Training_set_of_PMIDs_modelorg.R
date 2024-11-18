########## script aims ############
## selecting species and pmids to use for a train set
# we want ~ 500 papers in the set

# we want to have species that are well annotated and we have gene ID conversion
# tools for as a sort of positive controls/easy task to annotate

# we probably want to use pmids that are contributing a larger number of genes
# and/or GO terms to the annotations so that we get exposed to a wide array of 
# gene annotations, but also papers that contribute small and medium numbers as
# we want the model to be trained on a variety of methods, not just high 
# throughput studies
###################################

######## libs and read in #########
library(tidyverse)

library(rentrez) # for converting taxon id to species name
library(httr) # for getting 
library(jsonlite)
source("getting_pmids_from_other_id_types_functions.R")
get_species_name <- function(taxon_id) {
  record <- entrez_summary(db = "taxonomy", id = taxon_id)
  return(record$scientificname)
}

GAF_model_organisms <- read.delim("../output/model_species_merged_gaf.tsv", sep = "\t", header = TRUE)
GAF_headers <- read.csv("../data/GAF_headers.csv")
Evidence_Codes <- read.csv("../data/EvidenceCodesSelection.csv")
evidence_codes_for_training_set <- Evidence_Codes$Evidence_Code[Evidence_Codes$Include_Exclude == "include"]
##### tidying the merged GAF file for model organisms ####
# Step 1 - add meaningful headers from the GO consortium GAF headers description
head(GAF_model_organisms)
names(GAF_model_organisms) <- GAF_headers$clean_header

# Step 2 - filter data based on evidecne codes
GAF_model_organisms_filt <- GAF_model_organisms[GAF_model_organisms$Evidence_Code %in% evidence_codes_for_training_set, ]
# these are ~ 50% of entries

# Step 3 - convert references to PMIDs if possible
head(GAF_model_organisms_filt)
# get anything before the :, and extract unique list of options from other refs 
# to view all possibilities

# Split by "|" first, then split each component by ":", and extract unique prefixes

# Split by "|" to handle multiple references
refs_list <- GAF_model_organisms_filt$DB_Reference %>% 
  as.character() %>%
  strsplit("\\|") 
  
# if any of them start with PMID:, add them to a new column called PMIDs
# initialise the column
GAF_model_organisms_filt$PMIDs <- NA

# Process each reference in the list of refs
refs_list_processed <- lapply(refs_list, function(refs) {
  ref_prefixes <- lapply(refs, function(ref) {
    ref_parts <- strsplit(ref, ":")[[1]]
    
    # Check if reference starts with "PMID" and collect it if true
    if (ref_parts[1] == "PMID") {
      ref_id <- ref_parts[2]
      list(PMID = ref_id, prefix = NULL)
    } else {
      list(PMID = NULL, prefix = ref_parts[1])
    }
  })
  
  # Extract PMIDs (if any)
  PMIDs <- sapply(ref_prefixes, function(x) x$PMID)
  PMIDs <- PMIDs[!is.na(PMIDs)] # Filter out NAs
  
  # Extract prefixes
  prefixes <- unlist(sapply(ref_prefixes, function(x) x$prefix))
  prefixes <- prefixes[!is.na(prefixes)] # Filter out NAs
  
  # Return both PMIDs and prefixes
  list(PMID = if (length(PMIDs) > 0) paste(PMIDs, collapse = "|") else NA,
       prefixes = prefixes)
})


# Update dataframe with PMIDs and get unique reference prefixes for non-PMID references
GAF_model_organisms_filt$PMIDs <- sapply(refs_list_processed, function(x) x$PMID)

unique_prefixes <- unlist(lapply(refs_list_processed, function(x) x$prefixes)) %>%
  unique()

# > unique_prefixes
# "TAIR"     "GO_REF"   "PMID"     "Reactome" "WB_REF"   "CGD_REF"  "ZFIN"     
# "FB"       "DOI"      "MGI"      "GO"       "RGD"      "REACTOME" "SGN_ref"

# we want to convert some of the other ID types present to PMIDs where no PMID is 
# already available.

# first lets split the data into two subsets - with PMIDs and without. 


GAF_model_organisms_withpmid <- GAF_model_organisms_filt[!grepl("NULL", GAF_model_organisms_filt$PMIDs), ]
GAF_model_organisms_nopmid <- GAF_model_organisms_filt[grepl("NULL", GAF_model_organisms_filt$PMIDs), ]


# Extract the reference type for those without (prefix before ":")
GAF_model_organisms_nopmid <- GAF_model_organisms_nopmid %>%
  mutate(ref_type = sapply(strsplit(as.character(DB_Reference), ":"), `[`, 1))

# Calculate total and unique counts per reference type that don't have PMID
reference_summary <- GAF_model_organisms_nopmid %>%
  group_by(ref_type) %>%
  summarise(
    total_count = n(),
    unique_count = n_distinct(DB_Reference)
  )
# 
# > reference_summary
# # A tibble: 13 × 3
# ref_type total_count unique_count
# <chr>          <int>        <int>
#   1 CGD_REF            6            3
# 2 DOI               21            6
# 3 FB              8855           39
# 4 GO                 8            3
# 5 GO_REF        916927           26
# 6 MGI              157           70
# 7 REACTOME      100045        15117
# 8 RGD           202450          441
# 9 Reactome       92664        14225
# 10 SGN_ref           26            1
# 11 TAIR           42391          323
# 12 WB_REF             6            1
# 13 ZFIN            7632           21

# we should be able to convert DOI to PMID easily, for the rest it is a bit more 
# complicated, we can attempt a few and see how it goes. 

DOIs <- GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "DOI"] 
DOIs <- gsub("DOI:", "", DOIs) %>% unique()

PMIDs_from_DOI <- sapply(DOIs, get_pmid_from_doi)



reactome_IDs <- GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "REACTOME" |
                                                        GAF_model_organisms_nopmid$ref_type == "Reactome" ]
reactome_IDs <- gsub("Reactome:", "", reactome_IDs) 
reactome_IDs <- gsub("REACTOME:", "", reactome_IDs) %>% unique()

test_reactome_IDs <- head(reactome_IDs)
PMIDs_from_reactome <- sapply(test_reactome_IDs, get_pmids_from_reactome)
# this does not work, leave for now. 


MGIs <- GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "MGI"]
MGIs <- sub("MGI:", "", MGIs) %>% unique()
# MGI website actually expects the MGI: in the ref number, but we have 2 so we sub the first one out
PMIDs_from_MGI <- sapply(MGIs, get_pmids_from_db, db = "MGI")
# did not get any.

ZFINs <-  GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "ZFIN"]
ZFINs <- gsub("ZFIN:", "", ZFINs) %>% unique()
# test_ZFINs <- head(ZFINs)
PMIDs_from_ZFIN <- sapply(ZFINs, get_pmids_from_db, db = "ZFIN")
# worked for some! 

RGDs <-  GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "RGD"]
RGDs <- gsub("RGD:", "", RGDs) %>% 
  gsub("\\|.*", "", .) %>%   # Remove everything after and including "|"
  unique()
# strip anything after and including |
RGD_test <- head(RGDs)
PMIDs_from_RGD <- sapply(RGDs, get_pmids_from_db, db = "RGD")
# some work :) 


FBs <-  GAF_model_organisms_nopmid$DB_Reference[GAF_model_organisms_nopmid$ref_type == "FB"]
FBs <- RGDs <- sub("FB:", "", FBs) %>% 
  gsub("\\|.*", "", .) %>%   # Remove everything after and including "|"
  unique()

PMIDs_from_FB <- sapply(FBs, get_pmids_from_db, db = "FB")
# doesn't work atm. 

## it's proving difficult to  extract PMIDs from some of these. for the training set,
# I have decided to take forward only ones that were easily converted. 

# Step4: Filter the training set to only include papers that have been successfully converted to PMIDs.
# first re-include the new PMIDs from RGD and ZFIN in GAF_model_organisms_nopmid$PMID

# Convert the results of sapply to named vectors for easy matching
RGD_pmids <- PMIDs_from_RGD[!is.na(PMIDs_from_RGD)]
ZFIN_pmids <- PMIDs_from_ZFIN[!is.na(PMIDs_from_ZFIN)]
DOI_pmids <- PMIDs_from_DOI[!is.na(PMIDs_from_DOI)]
# add prefixes back in
names(RGD_pmids) <- paste0("RGD:", names(RGD_pmids))
names(ZFIN_pmids) <- paste0("ZFIN:", names(ZFIN_pmids))
names(DOI_pmids) <- paste0("DOI:", names(DOI_pmids))

# Update PMIDs for RGD and ZFIN references
# Step 2: Create a combined lookup list for RGD and ZFIN with prefixes
pmid_lookup <- c(RGD_pmids, ZFIN_pmids, DOI_pmids)

# Function to find and replace PMIDs
get_pmid_from_reference <- function(ref, pmid_lookup) {
  # Split each reference entry by "|"
  refs <- unlist(strsplit(ref, "\\|"))
  # Find matching PMID for any prefix in the lookup table
  matched_pmids <- pmid_lookup[refs]
  # Remove NAs and concatenate results if any matches are found
  pmid <- if (any(!is.na(matched_pmids))) {
    paste(unique(na.omit(matched_pmids)), collapse = "|")
  } else {
    NA
  }
  return(pmid)
}

# Apply the function to update PMIDs in one pass
GAF_model_organisms_nopmid$PMIDs <- mapply(
  function(ref, ref_type) {
    if (ref_type %in% c("RGD", "ZFIN")) {
      get_pmid_from_reference(ref, pmid_lookup)
    } else {
      NA  # Leave NA for other reference types
    }
  },
  GAF_model_organisms_nopmid$DB_Reference,
  GAF_model_organisms_nopmid$ref_type
)

# retain only ones with pmid now and get rid of last column
GAF_model_organisms_newpmids <- GAF_model_organisms_nopmid[!is.na(GAF_model_organisms_nopmid$PMIDs), 1:18]
# Combine the datasets with and without PMIDs
GAF_model_organisms_final <- bind_rows(GAF_model_organisms_withpmid, GAF_model_organisms_newpmids)


# Step 5: convert taxon IDs to species.

# convert taxon notation ti just IDs
head(GAF_model_organisms_final$taxon)
GAF_model_organisms_final$taxon_id <- gsub("taxon:", "", GAF_model_organisms_final$taxon)

# some have 2 ids, generate column for primary and secondary
GAF_model_organisms_final <- GAF_model_organisms_final %>%
  mutate(primary_taxon_id = sapply(strsplit(taxon_id, "\\|"), `[`, 1),
         secondary_taxon_id = sapply(strsplit(taxon_id, "\\|"), `[`, 2))

# Get unique taxon IDs for lookup to avoid duplicate API calls
unique_primary_taxon_ids <- unique(GAF_model_organisms_final$primary_taxon_id)
unique_secondary_taxon_ids <- unique(GAF_model_organisms_final$secondary_taxon_id[!is.na(GAF_model_organisms_final$secondary_taxon_id)])


# Create lookup tables for species names
primary_species_lookup <- setNames(sapply(unique_primary_taxon_ids, get_species_name), unique_primary_taxon_ids)

# safe function to go over warnign producing downstream errors
safe_get_species_name <- function(taxon_id) {
  tryCatch(
    get_species_name(taxon_id),
    error = function(e) return("")
  )
}
secondary_species_lookup <- setNames(sapply(unique_secondary_taxon_ids, safe_get_species_name), unique_secondary_taxon_ids)
# one null hit as warning produced: ID 489466 produced error 'cannot get document summary' 

# create that manually for the lookup
secondary_species_lookup["489466"] <- ""

secondary_species_lookup_vec <- unlist(secondary_species_lookup)


GAF_model_organisms_final$primary_taxon_id <- GAF_model_organisms_final$primary_taxon_id %>% as.character()
GAF_model_organisms_final$secondary_taxon_id <- GAF_model_organisms_final$secondary_taxon_id %>% as.character()
# Map species names back to the dataframe
GAF_model_organisms_final <- GAF_model_organisms_final %>%
  mutate(
    Primary_Species = primary_species_lookup[primary_taxon_id],
    Secondary_Species = ifelse(!is.na(secondary_taxon_id), secondary_species_lookup_vec[secondary_taxon_id], "")
  )


GAF_model_organisms_final$Secondary_Species <- GAF_model_organisms_final$Secondary_Species %>% as.factor()
summary(GAF_model_organisms_final$Secondary_Species)
# Check the resulting dataframe
head(GAF_model_organisms_final[c("taxon_id", "Primary_Species", "Secondary_Species")])


write.csv(GAF_model_organisms_final, "../output/model_species_with_pmids_and_species.csv")

GAF_model_organisms_final <- read.csv("../output/model_species_with_pmids_and_species.csv")
# Step 5: 

# for each primary species sort the PMIDs base don number of rows contributed and 
# pick 3 from the top 100 PMIDs, 3 from the top half (except top 100), and 4 from the bottom half 
# at random, then save this as training_set_model_species_with_pmids_and_species.csv

pmid_counts <- GAF_model_organisms_final %>%
  group_by(Primary_Species, PMIDs) %>%
  summarise(row_count = n(), .groups = "drop") %>%
  arrange(Primary_Species, desc(row_count))

# Only keep species with >100 PMIDs 
pmid_counts_filtered <- pmid_counts %>%
  group_by(Primary_Species) %>%
  filter(n() > 100) %>%  
  ungroup()
table(pmid_counts_filtered$Primary_Species)
# Function to sample PMIDs based on percentile splits
sample_pmids_for_species <- function(data) {
  # Rank PMIDs based on row counts
  data <- data %>% arrange(desc(row_count))
  n <- nrow(data)  # Total number of PMIDs for this species
  
  # Calculate indices for the percentile-based splits
  top_5_percent <- data %>% slice(1:ceiling(0.05 * n))
  top_50_to_95_percent <- data %>% slice(ceiling(0.05 * n) + 1:ceiling(0.95 * n))
  bottom_50_percent <- data %>% slice(ceiling(0.95 * n) + 1:n)
  
  # Randomly sample from each range
  sampled_pmids <- c(
    sample(top_5_percent$PMIDs, min(3, nrow(top_5_percent))),
    sample(top_50_to_95_percent$PMIDs, min(3, nrow(top_50_to_95_percent))),
    sample(bottom_50_percent$PMIDs, min(4, nrow(bottom_50_percent)))
  )
  
  return(sampled_pmids)
}

set.seed(99) # this seed was used for v2
# no seed was set for v1 :/ 
# Apply sampling to each species with sufficient PMIDs
training_set_pmids <- pmid_counts_filtered %>%
  group_by(Primary_Species) %>%
  summarise(selected_PMIDs = list(sample_pmids_for_species(pick(c(PMIDs, row_count)))), .groups = "drop") %>%
  unnest(cols = c(selected_PMIDs)) %>%
  rename(PMIDs = selected_PMIDs)  # Rename for joining

# Filter the main data based on the selected PMIDs for each species
training_set <- GAF_model_organisms_final %>%
  semi_join(training_set_pmids, by = c("Primary_Species", "PMIDs"))

# Check the resulting training set
head(training_set)
training_set$PMIDs %>% unique() %>% length()
training_set$Primary_Species %>% unique() %>% length()
# we have ended up with 229 papers spanning 23 species?


# Write to CSV
write.csv(training_set, "../output/training_set_model_species_with_pmids_and_species_v2.csv", row.names = FALSE)


