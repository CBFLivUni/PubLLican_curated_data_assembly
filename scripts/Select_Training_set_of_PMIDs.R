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

get_species_name <- function(taxon_id) {
  record <- entrez_summary(db = "taxonomy", id = taxon_id)
  return(record$scientificname)
}


GAF_merged <- read.delim("../output/merged_gaf.tsv", sep = "\t", header = TRUE)
GAF_headers <- read.csv("../data/GAF_headers.csv")
Evidence_Codes <- read.csv("../data/EvidenceCodesSelection.csv")
evidence_codes_for_training_set <- Evidence_Codes$Evidence_Code[Evidence_Codes$Include_Exclude == "include"]
##### tidying the merged GAF file
# Step 1 - add meaningful headers from the GO consortium GAF headers description
head(GAF_merged)

names(GAF_merged) <- c(GAF_headers$clean_header, "Species", "Strain")

# Step 2 - filter data based on evidecne codes
GAF_merged_filt <- GAF_merged[GAF_merged$Evidence_Code %in% evidence_codes_for_training_set, ]
# these are ~ 30% of entries

# Step 3 - convert references to PMIDs if possible
head(GAF_merged_filt)
# get anything before the :, and extract unique list of options from other refs 
# to view all possibilities

# Split by "|" first, then split each component by ":", and extract unique prefixes

# Split by "|" to handle multiple references
refs_list <- GAF_merged_filt$DB_Reference %>% 
  as.character() %>%
  strsplit("\\|") 

# if any of them start with PMID:, add them to a new column called PMIDs
# initialise the column
GAF_merged_filt$PMIDs <- NA

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
GAF_merged_filt$PMIDs <- sapply(refs_list_processed, function(x) x$PMID)

unique_prefixes <- unlist(lapply(refs_list_processed, function(x) x$prefixes)) %>%
  unique()

# > unique_prefixes
# "VEuPathDB" "AspGD_REF" "GO_REF"    "CGD_REF"   "NCBI"      "SGD_REF"   
# "DawsonLab" "MID"       "TrypTag" 


# GAF_merged_filt$DB_Reference[grepl("NCBI", GAF_merged_filt$DB_Reference)]
# we want to convert some of the other ID types present to PMIDs where no PMID is 
# already available.

# first lets split the data into two subsets - with PMIDs and without. 
# clean the NULLs to actually be "", even when separated by |
GAF_merged_filt$PMIDs <- gsub("NULL", "", GAF_merged_filt$PMIDs)
GAF_merged_filt$PMIDs <- gsub("\\|", "", GAF_merged_filt$PMIDs)
GAF_merged_filt$PMIDs <- gsub(" ", "", GAF_merged_filt$PMIDs)

GAF_merged_withpmid <- GAF_merged_filt[GAF_merged_filt$PMIDs != "", ]
GAF_merged_nopmid <- GAF_merged_filt[GAF_merged_filt$PMIDs == "", ]


# Extract the reference type for those without (prefix before ":")
GAF_merged_nopmid <- GAF_merged_nopmid %>%
  mutate(ref_type = sapply(strsplit(as.character(DB_Reference), ":"), `[`, 1))

# Calculate total and unique counts per reference type that don't have PMID
reference_summary <- GAF_merged_nopmid %>%
  group_by(ref_type) %>%
  summarise(
    total_count = n(),
    unique_count = n_distinct(DB_Reference)
  )
# 
# ref_type  total_count unique_count
# <chr>           <int>        <int>
# 1 AspGD_REF        7789         1635
# 2 DawsonLab         799          397
# 3 GO_REF         153936           17
# 4 MID                 1            1
# 5 NCBI              139            2
# 6 SGD_REF            19            1
# 7 TrypTag         14549         7079
# 8 VEuPathDB     1974857       113154

# we might be able to convert NCBI ones easily

NCBIs <-  GAF_merged_nopmid$DB_Reference[GAF_merged_nopmid$ref_type == "NCBI"]
# they are all from gene expresison omnibuss, no PMIDs, leave for now. 

# Step 4 - tidy up those with PMID to only include the PMID and nothing else
# remove species already in the consortium model species
model_species_training_set <- read.csv("../output/training_set_model_species_with_pmids_and_species.csv")

model_species_training_set$Primary_Species %>% unique()
GAF_merged_withpmid$Species %>%unique()

# Species to exclude
exclude_species <- c(
  "C.albicans",         # Candida albicans SC5314
  "N.glabratus",        # Nakaseomyces glabratus CBS 138
  "S.cerevisiae",       # Saccharomyces cerevisiae S288C
  "S.pombe",            # Schizosaccharomyces pombe
  "P.falciparum",       # Plasmodium falciparum 3D7
  "L.major",            # Leishmania major strain Friedlin
  "T.brucei"            # Trypanosoma brucei brucei TREU927
)

GAF_merged_withpmid_filtered <- GAF_merged_withpmid[!GAF_merged_withpmid$Species %in% exclude_species, ]

table(GAF_merged_withpmid_filtered$Species)
write.csv(GAF_merged_withpmid_filtered, "../output/rare_species_with_pmids_and_species_model_species_excluded.csv")



# Step 5: 

# for each species pick 10 papers where possible, don't filter species with < 100 PMIDs

pmid_counts <- GAF_merged_withpmid_filtered %>%
  group_by(Species, PMIDs) %>%
  summarise(row_count = n(), .groups = "drop") %>%
  arrange(Species, desc(row_count))

# Only keep species with >100 PMIDs 
pmid_counts_filtered <- pmid_counts %>%
  group_by(Species) %>%
  filter(n() > 1) %>%  
  ungroup()
table(pmid_counts_filtered$Species)


# Function to sample PMIDs based on percentile splits
sample_pmids_for_species <- function(data) {
  # Rank PMIDs based on row counts
  data <- data %>% arrange(desc(row_count))
  n <- nrow(data)  # Total number of PMIDs for this species
  if (n >= 100) {
    
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
    
    
  } else {
    
    # just pick 10 papers at random
    sampled_pmids <- sample(data$PMIDs, min(10, n))
  }

  return(sampled_pmids)
}

# Apply sampling to each species with sufficient PMIDs
training_set_pmids <- pmid_counts_filtered %>%
  group_by(Species) %>%
  summarise(selected_PMIDs = list(sample_pmids_for_species(pick(c(PMIDs, row_count)))), .groups = "drop") %>%
  unnest(cols = c(selected_PMIDs)) %>%
  rename(PMIDs = selected_PMIDs)  # Rename for joining

# Filter the main data based on the selected PMIDs for each species
training_set <- GAF_merged_withpmid_filtered %>%
  semi_join(training_set_pmids, by = c("Species", "PMIDs"))

# Check the resulting training set
head(training_set)
training_set$PMIDs %>% unique() %>% length()
training_set$Species %>% unique() %>% length()
# we have ended up with 215 papers spanning 25 species


# Write to CSV
write.csv(training_set, "../output/training_set_rare_species_with_pmids_and_species.csv", row.names = FALSE)



training_set_unique_pmids <- read.csv("../output/training_set_combined.csv")$PMIDs %>% unique()

writeLines(training_set_unique_pmids, "../output/training_set_unique_pmids.txt")





## repeat stpe 5 to increase number of papers - some fail at the retroeval step
## need to get workaround but doesn't have to impede training set development
# Step 5_repeated: 

GAF_rare_withpmid <- read.csv("../output/training_set_rare_species_with_pmids_and_species.csv")
# for each species pick 10 papers where possible, don't filter species with < 100 PMIDs

pmid_counts_rare <- GAF_rare_withpmid %>%
  group_by(Species, PMIDs) %>%
  summarise(row_count = n(), .groups = "drop") %>%
  arrange(Species, desc(row_count))

# Only keep species with >100 PMIDs 
pmid_counts_rare_filtered <- pmid_counts_rare %>%
  group_by(Species) %>%
  filter(n() > 1) %>%  
  ungroup()
table(pmid_counts_rare_filtered$Species)


# Function to sample PMIDs based on percentile splits
sample_pmids_for_species <- function(data) {
  # Rank PMIDs based on row counts
  data <- data %>% arrange(desc(row_count))
  n <- nrow(data)  # Total number of PMIDs for this species
  if (n >= 100) {
    
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
    
    
  } else {
    
    # just pick 10 papers at random
    sampled_pmids <- sample(data$PMIDs, min(10, n))
  }
  
  return(sampled_pmids)
}

set.seed(99)
# Apply sampling to each species with sufficient PMIDs
training_set_rare_pmids <- pmid_counts_rare_filtered %>%
  group_by(Species) %>%
  summarise(selected_PMIDs = list(sample_pmids_for_species(pick(c(PMIDs, row_count)))), .groups = "drop") %>%
  unnest(cols = c(selected_PMIDs)) %>%
  rename(PMIDs = selected_PMIDs)  # Rename for joining


# Filter the main data based on the selected PMIDs for each species
training_set <- GAF_rare_withpmid %>%
  semi_join(training_set_rare_pmids, by = c("Species", "PMIDs"))

# Check the resulting training set
head(training_set)
training_set$PMIDs %>% unique() %>% length()
training_set$Species %>% unique() %>% length()
# we have ended up with 215 papers spanning 25 species


# Write to CSV
write.csv(training_set, "../output/training_set_rare_species_with_pmids_and_species_v2.csv", row.names = FALSE)

#compare set 1 to set 2, only include the ones that arent in set 1 for running set 2 
training_set_v1_unique_pmids <- read.csv("../output/training_set_combined.csv")$PMIDs %>% unique()
training_set_v2_unique_pmids <- read.csv("../output/training_set_combined_v2.csv")$PMIDs %>% unique()

set2_not_in_set1_pmids <- setdiff(training_set_v2_unique_pmids, training_set_v1_unique_pmids)


# training_set_1_failed_pmids <- readLines("../output/failed_pmids.txt") # for reference, the ones that failed  - we can look into it further. 

writeLines(training_set_v1_unique_pmids, "../output/training_set_unique_pmids_v1.txt")
writeLines(training_set_v2_unique_pmids, "../output/training_set_unique_pmids_v2.txt")
writeLines(set2_not_in_set1_pmids, "../output/training_set_unique_pmids_set2_not_in_set1.txt") # only the ones in set 2 but not set 1
