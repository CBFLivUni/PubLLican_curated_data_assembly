######## libs and read in #########
library(tidyverse)
library(rentrez) # for converting taxon id to species name
library(httr) # for getting 
library(jsonlite)

## a large fraction of pmids cannot be easily obtained as JSON from pubmed,
## extend the training set so we get a larger number of papers actually available
## ~ 20% pmids succeed, so increasing total set to ~ 2500 would be good. 

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



# Function to sample PMIDs based on percentile splits
sample_pmids_for_species <- function(pmid_count_data, set_type = "model") {
  
  ########################## set params ###########
  #### change below number of papers and percentiles to select larger number. 
  #### perhaps ~ 1000 from model species and ~1500 from rest
  top_percentile = 0.1 # top 10 %
  
  mid_percentile = 0.5 # upper half of pmids when ranked by row number
  
  bot_percentile = 1 # (all other pmids)
  
  top_n_pmids = 12 # 4x more than original
  mid_n_pmids = 12
  bot_n_pmids = 16
  
  non_model_n = 60 # was 10 but we need ~ 6x more papers here
  
  
  # Rank PMIDs based on row counts
  data <- pmid_count_data %>% arrange(desc(row_count))
  n <- nrow(data)  # Total number of PMIDs for this species
  
  if (set_type == "model") {
    
    # Calculate indices for the percentile-based splits
    top <- data %>% slice(1:ceiling(top_percentile * n))
    mid <- data %>% slice((ceiling(top_percentile * n) + 1):ceiling(mid_percentile * n))
    bot <- data %>% slice((ceiling(mid_percentile * n) + 1):n)
    
    # Randomly sample from each range
    sampled_pmids <- c(
      sample(top$PMIDs, min(top_n_pmids, nrow(top))),
      sample(mid$PMIDs, min(mid_n_pmids, nrow(mid))),
      sample(bot$PMIDs, min(bot_n_pmids, nrow(bot)))
    )
    
    
  } else {
    
    # just pick up to 60 papers at random
    sampled_pmids <- sample(data$PMIDs, min(non_model_n, n))
  }
  
  return(sampled_pmids)
}


# set differnet seed to before just in case we were unlucky 
set.seed(43)
# Apply sampling to each species with sufficient PMIDs
training_set_rare_pmids <- pmid_counts_rare_filtered %>%
  group_by(Species) %>%
  summarise(selected_PMIDs = list(sample_pmids_for_species(pick(c(PMIDs, row_count)), 
                                                           set_type = "rare")),
            .groups = "drop") %>%
  unnest(cols = c(selected_PMIDs)) %>%
  rename(PMIDs = selected_PMIDs)  # Rename for joining


# Filter the main data based on the selected PMIDs for each species
training_set_rare <- GAF_rare_withpmid %>%
  semi_join(training_set_rare_pmids, by = c("Species", "PMIDs"))

# Check the resulting training set
head(training_set_rare)
training_set_rare$PMIDs %>% unique() %>% length()
training_set_rare$Species %>% unique() %>% length()
# we have ended up with 215 papers spanning 25 species
GAF_rare_withpmid$PMIDs %>% unique() %>% length()
# we actually only have 215 pmids of rare species! 
# Write to CSV
write.csv(training_set_rare, "../output/extended_training_set_rare_species.csv", row.names = FALSE)
rare_pmids_all <- training_set_rare$PMIDs %>% unique()
writeLines(rare_pmids_all, "../output/all_rare_pmids.txt")
# Step 5: 

GAF_model_organisms_final <- read.csv("../output/model_species_with_pmids_and_species.csv")
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
# we have ended up with 919 papers spanning 23 species?

model_pmids_extended <- training_set$PMIDs %>% unique()
# Write to CSV
write.csv(training_set, "../output/extended_training_set_model_species.csv", row.names = FALSE)
writeLines(model_pmids_extended, "../output/extended_model_pmids.txt")

writeLines(c(model_pmids_extended, rare_pmids_all) , "../output/all_pmids_to_process.txt")
