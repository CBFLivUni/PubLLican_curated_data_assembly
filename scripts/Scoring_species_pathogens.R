# Load required libraries
library(jsonlite)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyverse)


######## leftover todos ######
# leftover todos:
# TODO: if the species is Candida glabrata or Nakaseomyces glabratus 
# consider both as possibility

# TODO: gene level matching tidy up, not quite in the format I want it right now. 

##### step by step process ######
# Step 1: get the results from the status JSON file as a df with columns:

# PMID = file name (strip .json from end)
# species_latin_name
# gene_ID
# gene_name
# gene_symbol

# for all gene-species combinations of each pmid

# Define the folder containing the JSON files
pmid_status_folder <- "C:/Users/jtzve/Desktop/llm-testbed_JT/caches/status"


# Get a list of all JSON files in the folder
json_files <- list.files(pmid_status_folder, pattern = "\\.json$", full.names = TRUE)


# Initialise an empty df to store the results
results_df <- data.frame(
  PMID = character(),
  species_latin_name = character(),
  gene_ID = character(),
  gene_name = character(),
  gene_symbol = character(),
  stringsAsFactors = FALSE
)

# Loop through each JSON file
for (json_file in json_files) {
  # json_file_test_error <- "C:/Users/jtzve/Desktop/llm-testbed_JT/caches/status/19497128.json"
  # Extract the PMID from the filename (remove folder path and .json extension)
  pmid <- sub("\\.json$", "", basename(json_file))
  # print(pmid)
  # Read the JSON file
  json_data <- fromJSON(json_file, flatten = TRUE)
  
  # Check if "getPaperGenes" exists and succeeded
  if (!is.null(json_data$getPaperGenes) && json_data$getPaperGenes$success) {
    
    # Extract species and genes information
    species_df <- json_data$getPaperGenes$response$species
    # for each species
    for (i in 1:nrow(species_df)) {
      # print(i)
      # i <- 4
      # get the name of the species
      species_name <- species_df$name[i]
      # and the associated gene information
      genes_df <- species_df$genes[[i]]
      
      # check if any genes assigned to species, if not add the species
      # and assigne all gene IDs to "no genes"
      if (is.null(genes_df) || nrow(genes_df) == 0) {
        
        results_df <- rbind(
          results_df,
          data.frame(
            PMID = pmid,
            species_latin_name = species_name,
            gene_ID = "no genes",
            gene_name = "no genes",
            gene_symbol = "no genes",
            stringsAsFactors = FALSE
          )
        )
        
      } else {
        
        # if yes, add the gene-species pairs to the results
        # for each gene
        for (j in 1:nrow(genes_df)) {
          # extract the relevant identifiers and append to the df
          results_df <- rbind(
            results_df,
            data.frame(
              PMID = pmid,
              species_latin_name = species_name,
              gene_ID = genes_df$identifier[j],
              gene_name = genes_df$name[j],
              gene_symbol = genes_df$symbol[j],
              stringsAsFactors = FALSE
            )
          )
        }
        
        
      }
      
    }
  }
}

successful_pmids <- results_df$PMID %>% unique()
write_lines(successful_pmids, 
            file =  "../output/pahogen_training_set_pmids_data_exploration/successful_pmids.txt")

write.csv(results_df, "../output/pathogens_training_set/pathogens_single_genes_prompt_complete_outputs.csv", row.names = FALSE)


rm(list = c("genes_df", "json_data", "species_df", "i", "j", "json_file", "pmid", "species_name"))
gc()


# Step 2: get the original curated dataset and strip down to
# info related to the DB, gene names and symbols, and species; 
# we might consider evidence codes at another stage to see if certain codes are
# missed out more often
curated_df <- read.csv("../output/pathogens_training_set/pathogens_training_set.csv")

names(curated_df)  
head(curated_df)

# we want to keep the following columns
GAF_cols_of_interest <- c("DB",
                          "DB_ID",
                          "DB_Symbol",
                          # "Evidence_Code",
                          "Name",
                          "Synonym",
                          "Species",
                          "Category",
                          "PMIDs")


curated_df <- curated_df[, GAF_cols_of_interest] %>% distinct()

# step 3: assess model's recall on the species level

# get curated species data - keep only pmids the model has actually processed
# we only want the species and PMIDs at this stage

species_curated_df <- curated_df[curated_df$PMIDs %in% results_df$PMID, c("PMIDs", "Species")] %>% distinct()

# do reverse subsetting to ensure cases where PMID extraction was bugged won't get in the way
# here we include cases where no genes were assigned to the species, perhaps we should try scorin  with these
# cases excluded and compare 
species_test_results <- results_df[results_df$PMID %in% species_curated_df$PMIDs, c("PMID", "species_latin_name")] %>% distinct()


# we have the full latin name from the rest results
# in the curated set we either have full latin name + strain or First name initial.second name. 
# lets convert the curated names to the same format as the test set schema
curated_species <- species_curated_df$Species %>% unique()

# > curated_species
# [1] "A.fumigatus"                                          
# [2] "A.nidulans"                                           
# [3] "C.neoformans"                                         
# [4] "F.graminearum"                                        
# [5] "P.oryzae"                                             
# [6] "P.berghei"                                            
# [7] "Pchabaudichabaudi"                                    
# [8] "P.gallinaceum"                                        
# [9] "P.knowlesi"                                           
# [10] "P.malariae"                                           
# [11] "P.ovalecurtisi"                                       
# [12] "P.reichenowi"                                         
# [13] "P.relictum"                                           
# [14] "P.vivax"                                              
# [15] "P.yoeliiyoelii"                                       
# [16] "T.gondii"                                             
# [17] "L.braziliensis"                                       
# [18] "L.donovani"                                           
# [19] "L.infantum"                                           
# [20] "L.mexicana"                                           
# [21] "T.congolense"                                         
# [22] "T.cruzi"                                              
# [23] "T.vivax"                                              
# [24] "Candida albicans SC5314"                              
# [25] "Nakaseomyces glabratus CBS 138"                       
# [26] "Candidozyma auris"                                    
# [27] "Candida dubliniensis CD36"                            
# [28] "Candida parapsilosis CDC317"                          
# [29] "Leishmania major strain Friedlin"                     
# [30] "Plasmodium falciparum 3D7"                            
# [31] "Severe acute respiratory syndrome-related coronavirus"
# [32] "Trypanosoma brucei brucei TREU927"

# Function to expand and standardize species names from chatGPT
standardise_species_name <- function(species) {
  # Check for strain information (words after the first two parts)
  parts <- strsplit(species, " ")[[1]]
  if (length(parts) > 2) {
    species <- paste(parts[1], parts[2])
  }
  
  # Replace initials with full names if applicable
  species <- gsub("A\\.", "Aspergillus ", species)
  species <- gsub("C\\.", "Cryptococcus ", species)
  species <- gsub("F\\.", "Fusarium ", species)
  species <- gsub("P\\.", "Plasmodium ", species)
  species <- gsub("T\\.", "Toxoplasma ", species)
  species <- gsub("L\\.", "Leishmania ", species)
  species <- gsub("T\\.", "Trypanosoma ", species)
  
  # Remove extra information like strain names
  species <- gsub("strain .*", "", species)
  species <- gsub("CBS .*", "", species)
  species <- gsub("SC[0-9]+", "", species)
  species <- gsub("CD[0-9]+", "", species)
  species <- gsub("3D7", "", species)
  species <- gsub("T[REU927]+", "", species)
  
  # Trim leading/trailing spaces
  species <- trimws(species)
  
  return(species)
}

# Apply the function to the curated species
standardised_species <- sapply(curated_species, standardise_species_name)


# Add standardised species names to the curated and test dfs
species_curated_df$standardised_species <- sapply(species_curated_df$Species, standardise_species_name)

species_test_results$standardised_species <- sapply(species_test_results$species_latin_name, standardise_species_name)

# Merge test and curated data on standardised species and PMIDs
merged_results <- merge(
  species_curated_df,
  species_test_results,
  by.x = c("PMIDs", "standardised_species"),
  by.y = c("PMID", "standardised_species"),
  all = TRUE # Include all to identify missed and additional species
)

# Mark correctly recalled species
merged_results$retrieved <- !is.na(merged_results$species_latin_name) & !is.na(merged_results$Species)

# Group by PMID and calculate detailed recall information
recall_per_pmid <- merged_results %>%
  group_by(PMIDs) %>%
  summarise(
    total_curated = sum(!is.na(Species)),            # Total curated species
    total_retrieved = sum(!is.na(species_latin_name)), # Total retrieved species
    correctly_recalled = sum(retrieved),            # Correctly matched species
    correctly_recalled_species = paste(Species[retrieved], collapse = "; "), # Correctly recalled species as string
    missed_species = paste(Species[is.na(species_latin_name)], collapse = "; "), # Missed species as string
    additional_species = paste(species_latin_name[is.na(Species)], collapse = "; "), # Additional species as string
    recall = correctly_recalled / total_curated     # Recall score for this PMID
  ) %>%
  ungroup() # Ungroup to make it a regular dataframe

# Calculate overall recall (average of all PMIDs' recalls)
overall_recall <- mean(recall_per_pmid$recall, na.rm = TRUE)

# Save the recall per PMID dataframe to a CSV file
write.csv(recall_per_pmid, "../output/pathogens_training_set/species_recall_per_pmid.csv", row.names = FALSE)

####### lets save some plots #####
output_folder <- "../output/pathogens_training_set/"


# Distribution of Recall Scores
# Add a histogram with counts displayed above bars
ggplot(recall_per_pmid, aes(x = recall)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  stat_bin(binwidth = 0.1, geom = "text", aes(label = ..count..), vjust = -0.5, color = "black") +
  labs(
    title = "Distribution of Species Recall Scores",
    x = "Recall Score",
    y = "Count of PMIDs"
  ) +
  theme_classic()

ggsave(filename = paste0(output_folder, "distribution_species_recall_scores.png"), 
       dpi = 600, height = 12, width = 8, units = "cm")

# Species Missed the Most
missed_species_counts <- merged_results %>%
  filter(is.na(species_latin_name)) %>%
  count(Species, sort = TRUE)

ggplot(missed_species_counts, aes(x = reorder(Species, -n), y = n)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(
    title = "Species Missed the Most",
    x = "Species",
    y = "Number of PMIDs Missed"
  ) +
  coord_flip() +
  theme_classic()

ggsave(filename = paste0(output_folder, "species_missed_most_n.png"), 
       dpi = 600, height = 12, width = 16, units = "cm")

# Calculate species-level recall
species_total_counts <- species_curated_df %>%
  count(Species, name = "total_pmids")

species_recall <- missed_species_counts %>%
  full_join(species_total_counts, by = "Species") %>%
  mutate(
    missed_pmids = coalesce(n, 0),                     # Replace NA with 0 for missed counts
    percent_missed = (missed_pmids / total_pmids) * 100 # Calculate percentage missed
  )

# Plot Species Missed the Most with % Missed on X-axis and Number Missed as Text
ggplot(species_recall, aes(x = reorder(Species, percent_missed), y = percent_missed)) +
  geom_bar(stat = "identity", fill = "salmon") +
  geom_text(aes(label = missed_pmids), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    title = "Species Missed the Most (% Missed)",
    x = "Species",
    y = "% Missed (PMIDs)"
  ) +
  coord_flip() +
  theme_classic()
ggsave(filename = paste0(output_folder, "species_missed_most_percent.png"), 
       dpi = 600, height = 12, width = 16, units = "cm")

# Recall vs. Total Curated Species per PMID
# Calculate the frequency of each total_curated and recall combination
recall_freq <- recall_per_pmid %>%
  group_by(total_curated, recall) %>%
  summarise(frequency = n(), .groups = "drop")

ggplot(recall_freq, aes(x = total_curated, y = recall)) +
  geom_point(aes(size = frequency), color = "dodgerblue", alpha = 0.7) +
  
  scale_size_continuous(
    range = c(3, 20), # Set minimum and maximum point sizes
    breaks = c(1, 5, 50, 100, 500, 1000) # Explicitly define the breakpoints
  ) +
  geom_smooth(
    aes(weight = frequency), # Weight the LOESS line by frequency
    method = "loess", span = 1, color = "gold", linewidth = 1, se = FALSE
  ) +
  geom_text(
    aes(label = frequency), 
    hjust = -0.7, # Move text to the right of each dot
    color = "black", 
    size = 5, 
    fontface = "bold", # Make text bold
    family = "sans"
  ) + 
  labs(
    title = "Recall vs. Total Curated Species per PMID",
    x = "Total Curated Species per PMID",
    y = "Recall Score",
    size = "Frequency"
  ) +
  theme_classic()


# Save the plot
ggsave(filename = paste0(output_folder, "recall_vs_total_curated_species.png"), 
       dpi = 600, height = 12, width = 16, units = "cm")

# Top Additional Species - plot these appearing in more than 10 PMIDs
additional_species_counts <- merged_results %>%
  filter(is.na(Species)) %>%
  count(species_latin_name, sort = TRUE)

additional_species_counts <- additional_species_counts[additional_species_counts$n >=10,]

ggplot(additional_species_counts, aes(x = reorder(species_latin_name, -n), y = n)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  labs(
    title = "Top Additional Species Retrieved",
    x = "Species",
    y = "Number of PMIDs"
  ) +
  coord_flip() +
  theme_classic()

ggsave(filename = paste0(output_folder, "frequent_additional_species.png"), 
       dpi = 600, height = 12, width = 16, units = "cm")


############################### accuracy sanity check ################
plaintext_files_dir <- "C:/Users/jtzve/Desktop/llm-testbed_JT/caches/articles/plaintext"

# for each pmid look get {pmid}.txt from the plaintext folder and read it in 

# for each species in the test data for that PMID

# check if the species assigned by the AI was in the plaintext and count how many 
# times it appears


# check if the species assigned by curators was in the plaintext and copunt 
# how many times it appears
# handle common and scientific names
# handle cases where we have only the first or second name of a species
# handle the exception for the name change since 2022 for 
# Candida glabrata to Nakaseomyces glabratus

# calculate an accuracy score for both curators and the AI 


# step one create a map of scientific to common name: 
scientific_to_common <- c(
  "Homo sapiens" = "human",
  "Mus musculus" = "mouse",
  "Rattus norvegicus" = "rat",
  "Xenopus laevis" = "frog",
  "Gallus gallus" = "chicken",
  "Bos taurus" = "cow",
  "Saccharomyces cerevisiae" = "yeast",
  "Candida albicans" = "thrush",
  "Aspergillus fumigatus" = "aspergillus",
  "Drosophila melanogaster" = "fruit fly",
  "Arabidopsis thaliana" = "thale cress",
  "Plasmodium falciparum" = "malaria parasite",
  "Plasmodium vivax" = "malaria parasite",
  "Plasmodium knowlesi" = "malaria parasite",
  "Plasmodium berghei" = "malaria parasite",
  "Escherichia coli" = "e. coli",
  "Schizosaccharomyces pombe" = "fission yeast",
  "Leishmania major" = "leishmania",
  "Leishmania mexicana" = "leishmania",
  "Trypanosoma brucei" = "trypanosome",
  "Trypanosoma cruzi" = "chagas parasite",
  "Caenorhabditis elegans" = "roundworm",
  "Danio rerio" = "zebrafish",
  "Anopheles gambiae" = "mosquito",
  "Mycobacterium tuberculosis" = "tuberculosis bacterium",
  "Listeria monocytogenes" = "listeria",
  "Drosophila" = "fruit fly",
  "Anopheles stephensi" = "mosquito",
  "Giardia lamblia" = "giardia",
  "Cryptococcus neoformans" = "cryptococcus",
  "Eimeria tenella" = "eimeria",
  "Candida glabrata" = "glabrata yeast",
  "Trypanosoma congolense" = "trypanosome",
  "Babesia bovis" = "babesia",
  "Entamoeba histolytica" = "entamoeba",
  "Schistosoma mansoni" = "blood fluke",
  "Haemophilus influenzae" = "haemophilus",
  "Oryza sativa" = "rice",
  "Plasmodium relictum" = "malaria parasite",
  "Apicomplexa" = "apicomplexan",
  "Candida auris" = "auris yeast",
  "Aspergillus nidulans" = "nidulans aspergillus",
  "Fusarium graminearum" = "fusarium",
  "Magnaporthe oryzae" = "rice blast fungus",
  "Anopheles freeborni" = "mosquito",
  "Anopheles funestus" = "mosquito",
  "Eimeria sp." = "eimeria",
  "Eimeria spp." = "eimeria",
  "Theileria annulata" = "theileria",
  "Theileria parva" = "theileria",
  "Giardia intestinalis" = "giardia",
  "Cryptosporidium spp." = "cryptosporidium",
  "Cryptosporidium hominis" = "cryptosporidium",
  "Candida dubliniensis" = "dubliniensis yeast",
  "Candida tropicalis" = "tropicalis yeast",
  "Candida krusei" = "krusei yeast",
  "Candida lusitaniae" = "lusitaniae yeast",
  "Candida parapsilosis" = "parapsilosis yeast",
  "Candida orthopsilosis" = "orthopsilosis yeast",
  "Babesia divergens" = "babesia",
  "Babesia microti" = "babesia",
  "Plasmodium malariae" = "malaria parasite",
  "Plasmodium ovale" = "malaria parasite",
  "Plasmodium yoelii" = "malaria parasite",
  "Plasmodium spp." = "malaria parasite",
  "Leishmania tropica" = "leishmania",
  "Leishmania braziliensis" = "leishmania",
  "Leishmania infantum" = "leishmania",
  "Leishmania donovani" = "leishmania",
  "Leishmania guyanensis" = "leishmania",
  "Leishmania peruviana" = "leishmania",
  "Leishmania aethiopica" = "leishmania",
  "Leishmania amazonensis" = "leishmania",
  "Leishmania archibaldi" = "leishmania",
  "Leishmania brasiliensis" = "leishmania",
  "Leishmania tropica" = "leishmania",
  "Leishmania turanica" = "leishmania",
  "Trypanosoma congolense" = "trypanosome",
  "Trypanosoma evansi" = "trypanosome",
  "Trypanosoma grayi" = "trypanosome",
  "Trypanosoma rangeli" = "trypanosome",
  "Trypanosoma vivax" = "trypanosome"
)


# # Function to calculate accuracy and additional details for each PMID
# calculate_accuracy <- function(pmid, ai_species, curator_species) {
#   # Read the plaintext file
#   file_path <- file.path(plaintext_files_dir, paste0(pmid, ".txt"))
#   if (!file.exists(file_path)) {
#     warning(paste("Plaintext file not found for PMID:", pmid))
#     return(data.frame(
#       PMID = pmid,
#       ai_correct = 0,
#       curator_correct = 0,
#       ai_total = length(ai_species),
#       curator_total = length(curator_species),
#       ai_accuracy = NA,
#       curator_accuracy = NA,
#       correct_ai_guesses = NA,
#       incorrect_ai_guesses = NA,
#       correct_curator_guesses = NA,
#       incorrect_curator_guesses = NA
#     ))
#   }
#   text <- tolower(readLines(file_path)) %>% paste(collapse = " ")
#   
#   
#   # TODO: implement exceptions for cases where in text it would be mentioned 
#   # as a different name, just second name, First anme . last name, or common name
#   # as in the case of H. sapiens, human, or the name change since 2022 for 
#   # Candida glabrata to Nakaseomyces glabratus
#   
#   # Check occurrences for AI-assigned species
#   ai_matches <- sapply(ai_species, function(species) {
#     str_count(text, tolower(species))
#   })
#   
#   # Check occurrences for curator-assigned species
#   curator_matches <- sapply(curator_species, function(species) {
#     str_count(text, tolower(species))
#   })
#   
#   # TODO: change code to somehow capture the number of occurances of each variation
#   # of a species, can be useful to explore when engineering prompts. 
#   
#   # Identify correct and incorrect guesses
#   correct_ai <- ai_species[ai_matches > 0]
#   incorrect_ai <- ai_species[ai_matches == 0]
#   correct_curator <- curator_species[curator_matches > 0]
#   incorrect_curator <- curator_species[curator_matches == 0]
#   
#   # Calculate accuracies
#   ai_correct <- length(correct_ai)
#   curator_correct <- length(correct_curator)
#   ai_accuracy <- if (length(ai_species) > 0) ai_correct / length(ai_species) else NA
#   curator_accuracy <- if (length(curator_species) > 0) curator_correct / length(curator_species) else NA
#   
#   # Return results for this PMID
#   return(data.frame(
#     PMID = pmid,
#     ai_correct = ai_correct,
#     curator_correct = curator_correct,
#     ai_total = length(ai_species),
#     curator_total = length(curator_species),
#     ai_accuracy = ai_accuracy,
#     curator_accuracy = curator_accuracy,
#     correct_ai_guesses = paste(correct_ai, collapse = "; "),
#     incorrect_ai_guesses = paste(incorrect_ai, collapse = "; "),
#     correct_curator_guesses = paste(correct_curator, collapse = "; "),
#     incorrect_curator_guesses = paste(incorrect_curator, collapse = "; ")
#   ))
# }
# 


calculate_accuracy <- function(pmid, ai_species, curator_species, plaintext_files_dir, scientific_to_common, output_folder) {
  # Read the plaintext file
  file_path <- file.path(plaintext_files_dir, paste0(pmid, ".txt"))
  # print(pmid)
  if (!file.exists(file_path)) {
    warning(paste("Plaintext file not found for PMID:", pmid))
    return(data.frame(
      PMID = pmid,
      ai_correct = 0,
      curator_correct = 0,
      ai_total = length(ai_species),
      curator_total = length(curator_species),
      ai_accuracy = NA,
      curator_accuracy = NA,
      correct_ai_guesses = NA,
      incorrect_ai_guesses = NA,
      correct_curator_guesses = NA,
      incorrect_curator_guesses = NA
    ))
  }
  
  # Read and preprocess the text
  text <- tolower(readLines(file_path)) %>% paste(collapse = " ")
  
  # Helper function to generate variations for a species
  get_variations <- function(species) {
    
    # TODO: if the species is Candida glabrata or Nakaseomyces glabratus 
    # consider both as possibility
    full_name <- species
    first_name <- gsub(" .*", "", species)  # First word (genus)
    last_name <- gsub("^.* ", "", species)  # Last word (specific epithet or strain)
    
    f_dot_last <- paste0(substr(first_name, 1, 1), ". ", last_name)  # First initial. Last
    no_dots <- gsub("\\.", "", f_dot_last)     # Remove dots (e.g., "H. sapiens" -> "H sapiens")
    
    # Add common names if applicable
    common_name <- scientific_to_common[species]
    
    variations <- trimws(tolower(c(
      full_name,
      first_name,         # we will not match on that but want to keep it
      last_name,          # Allow matching on last name
      no_dots,            # Remove dots for matching
      f_dot_last,         # First letter dot Last name
      if (!is.na(common_name)) common_name else NULL
    )))
    
    types <- c(
      "full_name",
      "first_name",
      "last_name",
      "no_dots",
      "f_dot_last",
      if (!is.na(common_name)) "common_name" else NULL
    )
    
    return(list(variations = variations, types = types))
  }
  
  # Match species and calculate occurrences
  match_species_variations <- function(species_list) {
    results <- list()
    for (species in species_list) {
      # print(species)
      # species <- "Leishmania major" # test case
      variation_data <- get_variations(species)
      variations <- variation_data$variations
      types <- variation_data$types
      # first_name <- variation_data$first_name
      
      counts <- sapply(variations, function(variant) str_count(text, fixed(variant)))
      # first_name_count <- str_count(text, fixed(first_name))
      
      results[[species]] <- data.frame(
        Variation = variations,
        Count = counts,
        Type = types#,
        # FirstNameOnlyCount = first_name_count
      )
    }
    return(results)
  }
  
  # Match AI and curator species
  ai_results <- match_species_variations(ai_species)
  curator_results <- match_species_variations(curator_species)
  
  # Save variation counts as a CSV
  variation_counts_dir <- file.path(output_folder, "species_variation_counts")
  if (!dir.exists(variation_counts_dir)) {
    dir.create(variation_counts_dir, recursive = TRUE)
  }
  
  variation_counts_file <- file.path(variation_counts_dir, paste0(pmid, "_variation_counts.csv"))
  
  combined_results <- do.call(rbind, lapply(names(ai_results), function(species) {
    cbind(
      Species = species,
      Type = "AI",
      ai_results[[species]]
    )
  }))
  
  combined_results <- rbind(
    combined_results,
    do.call(rbind, lapply(names(curator_results), function(species) {
      cbind(
        Species = species,
        Type = "Curator",
        curator_results[[species]]
      )
    }))
  )
  
  # Save to CSV
  write.csv(combined_results, variation_counts_file, row.names = FALSE)
  
  # Identify correct and incorrect guesses
  ai_matches <- sapply(ai_species, function(species) {
    relevant_counts <- sum(ai_results[[species]]$Count[ai_results[[species]]$Type != "first_name"])
    return(relevant_counts > 0)
  })
  
  curator_matches <- sapply(curator_species, function(species) {
    relevant_counts <- sum(curator_results[[species]]$Count[curator_results[[species]]$Type != "first_name"])
    return(relevant_counts > 0)
  })
  
  correct_ai <- ai_species[as.logical(ai_matches)]
  incorrect_ai <- ai_species[!as.logical(ai_matches)]
  correct_curator <- curator_species[as.logical(curator_matches)]
  incorrect_curator <- curator_species[!as.logical(curator_matches)]
  
  # Calculate accuracies
  ai_correct <- length(correct_ai)
  curator_correct <- length(correct_curator)
  ai_accuracy <- if (length(ai_species) > 0) ai_correct / length(ai_species) else NA
  curator_accuracy <- if (length(curator_species) > 0) curator_correct / length(curator_species) else NA
  
  # Return results for this PMID
  return(data.frame(
    PMID = pmid,
    ai_correct = ai_correct,
    curator_correct = curator_correct,
    ai_total = length(ai_species),
    curator_total = length(curator_species),
    ai_accuracy = ai_accuracy,
    curator_accuracy = curator_accuracy,
    correct_ai_guesses = paste(correct_ai, collapse = "; "),
    incorrect_ai_guesses = paste(incorrect_ai, collapse = "; "),
    correct_curator_guesses = paste(correct_curator, collapse = "; "),
    incorrect_curator_guesses = paste(incorrect_curator, collapse = "; ")
  ))
}




# Extract unique PMIDs
pmid_list <- unique(merged_results$PMIDs)
# pmid_list <- "35972967" # test case
# pmid_list <- "14657221" # test case
# Function to get species for a given PMID and assign it to AI or curator groups
get_species_for_pmid <- function(pmid) {
  # Subset data for this PMID
  pmid_data <- merged_results[merged_results$PMIDs == pmid, ]
  
  # Extract AI-assigned species excluding "unknown species"
  ai_species <- pmid_data$standardised_species[!is.na(pmid_data$species_latin_name) & pmid_data$species_latin_name != "unknown species"]
  
  # Extract curator-assigned species
  curator_species <- pmid_data$standardised_species[!is.na(pmid_data$Species)]
  
  return(list(ai_species = ai_species, curator_species = curator_species))
}

# Run the accuracy calculation for all PMIDs
accuracy_results <- do.call(rbind, lapply(pmid_list, function(pmid) {
  # pmid <- pmid_list[1]
  species_data <- get_species_for_pmid(pmid)
  
  calculate_accuracy(
    pmid = pmid, 
    ai_species = species_data$ai_species, 
    curator_species = species_data$curator_species, 
    plaintext_files_dir = plaintext_files_dir, 
    scientific_to_common = scientific_to_common, 
    output_folder = output_folder
  )
}))


# Save results to CSV
write.csv(accuracy_results, "../output/pathogens_training_set/species_accuracy_results_updated_genus_discounted_as_match.csv", row.names = FALSE)


# Combine AI and curator accuracy into a long format for comparison
accuracy_long <- accuracy_results %>%
  select(PMID, ai_accuracy, curator_accuracy) %>%
  pivot_longer(cols = c(ai_accuracy, curator_accuracy), 
               names_to = "source", values_to = "accuracy") %>%
  mutate(source = recode(source, 
                         ai_accuracy = "AI", 
                         curator_accuracy = "Curator"))


# # Plot accuracy distribution with facet wrap
# ggplot(accuracy_long, aes(x = accuracy)) +
#   geom_histogram(fill = "skyblue", color = "black", bins = 20) +
#   facet_wrap(~source, ncol = 1) + # Separate plots for AI and Curator
#   labs(
#     title = "Accuracy Distribution: AI vs. Curators",
#     x = "Accuracy Score",
#     y = "Count of PMIDs"
#   ) +
#   theme_classic()


# Plot accuracy distribution with side-by-side bars
ggplot(accuracy_long, aes(x = accuracy, fill = source)) +
  geom_histogram(alpha = 0.8, position = "dodge", bins = 20, color = "black") +
  scale_fill_manual(values = c("AI" = "skyblue", "Curator" = "salmon")) +
  labs(
    title = "Accuracy Distribution based on plaintext: AI vs. Curators",
    x = "Accuracy Score",
    y = "Count of PMIDs",
    fill = "Source"
  ) +
  theme_classic()

ggsave(filename = paste0(output_folder, "species_accuracy_from_plaintext_byPMID_updated.png"), 
       dpi = 600, height = 12, width = 16, units = "cm")



# Calculate species-level accuracy:
# - for each species guessed by the AI, how many times
# was it in the accuracy_results$correct_ai_guesses versus incorrect_ai_guesses column

# perform the same for the curator-assigned species

# make separate plots for AI and curator data
# where species are ordered by total number of mentions by the AI/curator
# a bar is shown for the fraction of accurate assignments and a text next to the 
# bar is added showing label e.g paste0(ai_correct, "/" , ai_total)
# Split and aggregate accuracy results for AI and Curator

# Helper function to count species mentions
count_species_accuracy <- function(results_column) {
  results_column <- strsplit(results_column, "; ")
  results_column <- unlist(results_column)
  table(results_column)
}

# Calculate correct and incorrect counts for AI
ai_correct_counts <- count_species_accuracy(accuracy_results$correct_ai_guesses)
ai_incorrect_counts <- count_species_accuracy(accuracy_results$incorrect_ai_guesses)

# Combine AI counts into a single dataframe
ai_species_accuracy <- data.frame(
  Species = names(ai_correct_counts),
  ai_correct = as.numeric(ai_correct_counts),
  ai_incorrect = as.numeric(ai_incorrect_counts[names(ai_correct_counts)])
)

# Replace NA for missing incorrect counts with 0
ai_species_accuracy$ai_incorrect[is.na(ai_species_accuracy$ai_incorrect)] <- 0

# Calculate total assignments and accuracy
ai_species_accuracy <- ai_species_accuracy %>%
  mutate(
    ai_total = ai_correct + ai_incorrect,
    ai_accuracy = ai_correct / ai_total
  )

# Repeat for Curator data
curator_correct_counts <- count_species_accuracy(accuracy_results$correct_curator_guesses)
curator_incorrect_counts <- count_species_accuracy(accuracy_results$incorrect_curator_guesses)

curator_species_accuracy <- data.frame(
  Species = names(curator_correct_counts),
  curator_correct = as.numeric(curator_correct_counts),
  curator_incorrect = as.numeric(curator_incorrect_counts[names(curator_correct_counts)])
)

# Replace NA for missing incorrect counts with 0
curator_species_accuracy$curator_incorrect[is.na(curator_species_accuracy$curator_incorrect)] <- 0

# Calculate total assignments and accuracy
curator_species_accuracy <- curator_species_accuracy %>%
  mutate(
    curator_total = curator_correct + curator_incorrect,
    curator_accuracy = curator_correct / curator_total
  )

# Merge AI and Curator results into one dataframe for visualization
species_accuracy <- full_join(ai_species_accuracy, curator_species_accuracy, by = "Species")


write.csv(species_accuracy, "../output/pathogens_training_set/species_accuracy_aggregated_by_species.csv", row.names = FALSE)

# only include species wiht more than 10 guesses by the AI OR more than 1 by curator
species_accuracy_for_plotting <- species_accuracy[species_accuracy$ai_total >= 10 |
                                                  species_accuracy$curator_total >= 1,  ]

species_accuracy_for_plotting_ai <- species_accuracy_for_plotting[!is.na(species_accuracy_for_plotting$ai_total), ]
# Plotting Species Accuracy for AI
ggplot(species_accuracy_for_plotting_ai, aes(x = reorder(Species, ai_total), y = ai_accuracy)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(ai_correct, "/", ai_total)), 
            hjust = -0.05, size = 3, color = "black") +
  labs(
    title = "Species-Level Accuracy: AI",
    x = "Species",
    y = "Accuracy",
    fill = "Source"
  ) +
  coord_flip() +
  theme_classic()

ggsave(filename = paste0(output_folder, "species_accuracy_ai_morethan10guesses.png"), 
       dpi = 600, height = 20, width = 35, units = "cm")



# Plotting Species Accuracy for Curators
species_accuracy_for_plotting_curators <- species_accuracy_for_plotting[!is.na(species_accuracy_for_plotting$curator_total), ]
ggplot(species_accuracy_for_plotting_curators, aes(x = reorder(Species, curator_total), y = curator_accuracy)) +
  geom_bar(stat = "identity", fill = "salmon", color = "black", alpha = 0.8) +
  geom_text(aes(label = paste0(curator_correct, "/", curator_total)), 
            hjust = -0.2, size = 3, color = "black") +
  labs(
    title = "Species-Level Accuracy: Curators",
    x = "Species",
    y = "Accuracy",
    fill = "Source"
  ) +
  coord_flip() +
  theme_classic()

ggsave(filename = paste0(output_folder, "species_accuracy_curator.png"), 
       dpi = 600, height = 20, width = 35, units = "cm")



################### gene level scoring ###################
# this is where things get complicated. 
# in the model results we have:
# gene ID - often a database ID, sometimes it is the same as the gene name AND symbol,
# sometimes just as the gene symbol

# The symbol is the best place to start if converted to lower case and special 
# characters are removed we should have similar symbols accross different species. 
# we can try scoring here in R using this or using an LLM to score similarly to 
# how this was done for GO terms originally. 

# ensure pmids that have been processed to json are included but no others.
genes_curated_df <- curated_df[curated_df$PMIDs %in% results_df$PMID, ] %>% distinct()

# keep only gene symbol as starting point for symplicity
genes_test_df <- results_df[results_df$PMID %in% genes_curated_df$PMIDs, c("PMID", 
                                                                           "species_latin_name", 
                                                                           "gene_symbol")] %>% distinct()


# convert the gene symbols to lower case but do not remove any special characters at this stage

genes_test_df <- genes_test_df %>%
  mutate(
    gene_symbol_clean = tolower(gene_symbol) 
  )

# as a starting point ignore species and look at pmid-gene pairs. for each 
# cleaned gene symbol use regex to see if it is found as for the same 
# pmid in the curated df under any of the following columns:
cols_to_check <- c("DB_ID",
                   "DB_Symbol", 
                   "Name",      
                   "Synonym")

# keep in mind that $PMID in genes_test_df matches $PMIDs in genes_curated_df
# keep in mind that $Synonym can be empty "", include one sysnonym, or include 
# multiple synonyms sebarated by "|". All of the cols to check will have to be 
# converted to lowercase and stripped of special characters, for Synonym and 
# DB_ID the split will have to be done after separating at "|"

# prepare the curated df 
genes_curated_df_clean <- genes_curated_df %>%
  mutate(across(
    all_of(cols_to_check), 
    ~ tolower(.) 
  ))


# Expand `Synonym` into multiple rows
genes_curated_df_clean <- genes_curated_df_clean %>%
  separate_rows(Synonym, sep = "\\|") # Split Synonym into multiple rows


# Step 3: Perform column-specific matching

# make sure workign with PMIDs present in both Dfs
genes_curated_df_clean <- genes_curated_df_clean[genes_curated_df_clean$PMIDs %in% genes_test_df$PMID, ]
genes_test_df <- genes_test_df[genes_test_df$PMID %in% genes_curated_df_clean$PMIDs, ]

pmids_list <- genes_curated_df$PMIDs %>% unique()


# harmonise PM ID col name
names(genes_curated_df_clean)[8] <- "PMID"

### TODO: needs some final tidying, but works OK for now.
# store results for matches for each clumn that could contain a match
gene_match_hits <- data.frame(
  PMID = character(),
  ai_species = character(),
  curator_species = character(),
  ai_symbol = character(),
  curator_symbol = character(),
  matching_column = character(),
  gene_match_found = logical(),
  stringsAsFactors = FALSE
)

### helper function for regex error evasion 

library(stringr) # for str_replace_all

# Escape special characters in the pattern
safe_grepl <- function(pattern, text, ignore.case = TRUE) {
  # Escape regex characters
  safe_pattern <- str_replace_all(pattern, "([\\[\\]\\(\\)\\{\\}\\^\\$\\|\\.\\+\\?\\*\\\\])", "\\\\\\1")
  
  # Safely perform grepl inside a tryCatch
  tryCatch({
    grepl(safe_pattern, text, ignore.case = ignore.case)
  }, error = function(e) {
    # If regex compilation fails, return FALSE
    warning(paste("Regex error for pattern:", safe_pattern, "-", e$message))
    FALSE
  })
}

# pmids_list_test <- pmids_list[1:100]

pmids_to_process_again <- c()


for (j in 1:length(pmids_list)) {
  # print("pmid:")
  # print(pmid) # for torubleshooting
  # pmid <- "23617571"
  # Get records for the current PMID
  print("Processed pmids:")
  print(paste0(j,"/",length(pmids_list)))
  
  pmid <- pmids_list[j]
  pmid_curated_data <- genes_curated_df_clean[genes_curated_df_clean$PMID == pmid, ]
  pmid_test_gene_data <- genes_test_df[genes_test_df$PMID == pmid, ]
  
  # Extract unique gene-symbol and species pairs for the current PMID from the test dataset
  ai_gene_species_pairs <- pmid_test_gene_data %>%
    select(gene_symbol_clean, species_latin_name) %>%
    distinct()
  
  # Iterate through each AI gene-symbol for the current PMID
  for (i in seq_len(nrow(ai_gene_species_pairs))) {
   
    # i = 1
    # print("i:")
    # print(i) # for troubleshooting
  
    ai_gene_symbol <- ai_gene_species_pairs$gene_symbol_clean[i]
    if (length(ai_gene_symbol) > 100) {
      
      print("Too many genes annotated by AI. Leaving this PMID for later.")
      pmids_to_process_again <- c(pmids_to_process_again, pmid)
      
      next
      
      }
    ai_species <- ai_gene_species_pairs$species_latin_name[i]
 
    # Iterate through each column to check in the curated dataset
    for (col in cols_to_check) {
      # print("col:")
      # print(col) # for troubleshooting
      # col <- "Synonym"
      
      # Perform forward and reverse matching (vectorised using only valid patterns)
      valid_patterns <- pmid_curated_data[[col]][!is.na(pmid_curated_data[[col]]) & pmid_curated_data[[col]] != ""]
      
      # check if any valid patterns remain, if nor proceed on to the next column
      if (length(valid_patterns) == 0) {
        
        next  # Skip to the next column if no valid patterns
        
      } else if (length(valid_patterns) >= 100) {
        
        print("Too many genes provided by curator. Leaving this PMID for later.")
        pmids_to_process_again <- c(pmids_to_process_again, pmid)
        
        next
        
      }
      # Perform regex matching (forward) with valid patterns
      sym_to_col_match <- sapply(valid_patterns, function(pattern) {
        safe_grepl(ai_gene_symbol, pattern, ignore.case = TRUE)
      })
      
      # Perform reverse matching only on valid patterns
      reverse_match <- sapply(valid_patterns, function(pattern) {
        safe_grepl(pattern, ai_gene_symbol, ignore.case = TRUE)
      })
      
      # Combine results
      match_results <- sym_to_col_match | reverse_match
      
      # If a match is found, record the species and other details
      if (any(match_results, na.rm = TRUE)) {
        
        matching_rows <- which(match_results)
        
        for (row in matching_rows) {
          # row <- 1
          # print("row:")
          # print(row) # for troubleshooting
          gene_match_hits <- rbind(gene_match_hits, data.frame(
            PMID = pmid,
            ai_species = ai_species,
            curator_species = pmid_curated_data$Species[row],
            ai_symbol = ai_gene_symbol,
            curator_symbol = pmid_curated_data[[col]][row],
            matching_column = col,
            gene_match_found = TRUE,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        
        gene_match_hits <- rbind(gene_match_hits, data.frame(
          PMID = pmid,
          ai_species = ai_species,
          curator_species = "No match for gene",
          ai_symbol = ai_gene_symbol,
          curator_symbol = "No match",
          matching_column = "No match",
          gene_match_found = FALSE,
          stringsAsFactors = FALSE
        ))
        
      }
    }
  }
}


ai_genes_with_matches <- gene_match_hits[gene_match_hits$gene_match_found == TRUE, ]
names(ai_genes_with_matches)
# collapse matching_column by separating different entries with ; if all esle is the same

ai_genes_with_collapsed_matches <- ai_genes_with_matches %>%
  group_by(PMID, ai_species, curator_species, ai_symbol, curator_symbol) %>%
  summarize(
    matching_column = paste(unique(matching_column), collapse = ";"),
    .groups = "drop"  # Ungroup after summarization
  )

write.csv(ai_genes_with_collapsed_matches, 
          "../output/pathogens_training_set/ai_genes_with_matches.csv", 
          row.names = FALSE)
# calculate recall
# use the total number of genes in curated dataset for a pmid as per genes_curated_df_clean
# number matched as per ai_genes_with_collapsed_matches

head(genes_curated_df_clean)


head(ai_genes_with_collapsed_matches)


# Step 1: Count total curated genes per PMID
curated_gene_count <- genes_curated_df_clean %>%
  group_by(PMID) %>%
  summarize(total_curated_genes = n(), .groups = "drop")  # Number of rows per PMID

# Step 2: Count AI-matched genes per PMID
ai_gene_match_count <- ai_genes_with_collapsed_matches %>%
  group_by(PMID) %>%
  summarize(matched_genes = n_distinct(ai_symbol), .groups = "drop")

ai_genes_with_collapsed_matches[ai_genes_with_collapsed_matches$PMID == "24533860",]


# Step 3: Join and calculate recall
recall_df <- curated_gene_count %>%
  left_join(ai_gene_match_count, by = "PMID") %>%
  mutate(
    matched_genes = ifelse(is.na(matched_genes), 0, matched_genes), # Handle PMIDs with no matches
  ) 


# set anything greater than 1 to 1, it is because different gene notaitons for multiple 
# species were matched
# Calculate recall score
recall_df <- recall_df %>%
  mutate(
    recall_score_species_not_accounted_for = matched_genes / total_curated_genes,
    recall_score_species_not_accounted_for = pmin(recall_score_species_not_accounted_for, 1)  # Cap the recall score at 1
  )
write.csv(recall_df, "../output/pathogens_training_set/genes_recall_bypmid.csv", 
          row.names = FALSE)
# lets take species into account for the recall - we don't want to count towards recall genes
# assigned to incorrect species

ai_genes_with_collapsed_matches$ai_species_standardised <- sapply(ai_genes_with_collapsed_matches$ai_species, standardise_species_name)
ai_genes_with_collapsed_matches$curator_species_standardised <- sapply(ai_genes_with_collapsed_matches$curator_species, standardise_species_name)


ai_genes_with_collapsed_matches_filtered_for_species_match <- ai_genes_with_collapsed_matches[ai_genes_with_collapsed_matches$ai_species_standardised == ai_genes_with_collapsed_matches$curator_species_standardised, ]


# Step 2: Count AI-matched genes per PMID with matching species
ai_gene_match_count <- ai_genes_with_collapsed_matches_filtered_for_species_match %>%
  group_by(PMID) %>%
  summarize(matched_genes = n_distinct(ai_symbol), .groups = "drop")




# Step 3: Join and calculate recall
recall_df_with_species_accounted_for <- curated_gene_count %>%
  left_join(ai_gene_match_count, by = "PMID") %>%
  mutate(
    matched_genes = ifelse(is.na(matched_genes), 0, matched_genes), # Handle PMIDs with no matches
  ) 


# set anything greater than 1 to 1, it is because different gene notaitons for multiple 
# species were matched
# Calculate recall score
recall_df_with_species_accounted_for <- recall_df_with_species_accounted_for %>%
  mutate(
    recall_score_species_accounted_for = matched_genes / total_curated_genes,
    recall_score_species_accounted_for = pmin(recall_score_species_accounted_for, 1)  # Cap the recall score at 1
  )


write.csv(recall_df_with_species_accounted_for, "../output/pathogens_training_set/genes_recall_bypmid_correct_species_only.csv", 
          row.names = FALSE)
mean(recall_df_with_species_accounted_for$recall_score_species_accounted_for)
mean(recall_df$recall_score_species_not_accounted_for)


### NB this is overestomating scores as sometimes it is not the exact gene matching using our algo. 
## really need to work on matching using AI. 