######## libs and read in #########
library(tidyverse)
library(rentrez) # for converting taxon id to species name
library(httr) # for getting 
library(jsonlite)

## a large fraction of pmids cannot be easily obtained as JSON from pubmed, 
# so lets look at optimising the prompting on species thyat fiot our group's niche of
# working on DBs with parasites and other pathogens

GAF_rare_withpmid <- read.csv("../output/training_set_rare_species_with_pmids_and_species.csv")

GAF_model_organisms_final <- read.csv("../output/model_species_with_pmids_and_species.csv")

rare_species <- GAF_rare_withpmid$Species %>% unique()
# > rare_species
# [1] "A.fumigatus"       "A.nidulans"        "C.neoformans"      "F.graminearum"    
# [5] "N.crassa"          "P.oryzae"          "P.berghei"         "Pchabaudichabaudi"
# [9] "P.gallinaceum"     "P.knowlesi"        "P.malariae"        "P.ovalecurtisi"   
# [13] "P.reichenowi"      "P.relictum"        "P.vivax"           "P.yoeliiyoelii"   
# [17] "T.gondii"          "L.braziliensis"    "L.donovani"        "L.infantum"       
# [21] "L.mexicana"        "T.bruceigambiense" "T.congolense"      "T.cruzi"          
# [25] "T.vivax" 

model_species <- GAF_model_organisms_final$Primary_Species %>% unique()
# > model_species
# [1] "Arabidopsis thaliana"                                  "Bos taurus"                                            "Caenorhabditis elegans"                               
# [4] "Candida albicans SC5314"                               "Nakaseomyces glabratus CBS 138"                        "Candidozyma auris"                                    
# [7] "Candida dubliniensis CD36"                             "Candida parapsilosis CDC317"                           "Canis lupus familiaris"                               
# [10] "Danio rerio"                                           "Dictyostelium discoideum"                              "Drosophila melanogaster"                              
# [13] "Escherichia coli K-12"                                 "Gallus gallus"                                         "Homo sapiens"                                         
# [16] "Leishmania major strain Friedlin"                      "Mus musculus"                                          "Plasmodium falciparum 3D7"                            
# [19] "Pseudomonas aeruginosa PAO1"                           "Rattus norvegicus"                                     "Rotavirus A"                                          
# [22] "Corynephage beta"                                      "Clostridium tetani"                                    "Saccharomyces cerevisiae"                             
# [25] "Severe acute respiratory syndrome-related coronavirus" "Clostridium botulinum"                                 "Oryctolagus cuniculus"                                
# [28] "Sus scrofa"                                            "Bacillus anthracis"                                    "Saccharomyces cerevisiae S288C"                       
# [31] "Schizosaccharomyces japonicus"                         "Schizosaccharomyces pombe"                             "Solanum lycopersicum"                                 
# [34] "Trypanosoma brucei brucei TREU927"                     "Xenopus laevis"                                        "Xenopus tropicalis"




# pick out the ones that fall under: 
# eukaryotic pathogens
# fungal pathogens 
# insect vectors
# non-pathogen

# first classify each species by googling
rare_species_categories <- c(
  "A.fumigatus" = "fungal pathogen", 
  "A.nidulans" = "fungal pathogen", 
  "C.neoformans" = "fungal pathogen", 
  "F.graminearum" = "fungal pathogen", 
  "N.crassa" = "non-pathogen", 
  "P.oryzae" = "fungal pathogen", 
  "P.berghei" = "parasite", 
  "Pchabaudichabaudi" = "parasite", 
  "P.gallinaceum" = "parasite", 
  "P.knowlesi" = "parasite", 
  "P.malariae" = "parasite", 
  "P.ovalecurtisi" = "parasite", 
  "P.reichenowi" = "parasite", 
  "P.relictum" = "parasite", 
  "P.vivax" = "parasite", 
  "P.yoeliiyoelii" = "parasite", 
  "T.gondii" = "eukaryotic pathogen", # is also a parasite
  "L.braziliensis" = "eukaryotic pathogen", # parasite
  "L.donovani" = "eukaryotic pathogen", #parasite
  "L.infantum" = "eukaryotic pathogen", # parasite
  "L.mexicana" = "eukaryotic pathogen", # parasite
  "T.bruceigambiense" = "eukaryotic pathogen", # parasite
  "T.congolense" = "eukaryotic pathogen", # parasite
  "T.cruzi" = "eukaryotic pathogen", # parasite
  "T.vivax" = "eukaryotic pathogen" # parasite
)


model_species_categories <- c(
  "Arabidopsis thaliana" = "non-pathogen",
  "Bos taurus" = "non-pathogen",
  "Caenorhabditis elegans" = "non-pathogen",
  "Candida albicans SC5314" = "fungal pathogen",
  "Nakaseomyces glabratus CBS 138" = "fungal pathogen",
  "Candidozyma auris" = "fungal pathogen",
  "Candida dubliniensis CD36" = "fungal pathogen",
  "Candida parapsilosis CDC317" = "fungal pathogen",
  "Canis lupus familiaris" = "non-pathogen",
  "Danio rerio" = "non-pathogen",
  "Dictyostelium discoideum" = "non-pathogen",
  "Drosophila melanogaster" = "non-pathogen",
  "Escherichia coli K-12" = "bacterial pathogen", # although technically not pathogenic? 
  "Gallus gallus" = "non-pathogen",
  "Homo sapiens" = "non-pathogen",
  "Leishmania major strain Friedlin" = "eukaryotic pathogen",
  "Mus musculus" = "non-pathogen",
  "Plasmodium falciparum 3D7" = "eukaryotic pathogen",
  "Pseudomonas aeruginosa PAO1" = "bacterial pathogen",
  "Rattus norvegicus" = "non-pathogen",
  "Rotavirus A" = "viral pathogen",
  "Corynephage beta" = "bacterial pathogen",
  "Clostridium tetani" = "bacterial pathogen",
  "Saccharomyces cerevisiae" = "non-pathogen",
  "Severe acute respiratory syndrome-related coronavirus" = "viral pathogen",
  "Clostridium botulinum" = "bacterial pathogen",
  "Oryctolagus cuniculus" = "non-pathogen",
  "Sus scrofa" = "non-pathogen",
  "Bacillus anthracis" = "bacterial pathogen",
  "Saccharomyces cerevisiae S288C" = "non-pathogen",
  "Schizosaccharomyces japonicus" = "non-pathogen",
  "Schizosaccharomyces pombe" = "non-pathogen",
  "Solanum lycopersicum" = "non-pathogen",
  "Trypanosoma brucei brucei TREU927" = "eukaryotic pathogen",
  "Xenopus laevis" = "non-pathogen",
  "Xenopus tropicalis" = "non-pathogen"
)

# then filter the GAF dfs to only include ones that are 
# eucaryotic, viral, or fungal pathogens,parasites or insect vectors
categories_of_interest <- c(
  "eukaryotic pathogen",
  "fungal pathogen",
  "viral pathogen",
  "parasite",
  "insect vector"
)


# Add the 'Category' column to rare species
GAF_rare_pathogens <- GAF_rare_withpmid %>%
  mutate(Category = rare_species_categories[as.character(Species)]) %>%
  filter(Category %in% categories_of_interest)

# Add the 'Category' column to model species
GAF_model_organisms_pathogens <- GAF_model_organisms_final %>%
  mutate(Category = model_species_categories[as.character(Primary_Species)]) %>%
  filter(Category %in% categories_of_interest)


# include all papers for rare pathogens
rare_pathogen_pmids <- GAF_rare_pathogens$PMIDs %>% unique()

## for model ones apply a bit more of a selection  
# if more than 10 papers per speces, include 20% of papers per species chosen at random
# otherwise include all up to 10 papers
pmid_counts <- GAF_model_organisms_pathogens %>%
  group_by(Primary_Species, PMIDs) %>%
  summarise(row_count = n(), .groups = "drop") %>%
  arrange(Primary_Species, desc(row_count))

set.seed(99)

# Apply the selection logic
selected_model_pmids <- pmid_counts %>%
  group_by(Primary_Species) %>%
  mutate(
    total_pmids = n(),
    include_all = total_pmids <= 20,  # If <= 10 PMIDs, include all
    sampled_pmids = ifelse(
      include_all, 
      PMIDs, 
      sample(PMIDs, size = min(ceiling(0.2 * total_pmids), total_pmids), replace = FALSE) # Ensure valid sample size
    )
  ) %>%
  ungroup() %>%
  filter(include_all | PMIDs %in% sampled_pmids)

GAF_model_organisms_pathogens <- GAF_model_organisms_pathogens[, 2:25]
# merge the rare and model species data together
names(GAF_model_organisms_pathogens)
names(GAF_rare_pathogens)

GAF_model_organisms_pathogens <- GAF_model_organisms_pathogens %>%
  rename(Species = Primary_Species)

# Select and order columns for both datasets
GAF_rare_pathogens <- GAF_rare_pathogens %>%
  select(DB, DB_ID, DB_Symbol, Qualifier, GO_ID, DB_Reference, Evidence_Code, With_From, Aspect, Name, 
         Synonym, Type, taxon, Date, Assigned_By, Annotation_Extension, Gene_Product_Form_ID, 
         Species, Category, PMIDs)

GAF_model_organisms_pathogens <- GAF_model_organisms_pathogens %>%
  select(DB, DB_ID, DB_Symbol, Qualifier, GO_ID, DB_Reference, Evidence_Code, With_From, Aspect, Name, 
         Synonym, Type, taxon, Date, Assigned_By, Annotation_Extension, Gene_Product_Form_ID, 
         Species, Category, PMIDs)

# Merge the datasets
merged_pathogens <- bind_rows(GAF_rare_pathogens, GAF_model_organisms_pathogens)
# Save the merged dataset if needed
write.csv(merged_pathogens, "../output/pathogens_training_set.csv", row.names = FALSE)

pathogen_pmids <- merged_pathogens$PMIDs %>% unique()

# save as txt
write.table(pathogen_pmids, "../output/pathogen_training_set_pmids.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
