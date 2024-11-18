library(tidyverse)
## looking at the varierty of PMID substitutes we are getting

GAF_merged <- read.delim("../output/merged_gaf.tsv", sep = "\t", header = TRUE)
head(GAF_merged)

Total_IDs <- GAF_merged$X5 %>% length()
Unique_IDs <- GAF_merged$X5 %>% unique() %>% length()

PMIDs_total <- grepl("PMID:", GAF_merged$X5) %>% sum()
PMIDS_unique <- GAF_merged$X5[grepl("PMID:", GAF_merged$X5)] %>% unique() %>% length()





GO_REFs_total <- grepl("GO_REF:", GAF_merged$X5) %>% sum()
GO_REFS_unique <- GAF_merged$X5[grepl("GO_REF:", GAF_merged$X5)] %>% unique() %>% length()


ENSs_total <- grepl("VEuPathDB:ENS", GAF_merged$X5) %>% sum()
ENSS_unique <- GAF_merged$X5[grepl("VEuPathDB:ENS", GAF_merged$X5)] %>% unique() %>% length()

# there's other VEuPathDB: as well. not just ENS; break down wat types
# are there any cases I'm missing? ENS, GO_REF and PMID should cover all rows
PMIDs_total + GO_REFs_total + ENSs_total == Total_IDs

# get the ones that aren't of the types we have counted above

other_refs <- GAF_merged$X5[!grepl("TrypTag:|MID:|DawsonLab:|PMID:|GO_REF:|VEuPathDB:|AspGD_REF:|CGD_REF|NCBI:GEO:|SGD_REF:", GAF_merged$X5)]
# grepl("VEuPathDB:", other_refs) %>% sum()
head(other_refs[other_refs != ""])

norrefd <- GAF_merged[GAF_merged$X5=="",]]

species_with_PMIDs <- table(GAF_merged$Species[grepl("PMID:", GAF_merged$X5)])
GAF_merged$Species[!grepl("PMID:", GAF_merged$X5)] %>% unique()

species_without_PMIDs <- table(GAF_merged$Species[!grepl("PMID:", GAF_merged$X5)])
