######### read in #########
## read in the curated data
rare_for_scoring <- read.csv("../output/extended_training_set_rare_species.csv")
model_for_scoring <- read.csv("../output/extended_training_set_model_species.csv")

## subset to include only column of interest: pmid, taxon id, species, 
## gene symbol, gene name
rare_for_scoring <- rare_for_scoring[, c("PMIDs", "primary_taxon_id", "Species", "DB_Symbol", "Name", "Synonym", "Evidence_Code")]
names(rare_for_scoring) <- c("PMID", "taxon_id", "Species", "Gene_Symbol", "Gene_Name", "Gene_Synonym", "Evidence_Code")


model_for_scoring <- model_for_scoring[, c("PMIDs", "primary_taxon_id", "Primary_Species", "DB_Symbol", "Name", "Synonym", "Evidence_Code")]
names(model_for_scoring) <- c("PMID", "taxon_id", "Species", "Gene_Symbol", "Gene_Name", "Gene_Synonym", "Evidence_Code")


## read in list of selected pmids and list of failed

selected_for_training_set <- readLines("../data/all_pmids_to_process.txt")
failed <- readLines("../data/failed_pmids_v3.txt")
successfull <- setdiff(selected_for_training_set, failed)



# see how many failed/successful for each species
rare_failed <- rare_for_scoring[rare_for_scoring$PMID %in% failed,]
rare_succes <- rare_for_scoring[rare_for_scoring$PMID %in% successfull,]


model_failed <- model_for_scoring[model_for_scoring$PMID %in% failed,]
model_succes <- model_for_scoring[model_for_scoring$PMID %in% successfull,]

# have any species been affected disproportionately? 

# Combine rare and model data
all_data <- rbind(rare_for_scoring, model_for_scoring)
# retain only species and pmid cols
all_data <- all_data[, c("PMID", "Species")]# Subset to only include PMIDs from the selected training set

all_selected <- all_data[all_data$PMID %in% selected_for_training_set, ]  %>% distinct()

# Create subsets for failed and successful PMIDs by species
all_failed <- all_selected[all_selected$PMID %in% failed, ] %>% distinct()
all_success <- all_selected[all_selected$PMID %in% successfull, ] %>% distinct()

# Calculate total PMIDs per species
total_counts <- as.data.frame(table(all_selected$Species))
names(total_counts) <- c("Species", "Total_PMIDs")

# Calculate failed PMIDs per species
failed_counts <- as.data.frame(table(all_failed$Species))
names(failed_counts) <- c("Species", "Failed_PMIDs")

# Calculate successful PMIDs per species
success_counts <- as.data.frame(table(all_success$Species))
names(success_counts) <- c("Species", "Success_PMIDs")

# Merge the data frames
species_stats <- merge(total_counts, failed_counts, by = "Species", all = TRUE)
species_stats <- merge(species_stats, success_counts, by = "Species", all = TRUE)

# Replace NAs with zeros
species_stats[is.na(species_stats)] <- 0

# Calculate failure rate per species
species_stats$Failure_Rate <- species_stats$Failed_PMIDs / species_stats$Total_PMIDs

# View species with highest failure rates
species_stats <- species_stats[order(-species_stats$Failure_Rate), ]
print(species_stats)

# Optionally, visualize the failure rates
library(ggplot2)
ggplot(species_stats, aes(x = reorder(Species, -Failure_Rate), y = Failure_Rate)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Failure Rate by Species", x = "Species", y = "Failure Rate")


# Plot the number of successful PMIDs by species


# Sort species by the number of successful PMIDs
species_stats_success <- species_stats[order(species_stats$Success_PMIDs, decreasing = TRUE), ]

ggplot(species_stats_success, aes(x = reorder(Species, Success_PMIDs), y = Success_PMIDs)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Number of Successful PMIDs by Species",
       x = "Species",
       y = "Number of Successful PMIDs") +
  theme_minimal()





# Reshape the data for plotting
library(reshape2)

# Prepare data for plotting
species_melted <- melt(species_stats[, c("Species", "Success_PMIDs", "Failed_PMIDs")],
                       id.vars = "Species",
                       variable.name = "Status",
                       value.name = "PMID_Count")

# Rename levels for clarity
species_melted$Status <- factor(species_melted$Status,
                                levels = c("Success_PMIDs", "Failed_PMIDs"),
                                labels = c("Successful PMIDs", "Failed PMIDs"))

# Define color-blind friendly colors
cb_palette <- c("Successful PMIDs" = "dodgerblue",  # Bluish green
                "Failed PMIDs"     = "firebrick2")  # Vermilion

# Plot successful and failed PMIDs by species with color-blind friendly colors
ggplot(species_melted, aes(x = reorder(Species, PMID_Count), y = PMID_Count, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = cb_palette) +
  labs(title = "Successful vs. Failed PMIDs by Species",
       x = "Species",
       y = "PMID Count",
       fill = "Status") +
  theme_minimal()


library(ggrepel)


# Scatter plot of Failure Rate vs. Number of Successful PMIDs with jitter and alpha
ggplot(species_stats, aes(x = Success_PMIDs, y = Failure_Rate)) +
  geom_jitter(size = 2.5, color = "darkred", width = 0.2, height = 0.02, alpha = 0.7) +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = Inf) +
  labs(title = "Failure Rate vs. Number of Successful PMIDs",
       x = "Number of Successful PMIDs",
       y = "Failure Rate") +
  theme_minimal()

ggplot(species_stats, aes(x = Success_PMIDs, y = Failure_Rate)) +
  geom_jitter(size = 2.5, color = "darkred", width = 0.2, height = 0.02, alpha = 0.5) +
  geom_text_repel(aes(label = Species), size = 3) +
  labs(title = "Failure Rate vs. Number of Successful PMIDs",
       x = "Number of Successful PMIDs",
       y = "Failure Rate") +
  theme_minimal()


