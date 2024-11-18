library(tidyverse)
# read in extracted data for each pmid in extended set
pmid_data <- read.csv("../data/extended_set_pmid_data.csv")

# tidy publication date 
pmid_data$Publication.Date <- pmid_data$Publication.Date
pmid_data$publication_year <- sub(".*'Year': '([0-9]{4}).*", 
                                  "\\1", 
                                  pmid_data$Publication.Date) %>% as.numeric()


# read in list of pmids that failed
failed_pmids <- readLines("../data/failed_pmids_v3.txt")

# if pmid in pmid_data is in the list of failed, mark it as such
pmid_data$failed <- pmid_data$PMID %in% failed_pmids


# Get the top 10 journals by count
top_journals <- pmid_data %>%
  count(Journal) %>%
  top_n(10, wt = n) %>%
  pull(Journal)

# Filter data to include only top journals
filtered_data <- pmid_data %>%
  filter(Journal %in% top_journals)

# Plot the relationship between Journal and Failure Status
ggplot(filtered_data, aes(x = Journal, fill = failed)) +
  geom_bar(position = "dodge") +
  labs(title = "Failure Status by Journal",
       x = "Journal",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

names(pmid_data)


library(ggplot2)
# Plot the distribution of publication years by failure status
ggplot(pmid_data, aes(x = as.factor(failed), y = publication_year, fill = as.factor(failed))) +
  geom_violin(alpha = 0.6) + 
  geom_boxplot(width = 0.2, color = "black", alpha = 0.4) +  # Adjusted width for better overlay
  labs(title = "Publication Year vs. Failure Status",
       x = "Failed",
       y = "Publication Year") +
  scale_y_continuous(breaks = seq(min(pmid_data$publication_year, na.rm = TRUE),
                                  max(pmid_data$publication_year, na.rm = TRUE), by = 5)) +
  theme_minimal() +
  theme(legend.position = "none")  # Hide legend if failed status is self-evident



# Convert Free.Text.Availability and failed to factors 
pmid_data$Free.Text.Availability <- as.factor(pmid_data$Free.Text.Availability)
pmid_data$failed <- as.factor(pmid_data$failed)

# Plot the relationship between Free Text Availability and Failure Status
ggplot(pmid_data, aes(x = Free.Text.Availability, fill = failed)) +
  geom_bar(position = "dodge") +
  labs(title = "Free Text Availability vs. Failure Status",
       x = "Free Text Availability",
       y = "Count") +
  theme_minimal()


# Plot the relationship between Journal and Failure Status
library(dplyr)

# Summarize to find journals with only failed or only non-failed articles
journal_status_summary <- pmid_data %>%
  group_by(Journal) %>%
  summarize(all_failed = all(failed == TRUE),
            all_non_failed = all(failed == FALSE),
            total_articles = n())

# Filter for journals that are either only failed or only non-failed
only_failed_journals <- journal_status_summary %>% filter(all_failed == TRUE)
only_non_failed_journals <- journal_status_summary %>% filter(all_non_failed == TRUE)

print("Journals with only failed articles:")
print(only_failed_journals)

print("Journals with only non-failed articles:")
print(only_non_failed_journals)

table(pmid_data$Journal)
library(ggplot2)

# Plot the count of failed vs non-failed articles for each journal
ggplot(pmid_data, aes(x = Journal, fill = as.factor(failed))) +
  geom_bar() +
  labs(title = "Failure Status by Journal",
       x = "Journal",
       y = "Count of Articles",
       fill = "Failed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"))



library(dplyr)

# Ensure `failed` is logical (TRUE/FALSE) or convert it to numeric
filtered_data <- pmid_data %>%
  filter(publication_year > 2008) %>%
  mutate(failed = as.numeric(as.logical(failed)))

# Group by Journal and calculate counts
journal_summary <- filtered_data %>%
  group_by(Journal) %>%
  summarize(
    total_articles = n(),
    failed_count = sum(failed),
    non_failed_count = total_articles - failed_count
  ) %>%
  # Filter for journals with more than 5 articles
  filter(total_articles > 5)

# Identify journals with only failed or only non-failed publications
only_failed_journals <- journal_summary %>% filter(failed_count == total_articles)
only_non_failed_journals <- journal_summary %>% filter(non_failed_count == total_articles)

# Display results
print("Journals with only failed publications:")
print(only_failed_journals)

print("Journals with only non-failed publications:")
print(only_non_failed_journals)


library(dplyr)
library(ggplot2)

# Calculate the percentage of failed and successful publications by journal
journal_summary <- filtered_data %>%
  group_by(Journal) %>%
  summarize(
    total_articles = n(),
    failed_count = sum(failed),
    success_count = total_articles - failed_count,
    failed_percent = (failed_count / total_articles) * 100,
    success_percent = (success_count / total_articles) * 100
  ) %>%
  # Filter for journals with more than 5 articles
  filter(total_articles > 5)

# Reshape data for plotting
plot_data <- journal_summary %>%
  select(Journal, failed_percent, success_percent) %>%
  pivot_longer(cols = c(failed_percent, success_percent),
               names_to = "status",
               values_to = "percentage") %>%
  mutate(percentage = ifelse(status == "failed_percent", -percentage, percentage),
         status = ifelse(status == "failed_percent", "Failed", "Successful"))

# Plot
# Reshape data for plotting
plot_data <- journal_summary %>%
  select(Journal, failed_count, success_count, failed_percent, success_percent) %>%
  pivot_longer(cols = c(failed_count, success_count, failed_percent, success_percent),
               names_to = c("status", ".value"),
               names_pattern = "(failed|success)_(.*)") %>%
  mutate(status = ifelse(status == "failed", "Failed", "Successful"))

# Plot
ggplot(plot_data, aes(x = reorder(Journal, percent), y = status, size = count, fill = status)) +
  geom_point(shape = 21, alpha = 0.7) +
  geom_text(aes(label = count), color = "white", fontface = "bold", size = 3) +
  coord_flip() +
  scale_size(range = c(5, 20)) +  # Adjust bubble size range
  labs(title = "Percentage of Failed vs. Successful Publications by Journal",
       x = "Journal",
       y = "Publication Status",
       size = "Count",
       fill = "Publication Status") +
  scale_fill_manual(values = c("Failed" = "blue", "Successful" = "orange")) +
  theme_minimal()


library(dplyr)
library(ggplot2)

# Calculate the percentage of failed and successful publications by journal
journal_summary <- filtered_data %>%
  group_by(Journal) %>%
  summarize(
    total_articles = n(),
    failed_count = sum(failed),
    success_count = total_articles - failed_count,
    failed_percent = (failed_count / total_articles) * 100,
    success_percent = (success_count / total_articles) * 100
  ) %>%
  # Filter for journals with more than 5 articles
  filter(total_articles > 5)

# Reshape data for plotting
plot_data <- journal_summary %>%
  select(Journal, failed_percent, success_percent, total_articles) %>%
  pivot_longer(cols = c(failed_percent, success_percent),
               names_to = "status",
               values_to = "percentage") %>%
  mutate(percentage = ifelse(status == "failed_percent", -percentage, percentage),
         status = ifelse(status == "failed_percent", "Failed", "Successful"))

# Plot with total articles as text
ggplot(plot_data, aes(x = reorder(Journal, -percentage), y = percentage)) +
  geom_bar(aes(fill = status), stat = "identity", position = "identity") +
  geom_text(data = journal_summary, aes(x = Journal, y = 0, label = total_articles),
            vjust = -1, size = 3.5, color = "black") +
  scale_y_continuous(labels = abs) +  # Display y-axis labels as positive values
  labs(title = "Percentage of Failed vs. Successful Publications after 2008 by Journal",
       x = "Journal",
       y = "Percentage",
       fill = "Publication Status") +
  scale_fill_manual(values = c("Failed" = "blue", "Successful" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


