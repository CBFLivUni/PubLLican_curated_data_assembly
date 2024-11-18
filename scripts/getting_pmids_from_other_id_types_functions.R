library(rentrez)
library(httr)
library(jsonlite)


# DOI to PMID
get_pmid_from_doi <- function(doi) {
  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  response <- GET(base_url, query = list(db = "pubmed", term = doi, retmode = "json"))
  content <- content(response, as = "parsed", type = "application/json")
  if (length(content$esearchresult$idlist) > 0) {
    return(content$esearchresult$idlist[[1]])  # Return the PMID
  } else {
    return(NA)
  }
}


# Reactome to PMID
get_pmids_from_reactome <- function(reactome_id) {
  url <- paste0("https://reactome.org/content/detail/", reactome_id)
  page <- readLines(url, warn = FALSE)
  pmid_matches <- regmatches(page, gregexpr("PMID:[0-9]+", page))
  pmids <- unlist(pmid_matches)
  if (length(pmids) > 0) {
    return(unique(pmids))
  } else {
    return(NA)
  }
}



library(httr)
library(rvest)
library(xml2)
# install.packages("RSelenium")
library(RSelenium)

get_pmid_from_mgi <- function(mgi_id= "MGI:64874") {
  # Construct the search URL for MGI
  base_url <- "https://www.informatics.jax.org/quicksearch/summary"
  query_url <- paste0(base_url, "?queryType=exactPhrase&query=", mgi_id, "&submit=Quick+Search")
  
  # Start a Selenium server and browser instance
  rD <- rsDriver(browser = "firefox", port = 4545L)
  remDr <- rD[["client"]]
  remDr$navigate(query_url)
  # Try to read the page content
  page <- tryCatch(read_html(query_url), error = function(e) return(NA))
  
  
  # Print the entire page content to inspect if it includes reference text
  page_text <- page %>% html_text()
  print(page_text)
  
  # If page reading fails, return NA
  if (is.na(page)) return(NA)
  # TODO: this is not workign because we need to install java, the table is
  # dynamically loaded. don't fix for now as it's training data anyway, might need to 
  # for future applicaitons! 
  # Extract the table with id 'b3Table' and get the text in td elements
  ref_text <- page %>%
    xml_find_first("//table[@id='b3Table']/tbody/tr[2]/td[2]") %>%
    xml_text(trim = TRUE)
  
  print(ref_text)
  
  # Combine all text from the td elements to search for DOI or other identifiers
  ref_text_combined <- paste(ref_text, collapse = " ")
  
  # Define a pattern for DOI if it is present in the text
  doi_pattern <- "DOI:[0-9.]+/[a-zA-Z0-9.]+"
  doi_match <- regmatches(ref_text_combined, regexpr(doi_pattern, ref_text_combined))
  
  if (length(doi_match) > 0) {
    return(gsub("DOI:", "", doi_match))  # Return DOI without "DOI:" prefix
  } else {
    # No DOI found, return reference text or NA if empty
    if (nchar(ref_text_combined) > 0) {
      return(ref_text_combined)
    } else {
      return(NA)
    }
  }
}

test <- get_pmid_from_mgi("MGI%64874")

# DB to PMID
get_pmids_from_db <- function(db_id, db) {
  # Define base URLs for each database
  base_urls <- list(
    MGI = "http://www.informatics.jax.org/marker/",
    ZFIN = "https://zfin.org/",
    RGD = "https://rgd.mcw.edu/rgdweb/search/search.html?term=",
    FB = "https://flybase.org/reports/"
  )
  
  # Construct the URL for the database and ID
  url <- paste0(base_urls[[db]], db_id)
  
  # Attempt to read the page content
  page <- tryCatch(readLines(url, warn = FALSE), error = function(e) return(character(0)))
  
  # Check if page is empty
  if (length(page) == 0) return(NA)
  Sys.sleep(1)
  # Search for PMIDs in the page content
  pmid_matches <- regmatches(page, gregexpr("PMID:[0-9]+", page))
  pmids <- unlist(pmid_matches)
  
  # Return unique PMIDs if found, otherwise NA
  if (length(pmids) > 0) {
    return(unique(pmids))
  } else {
    return(NA)
  }
}


# wraper for above
convert_to_pmid <- function(ref, ref_type) {
  if (ref_type == "DOI") {
    return(get_pmid_from_doi(ref))
  } else if (ref_type == "Reactome" || ref_type == "REACTOME") {
    return(get_pmids_from_reactome(ref))
  } else if (ref_type %in% c("MGI", "ZFIN", "RGD", "FB")) {
    return(get_pmids_from_db(ref, ref_type))
  } else {
    return(NA)  # Return NA if no conversion is possible
  }
}

