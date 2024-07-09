import os


# from the set of curated data, get a list of all unique PMIDs that contributed to the curation
# this list will be used to run Noah's pipeline using different LLMs to scorre their performance
# and to test how fine-tuning impacts this (how similar are the GO terms annotations of the LLMs)

PMIDs_data_dir = os.path.join("..", "output")


# read in