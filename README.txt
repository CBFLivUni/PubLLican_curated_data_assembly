Initially GAF files were downloaded from VEuPAthDB and its associated databases (e.g. Fungi DB). 

Due to missing PMID annotations, GAF files were subsequently downloaded from the GO website.
Here are the links:

https://geneontology.org/docs/download-go-annotations/

https://current.geneontology.org/products/pages/downloads.html (species selected by the GO consortium)

https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/ (all species)


The pipeline was first applied to species selected by the GO consortium (as a sort of control to see how it operates 
on well annotated species.
 
A training dataset of papers was built using all species on the list (24) by selecting the top 10 papers
contributing the largest number of GO annotations for each species. 

The training dataset was subsequently expanded by adding the top 260 papers contributing the largest number 
of annotations to GAF files from VEuPAthDB that were not already included in the model species set (e.g. P.falciparum
is a model species). The number of papers per species was limited to a maximum of 25 to avoid overbiasing the training set to a certain species/gene symbol annotation. 

Gene symbols were subsequently converted to a common notation in R. 