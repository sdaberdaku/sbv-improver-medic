# prepare.datasets.R
# prepare the dataset and save into an RData file for easy loading

library(edgeR)

rm(list = ls())
class.labels.He <- read.delim("original_datasets/Class_labels_He.txt", header = TRUE, sep = "\t", colClasses = c("integer", "character", "factor"))
class.labels.Schirmer <- read.delim("original_datasets/Class_labels_Schirmer.txt", header = TRUE, sep = "\t", colClasses = c("integer", "character", "factor"))

taxonomy.abundance.He <- read.delim("original_datasets/TrainingHe_TaxonomyAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "TaxID")
taxonomy.abundance.Schirmer <- read.delim("original_datasets/TrainingSchirmer_TaxonomyAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "TaxID")
taxonomy.abundance.PMI <- read.delim("original_datasets/TestingDataset_TaxonomyAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "TaxID")

training.set.taxonomy.abundance <- as.data.frame(t(cbind(taxonomy.abundance.He, taxonomy.abundance.Schirmer)))
training.set.class.labels <- rbind(cbind(class.labels.He, dataset="He"), cbind(class.labels.Schirmer, dataset="Schirmer"))
test.set.taxonomy.abundance <- as.data.frame(t(taxonomy.abundance.PMI))

taxonomy.description <- read.delim("original_datasets/TaxID_Description.txt", header = TRUE, sep = "\t", colClasses = c("character", "character", "factor"))

log_transform <- function(x) log10(x+1)
training.set.taxonomy.abundance <- log_transform(training.set.taxonomy.abundance)
test.set.taxonomy.abundance <- log_transform(test.set.taxonomy.abundance)

save(training.set.taxonomy.abundance, 
     test.set.taxonomy.abundance,
     training.set.class.labels, 
     taxonomy.description, 
     file="./datasets/taxonomy.abundance.datasets.RData")

pathway.abundance.He <- read.delim("original_datasets/TrainingHe_PathwayAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "PathID")
pathway.abundance.Schirmer <- read.delim("original_datasets/TrainingSchirmer_PathwayAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "PathID")
pathway.abundance.PMI <- read.delim("original_datasets/TestingDataset_PathwayAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "PathID")

training.set.pathway.abundance <- as.data.frame(t(cpm(cbind(pathway.abundance.He, pathway.abundance.Schirmer))))
test.set.pathway.abundance <- as.data.frame(t(cpm(pathway.abundance.PMI)))

pathway.description <- read.delim("original_datasets/PathID_Description.txt", header = TRUE, sep = "\t", colClasses = c("character", "character"))
pathway.description <- cbind(pathway.description, Rank = as.factor(ifelse(grepl("\\|", pathway.description$Pathway), "species", "community")))

save(training.set.pathway.abundance, 
     test.set.pathway.abundance,
     training.set.class.labels, 
     pathway.description, 
     file="./datasets/pathway.abundance.datasets.RData")
