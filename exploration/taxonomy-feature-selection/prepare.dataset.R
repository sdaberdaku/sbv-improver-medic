# prepare.dataset.R
# prepare the dataset and save into an RData file for easy loading
rm(list = ls())
class.labels.He <- read.delim("./original_data/Class_labels_He.txt", header = TRUE, sep = "\t", colClasses = c("integer", "character", "factor"))
class.labels.Schirmer <- read.delim("./original_data/Class_labels_Schirmer.txt", header = TRUE, sep = "\t", colClasses = c("integer", "character", "factor"))

taxonomy.abundance.He <- read.delim("./original_data/TrainingHe_TaxonomyAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "TaxID")
taxonomy.abundance.Schirmer <- read.delim("./original_data/TrainingSchirmer_TaxonomyAbundance_matrix.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = "TaxID")

taxonomy.abbundance.all <- as.data.frame(t(cbind(taxonomy.abundance.He, taxonomy.abundance.Schirmer)))
class.labels <- rbind(cbind(class.labels.He, dataset="He"), cbind(class.labels.Schirmer, dataset="Schirmer"))

taxonomy.description <- read.delim("./original_data/TaxID_Description.txt", header = TRUE, sep = "\t", colClasses = c("character", "character", "factor"))

log_transform <- function(x) log(x+1)
dataset.taxonomy.abundance <- log_transform(taxonomy.abbundance.all)

output_dir = "./Data"
if (!dir.exists(output_dir)){dir.create(output_dir)}
save(dataset.taxonomy.abundance, class.labels, taxonomy.description, file="./Data/training.set.RData")
