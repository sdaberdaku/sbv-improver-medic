# prediction.R
rm(list=ls())
.libPaths( c( .libPaths(), "/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library") )
#### load required packages and sources ####
# install missing packages
list.of.packages <- c("caret", 
                      "e1071", 
                      "kernlab", 
                      "pROC", 
                      "mltools", 
                      "plyr", 
                      "zoo", 
                      "ggplot2", 
                      "tictoc", 
                      "doMC", 
                      "glmnet", 
                      "broom", 
                      "reshape2",
                      "doRNG",
                      "rngtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) { 
  install.packages(new.packages, dependencies = TRUE) 
}
# load required packages
for (package in list.of.packages) {
  eval(parse(text=(paste("library(", package, ")", sep="")))) 
}
rm(list.of.packages, new.packages, package)

source("./functions/randomised.logistic.regression.R")

#### Setup parallel backend to use many processors ####
cores <- 31 # setting all cores might end up using too much RAM!
registerDoMC(cores)  # change the number of CPU cores

#### various parameters ####
selection.threshold = 0.3

n.folds = 10
n.repeats = 10

C = 2^(-5:15)
sigma = 2^(-15:3)
grid <- expand.grid(C = C, sigma = sigma)

# Setting the seed for reproducible random numbers
seed <- 42
set.seed(seed)
seeds <- vector(mode = "list", length = n.folds * n.repeats + 1)
for(i in 1:(n.folds * n.repeats)) {
  seeds[[i]] <- sample.int(10000, dim(grid)[1] * dim(grid)[2])
}
## For the last model:
seeds[[n.folds * n.repeats + 1]] <- sample.int(10000, 1)
# trainControl object setup for 10 fold CV
control <- trainControl(method = "repeatedcv",
                        number = n.folds,
                        repeats = n.repeats,
                        seeds = seeds,
                        classProbs = TRUE,
                        summaryFunction = twoClassSummary,
                        search = "grid",
                        allowParallel = TRUE,
                        preProcOptions = c("center", "scale"),
                        savePredictions="final")

vars.and.functions <- c("vars.and.functions", ls())
dir.create("results", showWarnings = FALSE)

#### Taxonomy Abundance ####

#### 1. IBD vs. non-IBD ####
load("datasets/taxonomy.abundance.datasets.RData")

#keep only species level TaxIDs
taxIDs.to.keep <- taxonomy.description[which(taxonomy.description$Rank == "species"), "TaxID"]
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]
test.set.taxonomy.abundance <- test.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]

levels(training.set.class.labels$group) <- c("IBD", "nonIBD", "IBD")

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.taxonomy.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "IBD", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- taxonomy.description[which(taxonomy.description$TaxID %in% best.features), "Taxon"]

features.df = data.frame(TaxonomyID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6), 
                         Description = feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Taxonomy_IBD_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.taxonomy.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "IBD"), best.svm.model$pred[, "IBD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.taxonomy.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.taxonomy.abundance),	
                                        Confidence_Value_IBD = round(test.set.prediction$IBD, digits = 6),
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$IBD, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Taxonomy_IBD_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))

#### 2. UC vs. non-IBD ####
load("datasets/taxonomy.abundance.datasets.RData")

#keep only species level TaxIDs
taxIDs.to.keep <- taxonomy.description[which(taxonomy.description$Rank == "species"), "TaxID"]
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]
test.set.taxonomy.abundance <- test.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "CD")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.taxonomy.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "UC", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- taxonomy.description[which(taxonomy.description$TaxID %in% best.features), "Taxon"]

features.df = data.frame(TaxonomyID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6), 
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Taxonomy_UC_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.taxonomy.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "UC"), best.svm.model$pred[, "UC"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.taxonomy.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.taxonomy.abundance),	
                                        Confidence_Value_UC = round(test.set.prediction$UC, digits = 6),
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$UC, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Taxonomy_UC_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))
#### 3. CD vs. non-IBD ####
load("datasets/taxonomy.abundance.datasets.RData")

#keep only species level TaxIDs
taxIDs.to.keep <- taxonomy.description[which(taxonomy.description$Rank == "species"), "TaxID"]
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]
test.set.taxonomy.abundance <- test.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "UC")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.taxonomy.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "CD", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- taxonomy.description[which(taxonomy.description$TaxID %in% best.features), "Taxon"]

features.df = data.frame(TaxonomyID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6), 
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Taxonomy_CD_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.taxonomy.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "CD"), best.svm.model$pred[, "CD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.taxonomy.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.taxonomy.abundance),	
                                        Confidence_Value_CD = round(test.set.prediction$CD, digits = 6),
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$CD, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Taxonomy_CD_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))
#### 4. UC vs. CD ####
load("datasets/taxonomy.abundance.datasets.RData")

#keep only species level TaxIDs
taxIDs.to.keep <- taxonomy.description[which(taxonomy.description$Rank == "species"), "TaxID"]
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]
test.set.taxonomy.abundance <- test.set.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "nonIBD")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.taxonomy.abundance <- training.set.taxonomy.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.taxonomy.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "UC", 
                                              negative.label = "CD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- taxonomy.description[which(taxonomy.description$TaxID %in% best.features), "Taxon"]

features.df = data.frame(TaxonomyID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6),
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Taxonomy_CD_vs_UC_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.taxonomy.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "CD"), best.svm.model$pred[, "CD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.taxonomy.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.taxonomy.abundance),	
                                        Confidence_Value_CD = round(test.set.prediction$CD, digits = 6),
                                        Confidence_Value_UC = 1 - round(test.set.prediction$CD, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Taxonomy_CD_vs_UC_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))


#### Pathway Abundance ####

#### 1. IBD vs. non-IBD ####
load("datasets/pathway.abundance.datasets.RData")

#keep only species level PathIDs
pathIDs.to.keep <- pathway.description[which(pathway.description$Rank == "species"), "PathID"]
training.set.pathway.abundance <- training.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]
test.set.pathway.abundance <- test.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]

levels(training.set.class.labels$group) <- c("IBD", "nonIBD", "IBD")

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.pathway.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "IBD", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- pathway.description[which(pathway.description$PathID %in% best.features), "Pathway"]

features.df = data.frame(PathwayID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6),
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Pathways_IBD_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.pathway.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "IBD"), best.svm.model$pred[, "IBD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.pathway.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.pathway.abundance),	
                                        Confidence_Value_IBD = round(test.set.prediction$IBD, digits = 6),
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$IBD, digits = 6))
write.table(format(test.set.prediction.output, digits=5, scientific=F), 
            file=file.path(output.dir, "SC2-Processed_Pathways_IBD_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))

#### 2. UC vs. non-IBD ####
load("datasets/pathway.abundance.datasets.RData")

#keep only species level PathIDs
pathIDs.to.keep <- pathway.description[which(pathway.description$Rank == "species"), "PathID"]
training.set.pathway.abundance <- training.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]
test.set.pathway.abundance <- test.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "CD")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.pathway.abundance <- training.set.pathway.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.pathway.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "UC", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- pathway.description[which(pathway.description$PathID %in% best.features), "Pathway"]

features.df = data.frame(PathwayID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6), 
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Pathways_UC_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.pathway.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "UC"), best.svm.model$pred[, "UC"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.pathway.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.pathway.abundance),	
                                        Confidence_Value_UC = round(test.set.prediction$UC, digits = 6), 
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$UC, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Pathways_UC_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))
#### 3. CD vs. non-IBD ####
load("datasets/pathway.abundance.datasets.RData")

#keep only species level PathIDs
pathIDs.to.keep <- pathway.description[which(pathway.description$Rank == "species"), "PathID"]
training.set.pathway.abundance <- training.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]
test.set.pathway.abundance <- test.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "UC")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.pathway.abundance <- training.set.pathway.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.pathway.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "CD", 
                                              negative.label = "nonIBD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- pathway.description[which(pathway.description$PathID %in% best.features), "Pathway"]

features.df = data.frame(PathwayID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6),
                         Description= feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Pathways_CD_vs_nonIBD_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.pathway.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "CD"), best.svm.model$pred[, "CD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.pathway.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.pathway.abundance),	
                                        Confidence_Value_CD = round(test.set.prediction$CD, digits = 6),
                                        Confidence_Value_nonIBD = 1 - round(test.set.prediction$CD, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Pathways_CD_vs_nonIBD_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))
#### 4. UC vs. CD ####
load("datasets/pathway.abundance.datasets.RData")

#keep only species level PathIDs
pathIDs.to.keep <- pathway.description[which(pathway.description$Rank == "species"), "PathID"]
training.set.pathway.abundance <- training.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]
test.set.pathway.abundance <- test.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "nonIBD")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.pathway.abundance <- training.set.pathway.abundance[-subjects.to.remove, ]

# Perform feature selection with the randomised.logistic.regression algorithm
tic()
RLR.results <- randomised.logistic.regression(feat.select.samples = training.set.pathway.abundance,
                                              sample.labels = training.set.class.labels$group,
                                              positive.label = "UC", 
                                              negative.label = "CD",
                                              selection.threshold = selection.threshold,
                                              regularisation.type = 0.5,
                                              lambda.values = 2^seq(-10, 2, 2),
                                              seed = seed,
                                              parallel = TRUE)
toc()

best.features <- RLR.results$selected.features
feature.importance <- RLR.results$stability.scores[best.features]
feature.description <- pathway.description[which(pathway.description$PathID %in% best.features), "Pathway"]

features.df = data.frame(PathwayID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6),
                         Description = feature.description)
output.dir <- file.path("results", paste("SC2", selection.threshold, "species", sep="_"))
dir.create(output.dir, showWarnings = FALSE)
write.table(features.df, file=file.path(output.dir, "SC2-Processed_Pathways_CD_vs_UC_Features.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

## Add the labels to the training set
training.dataset <- cbind(training.set.pathway.abundance[, best.features], Label = training.set.class.labels$group)

tic()
set.seed(seed)
best.svm.model <- train(Label ~ ., 
                        data=training.dataset, 
                        method = "svmRadial",
                        preProcess = c("center", "scale"),
                        trControl=control, 
                        tuneGrid=grid,
                        maximize=TRUE,
                        metric = "ROC")
toc()

roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "CD"), best.svm.model$pred[, "CD"])
cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")

test.set.prediction <- predict(best.svm.model, newdata=test.set.pathway.abundance[, best.features, drop=FALSE], type = "prob")

test.set.prediction.output = data.frame(SampleID = rownames(test.set.pathway.abundance),	
                                        Confidence_Value_CD = round(test.set.prediction$CD, digits = 6), 
                                        Confidence_Value_UC = 1 - round(test.set.prediction$CD, digits = 6))
write.table(test.set.prediction.output, file=file.path(output.dir, "SC2-Processed_Pathways_CD_vs_UC_Prediction.txt"), sep="\t", quote = FALSE, row.names = F, col.names = T)

rm(list = setdiff(ls(), vars.and.functions))

