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

source("./functions/mccSummary.R")
source("./functions/adaptive.RFE.linearSVM.R")

#### Setup parallel backend to use many processors ####
cores <- 16 # setting all cores might end up using too much RAM!
registerDoMC(cores)  # change the number of CPU cores

#### various parameters ####
n.folds = 10
n.repeats = 10

C = 2^(-5:15)
sigma = 2^(-15:3)
grid <- expand.grid(C = C, sigma = sigma)

# Setting the seed for reproducible random numbers
seed <- 123
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

#### Pathway Abundance ####

#### 4. UC vs. CD ####
load("datasets/pathway.abundance.datasets.RData")

# #keep only species level PathIDs
# pathIDs.to.keep <- pathway.description[which(pathway.description$Rank == "species"), "PathID"]
# training.set.pathway.abundance <- training.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]
# test.set.pathway.abundance <- test.set.pathway.abundance[, pathIDs.to.keep, drop=FALSE]

subjects.to.remove <- which(training.set.class.labels$group == "nonIBD")
training.set.class.labels <- training.set.class.labels[-subjects.to.remove, ]
training.set.class.labels$group <- droplevels(training.set.class.labels$group)
training.set.pathway.abundance <- training.set.pathway.abundance[-subjects.to.remove, ]

tic()
RFE.results <- adaptive.RFE.linearSVM(feat.select.samples = training.set.pathway.abundance,
                                      class.labels = training.set.class.labels$group,
                                      positive.label = "UC", 
                                      negative.label = "CD",
                                      bootstrap.iterations = 100, # number of boostrapping iterations 
                                      percentage.to.remove = 0.05, 
                                      feature.number.threshold = 100,
                                      parallel.bootstrap = FALSE,
                                      parallel.inner.cv = TRUE,
                                      seed = seed)
toc()

best.features <- RFE.results$best.features
n.best.features <- RFE.results$optimal.number.of.features
feature.importance <- 1/RFE.results$sorted.average.feature.ranks[1:n.best.features]
feature.description <- pathway.description[which(pathway.description$PathID %in% best.features), "Pathway"]

features.df = data.frame(PathwayID = best.features, 
                         Importance_Optional = round(feature.importance, digits = 6), 
                         Description = feature.description)
output.dir <- file.path("results", paste("SC2", "adaptiveRFE", sep="_"))
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
