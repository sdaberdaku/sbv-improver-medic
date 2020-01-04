# taxonomy_based.feature.selection.R
rm(list=ls())

list.of.packages <- c("ape", "caret", "e1071", "kernlab", "pROC", "mltools", "plyr", "zoo", "ggplot2", "tictoc", "doMC", "glmnet", "broom", "infotheo", "reshape2", "plot3D")
# install missing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) { 
  install.packages(new.packages, dependencies = TRUE) 
}
# load required packages
for (package in list.of.packages) {
  eval(parse(text=(paste("library(", package, ")", sep="")))) 
}
rm(list.of.packages, new.packages, package)

# Setup parallel backend to use many processors
cores <- detectCores() # setting all cores might end up using too much RAM!
registerDoMC(cores)  # change the number of CPU cores

load("./Data/training.set.RData")

taxID.list <- taxonomy.description$TaxID

# load the phylogenetic tree in newick format created using the ETE toolkit for python
# the TaxIDs have been mapped on the complete NCBI plylogenetic tree
phylogenetic.tree <- read.tree(file = "taxonomy_tree.nw")

# compute the Cophenetic distances among all tree node pairs
feature.distances <- dist.nodes(phylogenetic.tree)
feature.distances.labels <- paste("TaxID_", c(phylogenetic.tree$tip.label, phylogenetic.tree$node.label), sep="")
dimnames(feature.distances) <- list(feature.distances.labels, feature.distances.labels)

#keep only genus and species level TaxIDs
taxIDs.to.keep <- taxonomy.description[which(taxonomy.description$Rank %in% c("genus", "species")), "TaxID"]
dataset.taxonomy.abundance <- dataset.taxonomy.abundance[, taxIDs.to.keep, drop=FALSE]

# remove constant features
constant.features <- which(apply(dataset.taxonomy.abundance, 2, sd) == 0)
if(length(constant.features) > 0) {
  dataset.taxonomy.abundance <- dataset.taxonomy.abundance[, -constant.features, drop=FALSE]
}
# update the list of features
taxID.list <- colnames(dataset.taxonomy.abundance)

# keep only the TaxIDs we are interested in
indices.to.keep <- which(feature.distances.labels %in% taxID.list)
feature.distances <- feature.distances[indices.to.keep, indices.to.keep]

cat("indices.to.keep has the same number of elements as taxID.list: ", length(indices.to.keep) == length(taxID.list), "\n", sep="")
unmatched <- taxID.list[which(!(taxID.list %in% feature.distances.labels))]
cat("there are no unmatched TaxIDs: ", identical(unmatched, character(0)), "\n", sep="")

# perform Hierarchical Clustering on the Cophenecit distances among features
hc <- hclust(as.dist(feature.distances))
n.clusters = 200

clusters <- cutree(hc, k = n.clusters)

levels(class.labels$group) <- c("IBD", "nonIBD", "IBD")

n.samples <- dim(class.labels)[1]

seed <- 123
set.seed(seed)
n.folds <- n.samples
folds <- createFolds(class.labels$group, k = n.folds, list = TRUE, returnTrain = FALSE)

# trainControl object setup for 10 fold CV
control <- trainControl(method = "cv",
                        number = 10,
                        classProbs = TRUE,
                        summaryFunction = twoClassSummary,
                        search = "grid",
                        allowParallel = TRUE,
                        preProcOptions = c("center", "scale"),
                        savePredictions="final")

C = 2^(-5:15)
sigma = 2^(-15:3)
grid <- expand.grid(C = C, sigma = sigma)

selected.features <- vector("list", n.folds)

sample.prediction <- foreach(k=1:n.folds, .combine = rbind) %do% {
  cat("Computing fold: ", k, "\n", sep="")
  test.set.idx <- folds[[k]]
  training.set.idx <- setdiff(1:n.samples, test.set.idx)
  
  training.set <- dataset.taxonomy.abundance[training.set.idx, ,drop=FALSE]
  training.set.class.labels <- class.labels[training.set.idx, ,drop=FALSE]
  
  test.set <- dataset.taxonomy.abundance[test.set.idx, ,drop=FALSE]
  test.set.class.label <- class.labels[test.set.idx, ,drop=FALSE]
  
  # Perform feature selection 
  tic()
  selected.features[[k]] <- data.frame(TaxID=character(n.clusters), importance=numeric(n.clusters), stringsAsFactors = FALSE)
  for (ci in 1:n.clusters) {
    current.cluster.features <- names(which(clusters==ci))
    Y <- discretize(training.set.class.labels$group, disc= "equalwidth")
    cluster.entropies <- vector("list", length(current.cluster.features))
    names(cluster.entropies) <- current.cluster.features
    for (f in current.cluster.features) {
      X <- discretize(training.set[, f], disc= "equalwidth")
      cluster.entropies[[f]] <- mutinformation(X, Y)
    }
    selected.features[[k]][ci, "TaxID"] <- names(which.max(unlist(cluster.entropies)))
    selected.features[[k]][ci, "importance"] <- cluster.entropies[[which.max(unlist(cluster.entropies))]]
  }
  toc()

  
  best.features <- selected.features[[k]]$TaxID
  
  ## Add the labels to the training set
  training.dataset <- cbind(training.set, training.set.class.labels$group)
  colnames(training.dataset)[dim(training.dataset)[2]] <- "Label"
  
  tic()
  best.svm.model <- train(Label~., 
                          data=training.dataset[, c(best.features, "Label")], 
                          method = "svmRadial",
                          preProcess = c("center", "scale"),
                          trControl=control, 
                          tuneGrid=grid,
                          maximize=TRUE,
                          metric = "ROC")
  toc()
  
  roc.training.set <- roc(as.integer(best.svm.model$pred[, "obs"] == "IBD"), best.svm.model$pred[, "IBD"])
  cat("ROC on the training set: ", roc.training.set$auc, ".\n", sep="")
  
  test.set.prediction <- predict(best.svm.model, newdata=test.set[, best.features, drop=FALSE], type = "prob")
  # roc.test.set <- roc(as.integer(test.set.class.label$group == "IBD"), test.set.prediction[, "IBD"])
  # cat("ROC on the test set: ", roc.test.set$auc, ".\n", sep="")
  
  rownames(test.set.prediction) <- test.set.class.label$sampleID
  return(test.set.prediction)
}

sample.prediction <- sample.prediction[order(match(rownames(sample.prediction), class.labels$sampleID)), ]

roc.sample.prediction <- roc(as.integer(class.labels$group == "IBD"), sample.prediction[, "IBD"])
print(roc.sample.prediction)
plot(roc.sample.prediction)

save(sample.prediction, roc.sample.prediction, selected.features, file="loocv.fs.results.RData")


