adaptive.RFE.linearSVM <- function(feat.select.samples, # the samples to use for the feature selection
                                   class.labels, # the class labels for each subject
                                   positive.label, 
                                   negative.label,
                                   bootstrap.iterations = 100, # number of boostrapping iterations 
                                   percentage.to.remove = 0.05, 
                                   feature.number.threshold = 100, # when to start single feature elimination
                                   parallel.bootstrap = FALSE,
                                   parallel.inner.cv = TRUE,
                                   seed = 42) {
  
  `%switchpar%` <- ifelse(parallel.bootstrap, `%dorng%`, `%do%`)
  
  positive.class.idx <- which(class.labels == positive.label)
  negative.class.idx <- which(class.labels == negative.label)
  
  set.seed(seed)
  feature.importance.results <- foreach (ii = 1:bootstrap.iterations) %switchpar% {
    cat("adaptive RFE iteration: ", ii, "\n", sep = "")
    ### bootstrap
    ### sample with replacement from the dataset the same number of subjects and repeat if
    ### there are no subjects in the bootstrap test set
    repeat.bootstrap <- TRUE
    while(repeat.bootstrap) {
      # Here we repeat the bootstrap sampling procedure for the positive and negative samples 
      # in order to keep the bootstrap training and test sets balanced
      bootstrap.positive.train.indices <- sample(positive.class.idx, size = length(positive.class.idx), replace = TRUE)
      bootstrap.positive.train.indices <- unique(bootstrap.positive.train.indices)
      bootstrap.positive.test.indices <- setdiff(positive.class.idx, bootstrap.positive.train.indices)
      
      bootstrap.negative.train.indices <- sample(negative.class.idx, size = length(negative.class.idx), replace = TRUE)
      bootstrap.negative.train.indices <- unique(bootstrap.negative.train.indices)
      bootstrap.negative.test.indices <- setdiff(negative.class.idx, bootstrap.negative.train.indices)
      
      bootstrap.train.indices <- c(bootstrap.positive.train.indices, bootstrap.negative.train.indices)
      bootstrap.test.indices <- c(bootstrap.positive.test.indices, bootstrap.negative.test.indices)
      
      if (length(bootstrap.test.indices)) {
        bootstrap.train <- feat.select.samples[bootstrap.train.indices, ]
        bootstrap.test <- feat.select.samples[bootstrap.test.indices, ]
        
        bootstrap.train.labels <- class.labels[bootstrap.train.indices]
        bootstrap.test.labels <- class.labels[bootstrap.test.indices]
        
        ### Remove features with zero standard deviation in bootstrap.train
        constant.features.train <- which(apply(bootstrap.train, 2, sd) == 0)
        if (length(constant.features.train) < dim(bootstrap.train)[2]) {
          if (length(constant.features.train)) {
            bootstrap.train <- bootstrap.train[ , -constant.features.train, drop=FALSE]
            bootstrap.test <- bootstrap.test[ , -constant.features.train, drop=FALSE]
          }
          repeat.bootstrap <- FALSE
        } 
      }
    }
    
    ### Add class labels to the training data
    training.dataset <- cbind(bootstrap.train, Label = bootstrap.train.labels)

    grid <- expand.grid(cost = 2^seq(-5, 15, 2))
    
    # Setting the seed for reproducible random numbers
    n.folds <- 5
    n.repeats <- 1
    seeds <- vector(mode = "list", length = n.folds * n.repeats + 1)
    for(i in 1:(n.folds * n.repeats)) {
      seeds[[i]] <- sample.int(10000, dim(grid)[1] * dim(grid)[2])
    }
    ## For the last model:
    seeds[[n.folds * n.repeats + 1]] <- sample.int(10000, 1)
    
    control <- trainControl(method = "cv", 
                            number = n.folds, 
                            summaryFunction = mccSummary, 
                            seeds = seeds,
                            allowParallel = parallel.inner.cv)

    ### 5%-RFE
    cat("5%-RFE\n")
    remaining.features <- colnames(bootstrap.train)
    ranked.features <- constant.features.train
    feature.performance <- data.frame(n.features=integer(), MCC=double())
    while(length(remaining.features) > feature.number.threshold) {
      # cat(length(remaining.features), "remaining features \n")
      ### train SVM 
      best.svm.model <- train(Label~., 
                              data=training.dataset[, c(remaining.features, "Label")], 
                              method="svmLinear2",
                              preProcess = c("center", "scale"),
                              trControl=control, 
                              tuneGrid=grid,
                              maximize=TRUE,
                              metric = "MCC")
      ### predict the test set
      svm.test.prediction <- predict(best.svm.model, newdata=bootstrap.test[, remaining.features, drop=FALSE])#, type = "prob")
      
      # compute delta.J
      w <- t(best.svm.model$finalModel$coefs) %*% best.svm.model$finalModel$SV
      delta.J <- w^2
      # sort them in descending order without loosing the column names
      delta.J <- delta.J[, order(delta.J[1,], decreasing=TRUE), drop=FALSE] 
      # ordered feature names
      if (length(remaining.features) > 1) {
        remaining.features <- colnames(delta.J)
      }
      
      MCC <- mcc(actuals=bootstrap.test.labels, preds=svm.test.prediction)
      # cat("MCC:", MCC, "\n")
      feature.performance[nrow(feature.performance) + 1, "n.features"] <- length(remaining.features)
      feature.performance[nrow(feature.performance), "MCC"] <- MCC 
      
      feature.idx.rm <- which(delta.J <= quantile(delta.J, probs = percentage.to.remove))
      ranked.features <- c(remaining.features[feature.idx.rm], ranked.features)
      remaining.features <- remaining.features[-feature.idx.rm]
    }
    cat("1-RFE\n")
    while(length(remaining.features)) {
      # cat(length(remaining.features), "remaining features \n")
      ### train SVM 
      best.svm.model <- train(Label~., 
                              data=training.dataset[, c(remaining.features, "Label")], 
                              method="svmLinear2",
                              preProcess = c("center", "scale"),
                              trControl=control, 
                              tuneGrid=grid,
                              maximize=TRUE,
                              metric = "MCC")
      ### predict the test set
      svm.test.prediction <- predict(best.svm.model, newdata=bootstrap.test[, remaining.features, drop=FALSE])
      
      # compute delta.J
      w <- t(best.svm.model$finalModel$coefs) %*% best.svm.model$finalModel$SV
      delta.J <- w^2
      delta.J <- delta.J[, order(delta.J[1,], decreasing=TRUE), drop=FALSE] # sort them in descending order without loosing the column names
      
      # ordered feature names
      if (length(remaining.features) > 1) {
        remaining.features <- colnames(delta.J)
      }
      
      # cat("predictions: ", paste(svm.test.prediction, collapse = " "), "\n", sep="")
      MCC <- mcc(actuals=bootstrap.test.labels, preds=svm.test.prediction)
      # cat("MCC:", MCC, "\n")
      feature.performance[nrow(feature.performance) + 1, "n.features"] <- length(remaining.features)
      feature.performance[nrow(feature.performance), "MCC"] <- MCC 
      
      feature.idx.rm <- which.max(delta.J == min(delta.J))
      ranked.features <- c(remaining.features[feature.idx.rm], ranked.features)
      remaining.features <- remaining.features[-feature.idx.rm]
    }
    return(list("ranked.features" = ranked.features, "feature.performance" = feature.performance))
  }
  
  # alphabetically sorted list of all features
  all.features <- sort(colnames(feat.select.samples))
  # initialize average feature ranks
  average.feature.ranks <- rep(0, length(all.features))
  names(average.feature.ranks) <- all.features
  
  # compute the average ranks of each features
  for (ii in 1:bootstrap.iterations) {
    average.feature.ranks <- average.feature.ranks + match(all.features, feature.importance.results[[ii]]$ranked.features)/bootstrap.iterations
  }
  sorted.average.feature.ranks <- sort(average.feature.ranks)
  sorted.features <- names(sorted.average.feature.ranks)
  
  
  all.MCC.values <- matrix(data = rep(NA, length(all.features)*bootstrap.iterations), 
                           nrow = bootstrap.iterations, ncol = length(all.features), 
                           dimnames = list(c(1:bootstrap.iterations), c(1:length(all.features))))
  for (ii in 1:bootstrap.iterations) {
    for (j in 1:dim(feature.importance.results[[ii]]$feature.performance)[1]) {
      n.features <- feature.importance.results[[ii]]$feature.performance[j, "n.features"]
      MCC <- feature.importance.results[[ii]]$feature.performance[j, "MCC"]
      all.MCC.values[ii, n.features] <- MCC
    }
    all.MCC.values[ii, ] <- na.locf(all.MCC.values[ii, ])
  }
  
  average.MCC <- colMeans(all.MCC.values)
  stdev.MCC <- apply(all.MCC.values, 2, sd)
  median.MCC <- apply(all.MCC.values, 2, median)
  
  # determine the optimal number of features
  optimal.number.of.features <- which.max(average.MCC == max(average.MCC))
  
  # Finally, extract the best features
  best.features <- sorted.features[1:optimal.number.of.features]
  
  return(list("average.MCC" = average.MCC,
              "median.MCC" = median.MCC,
              "stdev.MCC" = stdev.MCC,
              "optimal.number.of.features" = optimal.number.of.features, 
              "best.features" = best.features, 
              "sorted.average.feature.ranks" = sorted.average.feature.ranks))
}
