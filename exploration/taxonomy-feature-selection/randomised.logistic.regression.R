#' @title The function performs Stability Selection with Randomised Logistic Regression
#' 
#' @description Stability selection is an example of an embedded feature selection method that tends to work well
#' in high-dimensional, sparse, problems. Embedded methods are a catch-all group of techniques that perform 
#' feature selection as part of the model construction process. They typically take the interaction between 
#' feature subset search and the learning algorithm into account, at the cost of extra computational time.
#'
#' Stability selection wraps around an existing structure learning algorithm, and tries to enhance and 
#' improve it. The goal of these kind of algorithms is to learn the underlying relationships of the data. 
#' Structure learning can take different forms, but one example is the LASSO algorithm for supervised regression. 
#' Here, the algorithm tries to determine which features relate to the target variable, and the (covariance) 
#' structure is then defined as the set of features with non-zero values for their coefficients in the final 
#' linear model.
#'
#' What these structure learning algorithms typically have in common is a parameter *lambda* that controls the
#' amount of regularisation, mentioned above within the parentheses. For every value of *lambda*, we can obtain 
#' a structure estimate that indicates which features to select. 
#'
#' The rough idea behind stability selection is to inject more noise into the original problem by generating 
#' bootstrap samples of the data, and to use a base structure learning algorithm to find out which features 
#' are important in every sampled version of the data.
#'
#' For a feature to be considered stable (or important), it has to be selected in a high number of perturbed 
#' versions of the original problem. This tends to filter out features that are only weakly related to the 
#' target variables, because the additional noise introduced by the bootstrapping breaks that weak relationship.
#'
#' The algorithm takes as input a grid of regularisation parameters *lambda.values*, and the number of 
#' sampling iterations *n.bootstrap.iterations* that need to be generated. Stability selection returns a 
#' selection probability for every regularisation parameter value and for every feature, and the set of 
#' stable features.
#'
#' The base algorithm is run on every bootstrap sample for a grid of values of the penalisation parameter, 
#' and not just a single value. A-priori we don't know what the right level of regularisation for a problem 
#' is. If a variable is related to the target variable in a meaningful way, it should show up in most bootstrap 
#' samples for at least one value of the penalisation parameter.
#'
#' When the stability score for a variable exceeds the threshold *selection.threshold* for one value in 
#' *lambda.values*, it is deemed stable. This means that some regularisation is necessary. In specific cases 
#' one can obtain minimum bounds (as obtained in the Meinshausen and Buhlmann paper), but in general this 
#' requires some fiddling around with the mathematical details.
#'
#' In the original paper by Meinshausen and Buhlmann, the authors propose to inject even more noise into 
#' the problem, by using a learner called the randomised LASSO that essentially employs random regularisation 
#' scores for each feature. The randomised lasso is a new generalisation of the lasso. Whereas the lasso 
#' penalises the absolute value of every component with a penalty term proportional to *lambda*, the randomised 
#' lasso changes the penalty *lambda* to a randomly chosen value in {*lambda*, *lambda*/*alpha*} with 
#' probability *p.w*.
#'
#' References:
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society: 
#' Series B (Statistical Methodology), 72(4), 417-473. DOI:10.1111/j.1467-9868.2010.00740.x
#'
#' Shah, R. D., & Samworth, R. J. (2013). Variable selection with error control: another look at stability 
#' selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(1), 55-80.
#' DOI:10.1111/j.1467-9868.2011.01034.x
#'
#' @author Sebastian Daberdaku
#' 
#' @param feat.select.samples matrix or data frame containing the samples to be used for feature selection.
#' @param metadata the class labels for each subject in *feat.select.samples*.
#' @param positive.label the label of the positive class.
#' @param negative.label the label of the negative class.
#' @param n.bootstrap.iterations the number of bootstrapping iterations. Defaults to 100.
#' @param lambda.values list of regularisation parameters. By default, these values are selected according to the method described in Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. Journal of statistical software, 33(1), 1.
#' @param alpha scale value for the regularisation parameters used in the internal randomisation step. Defaults to 0.2.
#' @param p.w the probability for a feature to be penalised with a *lambda*/*alpha* regularisation factor. Defaults to 0.5.
#' @param sample.fraction fraction of the overall samples to be resampled at each bootstrapping iteration. Defaults to 0.5. 
#' @param regularisation.type "1": for LASSO regression; "0": for ridge regression; a value between 0 and 1 (i.e. 0.3) for elastic net regression. Defaults to 1 (LASSO).
#' @param epsilon threshold used to determine if the Logistic Regression model coeffients are non-null. Defaults to 1e-3.
#' @param selection.threshold The threshold value used to determine the selected features. This is a tuning parameter whose influence is very small. For sensible values in the range of [0.6, 0.9] results tend to be very similar. Defaults to 0.6.
#' @param parallel if TRUE, the iterations are performed in parallel. An appropriate parallel backed must be registered before hand, such as *doMC* or *doSNOW*. Defaults to TRUE.
#' 
#' @return a list containing the following objects:
#' \describe{
#'   \item{feature.selection.probabilities}{The probability that each feature is selected by the Randomised Logistic Regression algorithm for each regularisation parameter.}
#'   \item{stability.scores}{The stability scores for each feature.}
#'   \item{selected.features}{The list of selected features.}
#' }
#' @keywords stability selection, randomised logistic regression
#' @export
randomised.logistic.regression <- function(feat.select.samples, 
                                           metadata, 
                                           positive.label, 
                                           negative.label,
                                           n.bootstrap.iterations = 100,
                                           lambda.values = NULL,
                                           alpha = 0.2,
                                           p.w = 0.5,
                                           sample.fraction = 0.5,
                                           regularisation.type = 1, 
                                           epsilon = 1e-9,
                                           selection.threshold = 0.6,
                                           parallel = TRUE) {
  
  `%switchpar%` <- ifelse(parallel, `%dopar%`, `%do%`)
  
  positive.class.idx <- which(metadata$group == positive.label)
  negative.class.idx <- which(metadata$group == negative.label)
  
  all.features <- colnames(feat.select.samples)
  
  
  # Determine the regularisation parameter list. The default values depend on the sample size nobs
  # relative to the number of variables nvars. If nobs > nvars, the default is 0.0001, close to zero. 
  # If nobs < nvars, the default is 0.01. A very small value of lambda.min.ratio will lead to a 
  # saturated fit in the nobs < nvars case.
  if (is.null(lambda.values)) {
    # Remove features with zero standard deviation in bootstrap.train
    const.feats <- which(lapply(feat.select.samples, sd) == 0)
    if (length(const.feats)) {
      x <- scale(as.matrix(feat.select.samples[, -const.feats]))
    } else {
      x <- scale(as.matrix(feat.select.samples))
    }
    y <- ifelse(metadata$group == positive.label, 1, 0)
    e <- ifelse(dim(x)[1] > dim(x)[2], 0.0001, 0.01)
    a <- ifelse(regularisation.type > 0, regularisation.type, 1)
    lambda.max <- max(abs(t(y - mean(y) * (1 - mean(y))) %*% x)) / (a * dim(x)[1])
    lambda.values <- round(exp(seq(log(lambda.max), log(lambda.max * e), length.out = 100)), digits = 10)
  }

  feature.importance.results <- foreach (lambda = lambda.values, .final = function(x) setNames(x, lambda.values)) %:% 
    foreach (ii = 1:n.bootstrap.iterations, .combine = rbind) %switchpar% {
      ### bootstrap
      ### sample with replacement from the dataset the same number of subjects and repeat if
      ### there are no subjects in the bootstrap test set
      repeat.bootstrap <- TRUE
      while(repeat.bootstrap) {
        # Here we repeat the bootstrap sampling procedure for the positive and negative samples 
        # in order to keep the bootstrap training and test sets balanced
        bootstrap.positive.train.indices <- sample(positive.class.idx, size = floor(sample.fraction * length(positive.class.idx)), replace = FALSE)
        bootstrap.negative.train.indices <- sample(negative.class.idx, size = floor(sample.fraction * length(negative.class.idx)), replace = FALSE)
        
        bootstrap.train.indices <- c(bootstrap.positive.train.indices, bootstrap.negative.train.indices)
        bootstrap.train <- feat.select.samples[bootstrap.train.indices, ]
        bootstrap.train.labels <- metadata[bootstrap.train.indices, ]
        
        ### Remove features with zero standard deviation in bootstrap.train
        constant.features.train <- which(apply(bootstrap.train, 2, sd) == 0)
        if (length(constant.features.train) < dim(bootstrap.train)[2]) {
          if (length(constant.features.train)) {
            bootstrap.train <- bootstrap.train[, -constant.features.train, drop=FALSE]
          }
          repeat.bootstrap <- FALSE
        }
      }
      
      feature.list <- colnames(bootstrap.train)
      
      ### Add class labels to the training data
      bootstrap.train.classes <- cbind(bootstrap.train, bootstrap.train.labels$group)
      colnames(bootstrap.train.classes)[dim(bootstrap.train.classes)[2]] <- "Label"
      
      x.train <- model.matrix(Label ~ ., bootstrap.train.classes)[,-1]
      y.train <- ifelse(bootstrap.train.classes$Label == positive.label, 1, 0)
      
      w <- sample(c(1/alpha, 1), size = length(feature.list), replace = TRUE, prob = c(p.w, 1-p.w))
      
      lr.model <- glmnet(x = x.train, 
                         y = y.train, 
                         lambda = lambda,
                         penalty.factor = w,
                         alpha = regularisation.type, # "1": for lasso regression; "0": for ridge regression; a value between 0 and 1 (say 0.3) for elastic net regression.
                         family = "binomial", # Use "binomial" for a binary outcome variable
                         standardize = TRUE)
      
      model.coefs <- as.matrix(coef(lr.model, s = "lambda.min"))[-1, ]

      feature.weights <- as.numeric(abs(model.coefs) > epsilon)
      names(feature.weights) <- names(model.coefs)
      
      feature.importance <- matrix(0, nrow = 1, ncol = length(all.features), dimnames = list(ii, all.features))
      
      for (f in feature.list) {
        feature.importance[1, f] <- feature.importance[1, f] + feature.weights[[f]]
      }
      
      return(feature.importance)
    }
  
  feat.select.probs.list <- lapply(feature.importance.results, colMeans)
  feature.selection.probabilities <- do.call(rbind, feat.select.probs.list)
  
  stability.scores <- colMeans(feature.selection.probabilities)
  
  selected.features.idx <- which(apply(feature.selection.probabilities, 2, function(x) max(x) > selection.threshold))
  if (length(selected.features.idx) == 0) {
    cat("No features selected! Try lowering the selection.threshold parameter!\n")
  }
  selected.features <- colnames(feature.selection.probabilities)[selected.features.idx]
  
  feature.selection.probabilities <- cbind(lambda=lambda.values, feature.selection.probabilities)
  
  return(list("feature.selection.probabilities" = as.data.frame(feature.selection.probabilities), 
              "stability.scores" = stability.scores, 
              "selected.features" = selected.features))
}

