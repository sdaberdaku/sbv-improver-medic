# Taxonomy-based features (exploration, January 2020)

An approach I tried in the first days of the challenge, before the pipeline in the repository root. It was not used for the submissions.

The idea: taxa that are close in the NCBI taxonomy carry related information, so group them before learning, and give the classifier one feature per group instead of thousands of single taxa.

1. `prepare.dataset.R` reads the He and Schirmer training sets (taxonomy abundances and class labels) from `original_data/` and saves them to `Data/training.set.RData`.
2. The genus and species taxa are placed on the NCBI taxonomy tree, read from `taxonomy_tree.nw`, and clustered hierarchically on their distances in the tree. The tree was built with the ETE toolkit for Python and is not included.
3. Each script reduces the clusters to features in a different way:
   - `taxonomy_based.feature.extraction.R`: 2,000 clusters, each replaced by its first principal component; stability selection then keeps the informative ones.
   - `taxonomy_based.feature.extraction.and.selection.R`: the same, with the tree cut at a fixed height instead of a fixed number of clusters.
   - `taxonomy_based.feature.selection.R`: 200 clusters, keeping from each the taxon with the highest entropy-based score.
4. An SVM with an RBF kernel separates IBD from non-IBD samples, evaluated with leave-one-out cross-validation (10-fold for the extraction-and-selection variant).

`randomised.logistic.regression.R` is an earlier version of the stability selection in `functions/`.

The scripts run from this folder, with the challenge files in `original_data/` (the file names are in the main README). The challenge data is not included, for the reasons explained there; `original_data/`, `Data/`, the tree and the results are ignored by git.
