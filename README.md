# sbv IMPROVER MEDIC: team GiGi's code

Code for sub-challenge 2 of the sbv IMPROVER Metagenomics Diagnosis for Inflammatory Bowel Disease Challenge (MEDIC, 2019–2020): classifying IBD versus non-IBD, Crohn's disease versus ulcerative colitis, and the other class pairs from taxonomy and pathway abundance profiles of stool metagenomes.

I wrote this code in January 2020 as a member of team GiGi (University of Padua: Barbara Di Camillo, Sebastian Daberdaku, Ilaria Patuzzi, Mehdi Poursheikhali Asghari, Filippo Pietrobon). The team was one of the three winners of the challenge, the best-performing teams of sub-challenge 2. The results are described in:

> Khachatryan L. et al. **Results and lessons learned from the sbv IMPROVER metagenomics diagnostics for inflammatory bowel disease challenge.** *Scientific Reports* 13, 6303 (2023). https://doi.org/10.1038/s41598-023-33050-0

## Method

For each classification task, features are selected on the training sets with one of two methods, and an SVM with an RBF kernel is then trained on the selected features. Its cost and kernel width are tuned by 10-fold cross-validation repeated 10 times (`caret`, ROC as the metric). `Writeup.odt` is the method description submitted with the predictions.

- `functions/randomised.logistic.regression.R`: stability selection with randomised logistic regression (Meinshausen and Bühlmann, 2010)
- `functions/adaptive.RFE.linearSVM.R`: recursive feature elimination with a linear SVM
- `functions/mccSummary.R`: Matthews correlation coefficient as a `caret` summary function
- `utils/prepare.datasets.R`: converts the challenge text files into `datasets/*.RData`
- `utils/round_up.results.R`: formats the predictions for submission
- `complete_jobs/`: one R script and SGE job file per run:
  - `prediction_<threshold>_all` and `prediction_<threshold>_species`: stability selection with a selection threshold from 0.1 to 0.4, on taxa of all ranks or species only
  - `prediction_RFE_taxo` and `prediction_RFE_path_1` to `_4`: recursive feature elimination on taxa or pathways

## Data

The challenge data is not included. The challenge rules licensed it to participants only for preparing a submission and did not allow sharing it.

The scripts expect these files in `original_datasets/`:

- `TrainingHe_TaxonomyAbundance_matrix.txt`, `TrainingHe_PathwayAbundance_matrix.txt`, `Class_labels_He.txt`
- `TrainingSchirmer_TaxonomyAbundance_matrix.txt`, `TrainingSchirmer_PathwayAbundance_matrix.txt`, `Class_labels_Schirmer.txt`
- `TestingDataset_TaxonomyAbundance_matrix.txt`, `TestingDataset_PathwayAbundance_matrix.txt`
- `TaxID_Description.txt`, `PathID_Description.txt`

The organizers built these matrices (Kraken2/Bracken taxonomy, HUMAnN2 pathways) from three shotgun metagenomics datasets whose raw reads are public:

- He et al., GigaScience 2017: ENA PRJEB15371, GigaDB https://doi.org/10.5524/100317
- IBDMDB / Schirmer et al., Nature Microbiology 2018 (HMP2): https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/rawfiles
- Test set (105 samples): SRA PRJNA893901

The processing steps are in the Methods of the paper above. Reprocessing the raw reads will not give exactly the same matrices.

## Exploration

`exploration/taxonomy-feature-selection/` holds an approach I tried before this pipeline, building features from groups of taxa that are close in the NCBI taxonomy. It was not used for the submissions; its README explains it.

## Software

The jobs ran on an SGE cluster inside a Singularity image with R 3.6.1 (not included). The R package versions are listed in `Writeup.odt` (output of `sessionInfo()`): mainly caret, glmnet, kernlab, e1071, pROC, edgeR, limma, foreach and doMC.

## License

MIT, see [LICENSE](LICENSE).
