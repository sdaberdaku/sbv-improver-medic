# round_up.results.R

round_extremes <- function(x) ifelse(x < 0.1, 0, ifelse(x > 0.9, 1, x))

for (f in c("all")) {#, "species")) {
  for (t in c(0.3)){#seq(0.1, 0.35, 0.05)) {
    src_dir <- "results/SC2_adaptiveRFE" #paste("results/SC2_adaboost", t, f, sep="_")
    out_dir <- "results/SC2_adaptiveRFE_rounded" #paste("results/SC2_adaboost", t, f, "rounded2", sep="_")
    dir.create(out_dir, showWarnings = FALSE)
    
    file.copy(paste(src_dir, "Writeup.odt", sep="/"), out_dir, overwrite = TRUE)
    
    for (task in c("CD_vs_nonIBD", "CD_vs_UC", "IBD_vs_nonIBD", "UC_vs_nonIBD")) {
      for (p in c("Pathways", "Taxonomy")) {
        features.file <- paste("SC2-Processed", p, task, "Features.txt", sep="_")
        prediction.file <- paste("SC2-Processed", p, task, "Prediction.txt", sep="_")
        file.copy(file.path(src_dir, features.file), out_dir)
        
        pred <- read.csv(file.path(src_dir, prediction.file), header = TRUE, sep = "\t")
        pred.round <- pred
        pred.round[, 2] <- round_extremes(pred.round[, 2])
        pred.round[, 3] <- 1 - pred.round[, 2]
        
        write.table(pred.round, file=file.path(out_dir, prediction.file), sep="\t", quote = FALSE, row.names = F, col.names = T, append = FALSE)
      }
    }
  }
}
