# Hierarchical clustering with bootstrap (pvclust) for limb gene expression
# across Chicken, Emu, Ostrich, Turtle, and Crocodile (Si)

library(dendextend)
library(pvclust)
library(ggplot2)
library(reshape2)

setwd("./")

# Load expression data
chicken  <- read.table("06.Chicken_limbgene_TPM_AV_sqrt.txt",  header = TRUE, sep = "\t", row.names = 1)
si       <- read.table("06.Corpor_limbgene_TPM_AV_sqrt.txt",   header = TRUE, sep = "\t", row.names = 1)
turtle   <- read.table("06.Turtle_limbgene_TPM_AV_sqrt.txt",   header = TRUE, sep = "\t", row.names = 1)
emu      <- read.table("06.Emu_limbgene_TPM_AV_sqrt.txt",      header = TRUE, sep = "\t", row.names = 1)
ostrich  <- read.table("06.Ostirch_limbgene_TPM_AV_sqrt.txt",  header = TRUE, sep = "\t", row.names = 1)

# Run pvclust on combined expression matrix and plot result
run_pvclust <- function(mat, nboot = 1000) {
  zscore <- as.data.frame(scale(mat, center = FALSE, scale = TRUE))
  hc <- pvclust(zscore, method.hclust = "average", method.dist = "correlation", nboot = nboot)
  plot(hc)
  invisible(hc)
}

# Chicken: forelimb (col 1-3) vs hindlimb (col 8-10)
run_pvclust(cbind(chicken[, 1:3],  turtle[, 1:5],  si[, 1:5]))
run_pvclust(cbind(chicken[, 8:10], turtle[, 11:15], si[, 11:15]))

# Emu: forelimb (col 1-3) vs hindlimb (col 7-9)
run_pvclust(cbind(emu[, 1:3],  turtle[, 1:5],  si[, 1:5]))
run_pvclust(cbind(emu[, 7:9],  turtle[, 11:15], si[, 11:15]))

# Ostrich: forelimb (col 1-3) vs hindlimb (col 9-11)
run_pvclust(cbind(ostrich[, 1:3],  turtle[, 1:5],  si[, 1:5]))
run_pvclust(cbind(ostrich[, 9:11], turtle[, 11:15], si[, 11:15]))