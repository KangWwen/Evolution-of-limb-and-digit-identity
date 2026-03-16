library(TreeExp)
library(ape)
library(TreeTools)
library(ggplot2)

build_tree <- function(tpm_file, taxa, subtaxa, outgroup) {
  tetraExp <- TEconstruct(ExpValueFP = tpm_file, taxa = taxa, subtaxa = subtaxa)
  dismat   <- expdist(tetraExp, taxa = taxa, subtaxa = subtaxa, method = "pea")
  tr       <- NJ(dismat)
  expr_tab <- exptabTE(objects = tetraExp, taxa = "all", subtaxa = "all")
  tr       <- root(tr, outgroup, resolve.root = TRUE)
  bs       <- boot.exphy(phy = tr, x = expr_tab, method = "pea", outgroup = outgroup, B = 100)
  tr$node.label <- bs
  plot(tr, show.node.label = TRUE)
  edgelabels(round(tr$edge.length * 100, 2), col = "black", font = 1, frame = "none")
  return(list(tree = tr, tetraExp = tetraExp))
}

run_ztest <- function(tetraExp, subtaxa, x, y, outgroup) {
  exp_table <- exptabTE(tetraExp, taxa = "all", subtaxa = subtaxa)
  RelaRate.test(expTable = exp_table, x = x, y = y, outgroup = outgroup, alternative = "greater")
}

plot_ztest_bar <- function(labels, values) {
  dt <- data.frame(obj = factor(labels, levels = rev(labels)), val = values)
  ggplot(dt, aes(x = obj, y = val, fill = obj)) +
    coord_flip() +
    geom_bar(stat = "identity") +
    geom_text(aes(label = val)) +
    theme_minimal()
}

data    <- read.table("00.AllDigit_TPM.txt", header = TRUE, sep = "\t", row.names = 1)
limb    <- read.table("00_limb_genes.txt",   header = TRUE, sep = "\t", row.names = 1)
limb_TPM <- data[rownames(data) %in% rownames(limb), ]

taxa_eF    <- c("SiameseeF", "TurtleeF", "ChickeneF", "EmueF", "OstricheF")
taxa_lF    <- c("SiameselF", "TurtlelF", "ChickenlF", "EmulF", "OstrichlF")
subtaxa_D  <- c("D1", "D2", "D3", "D4", "D5")

res_eF <- build_tree("01.Limb_AllDigit_TPM.txt", taxa_eF, subtaxa_D, "TurtleeF_D1")
res_lF <- build_tree("01.Limb_AllDigit_TPM.txt", taxa_lF, subtaxa_D, "TurtlelF_D1")

for (x_taxa in c("ChickeneF", "EmueF", "OstricheF")) {
  print(run_ztest(res_eF$tetraExp, "D4", x_taxa, "SiameseeF", "TurtleeF"))
}

for (x_taxa in c("ChickenlF", "EmulF", "OstrichlF")) {
  print(run_ztest(res_lF$tetraExp, "D4", x_taxa, "SiameselF", "TurtlelF"))
}

plot_ztest_bar(ztest_labels, ztest_values)