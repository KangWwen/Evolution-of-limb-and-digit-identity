library(dplyr)
library(pheatmap)
library(RColorBrewer)

gene_list <- c("HAND2","SPON2","HOXD10","HOXD11",
               "HOXD12","CHRDL1","GSC","HOXD9","PTCH2","TBX15",
               "PAX9","MSX1","ZIC3","ALX1","ALX4","DLX6","FGF18","LHX9","MDGA2","SFRP2",
               "DLX5","FAM20A","GFRA4","FGF8","GREM1",
               "SMOC1","ACAN","BMPER","TFAP2B",
               "HOXC11",
               "HOXC10","LGR5","BMP3","CLSTN2",
               "PAPPA2","SHH","COL14A1","SIM1",
               "ENPP2","ISL1","SHOX",
               "SHOX2","TBX2","TBX3","HOXC9","TWIST2","WNT3","WNT7A")

sample_names <- c("eChiF","lChiF","eEmuF","lEmuF",
                  "eEmuH","lEmuH","eOstF","lOstF",
                  "eAxmF", "lAxmF", "eChiH","lChiH",   
                  "eAxmH","lAxmH","eCorF","eCorH",
                  "lCorF","lCorH","eAcF","eAmF","lAmF","eMmF",
                  "eTurF","lTurF","eTurH",
                  "lTurH","eOstH","lOstH")

processTauData <- function(df, filteredGenes, numColumns) {
  template_df <- data.frame(Gene = filteredGenes,
                            matrix(0, nrow = length(filteredGenes), ncol = numColumns))
  rownames(template_df) <- filteredGenes
  df_subset <- df[rownames(df) %in% filteredGenes, 1:numColumns]
  df_subset$Gene <- rownames(df_subset)
  result_df <- merge(template_df, df_subset, by = "Gene", all = TRUE)
  result_df <- result_df[, c(1, (numColumns + 2):(numColumns * 2 + 1))]
  result_df[is.na(result_df)] <- 0
  rownames(result_df) <- result_df$Gene
  result_df$Gene <- NULL
  return(result_df)
}

minmax_norm <- function(mat) {
  t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))
}

select_max_stage <- function(early_NormTau, late_NormTau) {
  tau_col <- ncol(early_NormTau) - 1
  result <- early_NormTau
  for (i in 1:nrow(early_NormTau)) {
    if (as.numeric(late_NormTau[i, tau_col]) > as.numeric(early_NormTau[i, tau_col])) {
      result[i, ] <- late_NormTau[i, ]
    }
  }
  return(result)
}

setwd(tau_path)
for (name in sample_names) {
  assign(name, read.table(paste0(name, "_AlltauValue.txt"), sep = "\t", row.names = 1, header = TRUE))
}

tau_list <- list()
for (name in sample_names) {
  df <- get(name)
  tau_list[[name]] <- processTauData(df, gene_list, length(df))
  tau_list[[name]] <- tau_list[[name]][match(gene_order, rownames(tau_list[[name]])), ]
}

setwd(tpm_path)
limb_names <- c("eChiF","lChiF","eEmuF","lEmuF",
                "eEmuH","lEmuH","eOstF","lOstF",
                "eAxmF", "lAxmF", "eChiH","lChiH",   
                "eAxmH","lAxmH","eCorF","eCorH",
                "lCorF","lCorH","eAcF","eAmF","lAmF","eMmF",
                "eTurF","lTurF","eTurH",
                "lTurH","eOstH","lOstH")
for (name in limb_names) {
  assign(name, read.table(paste0("03.", name, "_AV_TPM.txt"), sep = "\t", row.names = 1, header = TRUE))
}

limb_list <- list()
for (name in limb_names) {
  df <- get(name)
  limb_list[[name]] <- processTauData(df, gene_list, length(df))
  limb_list[[name]] <- limb_list[[name]][match(gene_order, rownames(limb_list[[name]])), ]
}

sample_col_idx <- list(
  sampleA_F = list(F = 1:4,  H = 5:9),
  sampleB_F = list(F = 1:5,  H = 6:10)
)

norm_tau_list <- list()
for (sp in names(sample_col_idx)) {
  for (phase in c("e", "l")) {
    limb_mat <- limb_list[[sp]]
    tau_mat  <- tau_list[[paste0(phase, sp)]]
    for (limb in c("F", "H")) {
      idx  <- sample_col_idx[[sp]][[limb]]
      key  <- paste0(phase, sp, limb)
      norm <- minmax_norm(limb_mat[, idx])
      norm_tau_list[[key]] <- cbind(norm * tau_mat$Tau, tau_mat$Tau, tau_mat$Tau_info)
    }
  }
}

all_list <- list()
for (sp in names(sample_col_idx)) {
  for (limb in c("F", "H")) {
    key_e <- paste0("e", sp, limb)
    key_l <- paste0("l", sp, limb)
    all_list[[paste0(sp, limb)]] <- select_max_stage(norm_tau_list[[key_e]], norm_tau_list[[key_l]])
  }
}

AllResults_forHeatmap <- do.call(cbind, lapply(all_list, function(x) {
  x[, 1:(ncol(x) - 2)]
}))

high_range_color <- colorRampPalette(c("#D73027", "#9F0C04"))(70)
middle           <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(120)
custom_colors    <- c(middle, high_range_color)

pheatmap(all_forHeatmap_norm,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         color = custom_colors)