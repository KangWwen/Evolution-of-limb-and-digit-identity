library(dplyr)
library(stringr)
library(tidyverse)
library(pheatmap)
library(grDevices)
library(RColorBrewer)

file_names <- c("Axm","Mm", "Ac","Tur","Am", "Cor","Chi", "Emu", "Ost")

for (name in file_names) {
  assign(name, read.table(paste0("03.", name, "_AV_TPM.txt"), sep = "\t", row.names = 1, header = TRUE))
}

processLimbData <- function(df, filteredGenes, numColumns) {
  template_df <- data.frame(Gene = filteredGenes, matrix(0, nrow = length(filteredGenes), ncol = numColumns))
  rownames(template_df) <- filteredGenes
  df_subset <- df[rownames(df) %in% filteredGenes, 1:numColumns]
  df_subset$Gene <- rownames(df_subset)
  result_df <- merge(template_df, df_subset, by = "Gene", all = TRUE)
  result_df <- result_df[, c(1, (numColumns+2):(numColumns*2+1))]
  result_df[is.na(result_df)] <- 0
  rownames(result_df) <- result_df$Gene
  result_df$Gene <- NULL
  return(result_df)
}
TurLimb <- processLimbData(Tur, filteredAll, length(Tur))
CorLimb <- processLimbData(Cor, filteredAll, length(Cor))
AxmLimb <- processLimbData(Axm, filteredAll, length(Axm))

EmuLimb <- Emu[rownames(Emu) %in% filteredAll,]
ChiLimb <- Chi[rownames(Chi) %in% filteredAll,]
OstLimb <- processLimbData(Ost, filteredAll, length(Ost))

# Convert from wide to long format
data_long <- all_forHeatmap %>%
  pivot_longer(cols = -gene, names_to = "condition", values_to = "value") %>%
  separate(condition, into = c("limb", "type"), sep = "_") %>%
  unite("gene_limb", gene, limb, sep = "_")

# Spread the 'type' (M, T) into separate columns
data_wide <- data_long %>%
  pivot_wider(names_from = type, values_from = value)
data_wide1 <- as.data.frame(data_wide[c(1:22,49:54,57:58),2:19])

rownames(data_wide1) <- data_wide[c(1:22,49:54,57:58),]$gene_limb
head(data_wide1)

high_range_color <- colorRampPalette(c("#D73027","#9F0C04"))(10)                           
low_range_color <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "RdYlBu")))(120)                         
custom_colors <- c(low_range_color,high_range_color) 

pheatmap(log2(data_wide1+1), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none", 
         color = custom_colors)
