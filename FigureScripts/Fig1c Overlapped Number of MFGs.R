library(RColorBrewer)
library(pheatmap)
file_names <- c("eChiF","lChiF","eEmuF","lEmuF",
                "eEmuH","lEmuH","eOstF","lOstF",
                "eAxmF", "lAxmF", "eChiH","lChiH",   
                "eAxmH","lAxmH","eCorF","eCorH",
                "lCorF","lCorH","eAcF","eAmF","lAmF","eMmF",
                "eTurF","lTurF","eTurH",
                "lTurH","eOstH","lOstH")
for (name in file_names) {
  assign(name, read.table(paste0(name, "_tau.txt"), sep = "\t", row.names = 1, header = TRUE))
}
# List of paired data frame names
pairs <- list(c("eChiF","lChiF"), c("eEmuF","lEmuF"), c("eEmuH","lEmuH"), 
              c("eOstF","lOstF"), c("eAxmF", "lAxmF"), c("eChiH","lChiH"),
              c("eAxmH","lAxmH"), c("eCorF","lCorF"), c("eCorH","lCorH"),  
              c("eAmF","lAmF"), c("eTurF","lTurF"), c("eTurH", "lTurH"),
              c("eOstH","lOstH"))

# Merge each pair
merged <- lapply(pairs, function(x) {
  df1 <- get(x[1]) 
  df2 <- get(x[2])
  cbind(gene = c(df1$gene, df2$gene),
        Tau_info = c(df1$Tau_info, df2$Tau_info))
})

names(merged) <- sapply(pairs, function(x) paste0(x, collapse = "_"))
for(name in names(merged)) {
  assign(name, merged[[name]])
}

samples <- c("eChiF_lChiF","eEmuF_lEmuF","eEmuH_lEmuH","eOstF_lOstF","eAxmF_lAxmF",
             "eChiH_lChiH","eAxmH_lAxmH","eCorF_lCorF","eCorH_lCorH","eAmF_lAmF", 
             "eTurF_lTurF","eTurH_lTurH","eOstH_lOstH",
             "eAcF","eMmF")

for(sample in samples){
  sample_df <- as.data.frame(get(sample))
  sample_split <- lapply(sample_df, function(x) x)
  max_matches <- sapply(sample_split, grepl, pattern="MAX_")
  max_data <- sample_df[rowSums(max_matches) > 0,]
  max_data <- max_data[!duplicated(max_data$gene),]
  min_matches <- sapply(sample_split, grepl, pattern="MIN_")
  min_data <- sample_df[rowSums(min_matches) > 0,]
  min_data <- min_data[!duplicated(min_data$gene),]
  assign(paste0(sample, "_Max"), max_data)
  assign(paste0(sample, "_Min"), min_data)
}

extract_data <- function(data, tau_info_prefix) {
  tau_info_values <- c("D1", "D3", "D4")
  result <- list()
  for (tau_info in tau_info_values) {
    result[[tau_info]] <- data[which(data$Tau_info == paste(tau_info_prefix, tau_info, sep = "")), ]
  }
  return(result)
}
calculate_num <- function(D_Max, D_Min) {
  intersections <- list()
  for (i in 1:5) {
    D <- paste0("D", i)
    intersections[[D]] <- c(D_Max[[D]]$gene,
                            D_Min[[D]]$gene)
  }
  return(sapply(intersections, length))
}

#################################
calculate_intersections <- function(Chi_Max, Chi_Min, D_Max, D_Min) {
  intersections <- list()
  for (i in 1:5) {
    D <- paste0("D", i)
    intersections[[D]] <- c(intersect(Chi_Max$gene, D_Max[[D]]$gene),
                            intersect(Chi_Min$gene, D_Min[[D]]$gene))
  }
  return(sapply(intersections, length))
}


high_range_color <- colorRampPalette(c("#D73027", "#9F0C04"))(10)
middle <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(120)
custom_colors <- c(middle, high_range_color)

MergeAllSp <- rbind.data.frame(ChiAllF_norm,
                               ChiAllH_norm,
                               EmuAllF_norm,
                               EmuAllH_norm,
                               OstAllF_norm,
                               OstAllH_norm,
                               AxmAllF_norm,
                               AxmAllH_norm)
pheatmap(MergeAllSp, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none", 
         color = custom_colors)