library(DESeq2)
run_DESeq2_auto <- function(input_file, group1_cols, group2_cols, 
                           group1_name, group2_name, output_prefix,
                           padj_cutoff = 0.05, log2fc_cutoff = 1) {
  
  data <- read.table(input_file, header = TRUE, sep = "\t")
  
  Count <- as.data.frame(cbind(data[, group1_cols], data[, group2_cols]))
  rownames(Count) <- data$Geneid
  colnames(Count) <- c(paste0(group1_name, seq_along(group1_cols)), 
                       paste0(group2_name, seq_along(group2_cols)))
  
  Count_condition <- factor(c(rep(group1_name, length(group1_cols)), 
                              rep(group2_name, length(group2_cols))))
  coldata <- data.frame(row.names = colnames(Count), Count_condition)
  
  dds <- DESeqDataSetFromMatrix(countData = Count, colData = coldata, 
                                design = ~Count_condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  resOrdered <- res[order(res$padj), ]
  resdata <- merge(as.data.frame(resOrdered), 
                   as.data.frame(counts(dds, normalized = TRUE)), 
                   by = "row.names", sort = FALSE)
  names(resdata)[1] <- "Gene"
  
  resdata_filter <- subset(resdata, padj < padj_cutoff & 
                          abs(log2FoldChange) > log2fc_cutoff)
  
  write.table(resdata, file = paste0(output_prefix, "_info.txt"), 
              sep = "\t", row.names = FALSE)
  write.table(resdata_filter, file = paste0(output_prefix, "_DE.txt"), 
              sep = "\t", row.names = FALSE)
  
  cat("Total DE genes:", nrow(resdata_filter), "\n")
  return(list(all = resdata, DE = resdata_filter, dds = dds))
}
