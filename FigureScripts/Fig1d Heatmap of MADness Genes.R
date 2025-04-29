library(dplyr)
library(stringr)
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

# Apply the function to different data frames
TurLimb <- processLimbData(Tur, filteredAll, length(Tur))
CorLimb <- processLimbData(Cor, filteredAll, length(Cor))
AxmLimb <- processLimbData(Axm, filteredAll, length(Axm))

EmuLimb <- processLimbData(Emu, filteredAll, length(Emu))
ChiLimb <- processLimbData(Chi, filteredAll, length(Chi))
OstLimb <- processLimbData(Ost, filteredAll, length(Ost))
MmLimb <- processLimbData(Mm, filteredAll, length(Mm))

for (name in elements_names) {
  assign(name, as.data.frame(t(get(name))))
}

library(pheatmap)
elements <- c("AxmF","AxmH","TurF","TurH","CorF","CorH","ChiF","ChiH","EmuF","EmuH","OstF","OstH")

# Function to process each string
process_string <- function(s) {
  # Split the string into individual elements
  split_elements <- strsplit(s, ";")[[1]]
  # Initialize a list to keep track of the processed elements
  processed_elements <- setNames(rep(NA, length(elements)), paste0("e", elements))
  # Loop over each split element
  for (elem in split_elements) {
    # Extract the key part of the element (i.e., without the 'e' or 'l' prefix)
    key <- substr(elem, 2, nchar(elem))
    # Check if this element is one of the 12 we're interested in
    if (key %in% elements) {
      # Determine the current and new statuses
      current_status <- processed_elements[[paste0("e", key)]]
      new_status <- substr(elem, 1, 1)
      # Update the status, preferring 'e' over 'l'
      if (is.na(current_status) || new_status == "e") {
        processed_elements[[paste0("e", key)]] <- new_status
      }
    }
  }
  
  # Replace NA values with "e", indicating missing elements should be "e"
  processed_elements[is.na(processed_elements)] <- "e"
  # Combine the status with the elements to get the final transformed elements
  transformed_elements <- paste0(processed_elements, elements)
  return(transformed_elements)
}

# Apply the function to each element in your data
transformed_data <- sapply(stage_info_Max$names, process_string)

colnames(transformed_data) <- stage_info_Max$gene
combined_data <-as.data.frame(t(transformed_data))

high_range_color <- colorRampPalette(c("#D73027", "#9F0C04"))(70)
middle <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(120)
custom_colors <- c(middle, high_range_color)
# Draw the heatmap
pheatmap(combined_data, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none", 
         color = custom_colors)