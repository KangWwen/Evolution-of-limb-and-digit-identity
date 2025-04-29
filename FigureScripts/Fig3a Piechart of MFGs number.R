library(dplyr)
library(tidyr)
library(stringr)

aggregate_counts <- function(df, counts_column_name) {
  if (!counts_column_name %in% names(df)) {
    stop("The specified counts column does not exist in the dataframe.")
  }
  df_processed <- df %>%
    separate_rows(!!sym(counts_column_name), sep = "; ") %>%
    mutate(Identifier = str_extract(!!sym(counts_column_name), "^[A-Za-z]+\\d+"),
           Count = as.integer(str_extract(!!sym(counts_column_name), "\\d+$"))) %>%
    select(-!!sym(counts_column_name)) %>%
    group_by(Identifier) %>%
    summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = 'drop')
  return(df_processed)
}

process_cell <- function(cell) {
  elements <- unlist(strsplit(cell, ";"))
  last_three_chars <- substr(elements, nchar(elements)-2, nchar(elements))
  counts <- table(last_three_chars)
  result <- paste(names(counts), counts, sep=":", collapse="; ")
  return(result)
}

MAXgenes_processed <- MAXgenes %>%
  rowwise() %>%
  mutate(Counts = process_cell(Tau_names))
MINgenes_processed <- MINgenes %>%
  rowwise() %>%
  mutate(Counts = process_cell(Tau_names))

orth_genes <- orth$AllGenes

MAXgenes_orth <- MAXgenes_processed[MAXgenes_processed$gene %in% orth_genes,]
MINgenes_orth <- MINgenes_processed[MINgenes_processed$gene %in% orth_genes,]
MAXgenes_orth_single <- MAXgenes_orth[which(MAXgenes_orth$species_num==1),]
MINgenes_orth_single <- MINgenes_orth[which(MINgenes_orth$species_num==1),]

species_list <- c("Axm", "Tur", "Cor", "Chi", "Emu", "Ost")
results <- list()
for (spec in species_list) {
  max_subset <- MAXgenes_orth_single[MAXgenes_orth_single$species == spec, ]
  min_subset <- MINgenes_orth_single[MINgenes_orth_single$species == spec, ]
  combined <- rbind.data.frame(max_subset, min_subset)
  results[[spec]] <- combined
}

library(stringr)  
count_elements <- function(df) {
  counts_names <- sapply(str_split(df$Counts, ":"), `[`, 1)
  counts_freq <- table(counts_names)
  return(counts_freq)
}
datasets <- c("Axm", "Tur", "Cor", "Chi", "Emu", "Ost")
counts_list <- lapply(results[names(results) %in% datasets], count_elements)
names(counts_list) <- datasets
head(counts_list)

library(ggplot2)
library(gridExtra)
custom_colors <- c("#FFCE7E", "#F57E1F", "#57C9FB","#4575B4","#00870E")
draw_pie <- function(counts, title) {
  df <- as.data.frame(counts)
  names(df) <- c("name", "value")
  p <- ggplot(df, aes(x = "", y = value, fill = name)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y") +
    theme_void() +
    geom_text(aes(label = value), position = position_stack(vjust = 0.5), color = "white") +
    scale_fill_manual(values = custom_colors) +
    labs(title = title, fill = "Name")+theme(legend.position = "none") 
  
  return(p)
}

plots <- list()

for (species in names(counts_list)) {
  fd_counts <- counts_list[[species]][grepl("^FD", names(counts_list[[species]]))]
  hd_counts <- counts_list[[species]][grepl("^HD", names(counts_list[[species]]))]
  fd_plot <- draw_pie(fd_counts, paste(species, "FD"))
  hd_plot <- draw_pie(hd_counts, paste(species, "HD"))
  plots[[paste(species, "FD")]] <- fd_plot
  plots[[paste(species, "HD")]] <- hd_plot
}
