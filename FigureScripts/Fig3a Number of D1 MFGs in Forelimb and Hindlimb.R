library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cowplot)

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

orth <- readRDS("orth.rds")
orth_genes <- orth$AllGenes

MAXgenes_orth <- MAXgenes_processed[MAXgenes_processed$gene %in% orth_genes,]
MINgenes_orth <- MINgenes_processed[MINgenes_processed$gene %in% orth_genes,]
MAXgenes_orth_single <- MAXgenes_orth[which(MAXgenes_orth$species_num==1),]
MINgenes_orth_single <- MINgenes_orth[which(MINgenes_orth$species_num==1),]

species_list <- c("Tur", "Cor", "Chi", "Emu", "Ost")
results <- list()
for (spec in species_list) {
  max_subset <- MAXgenes_orth_single[MAXgenes_orth_single$species == spec, ]
  min_subset <- MINgenes_orth_single[MINgenes_orth_single$species == spec, ]
  combined <- rbind.data.frame(max_subset, min_subset)
  results[[spec]] <- combined
}
LBGs <- read.csv("00.All_LBGs.csv")

matrix_list <- list(Tur_filtered, Cor_filtered, Chi_filtered, Emu_filtered, Ost_filtered)
matrix_names <- c( "Tur", "Cor", "Chi", "Emu", "Ost")

overlap_list <- list()
overlap_un_list <- list()
frq_list <- list()
for (i in seq_along(matrix_list)) {
  overlap <- merge(matrix_list[[i]], LBGs, by = "gene")
  overlap_un <- overlap[overlap$Unique_Sp == matrix_names[i], ]
  overlap_un_unique <- overlap_un %>%
    distinct(gene, .keep_all = TRUE)
  overlap_list[[matrix_names[i]]] <- overlap
  overlap_un_list[[matrix_names[i]]] <- overlap_un_unique
  frq <- as.data.frame(table(overlap_un_unique$Tau_info))
  frq_list[[matrix_names[i]]] <- frq
}


ggplot(frq_list[Tur], aes(fill=Part, y=Value, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Value", x = "Group", title = "")
ggplot(frq_list[Chi], aes(fill=Part, y=Value, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Value", x = "Group", title = "")
ggplot(frq_list[Emu], aes(fill=Part, y=Value, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Value", x = "Group", title = "")
ggplot(frq_list[Ost], aes(fill=Part, y=Value, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Value", x = "Group", title = "")
ggplot(frq_list[Cor], aes(fill=Part, y=Value, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  labs(y = "Value", x = "Group", title = "")