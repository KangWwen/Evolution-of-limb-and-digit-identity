library(pheatmap)
raw_names <- c("eChiF","lChiF","eChiH","lChiH", "eEmuF","lEmuF","eEmuH","lEmuH",
               "eOstF","lOstF","eOstH","lOstH", "eCorF","lCorF","eCorH","lCorH",
               "eTurF","lTurF","eTurH","lTurH", "eAxmH","lAxmH", "eAxmF","lAxmF")

tau_data_list <- list()
for (n in raw_names) {
  file_name <- paste0(n, "_tau.txt")
  if (file.exists(file_name)) {
    tau_data_list[[n]] <- read.table(file_name, sep = "\t", row.names = 1, header = TRUE)
  }
}

get_d_genes <- function(df, prefix) {
  lapply(1:5, function(i) {
    pattern <- paste0(prefix, "_ D", i)
    unique(df$gene[df$Tau_info == pattern])
  })
}

samples <- c("ChiF","ChiH","EmuF","EmuH","OstF","OstH","CorF","CorH","TurF","TurH","AxmH","AxmF")
data_tree <- list()

for (s in samples) {
  for (phase in c("e", "l")) {
    raw_key <- paste0(phase, s)
    if (!is.null(tau_data_list[[raw_key]])) {
      df <- tau_data_list[[raw_key]]
      data_tree[[s]][[phase]][["MAX"]] <- get_d_genes(df, "MAX")
      data_tree[[s]][[phase]][["MIN"]] <- get_d_genes(df, "MIN")
    }
  }
}

calc_integrated_overlap <- function(src_name, src_d_idx, tgt_name) {
  sapply(1:5, function(tgt_d_idx) {
    all_conserved_genes <- c()
    for (phase in c("e", "l")) {
      if (!is.null(data_tree[[src_name]][[phase]]) && !is.null(data_tree[[tgt_name]][[phase]])) {
        min_inter <- intersect(data_tree[[src_name]][[phase]][["MIN"]][[src_d_idx]],
                               data_tree[[tgt_name]][[phase]][["MIN"]][[tgt_d_idx]])
        max_inter <- intersect(data_tree[[src_name]][[phase]][["MAX"]][[src_d_idx]],
                               data_tree[[tgt_name]][[phase]][["MAX"]][[tgt_d_idx]])
        all_conserved_genes <- c(all_conserved_genes, min_inter, max_inter)
      }
    }
    return(length(unique(all_conserved_genes)))
  })
}

comparison_list <- list(
  ChiFC2 = c("ChiF", 1, "CorF"), ChiFC3 = c("ChiF", 3, "CorF"), ChiFC4 = c("ChiF", 4, "CorF"),
  EmuFC2 = c("EmuF", 1, "CorF"), EmuFC3 = c("EmuF", 3, "CorF"), EmuFC4 = c("EmuF", 4, "CorF"),
  OstFC2 = c("OstF", 1, "CorF"), OstFC3 = c("OstF", 3, "CorF"), OstFC4 = c("OstF", 4, "CorF"),
  AxmFC1 = c("AxmF", 1, "CorF"), AxmFC2 = c("AxmF", 2, "CorF"), AxmFC3 = c("AxmF", 3, "CorF"), AxmFC4 = c("AxmF", 4, "CorF"),
  ChiHC1 = c("ChiH", 1, "CorH"), ChiHC2 = c("ChiH", 2, "CorH"), ChiHC3 = c("ChiH", 3, "CorH"), ChiHC4 = c("ChiH", 4, "CorH"),
  EmuHC2 = c("EmuH", 1, "CorH"), EmuHC3 = c("EmuH", 3, "CorH"), EmuHC4 = c("EmuH", 4, "CorH"),
  OstHC1 = c("OstH", 1, "CorH"), OstHC2 = c("OstH", 2, "CorH"), OstHC3 = c("OstH", 3, "CorH"), OstHC4 = c("OstH", 4, "CorH"), OstHC5 = c("OstH", 5, "CorH"),
  AxmHC1 = c("AxmH", 1, "CorH"), AxmHC2 = c("AxmH", 2, "CorH"), AxmHC3 = c("AxmH", 3, "CorH"), AxmHC4 = c("AxmH", 4, "CorH"), AxmHC5 = c("AxmH", 5, "CorH")
)

tur_comparison_list <- list()
for (n in names(comparison_list)) {
  params <- comparison_list[[n]]
  target_tur <- ifelse(grepl("F$", params[1]), "TurF", "TurH")
  tur_comparison_list[[paste0(n, "_vs_Tur")]] <- c(params[1], params[2], target_tur)
}

final_task_list <- c(comparison_list, tur_comparison_list)

results_matrix <- t(sapply(final_task_list, function(x) {
  calc_integrated_overlap(x[1], as.numeric(x[2]), x[3])
}))
colnames(results_matrix) <- paste0("D", 1:5)

calc_integrated_jaccard <- function(src_name, src_d_idx, tgt_name) {
  sapply(1:5, function(tgt_d_idx) {
    all_src_genes <- c()
    all_tgt_genes <- c()
    all_intersect_genes <- c()
    for (phase in c("e", "l")) {
      if (!is.null(data_tree[[src_name]][[phase]]) && !is.null(data_tree[[tgt_name]][[phase]])) {
        src_set <- unique(c(data_tree[[src_name]][[phase]][["MIN"]][[src_d_idx]],
                            data_tree[[src_name]][[phase]][["MAX"]][[src_d_idx]]))
        tgt_set <- unique(c(data_tree[[tgt_name]][[phase]][["MIN"]][[tgt_d_idx]],
                            data_tree[[tgt_name]][[phase]][["MAX"]][[tgt_d_idx]]))
        all_src_genes <- unique(c(all_src_genes, src_set))
        all_tgt_genes <- unique(c(all_tgt_genes, tgt_set))
        min_inter <- intersect(data_tree[[src_name]][[phase]][["MIN"]][[src_d_idx]],
                               data_tree[[tgt_name]][[phase]][["MIN"]][[tgt_d_idx]])
        max_inter <- intersect(data_tree[[src_name]][[phase]][["MAX"]][[src_d_idx]],
                               data_tree[[tgt_name]][[phase]][["MAX"]][[tgt_d_idx]])
        all_intersect_genes <- unique(c(all_intersect_genes, min_inter, max_inter))
      }
    }
    union_set <- unique(c(all_src_genes, all_tgt_genes))
    if (length(union_set) == 0) return(0)
    return(length(all_intersect_genes) / length(union_set))
  })
}

results_matrix_jaccard <- t(sapply(final_task_list, function(x) {
  calc_integrated_jaccard(x[1], as.numeric(x[2]), x[3])
}))
colnames(results_matrix_jaccard) <- paste0("D", 1:5)

my_colors_j <- colorRampPalette(c("white", brewer.pal(n = 7, name = "OrRd")))(100)
my_breaks_j <- seq(0, max(results_matrix_jaccard), length.out = 100)

pheatmap(results_matrix_jaccard,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = my_colors_j,
         main = "Integrated Jaccard Similarity (Normalized)",
         display_numbers = FALSE,
         number_format = "%.3f",
         fontsize_number = 7)