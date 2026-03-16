library(matrixStats)
library(ggplot2)
library(gridExtra)
extract_stage <- function(raw_list, sp, cols, new_names = NULL) {
  df <- raw_list[[sp]][, cols, drop = FALSE]
  colnames(df) <- if (is.null(new_names)) paste0("D", 1:ncol(df)) else new_names
  return(df)
}

file_names <- c("Axm", "Mm", "Tur", "Cor", "Chi", "Emu", "Ost")
raw_data <- setNames(lapply(file_names, function(x) {
  read.table(paste0("03.", x, "_AV_TPM.txt"), sep = "\t", row.names = 1, header = TRUE)
}), file_names)

common_genes <- Reduce(intersect, lapply(raw_data, rownames))
raw_data <- lapply(raw_data, function(df) df[common_genes, , drop = FALSE])

get_MFG_genes <- function(df) {
  mat <- as.matrix(df)
  row_maxs <- rowMaxs(mat)
  mat <- mat[row_maxs >= 3, , drop = FALSE]
  if (nrow(mat) == 0) return(character(0))
  n <- ncol(mat)
  gene_names <- rownames(mat)
  row_maxs_filt <- rowMaxs(mat)
  tau_max <- rowSums(1 - (mat / row_maxs_filt)) / (n - 1)
  genes_max <- gene_names[tau_max > 0.5 & !is.na(tau_max)]
  row_mins <- rowMins(mat)
  row_mins[row_mins == 0] <- 1e-6
  tau_min <- 1 - (1 / (rowSums(mat / row_mins) / (n - 1)))
  genes_min <- gene_names[tau_min > 0.7 & !is.na(tau_min)]
  return(unique(c(genes_max, genes_min)))
}

final_data_fixed <- list(
  eChiF = extract_stage(raw_data, "Chi", c("eCFC2","eCFC3","eCFC4"),c("D1","D3","D4")),
  lChiF = extract_stage(raw_data, "Chi", c("lCFC2","lCFC3","lCFC4"),c("D1","D3","D4")),
  eChiH = extract_stage(raw_data, "Chi", c("eCHC1","eCHC3","eCHC4"),c("D1","D3","D4")),
  lChiH = extract_stage(raw_data, "Chi", c("lCHC1","lCHC3","lCHC4"),c("D1","D3","D4")),
  eEmuF = extract_stage(raw_data, "Emu", c("eEFC2","eEFC3","eEFC4"),c("D1","D3","D4")),
  lEmuF = extract_stage(raw_data, "Emu", c("lEFC2","lEFC3","lEFC4"),c("D1","D3","D4")),
  eEmuH = extract_stage(raw_data, "Emu", c("eEHC2","eEHC3","eEHC4"),c("D1","D3","D4")),
  lEmuH = extract_stage(raw_data, "Emu", c("lEHC2","lEHC3","lEHC4"),c("D1","D3","D4")),
  eOstF = extract_stage(raw_data, "Ost", c("eOFC2","eOFC3","eOFC4"),c("D1","D3","D4")),
  lOstF = extract_stage(raw_data, "Ost", c("lOFC2","lOFC3","lOFC4"),c("D1","D3","D4")),
  eOstH = extract_stage(raw_data, "Ost", c("eOHC1","eOHC3","eOHC4"),c("D1","D3","D4")),
  lOstH = extract_stage(raw_data, "Ost", c("lOHC1","lOHC3","lOHC4"),c("D1","D3","D4")),
  eAxmF = extract_stage(raw_data, "Axm", c("Axm_FD1_52","Axm_FD3_52","Axm_FD4_52"), c("D1","D3","D4")),
  lAxmF = extract_stage(raw_data, "Axm", c("Axm_FD1_55","Axm_FD3_55","Axm_FD4_55"), c("D1","D3","D4")),
  eAxmH = extract_stage(raw_data, "Axm", c("Axm_HD1_55","Axm_HD3_55","Axm_HD4_55"), c("D1","D3","D4")),
  lAxmH = extract_stage(raw_data, "Axm", c("Axm_HD1_57","Axm_HD3_57","Axm_HD4_57"), c("D1","D3","D4")),
  eCorF = extract_stage(raw_data, "Cor", c("eFCor1","eFCor3","eFCor4"),c("D1","D3","D4")),
  lCorF = extract_stage(raw_data, "Cor", c("lFCor1","lFCor3","lFCor4"),c("D1","D3","D4")),
  eCorH = extract_stage(raw_data, "Cor", c("eHCor1","eHCor3","eHCor4"),c("D1","D3","D4")),
  lCorH = extract_stage(raw_data, "Cor", c("lHCor1","lHCor3","lHCor4"),c("D1","D3","D4")),
  eTurF = extract_stage(raw_data, "Tur", c("eTFC1","eTFC3","eTFC4"),c("D1","D3","D4")),
  lTurF = extract_stage(raw_data, "Tur", c("lTFC1","lTFC3","lTFC4"),c("D1","D3","D4")),
  eTurH = extract_stage(raw_data, "Tur", c("eTHC1","eTHC3","eTHC4"),c("D1","D3","D4")),
  lTurH = extract_stage(raw_data, "Tur", c("lTHC1","lTHC3","lTHC4"),c("D1","D3","D4")),
  eMmF  = extract_stage(raw_data, "Mm",  c("MmFC1","MmFC3","MmFC4"),c("D1","D3","D4"))
)

species_groups <- list(
  Chi = c("eChiF","lChiF","eChiH","lChiH"),
  Emu = c("eEmuF","lEmuF","eEmuH","lEmuH"),
  Ost = c("eOstF","lOstF","eOstH","lOstH"),
  Axm = c("eAxmF","lAxmF","eAxmH","lAxmH"),
  Cor = c("eCorF","lCorF","eCorH","lCorH"),
  Tur = c("eTurF","lTurF","eTurH","lTurH"),
  Mm  = c("eMmF")
)

species_type <- c(Chi="Avian", Emu="Avian", Ost="Avian",
                  Axm="Reptile", Cor="Reptile", Tur="Reptile", Mm="Reptile")

results_fixed <- do.call(rbind, lapply(names(species_groups), function(group_name) {
  all_genes <- unique(unlist(lapply(species_groups[[group_name]], function(stage) {
    if (stage %in% names(final_data_fixed)) get_MFG_genes(final_data_fixed[[stage]])
  })))
  data.frame(Species = group_name, MFG_Count = length(all_genes),
             Type = species_type[group_name], Method = "Fixed (D1,D3,D4)",
             stringsAsFactors = FALSE)
}))

digit_columns <- list(
  eChiF = c("eCFC2","eCFC3","eCFC4"),lChiF = c("lCFC2","lCFC3","lCFC4"),
  eChiH = c("eCHC1","eCHC3","eCHC4"),lChiH = c("lCHC1","lCHC3","lCHC4"),
  eEmuF = c("eEFC2","eEFC3","eEFC4"),lEmuF = c("lEFC2","lEFC3","lEFC4"),
  eEmuH = c("eEHC2","eEHC3","eEHC4"),lEmuH = c("lEHC2","lEHC3","lEHC4"),
  eOstF = c("eOFC2","eOFC3","eOFC4"), lOstF = c("lOFC2","lOFC3","lOFC4"),
  eOstH = c("eOHC1","eOHC3","eOHC4"), lOstH = c("lOHC1","lOHC3","lOHC4"),
  eAxmF = c("Axm_FD1_52","Axm_FD2_52","Axm_FD3_52","Axm_FD4_52"),
  lAxmF = c("Axm_FD1_55","Axm_FD2_55","Axm_FD3_55","Axm_FD4_55"),
  eAxmH = c("Axm_HD1_55","Axm_HD2_55","Axm_HD3_55","Axm_HD4_55","Axm_HD5_55"),
  lAxmH = c("Axm_HD1_57","Axm_HD2_57","Axm_HD3_57","Axm_HD4_57","Axm_HD5_57"),
  eCorF = c("eFCor1","eFCor2","eFCor3","eFCor4","eFCor5"),
  lCorF = c("lFCor1","lFCor2","lFCor3","lFCor4","lFCor5"),
  eCorH = c("eHCor1","eHCor2","eHCor3","eHCor4","eHCor5"),
  lCorH = c("lHCor1","lHCor2","lHCor3","lHCor4","lHCor5"),
  eMmF  = c("MmFC1","MmFC2","MmFC3","MmFC4","MmFC5"),
  eTurF = c("eTFC1","eTFC2","eTFC3","eTFC4","eTFC5"),
  lTurF = c("lTFC1","lTFC2","lTFC3","lTFC4","lTFC5"),
  eTurH = c("eTHC1","eTHC2","eTHC3","eTHC4","eTHC5"),
  lTurH = c("lTHC1","lTHC2","lTHC3","lTHC4","lTHC5")
)

need_sampling <- sapply(digit_columns, length) > 3

set.seed(123)
n_sim <- 1000

results_random <- do.call(rbind, lapply(names(species_groups), function(group_name) {
  stages <- species_groups[[group_name]]
  sim_counts <- sapply(1:n_sim, function(i) {
    all_genes <- unique(unlist(lapply(stages, function(stage) {
      if (!stage %in% names(digit_columns)) return(character(0))
      cols <- digit_columns[[stage]]
      selected_cols <- if (need_sampling[stage]) sample(cols, 3) else cols
      sp <- sub("[FH]$", "", sub("^[el]", "", stage))
      temp_df <- extract_stage(raw_data, sp, selected_cols, c("D1","D2","D3"))
      get_MFG_genes(temp_df)
    })))
    length(all_genes)
  })
  data.frame(Species = group_name, MFG_Count = sim_counts,
             Type = species_type[group_name], Method = "Random (3 digits)",
             stringsAsFactors = FALSE)
}))

plot_data_fixed  <- results_fixed
plot_data_random <- results_random