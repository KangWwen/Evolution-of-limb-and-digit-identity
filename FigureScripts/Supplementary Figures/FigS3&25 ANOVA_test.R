library(edgeR)
library(ggplot2)

setwd("path/to/data")
data <- read.table("02.Read_Counts.txt", header = TRUE, row.names = 1)

run_edgeR_anova <- function(counts, group, stage = NULL, coef, fc_cols, out_file, title) {
  y <- DGEList(counts = counts, group = group)
  y <- y[rowSums(cpm(y) > 1) >= 2, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- if (!is.null(stage)) model.matrix(~group + stage) else model.matrix(~group)
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = coef)
  fc_filter <- Reduce(`|`, lapply(fc_cols, function(col) abs(lrt$table[[col]]) > 1))
  filtered <- lrt$table[fc_filter, ]
  genes <- rownames(filtered[filtered$PValue < 0.05, ])
  print(
    qplot(filtered$PValue, geom = "histogram", binwidth = .02, xlab = "p value") +
      theme_bw() + labs(title = title) + ylim(0, 300)
  )
  write.table(genes, out_file, sep = "\t", row.names = TRUE, col.names = TRUE)
  return(genes)
}

run_edgeR_anova(
  counts = data[, c("eCFLC1","eCFLC2","eCFLC3","eCFRC1","eCFRC2","eCFRC3",
                    "lCFLC1","lCFLC2","lCFLC3","lCFRC1","lCFRC2","lCFRC3")],
  group = rep(c("D1","D2","D3"), 4),
  stage = c(rep("early", 6), rep("late", 6)),
  coef = 2:3,
  fc_cols = c("logFC.groupD2","logFC.groupD3"),
  out_file = "01.filter_CF_limbgene_ANOVA.txt", title = "CF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eCHLC1","eCHLC2","eCHLC3","eCHLC4","eCHRC1","eCHRC2","eCHRC3","eCHRC4",
                    "lCHLC1","lCHLC2","lCHLC3","lCHLC4","lCHRC1","lCHRC2","lCHRC3","lCHRC4")],
  group = rep(c("D1","D2","D3","D4"), 4),
  stage = c(rep("early", 8), rep("late", 8)),
  coef = 2:4,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4"),
  out_file = "01.filter_CH_limbgene_ANOVA.txt", title = "CH ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eTFLD1","eTFLD2","eTFLD3","eTFLD4","eTFLD5",
                    "eTFRD1","eTFRD2","eTFRD3","eTFRD4","eTFRD5",
                    "lTFLD1","lTFLD2","lTFLD3","lTFLD4","lTFLD5",
                    "lTFRD1","lTFRD2","lTFRD3","lTFRD4","lTFRD5")],
  group = rep(c("D1","D2","D3","D4","D5"), 4),
  stage = c(rep("early", 10), rep("late", 10)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_TF_limbgene_ANOVA.txt", title = "TF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eTHLD1","eTHLD2","eTHLD3","eTHLD4","eTHLD5",
                    "eTHRD1","eTHRD2","eTHRD3","eTHRD4","eTHRD5",
                    "lTHLD1","lTHLD2","lTHLD3","lTHLD4","lTHLD5",
                    "lTHRD1","lTHRD2","lTHRD3","lTHRD4","lTHRD5")],
  group = rep(c("D1","D2","D3","D4","D5"), 4),
  stage = c(rep("early", 10), rep("late", 10)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_TH_limbgene_ANOVA.txt", title = "TH ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eEFLC2","eEFRC2","eEFLC3","eEFRC3","eEFLC4","eEFRC4",
                    "lEFC2","lEFRC2","lEFLC3","lEFRC3","lEFLC4","lEFRC4")],
  group = rep(c("D2","D2","D3","D3","D4","D4"), 2),
  stage = c(rep("early", 6), rep("late", 6)),
  coef = 2:3,
  fc_cols = c("logFC.groupD3","logFC.groupD4"),
  out_file = "01.filter_EF_limbgene_ANOVA.txt", title = "EF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eEHLC2","eEHRC2","EHLC2_28","eEHLC3","EHLC3_28","eEHRC3","eEHRC4","EHLC4_28",
                    "lEHLC2","NE32HLC2","lEHLC3","NE32HLC3","lEHLC4","NE32HLC4")],
  group = c("D2","D2","D2","D3","D3","D3","D4","D4","D2","D2","D2","D2","D4","D4"),
  stage = c(rep("early", 8), rep("late", 6)),
  coef = 2:3,
  fc_cols = c("logFC.groupD3","logFC.groupD4"),
  out_file = "01.filter_EH_limbgene_ANOVA.txt", title = "EH ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eOFLC2","eOFLC3","eOFLC4","eOFRC2","eOFRC3","eOFRC4",
                    "lOFLC2","lOFLC3","lOFLC4","lOFRC2","lOFRC3","lOFRC4")],
  group = rep(c("D2","D3","D4"), 4),
  stage = c(rep("early", 6), rep("late", 6)),
  coef = 2:3,
  fc_cols = c("logFC.groupD3","logFC.groupD4"),
  out_file = "01.filter_OF_limbgene_ANOVA.txt", title = "OF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eOHLC1","eOHLC2","eOHLC3","eOHLC4","eOHLC5",
                    "eOHRC1","eOHRC2","eOHRC3","eOHRC4","eOHRC5",
                    "lOHLC1","lOHLC2","lOHLC3","lOHLC4","lOHLC5",
                    "lOHRC1","lOHRC2","lOHRC3","lOHRC4","lOHRC5")],
  group = rep(c("D1","D2","D3","D4","D5"), 4),
  stage = c(rep("early", 10), rep("late", 10)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_OH_limbgene_ANOVA.txt", title = "OH ANOVA"
)

run_edgeR_anova(
  counts = data[, c("eAmFD1_1","eAmFD1_2","eAmFD1_3","eAmFD2_1","eAmFD2_2","eAmFD2_3",
                    "eAmFD3_1","eAmFD3_2","eAmFD4_1","eAmFD4_2","eAmFD4_3",
                    "eAmFD5_1","eAmFD5_2","eAmFD5_3","lAmFD1_1","lAmFD1_2","lAmFD1_3",
                    "lAmFD2_1","lAmFD2_2","lAmFD3_1","lAmFD3_2","lAmFD3_3",
                    "lAmFD4_1","lAmFD4_2","lAmFD4_3","lAmFD5_1","lAmFD5_2","lAmFD5_3")],
  group = c("D1","D1","D1","D2","D2","D2","D3","D3","D4","D4","D4","D5","D5","D5",
            "D1","D1","D1","D2","D2","D3","D3","D3","D4","D4","D4","D5","D5","D5"),
  stage = c(rep("early", 14), rep("late", 14)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_AmF_limbgene_ANOVA.txt", title = "AmF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("ac_FL_D1_1","ac_FL_D1_2","ac_FL_D1_3","ac_FL_D2_1","ac_FL_D2_2","ac_FL_D2_3",
                    "ac_FL_D3_1","ac_FL_D3_2","ac_FL_D3_3","ac_FL_D4_1","ac_FL_D4_2","ac_FL_D4_3",
                    "ac_FL_D5_1","ac_FL_D5_2","ac_FL_D5_3")],
  group = c("D1","D1","D1","D2","D2","D2","D3","D3","D3","D4","D4","D3","D5","D5","D5"),
  stage = NULL,
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_AcF_limbgene_ANOVA.txt", title = "AcF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("mm_FL_D1_1","mm_FL_D1_2","mm_FL_D1_3","mm_FL_D2_1","mm_FL_D2_2","mm_FL_D2_3",
                    "mm_FL_D3_1","mm_FL_D3_2","mm_FL_D3_3","mm_FL_D4_1","mm_FL_D4_2","mm_FL_D4_3",
                    "mm_FL_D5_1","mm_FL_D5_2","mm_FL_D5_3")],
  group = c("D1","D1","D1","D2","D2","D2","D3","D3","D3","D4","D4","D3","D5","D5","D5"),
  stage = NULL,
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_MMF_limbgene_ANOVA.txt", title = "MMF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("CorPor_1_A_M18_LFD1","CorPor_6_A_M18_RFD1","CorPor_2_A_M18_LFD2","CorPor_7_A_M18_RFD2",
                    "CorPor_3_A_M18_LFD3","CorPor_8_A_M18_RFD3","CorPor_4_A_M18_LFD4","CorPor_9_A_M18_RFD4",
                    "CorPor_5_A_M18_LFD5","CorPor_10_A_M18_RFD5","X26_M20_FLD1","X31_M20_FRD1",
                    "X27_M20_FLD2","X32_M20_FRD2","X28_M20_FLD3","X33_A_MS20_FRD3",
                    "X29_M20_FLD4","X34_A_MS20_FRD4","X30_M20_FLD5","X35_A_MS20_FRD5")],
  group = rep(c("D1","D1","D2","D2","D3","D3","D4","D4","D5","D5"), 2),
  stage = c(rep("early", 10), rep("late", 10)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_SiF_limbgene_ANOVA.txt", title = "SiF ANOVA"
)

run_edgeR_anova(
  counts = data[, c("CorPor_11_A_M18_LHD1","CorPor_16_A_M18_RHD1","CorPor_12_A_M18_LHD2","CorPor_17_A_M18_RHD2",
                    "CorPor_13_A_M18_LHD3","CorPor_18_A_M18_RHD3","CorPor_14_A_M18_LHD4","CorPor_19_A_M18_RHD4",
                    "CorPor_15_A_M18_LHD5","CorPor_20_A_M18_RHD5","X36_A_MS20_HLD1","X41_A_MS20_HRD1",
                    "X37_A_MS20_HLD2","X42_A_MS20_HRD2","X38_A_MS20_HLD3","X43_A_MS20_HRD3",
                    "X39_A_MS20_HLD4","X44_A_MS20_HRD4","X40_A_MS20_HLD5","X45_A_MS20_HRD5")],
  group = rep(c("D1","D1","D2","D2","D3","D3","D4","D4","D5","D5"), 2),
  stage = c(rep("early", 10), rep("late", 10)),
  coef = 2:5,
  fc_cols = c("logFC.groupD2","logFC.groupD3","logFC.groupD4","logFC.groupD5"),
  out_file = "01.filter_SiH_limbgene_ANOVA.txt", title = "SiH ANOVA"
)