library(Seurat)
library(patchwork)
library(ggplot2)

make_first_uppercase <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

run_seurat_pipeline <- function(obj, dims_pca = 1:30, dims_nn = 1:15, resolution = 0.8, nfeatures = 5000) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  obj <- FindNeighbors(obj, dims = dims_nn)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = dims_pca)
  return(obj)
}

set.seed(123)
setwd("path/to/MFG_scRNA")

pbmc.data <- Read10X(data.dir = "path/to/M_limbE10.5/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "M_E10.5", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
pbmc <- run_seurat_pipeline(pbmc)
DimPlot(pbmc, label = TRUE, label.box = TRUE) + NoLegend()

clusters_of_interest <- c(0, 1, 2, 3, 4, 6, 7, 9, 10, 13)
pbmc_subset <- subset(pbmc, subset = seurat_clusters %in% clusters_of_interest)
pbmc_subset <- run_seurat_pipeline(pbmc_subset)
saveRDS(pbmc_subset, "E10.5_mesen.rds")

pbmc_subset <- readRDS("E10.5_mesen.rds")
pbmc_subset@reductions[["umap"]]@cell.embeddings <- -pbmc_subset@reductions[["umap"]]@cell.embeddings
DimPlot(pbmc_subset, label = TRUE, label.box = TRUE) + NoLegend()

pbmc.markers <- FindAllMarkers(pbmc_subset, only.pos = TRUE)
write.csv(pbmc.markers, "01.mesenchymal_E10.5_markers.csv")

top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

DotPlot(pbmc_subset, features = unique(top_markers$gene)) + RotatedAxis()

GeneFeatures <- c("Cyp26b1","Lmo1","Sall3","Lmo2","Hoxd13","Hoxa13","Socs2","Jag2",
                  "Scx","Gsc","Meis2","Pkdcc","Irx3",
                  "Pax1","Alx4","Asb4","Rspo4","Msx2","Lhx2","Prrx2","Msx1",
                  "Shh","Osr1","Hand2","Bmp4","Ptch1","Rspo3")

DotInfo <- DotPlot(pbmc_subset, features = GeneFeatures, dot.scale = 10) + RotatedAxis()
DotInfo$data$id <- factor(DotInfo$data$id, levels = rev(c(12,3,1,5,11,10,9,2,8,7,4,6,0)))

ggplot(DotInfo$data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "#9F0C04", na.value = "transparent") +
  scale_size(name = "Size", range = c(0, 7)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14))

FeaturePlot(pbmc_subset, features = c("Shh","Osr1","Hand2","Bmp4","Ptch1","Rspo3"),
            cols = c("grey","red"), reduction = "umap")

FeaturePlot(pbmc_subset, features = c("Pax1","Alx4","Asb4","Rspo4","Msx2","Lhx2","Prrx2","Msx1"),
            reduction = "umap")

mfg_genes <- c("ALX1","CRLF1","GREM1","GSC","HAND2",
               "HOXD11","HOXD12","HOXD4","HOXD9","LHX9",
               "MAB21L2","MSX1","RSPO4","SHOX2","TBX2","TBX3")
FeaturePlot(pbmc_subset, features = make_first_uppercase(mfg_genes), reduction = "umap")