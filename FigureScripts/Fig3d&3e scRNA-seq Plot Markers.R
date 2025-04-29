library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
####anteiror
high_range_color <- colorRampPalette(c("#FFA000"))(10)
#middle <- colorRampPalette(brewer.pal(n = 7, name = "Oranges"))(50)
middle <- colorRampPalette(rev(c("#FFA000","#FEEDDE")))(50)
low <- colorRampPalette(c("#FEEDDE"))(20)
custom_colors <- c(low,middle, high_range_color)
####posterior
high_range_color <- colorRampPalette(c("#4575B4"))(10)
#middle <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50)
middle <- colorRampPalette(rev(c("#4575B4","#EFF3FF")))(50)
low <- colorRampPalette(c("#EFF3FF"))(20)
custom_colors <- c(low,middle, high_range_color)

###distal
high_range_color <- colorRampPalette(c("#006A00"))(10)
#middle <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(50)
middle <- colorRampPalette(rev(c("#006A00","#EFF3FF")))(50)
low <- colorRampPalette(c("#EFF3FF"))(20)
custom_colors <- c(low,middle, high_range_color)

###prorixmal
high_range_color <- colorRampPalette(c("#FF7600"))(10)
#middle <- colorRampPalette(brewer.pal(n = 7, name = "Oranges"))(50)
middle <- colorRampPalette(rev(c("#FF7600","#FEEDDE")))(50)
low <- colorRampPalette(c("#FEEDDE"))(20)
custom_colors <- c(low,middle, high_range_color)

#brewer.pal(n = 7, name = "Blues")
genes_to_plot <- c("Pax1","Alx4","Rspo4","Msx2",
                   "Lhx2","Msx1")
genes_to_plot <- c("Shh","Osr1","Hand2","Bmp4","Ptch1","Rspo3")
genes_to_plot <- c("Cyp26b1","Lmo1","Sall3","Lmo2",
                   "Hoxd13","Hoxa13","Socs2","Jag2")
genes_to_plot <- c("Scx","Gsc","Meis2","Pkdcc","Irx3")

expression_data <- GetAssayData(pbmc_subset, slot = "data")[genes_to_plot, ]
average_expression <- colMeans(expression_data)
names(average_expression) <- colnames(pbmc_subset)
pbmc_subset$average_expression <- average_expression
feature_plot <- FeaturePlot(pbmc_subset, features = "average_expression") +
  scale_color_gradientn(colors = custom_colors)
print(feature_plot)
brewer.pal(n = 7, name = "Oranges")