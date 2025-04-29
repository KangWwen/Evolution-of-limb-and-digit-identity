library(reshape2)
GOterms <- read.csv("02.Merge_filtered_pvalue.csv")
GOterm_melt <- melt(GOterms[,1:8],id = c("Description","X.GeneInGOAndHitList"))
GOterm_melt$Description <- factor(GOterms$Description, 
                                  levels=unique(rev(GOterms$Description)))

GOterm_melt$value[GOterm_melt$value == 0] <- NA 

ggplot(GOterm_melt, aes(x =variable, y = Description, 
                        size = X.GeneInGOAndHitList, color = -value)) +
  geom_point() + scale_color_gradient(low = "#FFCE7E", 
                                      high = "#BD0026", na.value = NA)+
  scale_size(name = "Size", range = c(3, 7))

