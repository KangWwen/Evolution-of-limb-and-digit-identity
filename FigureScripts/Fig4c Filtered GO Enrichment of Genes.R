library(ggplot2)
GOterms <- read.csv("01.PickedGOs_Divided.csv")
GOterms <- GOterms[,c(1,3:14)]
GOterms_divied <- read.csv("./DividedGOs/Enrichment_GO/GO_AllLists.csv")
GOterms_divied_picked <- GOterms_divied[GOterms_divied$Description %in% GOterms$Description,]

GOterm_melt <- GOterms_divied_picked[,c(4,6,12,18)]
head(GOterm_melt)
colnames(GOterm_melt) <- c("Description","logP","Number","Clades")
GOterm_melt$Description <- factor(GOterm_melt$Description, 
                                  levels=unique(rev(GOterms$Description)))

GOterm_melt$Clades <- factor(GOterm_melt$Clades, 
                             levels=unique(colnames(GOterms)))

ggplot(GOterm_melt, aes(x =Clades, y = Description, size = Number, color = -logP)) +
  geom_point() + scale_color_gradient(low = "#FFD9D7", high = "#9F0C04", na.value = NA)+
  scale_size(name = "#Gene in GO", range = c(2, 12))

ggplot(GOterm_melt, aes(x =Clades, y = Description, size = Number, color = -logP)) +
  geom_point() + scale_color_gradient(low = "#D6DDE9", high = "#4575B4", na.value = NA)+
  scale_size(name = "#Gene in GO", range = c(2, 12))

