
library(ggplot2)
library(gridExtra)
library(reshape2)

df <- read.csv("./heatmaps/05.PickStage_forHeatmap_allLimbGenes_Orders.csv",row.names = 1)
df_forelimb <- df[, grepl("^...F", names(df))]
df_hindlimb <- df[, grepl("^...H", names(df))]

df_numeric <- apply(df, 2, as.numeric)     
rownames(df_numeric) <- rownames(df)

df_numeric[is.na(df_numeric)] <- 0

############Anterior genes 
genes <- c("PAX9","ZIC3","ALX4")
plots <- list()

df_numeric_picked <- df_numeric[rownames(df_numeric)%in% genes,]
df_melt <- melt(df_numeric_picked)

data1 <- df_melt[grep("D1$", df_melt$Var2), ]$value
data2 <- df_melt[grep("D[2-5]$", df_melt$Var2), ]$value

dens1 <- density(data1, bw = 0.05)
dens2 <- density(data2, bw = 0.05)

df_dens <- data.frame(x = c(dens1$x, dens2$x),
                      y = c(dens1$y, dens2$y),
                      group = c(rep("Group 1", length(dens1$x)), rep("Group 2", length(dens2$x))))

p <- ggplot(df_dens, aes(x = x, y = y, fill = group)) +
  geom_density(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("#FFCC70", "#4575B4")) +
  theme_minimal() +
  labs(x = "", y = "", title = "anterior") +
  scale_x_continuous(limits = c(0, 1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none") 
p

############Posterior genes 
genes <- c("TBX2","TBX3","HAND2","HOXD10","HOXD11","HOXD12","SHH","PTCH2")
plots <- list()

df_numeric_picked <- df_numeric[rownames(df_numeric)%in% genes,]
df_melt <- melt(df_numeric_picked)

data1 <- df_melt[grep("D[123]$", df_melt$Var2), ]$value
data2 <- df_melt[grep("D[45]$", df_melt$Var2), ]$value

dens1 <- density(data1, bw = 0.05)
dens2 <- density(data2, bw = 0.05)

df_dens <- data.frame(x = c(dens1$x, dens2$x),
                      y = c(dens1$y, dens2$y),
                      group = c(rep("Group 1", length(dens1$x)), rep("Group 2", length(dens2$x))))

p <- ggplot(df_dens, aes(x = x, y = y, fill = group)) +
  geom_density(stat = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("#FFCC70", "#4575B4")) +
  theme_minimal() +
  labs(x = "", y = "", title = "Posterior") +
  scale_x_continuous(limits = c(0, 1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none") 
p