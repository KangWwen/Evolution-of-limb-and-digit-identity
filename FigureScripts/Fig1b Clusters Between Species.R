library(pvclust)
library(dendextend)
library(magrittr)

Chicken=read.table("04.Axm_TPM_AV_sqrt.txt",header = TRUE,sep = "\t",row.names = 1)
Si=read.table("04.Si_TPM_AV_sqrt.txt",header = TRUE,sep = "\t",row.names = 1)
Turtle=read.table("04.Turtle_TPM_AV_sqrt.txt",header = TRUE,sep = "\t",row.names = 1)

early <- cbind(Chicken[,1:4],Turtle[,1:5],Si[,1:5])
eZscore <- as.data.frame(scale(early, center=F, scale=TRUE))
hc <- pvclust(eZscore, method.hclust="average", method.dist="correlation", nboot=1000)

# with a dendrogram of pvrect
dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

late <- cbind(Chicken[,5:9],Turtle[,6:10],Si[,6:10])
lZscore <- as.data.frame(scale(late, center=F, scale=TRUE))
hc <- pvclust(lZscore, method.hclust="average", method.dist="correlation", nboot=1000)

dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

####################

early <- cbind(Chicken[,1:4],Turtle[,1:5])
eZscore <- as.data.frame(scale(early, center=F, scale=TRUE))
hc <- pvclust(eZscore, method.hclust="average", method.dist="correlation", nboot=1000)

# with a dendrogram of pvrect
dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

late <- cbind(Chicken[,5:9],Turtle[,6:10])
lZscore <- as.data.frame(scale(late, center=F, scale=TRUE))
hc <- pvclust(lZscore, method.hclust="average", method.dist="correlation", nboot=1000)

dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

#######################
early <- cbind(Chicken[,1:4],Si[,1:5])
eZscore <- as.data.frame(scale(early, center=F, scale=TRUE))
hc <- pvclust(eZscore, method.hclust="average", method.dist="correlation", nboot=1000)

# with a dendrogram of pvrect
dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

late <- cbind(Chicken[,5:9],Si[,6:10])
lZscore <- as.data.frame(scale(late, center=F, scale=TRUE))
hc <- pvclust(lZscore, method.hclust="average", method.dist="correlation", nboot=1000)

dend <- as.dendrogram(hc)
hc %>% as.dendrogram %>% 
  plot(main = "")
hc %>% text

library(ggplot2)

all <- cbind(Chicken[,1:9],Turtle[,1:10],Si[,1:10])
allZscore <- as.data.frame(scale(all, center=F, scale=TRUE))

pca <- prcomp(t(allZscore),scale=F)
summary(pca)
plot(pca$x[,1:2])
group <- factor(c("C1","C2","C3","C4","C1","C2","C3","C4","C5",
                  "C1","C2","C3","C4","C5","C1","C2","C3","C4","C5",
                  "C1","C2","C3","C4","C5","C1","C2","C3","C4","C5"))
group2 <- data.frame(group)
group3 <- factor(c(rep("CF",4),rep("CH",5),
                   rep("TF",5),rep("TH",5),
                   rep("TF",5),rep("TH",5)))
group4 <- data.frame(group3)
pca_result <- as.data.frame(pca$x)
pca_result <- cbind(pca_result,group2,group4)
pca_result$group3
p <- ggplot(pca_result)+geom_point(aes(x=pca_result[,1],y=pca_result[,2],
                                       color=group,shape=group3),size=5)
p <- p+theme(legend.title = element_blank())+labs(x="PCA1",y="PCA2")+ 
  scale_color_manual(values=c("#FECA73", "#FB7F04", "#57C9FB","#2C7DB4", "#00870E",
                              "#FECA73", "#FB7F04", "#57C9FB","#2C7DB4", "#00870E"))+
  scale_shape_manual(values=c(16,15,10,7))
p
