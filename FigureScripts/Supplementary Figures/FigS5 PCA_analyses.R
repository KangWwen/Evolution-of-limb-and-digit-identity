Chicken=read.table("06.Chicken_TPM.txt",header = TRUE,sep = "\t",row.names = 1)
Si=read.table("06.Corpor_TPM.txt",header = TRUE,sep = "\t",row.names = 1)
Turtle=read.table("06.Turtle_TPM.txt",header = TRUE,sep = "\t",row.names = 1)

library(ggplot2)

all <- cbind(Chicken[,1:7],Turtle[,1:10],Si[,1:10])
allZscore <- as.data.frame(scale(all, center=F, scale=TRUE))
pca <- prcomp(t(allZscore),scale=F)
summary(pca)
plot(pca$x[,1:2])
group <- factor(c("C2","C3","C4","C1","C2","C3","C4",
                  "C1","C2","C3","C4","C5","C1","C2","C3","C4","C5",
                  "C1","C2","C3","C4","C5","C1","C2","C3","C4","C5"))
group2 <- data.frame(group)
group3 <- factor(c(rep("CF",3),rep("CH",4),
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
