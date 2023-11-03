library(FactoMineR)
library(factoextra)
library(corrplot)
library(cluster)
library(pheatmap)
library(stats)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(dendextend)
library(zoom)
library(ape)

Corals<-read.csv('C:/Users/Adam/Desktop/Fall 2022/EVR 715 Machine Learning and Predictive Analytics/Project/Data/Strain/Coral_Strain_no_Ch.csv',header=TRUE,sep =",")

Corals<-Corals[, colSums(Corals != 0) > 0]

coral<-Corals[,-1]
rownames(coral)<-coral[,1] 
coral<-coral[,-1]
coral<-coral[,-1]

c_dist<-dist(coral[,1:4173])

matrix<-as.matrix(c_dist)
write.csv(matrix,'C:/Users/Adam/Desktop/Fall 2022/EVR 715 Machine Learning and Predictive Analytics/Project/Data/Strain/matrix.csv', row.names = TRUE)

#Single
SLC<-hclust(c_dist,method="single") 
SLC
plot(SLC) 

#Average
ALC<-hclust(c_dist,method="average")
ALC
plot(ALC)

#Complete
CLC<-hclust(c_dist,method="complete")
CLC
plot(CLC)

plot(SLC,main="Single Linkage Cluster")
plot(ALC,main="Average Linkage Cluster")
plot(CLC,main="Complete Linkage Cluster")

#Deciding Where to Cut the Trees

#Cluster Analysis
labels<-paste(coral$Species,sep="")
coralTib<-as_tibble(Corals)

#DEFINE CLUSTERING METRIC FUNCTION 
cluster_metrics <- function(data, clusters, dist_matrix) {list(db = clusterSim::index.DB(data, clusters)$DB, G1 = clusterSim::index.G1(data, clusters), dunn = clValid::dunn(dist_matrix, clusters), clusters = length(unique(clusters)))} 

#GENERATE BOOTSTRAP SAMPLES OF DATA
coralBoot <- map(1:100, ~ {coral %>% as_tibble() %>% sample_n(size = nrow(.), replace = TRUE)})


#CALCULATE PERFORMANCE METRICS

metricsTib <- map_df(coralBoot, function(boot) {
  d <- dist(boot, method = "euclidean")
  cl <- hclust(d, method = "ward.D2")
  
  map_df(3:8, function(k) {
    cut <- cutree(cl, k = k)
    cluster_metrics(boot, clusters = cut, dist_matrix = d)
  })
})
metricsTib

# MUTATE A BOOTSTRAP COLUMN AND GATHER ----
metricsTib <- metricsTib %>%
  mutate(bootstrap = factor(rep(1:100, each = 6))) %>%
  gather(key = "Metric", value = "Value", -clusters, -bootstrap)

# PLOT THE BOOTSTRAP RESULTS ----
ggplot(metricsTib, aes(as.factor(clusters), Value)) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_line(size = 0.1, aes(group = bootstrap)) +
  geom_line(stat = "summary", fun = "mean", aes(group = 1)) +
  stat_summary(fun.data="mean_cl_boot", 
               geom="crossbar", width = 0.5, fill = "white") +
  theme_bw()

#ggsave("hclustBootstrap.pdf", width = 10, height = 5)
# CUTTING THE DENDROGRAM TO A FIXED NUMBER OF CLUSTERS ----
SingleCut <- cutree(SLC, k = 5)
AverageCut <- cutree(ALC, k=5)
CompleteCut <- cutree(CLC, k=5)
plot(SingleCut)
plot(AverageCut)
plot(CompleteCut)

par(mfrow=c(1,1))

groupcodes = c(rep("Ctenactic", 8), rep("Herpolithia", 7), rep("Seriatopora", 1), rep("Acropera", 1), rep("Montipora", 1), rep("PoritesC", 1), rep("Monastraea", 5), rep("PoritesA", 6), rep("AstrangiaB", 36), rep("AnstrangiaW", 36), rep("Eunicella", 12), rep("Stylophora", 5))
man.codes = c(Ctenactis = "blue", Herpolithia = "green", Seriatopora = "aquamarine", Acropera = "chartreuse", Montipora = "darkgoldenrod", PoritesC = "darkblue", Monastraea = "deeppink", PoritesA = "purple", AstrangiaB = "brown", AnstrangiaW = "red", Eunicella = "black", Stylophora = "pink")

# Dendrogram Construction

SLCDend<-as.dendrogram(SLC)
plot(SLC,main="Single Linkage")
labels_colors(SLCDend) = man.codes[groupcodes][order.dendrogram(SLCDend)]
plot(SLCDend,main="Single Linkage Dendrogram")
rect.hclust(SLC, k = 5)

ALCDend<-as.dendrogram(ALC)
plot(ALC,main="Average Linkage")
labels_colors(ALCDend) = man.codes[groupcodes][order.dendrogram(ALCDend)]
plot(ALCDend,main="Complete Linkage Dendrogram")
rect.hclust(ALC, k = 5)

CLCDend<-as.dendrogram(CLC)
plot(CLC,main="Complete Linkage")
labels_colors(CLCDend) = man.codes[groupcodes][order.dendrogram(CLCDend)]
plot(CLCDend,main="Complete Linkage Dendrogram")
rect.hclust(CLC, k = 4)

colors = c("red", "blue", "green", "black", "purple")
clus = cutree(CLCDend,3)
plot(as.phylo(CLCDend), type = "fan", tip.color = colors[clus], cex = 0.6, label.offset = 0.5)

pheatmap(t(coral), clustering_method = "complete", cutree_cols=3)

zm()

cov.coral<-round(cov(coral),4)
cov.coral
eigen(cov.coral)

pca<-PCA(coral,graph=FALSE) 
summary(pca) 
names(pca)
pca$eig

pc.out<-prcomp(coral, scale=FALSE) 
names(pc.out) 
summary(pc.out)
pc.out$center
pc.out$scale
pc.out$rotation

a1<-pc.out$rotation[,1]

PC1rot<-as.matrix(round(a1,4))
write.csv(PC1rot,'C:/Users/Adam/Desktop/Fall 2022/EVR 715 Machine Learning and Predictive Analytics/Project/Data/Strain/PC1rot.csv', row.names = TRUE)

a2<-pc.out$rotation[,2]
head(round(a2,4))

PC2rot<-as.matrix(round(a2,4))
write.csv(PC2rot,'C:/Users/Adam/Desktop/Fall 2022/EVR 715 Machine Learning and Predictive Analytics/Project/Data/Strain/PC2rot.csv', row.names = TRUE)

a3<-pc.out$rotation[,3]
head(round(a3,4))

PC3rot<-as.matrix(round(a3,4))
write.csv(PC3rot,'C:/Users/Adam/Desktop/Fall 2022/EVR 715 Machine Learning and Predictive Analytics/Project/Data/Strain/PC3rot.csv', row.names = TRUE)


head(predict(pc.out))
plot(pc.out)
biplot(pc.out,main="PCA")+
  abline(h=0,lty=2,lwd=2,col="blue")+ 
  abline(v=0,lty=2,lwd=2,col="blue") 

apply(l.turtles,2,mean)
apply(l.turtles,2,var)

fviz_eig(pc.out)

cor.coral<-round(cor(coral),4)

s.coral<-scale(coral)

cov.s.turtles<-round(cov(s.l.turtles),4)

cov.s.turtles
cor.turtles
eigen(cor.turtles)

pcr.out<-prcomp(coral, scale=TRUE) 

names(pcr.out) 
summary(pcr.out)
pcr.out$center
pcr.out$scale
pcr.out$rotation
head(predict(pcr.out))
plot(pcr.out)
biplot(pcr.out,main="Scaled PCA")+
  abline(h=0,lty=2,lwd=2,col="blue")+ 
  abline(v=0,lty=2,lwd=2,col="blue") 

apply(l.turtles,2,mean)
apply(l.turtles,2,var)

fviz_eig(pcr.out)

summary(pca)
fviz_pca_biplot(pca,repel=TRUE)

summary(T.pca)
fviz_pca_biplot(T.pca,repel=TRUE)

options(ggrepel.max.overlaps = Inf)
fviz_pca_biplot(pc.out,repel=TRUE)