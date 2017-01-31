library (gplots)
library(Heatplus) 
library(vegan)
library(RColorBrewer)
library(gdata)

data.prop  <- read.csv("/home/astride/MEGA/Microbiota_CX/R/silvaGG_otu_table_forR .csv")
row.names(data.prop) <- data.prop$OTU.ID
data.prop <- data.prop[, -1]
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
#mi primer heatmap
heatmap(as.matrix(data.prop), Rowv = NA, Colv = NA, col = scaleyellowred)

# determine the maximum relative abundance for each column
maxab <- apply(data.prop, 2, max)
head(maxab)
Acetivibrio     Acetobacter   Achromobacter Acidaminococcus
1.088e-04       5.972e-04       5.441e-04       8.578e-05

# remove the genera with less than 1% as their maximum relative abundance
n1 <- names(which(maxab < 0.01))
data.prop.1 <- data.prop[, -which(names(data.prop) %in% n1)]

# the margins command sets the width of the white space around the plot. The first element is the bottom margin and the second is the right margin
heatmap(as.matrix(data.prop.1), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))
# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
data.dist <- vegdist(data.prop, method = "bray")
# Do average linkage hierarchical clustering. Other options are 'complete' or 'single'. You'll need to choose the one that best fits the needs of your situation and your data.
row.clus <- hclust(data.dist, "aver")
rows.clus<- hclust(data.dist, "single")
# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = NA, col = scaleyellowred, margins = c(10, 3))
# you have to transpose the dataset to get the genera as rows
data.dist.g <- vegdist(t(data.prop.1), method = "bray")
col.clus <- hclust(data.dist.g, "aver")

# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(10, 3))

#diagnóstico
var1<- read.csv("/home/astride/MEGA/Microbiota_CX/R/dx.csv")
row.names(var1) <- var1$OTU.ID
var1 <- var1[, -1]
var1 <- replace(var1, which(var1 == 0), "red")
var1 <- replace(var1, which(var1 == 1), "blue")
var1 <- replace(var1, which(var1 == 2), "orange")
var1 <- replace(var1, which(var1 == 3), "green")

var2<- read.csv("/home/astride/MEGA/Microbiota_CX/R/dx2.csv")
row.names(var2) <- var2$OTU.ID
var2 <- var2[, -1]
var2 <- replace(var2, which(var2 == 0), "blue")
var2 <- replace(var2, which(var2 == 1), "orange")

# Heatmap con barrita =)
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var2, margins = c(10, 3))
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var2, margins = c(11, 6), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", main = "Cervical Microbiome Cancer", lhei = c(2, 8))

# make a data frame of variables for annotation
ann.dat <- data.frame(var1 = c(rep("cat1", 4), rep("cat2", 8)), var2 = rnorm(12,  mean = 50, sd = 20))

# the annHeatmap2 function needs to be wrapped in the plot function in order to display the results
plot(annHeatmap2(as.matrix(data.prop.1),
col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(43), breaks = 42,
dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))), legend = 3,
labels = list(Col = list(nrow =29))))
# make a data frame of variables for annotation
ann.dat <- data.frame(var1 = c(rep("cat1", 4), rep("cat2", 8)), var2 = rnorm(12, mean = 50, sd = 20))

heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, 
          RowSideColors = var1, margins = c(11, 6), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", 
          main = "Cervical Microbiome Cancer", lhei = c(2, 8), cex.axis=4)



plot(annHeatmap2(as.matrix(data.prop.1),
col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(51),
breaks = 50,
dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))),
legend = 3,
labels = list(Col = list(nrow = 29)),
ann = list(Row = list(data = ann.dat)),
cluster = list(Row = list(cuth = 0.25, col = brewer.pal(3, "Set2")))) 
# cuth gives the height at which the dedrogram should be cut to form clusters, and col specifies the colours for the clusters
                 ))

# PCA
pca1 = prcomp(data.prop.1, scale. = TRUE)
head(pca1$x)
JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

data.dist=dist.JSD(data)
pam(as.dist(data.dist.g), 6, diss=TRUE) # x is a distance matrix and k the number of clusters
nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist.g, k)
    nclusters[k]=index.G1(t(data.prop.1),data.cluster_temp,  d = data.dist.g,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

require(clusterSim)
nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids") 