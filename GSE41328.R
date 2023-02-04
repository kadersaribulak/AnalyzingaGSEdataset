library(GEOquery) #calling library
data <- getGEO("GSE41328") #to download GSE41328 dataset
data <- getGEO("GSE41328",  AnnotGPL = TRUE)
data = data[[1]] #unzziped

numerical_data =  data@assayData$exprs #numerical data
#hist(numerical_data)
numerical_data = log2(numerical_data) #log transformation is done
#hist(numerical_data)
sample =data@phenoData@data # pData(data) sample data
probe = data@featureData@data #probes info ; fData(data)



normal = sample[c(1:5,11:15),] #samples of normal tissue

tumor = sample[c(6:10,16:20),] #samples of tumor tissue


normal_num = numerical_data[, c(1:5,11:15)] #gene expression levels of normal tissue

tumor_num = numerical_data[,c(6:10,16:20)] #gene expression levels of tumor tissue

#controlling the code gives the right result
which(row.names(normal) == colnames(normal_num )) 
which(row.names(tumor) == colnames(tumor_num ))


# applying significant test by using t.test()
p = NULL
for (i in 1:nrow(probe) ) {
  p[i] = t.test(tumor_num[i,], normal_num[i,])$p.value
}

#finding indices genes with cut off 0.05 to BH correction method
FDR = p.adjust(p,method="BH") 
FDRs =  which (FDR < 0.05) #index

length(FDRs) #7064

BH = probe$`Gene symbol`[FDRs] #gene names of significantly changed according to BH.

unique_BH = unique(BH) 

length(unique_BH) #4954 significantly changed genes of 7064 have unique gene names


# forming data frame with gene symbols and p values
df = data.frame(probe$`Gene symbol`, p)
#to find three most significant genes according to BH by examining p value.
#sort the data frame by increasing
sorting_data = df[order(df$p , decreasing = FALSE), ]  

#gene symbols of most significant three genes

sign = sorting_data[1:3,]$probe..Gene.symbol.

sorting_data[1:3,] #12704, 31829,36731

#controlling

probe[12704,3]
# "CDH3"
probe[31829,3]
# "CLDN1"
probe[36731,3]
# "FOXQ1"

# to apply pearson and spearman correlation for both condition, expression level of them are determined
#for tumor condition
CDH3 = numerical_data[12704, c(6:10,16:20)]

CLDN1 = numerical_data[31829, c(6:10,16:20)]

FOXQ1 = numerical_data[36731, c(6:10,16:20)]
#data frame for tumor to see correlation of them easily
df_tumor = data.frame(CDH3,CLDN1,FOXQ1)

#correlation of each two gene for pearson correlation
cor.test(CDH3, CLDN1, method = "pearson")
cor.test(CDH3, FOXQ1, method = "pearson")
cor.test(FOXQ1, CLDN1, method = "pearson")

cor(df_tumor, method = "pearson")


#correlation of each two gene for spearman correlation
cor.test(CDH3, CLDN1, method = "spearman")
cor.test(CDH3, FOXQ1, method = "spearman")
cor.test(FOXQ1, CLDN1, method = "spearman")

cor(df_tumor, method = "spearman")

#for healthy condition
CDH3_ = numerical_data[12704,  c(1:5,11:15)]

CLDN1_ = numerical_data[31829,  c(1:5,11:15)]

FOXQ1_ = numerical_data[36731,  c(1:5,11:15)]

#data frame for healthy to see correlation of them easily
df_normal = data.frame(CDH3_,CLDN1_,FOXQ1_)

#correlation of each two gene for pearson correlation
cor.test(CDH3_, CLDN1_, method = "pearson")
cor.test(CDH3_, FOXQ1_, method = "pearson")
cor.test(FOXQ1_, CLDN1_, method = "pearson")

#correlation from data frame
cor(df_normal, method = "pearson")

#correlation of each two gene for spearman correlation
cor.test(CDH3_, CLDN1_, method = "spearman")
cor.test(CDH3_, FOXQ1_, method = "spearman")
cor.test(FOXQ1_, CLDN1_, method = "spearman")
#correlation from data frame
cor(df_normal, method = "spearman")




# Clustering

Clustdata = numerical_data[FDRs,] #for samples
genes_clust = t(Clustdata) # for genes


dist_samples = as.dist(1-cor(Clustdata, method = "pearson")) #distance measure for pearson
Hist = hclust(dist_samples, method = "centroid") #centroid method for clustering

plot(as.dendrogram(Hist))  #clustering samples
  

dist_genes = as.dist(1-cor(genes_clust, method = "pearson"))
Hist_genes = hclust(dist_genes, method = "complete") 

plot(as.dendrogram(Hist_genes))  #clustering genes



fit <- cutree(Hist_genes, k = 3) #cutting three cluster with k parameter, desired number of groups
table(fit)

a = cutree(Hist_genes, h = 1.2)  #cutting three cluster with h parameter, height of clustering, about 1.2
table(a)

three = which(fit == 3) #indices of genes in third cluster
third_cluster_genes = BH[three] #1309 #gene symbols of them

#unique genes belonging to third cluster
unique_third= unique(third_cluster_genes) #1018

#to apply GO:BP enrichment anaysis
write.table(unique_third, file = "third_cluster.txt", row.names = FALSE) 


#heatmap to cluster both rows and colunms 
heatmap(Clustdata, Rowv = as.dendrogram(Hist_genes), Colv = as.dendrogram(Hist), 
        scale = "row")


DataPCA = t(numerical_data) #transpose of numerical data to analyze samples in PCA

PCA = prcomp(DataPCA, center =  T, scale. = T, retx = T) #PCA analysis

summary(PCA)

str(PCA) #information of output PCA 

screeplot(PCA) #scree plot

#plotting with advanced graphical parameter
#plotting PC1 (most variation), PC2 (second most variation)

plot(PCA$x[,1], PCA$x[,2],main="PC1 vs. PC2",
     xlab="PC1", ylab="PC2", pch=18, col="blue") 
text(PCA$x[,1], PCA$x[,2], row.names(DataPCA), cex=0.6, pos=3, col="red")




