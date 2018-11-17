#!/usr/bin/Rscript

################ WGCNA RNA-seq Analysis project snapper ############
####################################################
ls()
rm(list=ls())
ls()
#source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#install.packages("WGCNA")

#install.packages("scales")
#install.packages("assertthat")
library("assertthat")
library("scales")
library("WGCNA")

#Global variables
setwd("01_projects/gamma/")
directory="01_projects/gamma/"

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the expr data set
logcpm_test<-read.table("log_matrix_gamma.txt", header=T,check.names = FALSE)
head(logcpm_test)


temp<-read.table("01_info_files/data_trait.txt",sep="\t", header=T,na.strings="NA")
vector<-temp$Ind

logcpm_test<-logcpm_test[, names(logcpm_test) %in% vector]

# Take a quick look at what is in the data set:
dim(logcpm_test)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
datExpr0 = as.data.frame(t(logcpm_test))

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#
#  Code chunk 6: check outliers
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 220, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 220, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)

# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0.2
table(KeepGenes)
datExpr=datExpr0[, KeepGenes]
name_datExpr <-colnames(datExpr)
#=====================================================================================
#
#  Code chunk 7: add traits
#
#=====================================================================================


allTraits = read.table("01_info_files/data_trait.txt",sep="\t", header=T,na.strings="NA");
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Ind);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
str(datTraits)

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits,signed= FALSE);
# Plot the sample dendrogram and the colors underneath.
#pdf("dendo_heatmap_subset.pdf",width=12,height=9)
par(mar=c(1, 10, 1, 1))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

dev.off()

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "wgcna/dataInput_gamma.Rda")



########################## Module construction step-by-step #################################
################################################################################
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "wgcna/dataInput_gamma.Rda");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 20, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")
#View(sft$fitIndices)

# Plot the results:
#load("sft_signed_cell_culture_clam.Rda")
#tiff(file = "wgcna_mean_connectivity_cell.tiff", width = 20, height = 20, units = "cm", res = 300)
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
 #    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
 #    main = paste("Scale independence"));
#text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
#abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
#plot(sft$fitIndices[,1], sft$fitIndices[,5],
 #    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  #   main = paste("Mean connectivity"))
#text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#dev.off()

save(sft,file="wgcna/sft_signed_gamma.Rda")


########################## Module construction step-by-step #################################
################################################################################
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "wgcna/dataInput_gamma.Rda");
#The variable lnames contains the names of loaded variables.
lnames
head(datExpr)

load("wgcna/dataInput_gamma.Rda")
load("wgcna/sft_signed_gamma.Rda")
message("data loaded")

softPower = 26; #85 %
adjacency = adjacency(datExpr, power = softPower,type="signed");
save(adjacency,file="wgcna/adjacency_gamma.Rda")

TOM = TOMsimilarity(adjacency,TOMType = "signed");
save(TOM,file="wgcna/TOM.Rda")
load("wgcna/TOM.Rda")
dissTOM = 1-TOM

#save(dissTOM,file="disTOM.Rda")
message("dissimilarity done")
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf("genetree_subset.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf("genetree_dyn_subset.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf("cluster_modules_subset.pdf",width=8,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
dev.off()

pdf(file = "geneDendro-3_subset.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna/network_gamma.Rda")
