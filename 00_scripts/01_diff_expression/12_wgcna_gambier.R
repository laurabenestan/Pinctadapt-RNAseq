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

########################## Module construction step-by-step #################################
################################################################################
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

#setting is important, do not omit.
# Load the data saved in the first part
lnames = load(file = "wgcna/dataInput_gambier.Rda");
#The variable lnames contains the names of loaded variables.
lnames
head(datExprGambier)

load("wgcna/dataInput_gambier.Rda")
load("wgcna/sft_signed_gambier.Rda")
message("data loaded")

softPower = 21; #90 %
adjacency = adjacency(datExprGambier, power = softPower,type="signed");
save(adjacency,file="wgcna/adjacency_gambier.Rda")

TOM = TOMsimilarity(adjacency,TOMType = "signed");
save(TOM,file="wgcna/TOM_gambier.Rda")
load("wgcna/TOM_gambier.Rda")
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
MEList = moduleEigengenes(datExprGambier, colors = dynamicColors)
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
merge = mergeCloseModules(datExprGambier, dynamicColors, cutHeight = MEDissThres, verbose = 3)
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
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna/network_gambier.Rda")
