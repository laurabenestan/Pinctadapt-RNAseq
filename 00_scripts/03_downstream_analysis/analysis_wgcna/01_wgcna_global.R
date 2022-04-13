#!/usr/bin/Rscript

################ WGCNA RNA-seq Analysis project snapper ############
####################################################
ls()
rm(list=ls())
ls()
#source("http://bioconductor.org/biocLite.R") 
#biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 


library("assertthat")
library("scales")
library("WGCNA")

#Global variables
setwd("analysis_wgcna/")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in the expr data set
logcpm_test<-read.table("vst_matrix_genome_110621.txt", header=T,check.names = FALSE)
head(logcpm_test)

temp<-read.table("../../01_info_files/data_trait.txt",sep="\t", header=T,na.strings="NA")
temp<- temp[ which(temp$Days != "0"),]
vector<-temp$Ind
vector
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
nSamples
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)


# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0.05 # TEST for 0.05
table(KeepGenes)
datExpr=datExpr[, KeepGenes]

name_datExpr <-colnames(datExpr)
#=====================================================================================
#
#  Code chunk 7: add traits
#
#=====================================================================================


allTraits = read.table("data_trait.txt",sep="\t", header=T,na.strings="NA");
names(allTraits)

#allTraits<- allTraits[ which(!allTraits$Ind =="62"),]
allTraits$Ind
summary(allTraits)
#allTraits$Genotype<-NULL
allTraits$Sex <-NULL
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


save(datExpr, datTraits, file = "dataInput_wgcna_global.Rda")



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
lnames = load(file = "dataInput_wgcna_global.Rda");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 15, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")
View(sft$fitIndices)
# Plot the results:
#load("sft_signed_cell_culture_clam.Rda")
#tiff(file = "wgcna_mean_connectivity_cell.tiff", width = 20, height = 20, units = "cm", res = 300)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


save(sft,file="sft_signed_global.Rda")







####### network  #######

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "dataInput_wgcna_global.Rda");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "network_global.Rda");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# load SFT  ### not made in the proper package
lnames=load(file="sft_signed_global.Rda")
lnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
head(moduleTraitCor)
head(moduleTraitPvalue)

sizeGrWindow(25,15)


# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");



dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textAdj = c(0.5, 0.5),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# reformat matrix
par(mar = c(1,1, 1 ,1));
library(corrplot)
moduleTraitPvalue[moduleTraitPvalue > 0.01]<-0
moduleTraitCor[moduleTraitCor > -0.32 & moduleTraitCor < 0.32]<-0
df_module<-as.data.frame(moduleTraitPvalue)
df_module<-as.data.frame(moduleTraitCor)
new_mat <- as.data.frame(moduleTraitCor[ rownames(moduleTraitCor) %in% c("MEcyan","MEdarkgrey","MElightgreen","MEmidnightblue","MEturquoise",
                                                                         "MEred","MEblue","MEbrown","MEtan","MEdarkgreen",
                                                                         "MEblack","MEdarkred","MEpink","MEgrey60","MEyellow"),])
row.names(new_mat)<-gsub("ME", "", row.names(new_mat))
row.names(new_mat)<-c("cyan (341)","darkgrey (85)","lightgreen (159)","midnightblue (291)","turquoise (9,863)",
                      "red (1,239)","blue (6,183)","brown (3,659)","tan (493)","darkgreen (105)",
                      "black (1,083)","darkred (113)","pink (770)","grey60 (181)","yellow (2,199)")
mat<-as.matrix(new_mat)
mat<-mat[,c(1,2,3)]

par(mar=c(2,5,6,5))

corrplothost<-corrplot(mat, method = "number",tl.col = "black",cl.pos = "n", col= colorRampPalette(c("darkblue","white", "darkred"))(100))

corrplothost


#### Get info
temperature= as.data.frame(datTraits$Temperature);
names(temperature) = "temperature"
str(temperature)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
str(geneModuleMembership)


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
str(MMPvalue)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, temperature, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(temperature), sep="");
names(GSPvalue) = paste("p.GS.", names(temperature), sep="");

## plot correlation
### significance by module

module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(file = "module_black.pdf", width = 6, height = 6)
plot<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in host", module, "module"),
                         ylab = "Gene significance for temperature",
                         main = NULL, cex.lab = 1.2, cex.axis = 1.2, col = module)
plot + text(0.3, 0.78, "abs(cor) = 0.82\np-value < 0.001",cex = 1.5,font=2)
dev.off()



# export geneinfo
annot = read.csv(file = "annotation_genome_pmarg.110621.csv",sep=";",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$Sp[probes2annot],
                       LocusLinkID = annot$ID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for temperature
modOrder = order(-abs(cor(MEs, temperature, use = "p")));

# get relevant modules
# Get the corresponding Locuis Link IDs
allLLIDs = annot$name[probes2annot];
# $ Choose interesting modules
intModules = c("cyan","darkgrey","lightgreen","midnightblue","turquoise",
               "red","blue","brown","tan","darkgreen","black","darkred",
               "pink","grey60","yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("module_global_", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote=FALSE)
}

# A
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.temperature));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneinfo_globaltemperature.csv")


###################################
 ## for genotype
###################################


#### Get info
genotype= as.data.frame(datTraits$Genotype);
names(genotype) = "genotype"
str(genotype)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
str(geneModuleMembership)


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
str(MMPvalue)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, genotype, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(genotype), sep="");
names(GSPvalue) = paste("p.GS.", names(genotype), sep="");


# export geneinfo
annot = read.csv(file = "annotation_genome_pmarg.110621.csv",sep=";",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$Sp[probes2annot],
                       LocusLinkID = annot$ID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for temperature
modOrder = order(-abs(cor(MEs, genotype, use = "p")));

# get relevant modules
# Get the corresponding Locuis Link IDs
allLLIDs = annot$name[probes2annot];
# $ Choose interesting modules
intModules = c("cyan","darkgrey","lightgreen","midnightblue","turquoise",
               "red","blue","brown","tan","darkgreen","black","darkred",
               "pink","grey60","yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("module_global_", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote=FALSE)
}

# A
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.genotype));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneinfo_globalgenotype.csv")



###################################
## for days
###################################


#### Get info
days= as.data.frame(datTraits$Days);
names(days) = "days"
str(days)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
str(geneModuleMembership)


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
str(MMPvalue)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, days, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(days), sep="");
names(GSPvalue) = paste("p.GS.", names(days), sep="");


# export geneinfo
annot = read.csv(file = "annotation_genome_pmarg.110621.csv",sep=";",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$Sp[probes2annot],
                       LocusLinkID = annot$ID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for temperature
modOrder = order(-abs(cor(MEs, days, use = "p")));

# get relevant modules
# Get the corresponding Locuis Link IDs
allLLIDs = annot$name[probes2annot];
# $ Choose interesting modules
intModules = c("cyan","darkgrey","lightgreen","midnightblue","turquoise",
               "red","blue","brown","tan","darkgreen","black","darkred",
               "pink","grey60","yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("module_global_", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote=FALSE)
}

# A
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.days));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneinfo_globaldays.csv")
