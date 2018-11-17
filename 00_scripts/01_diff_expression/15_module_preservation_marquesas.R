
ervation module Tuot I####
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
setwd("~/Documents/Projets/Gamma/")
directory="~/Documents/Projets/Gamma/"
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

### load existing data on symbiodinium culture
#### load symbiont clam

# load symbiodinium culture
lnames = load(file = "03_results/dataInput_wgcna_gambier.Rda");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "03_results/network_gambier.Rda");
lnames
nGenes = ncol(datExprGambier)
nSamples = nrow(datExprGambier)
# load SFT  ### not made in the proper package
lnames=load(file="03_results/sft_signed_Gambier.Rda")
lnames


#### load symbiont culture
# load symbiodinium culture
lnames = load(file = "03_results/dataInput_wgcna_marquesas.Rda");
#The variable lnames contains the names of loaded variables.
lnames
nGenes = ncol(datExprMarquesas)
nSamples = nrow(datExprMarquesas)


# Calculation module preservation
setLabels = c("Gambier","Marquesas");
multiExpr = list(Gambier = list(data = datExprGambier),Marquesas = list(data = datExprMarquesas)); #remove best connectivity and kept all
multiColor = list(Gambier = moduleColors);

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "03_results/modulePreservation_marquesas.RData");

load("03_results/modulePreservation_marquesas.RData")
#pdf(file="03_results/figures/wgcna/Zsummary_clam_culture.pdf")
# analyse an display module preservation
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# TEst
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

## plot the results
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold","grey60","royalblue","black","purple"));
plotMods = (modColors %in% c("salmon","brown","turquoise","pink","lightgreen","greenyellow","magenta","cyan","tan","blue","green"));
# Text labels for points
text = modColors[plotMods];
text
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
plotData

# Main titles for the plot
mains =  "Preservation Zsummary in Marquesas"
# Start the plot
pdf("03_results/module_preservation_marquesas.pdf",width=6,height=6)
par(mar=c(5,5,3,2))
zsumplot<-plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
               cex = 3,
               ylab = "Preservation Zsummary in Marquesas" , xlab = "Module size", log = "x",
               ylim = c(-2,20),
               font.lab=2,
               xlim = c(100, 2000), cex.lab = 1.8, cex.axis = 1.8, cex.main =2)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1.4, offs = 0.1) 
abline(h=0) +
  abline(h=2, col = "red", lty = 2) +
  abline(h=10, col = "darkgreen", lty = 2) 
dev.off()

