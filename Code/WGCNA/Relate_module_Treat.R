################################################################
# This script is for relating module and traits


setwd('~/Desktop/Projects/WQE/Data/Intensity/')
# load data
rm(list = ls())
strains <- c('c57', 'balbc', 'dba', 'fvb', 'aj', 'cej')
treat <- c('ctrl', 'iso')
timepoints <- c(0,1,  3,  5,  7, 10, 14)


lnames = load(file = "Expression_data_outlier_removed.RData");
lnames
lnames = load(file = 'networkConstruction-auto.RData')
lnames


mouse_strain = unlist(lapply(strsplit(rownames(dat), '_'), `[[`, 1))
treatment = unlist(lapply(strsplit(rownames(dat), '_'), `[[`, 2))
days = unlist(lapply(strsplit(rownames(dat), '_'), `[[`, 3))

traits = data.frame( strain = mouse_strain, 
                    treat = treatment, days = days)
rownames(traits) = rownames(dat)

### correlate module with traits
# Define numbers of genes and samples
nPros = ncol(dat);
nSamples = nrow(dat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# check the data before anova
library(ggplot2)
library(dplyr)
library(easyGgplot2)
temp <- data.frame(group = traits$strain, value = MEs$MEgrey)
ggplot2.histogram(data = temp, xName = 'value', groupName = 'group', faceting = T, 
                  facetingVarNames = 'group')

# anova test 
for(i in 1:3){
  traits[,i] = as.factor(traits[,i])
}


temp = cbind(traits, MEs)
p_mat = matrix(0, nrow = 3, ncol = dim(temp)[2]-3)
for(i in 1:3){
  for(j in 4:dim(temp)[2]){
    p_mat[i,j-3] = summary(aov(temp[,j] ~ temp[,i]))[[1]][['Pr(>F)']][1]
  }
}

library(RColorBrewer)
my_color = colorRampPalette(brewer.pal(9, 'OrRd'))
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  as.character(p_mat)
par(mar = c(6, 8.5, 3, 3));
p_mat_col = -log10(p_mat)
p_mat_col[p_mat_col< -log10(0.05)] = 0
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(p_mat_col),
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = my_color(50),
               textMatrix = t(trunc(-log10(p_mat)*10^2)/10^2),#t(formatC(log10(p_mat), format = "e", digits = 0)),#t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(0,3),
               main = "Module-trait relationships\n-log10(p-value)")

# biserial correlation
library(dummies)
traits_d = dummy.data.frame(traits, sep = '.')
traits_d = traits_d[,c(1:6,8)]


temp = cbind(traits_d, MEs)
p_mat = matrix(0, nrow = 15, ncol = dim(temp)[2]-15)
c_mat = matrix(0, nrow = 15, ncol = dim(temp)[2]-15)
for(i in 1:15){
  for(j in 16:dim(temp)[2]){
    a = cor.test(temp[,j], temp[,i])
    p_mat[i,j-15] = a$p.value
    c_mat[i,j-15] = a$estimate
  }
}

library(RColorBrewer)
my_color = colorRampPalette(brewer.pal(9, 'OrRd'))
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  as.character(p_mat)
par(mar = c(6, 8.5, 3, 3));

textMatrix =  paste(formatC(c_mat, format = "e", digits = 2), "\n(",
                    formatC(p_mat, format = "e", digits = 1), ")", sep = "");
dim(textMatrix) = c(15,17)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(c_mat),
               xLabels = names(traits_d),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = t(textMatrix),#t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.27,
               zlim = c(-1,1),
               main = "Module-trait relationships")

# module yellow is considered significant
modules <- data.frame(module = moduleColors, protein = names(moduleLabels))
write.csv(modules, './results/modules.csv', row.names = F, quote =F)
# TODO: KW test


# Define variable weight containing the weight column of datTrait
iso = as.data.frame(traits_d$treat.iso);
names(iso) = "iso"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(dat, iso, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(iso), sep="");
names(GSPvalue) = paste("p.GS.", names(iso), sep="");


module = "yellow" #'grey'
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for iso-treated group",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'red')


