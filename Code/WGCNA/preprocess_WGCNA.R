#######################################################################
# This script is for preprocessing data for WGCNA 

setwd('~/Desktop/Projects/WQE/Data/Intensity/')
library(WGCNA)
# load data
rm(list = ls())
strains <- c('c57', 'balbc', 'dba', 'fvb', 'aj', 'cej')
treat <- c('ctrl', 'iso')
timepoints <- c(0,1,  3,  5,  7, 10, 14)
temp <- read.csv('c57_ctrl_imputed.csv',header = T)
dat <- data.frame(pro = temp[,1])
for(strain in strains){
  for(tre in treat){
    temp <- read.csv(paste0(strain,'_',tre,'_imputed.csv'),header = T)
    names(temp) <- sapply(c('pro', 'd0', 'd1', 'd3','d5','d7','d10','d14'),
                          function(x){paste(strain, tre, x, sep = '_')})
    temp <- temp[,-1]
    temp <- exp(temp)
    dat <- cbind(dat, temp)
  }}
rownames(dat) = dat$pro
dat = dat[,-1]
dat = t(dat)

# discard extreme values
dat[dat>0.0001] = NA
dat[is.na(dat)] = 0.0001

# check missing values
gsg = goodSamplesGenes(dat, verbose = 3);
gsg$allOK

# clustering to detect outlier
sampleTree = hclust(dist(dat), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# # Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.01, minSize = 10)
keepSamples = (clust==1)
dat = dat[keepSamples, ]

sampleTree2 = hclust(dist(dat), method = "average")
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

save(dat, file = 'Expression_data_outlier_removed.RData')
