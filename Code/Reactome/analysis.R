############################################################
# reactome finding associated pathways

setwd('~/Desktop/Projects/WQE/FunctionalAnnotation/Reactome/')

rel = read.table('ReactomePathwaysRelation.txt')
names(rel) = c('parent', 'child')
rel$spe = unlist(lapply(strsplit(rel$parent, '-'), `[[`, 2))
rel = rel[rel$spe=='MMU',]

# judge if given pathway are parent pathway of child pathway
findparent <- function(rel, child, parent){
  if(parent == child){return(1)}
  if(parent != child){
    d_parent = rel$parent[which(rel$child==child)]
    if(length(d_parent)==0){return(0)}
    for(d_p in d_parent){
      found = findparent(rel, d_p, parent)
      if(found==1){
        return(1)
      }
    }
    return(0)
  }
}
test_rel = data.frame(parent = c(1,1,1,2,2,3,4), child = 2:8)
findparent(test_rel, 5,1)


diff = read.csv('diff_1.5_mouse.csv')
diff = diff[diff$Entities.pValue<0.05,]
yellow = read.csv('yellow_mouse.csv')
yellow = yellow[yellow$Entities.pValue<0.05,]
case = read.csv('top_mouse.csv')
case = case[case$Entities.pValue<0.05,]

pathways = c('R-MMU-390522','R-MMU-70171','R-MMU-1445148','R-MMU-75109',
             'R-MMU-70614','R-MMU-70268','R-MMU-1650814','R-MMU-446652',
             'R-MMU-975634','R-MMU-202733','R-MMU-5673001','R-MMU-611005')

diff_dist = matrix(NA, nrow = dim(diff)[1], ncol = length(pathways))
for(i in 1:dim(diff)[1]){
  for(j in 1:length(pathways)){
    diff_dist[i,j] = findparent(rel, diff$Pathway.identifier[i], pathways[j]) | findparent(rel, pathways[j], diff$Pathway.identifier[i])
  }
}
table(diff_dist)

yellow_dist = matrix(NA, nrow = dim(yellow)[1], ncol = length(pathways))
for(i in 1:dim(yellow)[1]){
  for(j in 1:length(pathways)){
    yellow_dist[i,j] = findparent(rel, yellow$Pathway.identifier[i], pathways[j]) | findparent(rel, pathways[j], yellow$Pathway.identifier[i])
  }
}
table(yellow_dist)

case_dist = matrix(NA, nrow = dim(case)[1], ncol = length(pathways))
for(i in 1:dim(case)[1]){
  for(j in 1:length(pathways)){
    case_dist[i,j] = findparent(rel, case$Pathway.identifier[i], pathways[j]) | findparent(rel,  pathways[j], case$Pathway.identifier[i])
  }
}
table(case_dist)

library(VennDiagram)
library(RColorBrewer)

myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(diff$Pathway.identifier, yellow$Pathway.identifier, case$Pathway.identifier),
  category.names = c('differential', 'yellow', 'text mining'),
  filename = 'venndiagram.png',
  output = T,
  fill = myCol,
  scaled = T,
  euler.d = T, 
  sep.dist = 5,
  lwd = 2,
  lty = 'blank',
  cex = 1.6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.default.pos = "outer",
  rotation = 1
  #alpha = .8
  )
      
