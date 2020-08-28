################################################################
# This script is for module analysis, including connectivity analysis
# p-value distribution, p-value connectivity association, network 
# visualization


library(WGCNA)
library(ggplot2)

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

# module mean p-value
modules <- read.csv('./results/modules.csv')
p_val = read.csv('./results/p_value_wilcox.csv')
p_val$module = modules$module[match(p_val$protein,modules$protein)]
p_mean = aggregate(p_val$p, list(p_val$module), mean)
ggplot(data = p_mean, aes(x=Group.1, y = -log10(x), fill = Group.1)) + 
  geom_bar(stat= 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(x = 'module', y = 'mean -log10(p-value)')

ggplot(data = p_val, aes(x = module, y = -log10(p), fill = module)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(x = 'module', y = '-log10(p-value)')


# p-val dist
ggplot(data = p_val[p_val$module=='yellow',], aes(x=-log10(p))) + 
  geom_histogram() + ggtitle('module yellow')
ggplot(data = p_val[p_val$module=='salmon',], aes(x=-log10(p))) + 
  geom_histogram() + ggtitle('module salmon')
ggplot(data = p_val[p_val$module=='grey',], aes(x=-log10(p))) + 
  geom_histogram() + ggtitle('module grey')

# calculate connectivity
Tom = TOMsimilarityFromExpr(dat, power=3)
diag(Tom) = 0
p_val$intraC = 0
p_val$wholeC = 0
p_val$interC = 0

for(m in names(table(p_val$module))){
  indexes = which(p_val$module==m)
  for(i in indexes){
    p_val$intraC[i] = sum(Tom[i,indexes])/(length(indexes)-1)
    p_val$wholeC[i] = sum(Tom[i,])/(dim(Tom)[1]-1)
    p_val$interC[i] = sum(Tom[i,-indexes])/(dim(Tom)[1]-1)
  }
}
# p_val$interC = p_val$wholeC - p_val$intraC
ggplot(data = p_val, aes(x=reorder(module, wholeC, mean), y=wholeC, color=module)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ggtitle('whole connectivity boxplot') + 
  labs(x = 'module')
ggplot(data = p_val, aes(x=reorder(module,intraC,mean), y=intraC, color=module)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ggtitle('intra connectivity boxplot') + 
  labs(x = 'module')
ggplot(data = p_val, aes(x=reorder(module,interC, mean), y=interC, color=module)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ggtitle('inter connectivity boxplot') + 
  labs(x = 'module')

p_val$traVter = p_val$intraC/p_val$interC
ggplot(data = p_val, aes(x=reorder(module,traVter, mean), y=traVter, color=module)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ggtitle('intra connectivity/ inter connectivity boxplot') + 
  labs(x = 'module')
# turquoise standout
# grey are highly not intra-connected 
# check connectivity association with p-value
ggplot(data=p_val, aes(x=wholeC, y=-log10(p), color=module)) +
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  facet_wrap(~ module, ncol = 4 )

library(data.table)
p_temp <- data.table(p_val)
p_temp[,list(intercept=coef(lm(p~wholeC))[1], coef=coef(lm(p~wholeC))[2]),by=module]
p_temp[,cor.test(-log10(p),wholeC, method='pearson')[3:4], by=module]
# salmon, lightgreen has sig positive association

ggplot(data=p_val, aes(x=intraC, y=-log10(p), color=module)) +
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  facet_wrap(~ module, ncol = 4 )
p_temp[,cor.test(-log10(p),intraC, method='pearson')[3:4], by=module]
# grey and salmonhas sig positive association

p_temp[,cor.test(-log10(p),traVter, method='pearson')[3:4], by=module]


p_temp[,list(intercept=coef(lm(p~intraC))[1], coef=coef(lm(p~intraC))[2]),by=module]
# no strong association shown between connectivity and p-value

write.csv(p_val, 'all_info.csv', row.names = F)
write.csv(p_val[p_val$module=='grey' | p_val$module=='pink' |p_val$module=='salmon',],
          'select_module.csv', row.names = F)

# get hub proteins, top5
a <- p_val[p_val$module=='grey' | p_val$module=='pink' |p_val$module=='salmon',]
a = a[order(a$module, a$intraC, decreasing = T),]


# see module distribution of sigificant proteins 
sig_pro = read.csv("./results/diff_pro_wilcox.csv")
dis = modules$module[match(sig_pro$protein, modules$protein)]
dis_dat = data.frame(module = names(table(dis)), count = table(dis))
names(dis_dat)[3] = 'count'
# normalized the count
module_tab <- table(modules$module)
dis_dat$count = dis_dat$count/module_tab[match(dis_dat$module, names(module_tab))]
ggplot(data = dis_dat, aes(x= module, y=count, fill = module)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  labs(y = 'normalized count')
  







p_tt = p_val[p_val$module=='yellow',]
p_tt = p_tt[order(p_tt$interC, decreasing = T),]
p_tt = p_tt[order(p_tt$intraC, decreasing = T),]
ggplot(data=p_tt, aes(x=intraC, y=-log10(p))) +
  geom_point() +
  stat_smooth(method = "lm", col = "black")


##################################################################
# network visualization
d = 1-Tom
nodes = colnames(dat)
thres = 0.5
ind = which(Tom>0.5, arr.ind = T)

link = data.frame(from = rep(0, length(ind)), to = 0, weight = 0)
for(i in 1:length(ind)){
  link[i,] = c(nodes[ind[i,1]], nodes[ind[i,2]], Tom[ind[i,1], ind[i,2]])
}
temp = apply(ind, 1, function(x){return(c(nodes[x[1]], nodes[x[2]], Tom[x[1], x[2]]))})







