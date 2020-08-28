####################################################################
# This script is for differential expression analysis

### Differential expression analysis
# determine whether the difference is normal distribution
# data size small to perform normality test
c57.d = c57_ctrl.dat - c57_iso.dat
hist(c57.d[1,])

# Does not seem normally distributed, or data too less to tell
# wilcoxon test
ctrl = Reduce(cbind, list(c57_ctrl.dat, balbc_ctrl.dat, dba_ctrl.dat,
                          fvb_ctrl.dat, aj_ctrl.dat, cej_ctrl.dat))
iso = Reduce(cbind, list(c57_iso.dat, balbc_iso.dat, dba_iso.dat,
                          fvb_iso.dat, aj_iso.dat, cej_iso.dat))
dif = iso - ctrl

res = data.frame(protein = rownames(ctrl), p = 0, median = 0)

res$p <- apply(dif, 1, function(x){wilcox.test(x)$p.value})
res$median <- apply(dif, 1, median)
res$median <- log2(exp(res$median))



# volcano plot
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = res$protein,
                x = 'median',
                y = 'p',
                xlim = c(-2,2),
                ylim = c(0, 7),
                xlab = 'ctrl/iso expression fold change',
                FCcutoff = 1,
                transcriptPointSize = 2
                )


library(ggplot2)
library(ggrepel)
res$diff <- NA
res$diff[res$median>log2(1.5) & res$p<0.05] <- 'up'
res$diff[res$median< -log2(1.5) & res$p<0.05] <- 'down'
res$label <- NA
res$label[!is.na(res$diff)] <- res$protein[!is.na(res$diff)]

ggplot(data = res, aes(x = median, y = -log10(p), col = diff, label = label)) +
  geom_point() +
  theme_minimal() + 
  geom_text_repel() +
  labs(y = '-log10(P-value)', x = 'log2foldchange')
  
diff_pro <- res[res$p<0.05,]
diff_pro <- diff_pro[abs(diff_pro$median)>log2(1.5),]
diff_pro <- diff_pro[,c(-4,-5)]
# store the results
write.csv(diff_pro,'./results/diff_pro_wilcox.csv', quote = F, row.names = F)
write.csv(res[,1:3], './results/p_value_wilcox.csv', quote = F, row.names = F)
write.csv(diff_pro$protein, './results/diff_pro_names.txt',
          quote = F, row.names = F, col.names = F)
