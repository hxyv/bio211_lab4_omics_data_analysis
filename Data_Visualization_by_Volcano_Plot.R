# [1] Data preparation
library(stringr)
rownames(DEG2) <- str_sub(rownames(DEG2), start=1, end=15)
DEG2$ENSEMBL <- rownames(DEG2)
diff_gene <- DEG2[DEG2$change!='NOT',]
write.table(diff_gene,
            file = 'Data/LIHC_diff_gene.txt',
            sep = '\t', quote = F)
# Change the path accordingly

# [2] Construct a volcano plot
library(cowplot)
library(patchwork)
library(ggplotify)
library(ggplot2)

loc_up <- intersect(which(DEG2$FDR < 0.05), which(DEG2$logFC >= 1))
loc_down <- intersect(which(DEG2$FDR < 0.05), which(DEG2$logFC < (-1)))
significant <- rep('normal', times = nrow(DEG2))
significant[loc_up] <- 'up'
significant[loc_down] <- 'down'
significant <- factor(significant, levels = c('up', 'down', 'normal'))
p <- qplot(x = DEG2$logFC, y = -log10(DEG2$FDR), 
           xlab='log2(FC)', ylab = '-log10(FDR)',
           size = I(0.7), colour = significant)
p <- p + scale_color_manual(values = c('up' = 'red', 
                                       'normal' = 'grey',
                                         'down' = 'blue'))

# Add cut-off lines & export image
xline <- c(-log2(2), log2(2))
p <- p + geom_vline(xintercept = xline, lty = 2, size = I(0.2), color = 'grey11')

yline <- -log(0.05, 10)
p <- p + geom_hline(yintercept = yline, lty = 2, size = I(0.2), color = 'grey11')

p <- p + theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          size = 1,
                                          fill = 'white'),
                            panel.grid = element_blank())
pdf('Graph/deg2_volcano.pdf')
print(p)
dev.off()
