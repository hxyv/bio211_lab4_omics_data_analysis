# [1] Data preparation
prot <- read.csv('Data/proteomics data analysis.csv')
dim(prot)
prot$logFC <- log2(prot$Fold.Change)
prot$logp <- -log10(prot$p.value)

# [2] Construct a volcano plot
library(cowplot)
library(patchwork)
library(ggplotify)
library(ggplot2)

up <- intersect(which(prot$p.value < 0.05), which(prot$Fold.Change >= 2))
down <- intersect(which(prot$p.value < 0.05), which(prot$Fold.Change <= 1 / 2))
significant <- rep('normal', times = nrow(prot))
significant[up] <- 'up'
significant[down] <- 'down'
significant <- factor(significant, levels=c('up', 'down', 'normal'))
p <- qplot(x = prot$logFC, y = prot$logp, xlab = 'log2(FC)', ylab = '-log10(p-value)',
      size = I(0.7), colour = significant)
p <- p + scale_color_manual(values = c('up' = 'red',
                                       'normal' = 'grey',
                                       'down' = 'blue'))

## Add cut-off lines & export image
xline <- c(-log2(2), log2(2))
p <- p + geom_vline(xintercept = xline, lty = 2, size = I(0.2), color = 'grey11')

yline <- -log(0.05, 10)
p <- p + geom_hline(yintercept = yline, lty = 2, size = I(0.2), color = 'grey11')

p <- p + theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          size = 1,
                                          fill = 'white'),
          panel.grid = element_blank())
pdf('Graph/proteomics.pdf')
print(p)
dev.off()
