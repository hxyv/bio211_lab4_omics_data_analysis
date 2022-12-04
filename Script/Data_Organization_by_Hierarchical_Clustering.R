install.packages('pheatmap')

cg1 <- rownames(edgeR_DEG)[edgeR_DEG$change != "NOT"]
cg2 <- rownames(limma_voom_DEG)[limma_voom_DEG$change != "NOT"]
library(pheatmap)
library(RColorBrewer)
color <- colorRampPalette(c("#436EEE", "white", "#EE0000"))(100)

# Hierarchical clustering on edgeR produced significant genes
mat1 <- a[cg1,]
n1 <- t(scale(t(mat1)))
n1[n1 > 1] <- 1
n1[n1 < -1] <- -1
ac <- data.frame(group = group_list)
rownames(ac) <- colnames(mat1)
ht1 <- pheatmap(n1, show_rownames = F, show_colnames = F,
                cluster_rows = F, cluster_cols = T,
                annotation_col = ac, color = color)
ggsave("Graph/ht_edgeR.pdf")
dev.off()

# Hierarchical clustering on limma-voom produced significant genes
mat2 <- a[cg2,]
n2 <- t(scale(t(mat2)))
n2[n2 > 1] <- 1
n2[n2 < -1] <- -1
ac <- data.frame(group = group_list)
rownames(ac) <- colnames(mat2)
ht2 <- pheatmap(n2, show_rownames = F, show_colnames = F,
                cluster_rows = F, cluster_cols = T,
                annotation_col = ac, color = color)
ggsave("Graph/ht_limma.pdf")
dev.off()

# Hierarchical clustering on edgeR & limma-voom overlapped significant genes
UP <- function(df){
    rownames(df)[df$change == "UP"]
}
DOWN <- function(df){
    rownames(df)[df$change == "DOWN"]
}
up <- intersect(UP(edgeR_DEG), UP(limma_voom_DEG))
down <- intersect(DOWN(edgeR_DEG), DOWN(limma_voom_DEG))
mat_total <- a[c(up, down),]
n4 <- t(scale(t(mat_total)))
n4[n4 > 1] <- 1
n4[n4 < -1] <- -1
ac <- data.frame(group = group_list)
rownames(ac) <- colnames(mat_total)
ht_combine <- pheatmap(n4, show_rownames = F, show_colnames = F,
                       cluster_rows = F, cluster_cols = T,
                       annotation_col = ac, color = color)
ggsave("Graph/ht_combine.pdf")
dev.off()
