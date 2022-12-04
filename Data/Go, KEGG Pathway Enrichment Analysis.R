# KEGG pathway
pathway <- read.table("Data/KEGG_pathway.txt", header = F, sep = "\t",
                      stringsAsFactors = F, fill = T)
colnames(pathway) <- c("Term", "Database","ID", "Input.number", "Background.number",
                       "P-value", "Corrected.P.Value")
kegg_pathway = pathway[1:264,]

## Calculate richFactor
kegg_pathway$richFactor <- kegg_pathway$Input.number / kegg_pathway$Background.number
kegg_top <- kegg_pathway[1:20,]
library(ggplot2)
colnames(kegg_top)
p <- ggplot(kegg_top, aes(x = richFactor, y = Term))

## Draw a scatter plot
p + geom_point()

## Map the size of the point to Input.number and map the color to Corrected.P.vlaue
p + geom_point(aes(size = Input.number, color = Corrected.P.Value))
kegg <- p + 
    geom_point(aes(size = Input.number, color = Corrected.P.Value)) +
    scale_color_gradient(low = "red", high = "green") +
    labs(title = "Statistics of Pathway Enrichment", 
         x = "Rich Factor", y = "",
         color = "qvalue",
         size = "Gene_number")+
    theme_bw()
ggsave("Graph/kegg.pdf")

# Go
go_enrich <- read.table("Data/Go_enrich.txt", header = F, sep = "\t",
                        stringsAsFactors = F, fill = T, quote = "")
go_enrich_clear <- go_enrich[3:3096,]
colnames(go_enrich_clear) <- c("Term", "Database","ID", "Input.number", "Background.number",
                               "P-value", "Corrected.P.Value", "Input", "Hyperlink")
go_enrich_top <- go_enrich_clear[1:10,]
p <- ggplot(go_enrich_top, aes(x = -log10(Corrected.P.Value), y = Term)) +
    geom_bar(stat = "identity", width = 0.7, fill = "#8DA1CB") +
    theme_test() +
    theme(axis.text = element_text(face = "bold", color = "gray50")) +
    labs(title = "The Most Enriched GO Terms")
ggsave("Graph/go.pdf")
dev.off()
