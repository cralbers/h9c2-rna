library(pheatmap)
library(org.Rn.eg.db)
library(extrafont)
library(ComplexHeatmap)
library(circlize)

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)

font_import()
loadfonts(device="all") 


load("lcpm.RData")
# lcpm

lcpm <- lcpm[,-c(13:18)]


# annots <- select(org.Rn.eg.db, keys=rownames(lcpm), columns=c("ENTREZID", "ENSEMBL", "SYMBOL"), keytype="ENSEMBL")
# m <- match(rownames(lcpm), annots$ENSEMBL)
# lcpm$EntrezGene <- annots$ENTREZID[m]

colnames(lcpm) <- c("DMSO1", "DMSO2", "DMSO3", 'DMSO4', 'DMSO5', 'DMSO6', 'CAB1', 'CAB2', 'CAB3', 'CAB4', 'CAB5', 'CAB6')

load('/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/gprofiler_all.RData')

lcpm_top20 <- rownames(cabtop20)

lcpm_top20 <- lcpm[c(cabtop20_ens),]





symbol_names <- mapIds(org.Rn.eg.db, rownames(lcpm_top20), keytype="ENSEMBL", column="SYMBOL")

rownames(lcpm_top20) <- symbol_names

row.names(lcpm_top20)[20] <- "Aspn"

rownames(lcpm_top20) <- toupper(rownames(lcpm_top20))


pheatmap(lcpm_top20, name= "lcpm", cluster_cols = FALSE, color = hcl.colors(50, 'Blue-Red 2'), fontfamily="Times") 


# exp_heatmap <- ggplot(data = lcpm_top20, mapping = aes(x = colnames(lcpm_top20),
#                                                      y = rownames(lcpm_top20),
#                                                      fill = log.expression)) +
#   geom_tile() +
#   xlab(label = "Subject") + # Add a nicer x-axis title
#   theme(axis.title.y = element_blank(), # Remove the y-axis title
#         axis.text.x = element_text(angle = 45, vjust = 0.5)) # Rotate the x-axis labels
# 
# exp_heatmap

col_fun = colorRamp2(seq(min(lcpm_top20), max(lcpm_top20), length = 10), hcl_palette = "Blue-Red 2")

pushViewport(viewport(gp = gpar(fontfamily = "serif")))
hm <- Heatmap(lcpm_top20, name= "lcpm", col = col_fun, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd =1))

draw(hm, newpage = FALSE)
popViewport()



