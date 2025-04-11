library(pheatmap)
library(org.Rn.eg.db)
library(extrafont)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggplot2)
library(tidyr)
# library(openxlsx)
# library(expss)
# library(xlsx)
# library(do)
library(XLConnect)

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)

# font_import()
# loadfonts(device="all") 

# load in files ----

load('/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/gprofiler_all.RData')

load("lcpm.RData")
lcpm <- lcpm[,-c(13:18)]
colnames(lcpm) <- c("DMSO1", "DMSO2", "DMSO3", 'DMSO4', 'DMSO5', 'DMSO6', 'CAB1', 'CAB2', 'CAB3', 'CAB4', 'CAB5', 'CAB6')

sig_de <- read.csv("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/sig_de.csv")
# annots <- select(org.Rn.eg.db, keys=rownames(lcpm), columns=c("ENTREZID", "ENSEMBL", "SYMBOL"), keytype="ENSEMBL")
# sig_de <- sig_de %>% drop_na(Symbol)

wb <- loadWorkbook("heatmap_DE_info.xlsx")

# make heatmap ----

go_id_numbers <- '0030872'
go_id <- paste0('HP:', go_id_numbers)
go_num <- paste0('HP', go_id_numbers)
go_name <- 'AbnormalCardiacVentFunction'
top20 <- FALSE
filename <- paste0("heatmap_", go_num, "_", go_name, ".tiff")
sheetname <- str_sub(paste0("heatmap_", go_num, "_", go_name), 1, 30)

# str_sub(string1, 1, 3)  # prints "Pro"





## get go term gprofiler info ----
df <- gprofiler_results_oe$result[gprofiler_results_oe$result$term_id == go_id,]

go_term_genes <- unlist(str_split(df$intersection, ","), use.names = FALSE)

## match symbol names ----
m <- match(go_term_genes, sig_de$Symbol)

sig_de_spec_genes <- sig_de$X[m]

sig_de_spec_genes_df <- sig_de[sig_de$X %in% sig_de_spec_genes,]
createSheet(wb, name = sheetname)
appendWorksheet(wb, sig_de_spec_genes_df, sheet = sheetname)
saveWorkbook(wb)


# # get rid of NAs
# new_genes_list <- list()
# for (i in 1:length(sig_de_spec_genes)) {
#   gene <- sig_de_spec_genes[i]
#   if (gene %in% rownames(lcpm)) {
#     new_genes_list[[length(new_genes_list)+1]] = gene
#   }
# }
# 
# new_genes_list <- unlist(new_genes_list, use.names=FALSE)


# for (i in 1:length(sig_de_spec_genes)) {
#   gene <- sig_de_spec_genes[i]
#   if (gene %in% new_genes_list) {
#     next} else {
#       print(gene)
#     }
# }

ens_genes_list <- sig_de_spec_genes_df$X
sym_genes_list <- sig_de_spec_genes_df$Symbol

lcpm_spec_genes <- lcpm[ens_genes_list,]





rownames(lcpm_spec_genes) <- sym_genes_list

# row.names(lcpm_top20)[20] <- "Aspn"

rownames(lcpm_spec_genes) <- toupper(rownames(lcpm_spec_genes))

ifelse(top20 == TRUE, 
       (heatmap <- lcpm_spec_genes[c(1:20),]), 
       (heatmap <- lcpm_spec_genes)
       ) 
  

# top20 <- lcpm_spec_genes[c(1:20),]



# pushViewport(viewport(gp = gpar(fontfamily = "serif")))
p <- pheatmap(heatmap, name= "lcpm", cluster_cols = FALSE, color = hcl.colors(50, 'Blue-Red 2'), fontfamily="Times") #, rect_gp = gpar(col = "lightgray", lwd =1)) #, fontfamily="Times") #, filename = filename)
# draw(p, newpage = FALSE)
# ggsave(filename, plot=p$gtable)
# popViewport()
p
print(filename)

# setwd('/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/heatmaps')
# ggsave(plot=p,
#        filename = "irma_fatalities.pdf",
#        device = "pdf",
#        height = 6, width = 5, units = "in")
# tiff(filename, width=1000, height = 918, units = "px")
# print(p)
# dev.off()
# setwd(workdir)



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



