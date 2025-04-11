library(gprofiler2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(docstring)

setwd("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster")

cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs"
cluster_files <- list.files(cluster_path)

filename <- "h9c2_cab_filt35sum_kmeanscluster2.csv"

load("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/et.CABvsDMSO_toptags.RData")

all_de_genes <- row.names(et.CABvsDMSO_toptags)


# write.table(all_de_genes$input, file = "all_de_genes.txt", sep = " ",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)





all_de_genes <- gconvert(all_de_genes, organism = "rnorvegicus", target = "ENTREZGENE")

get_cluster_genes <- function(filename){
  
  cluster <- read.csv(file.path(cluster_path, filename))
  
  # get cluster number 
  
  string <- strsplit(filename, "\\.")[[1]]
  number <- strsplit(string, "cluster")[[1]][2]
  
  genes <- cluster$X
  write.table(genes, file = paste0("cluster", number, "genes.txt"), sep = "",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  ## get list of ens genes from cluster
  
  cluster_genes <- cluster$X
  return(cluster_genes)
}

cluster <- get_cluster_genes(filename)

cluster <- gconvert(cluster, 
                     organism = "rnorvegicus",
                     target = "ENTREZGENE"
                     )


gprofiler_results_oe <- gost(query = cluster$target, 
                                  organism = "rnorvegicus",
                                  ordered_query = F, 
                                  exclude_iea = F, 
                                  user_threshold = 0.05, 
                                  # max_set_size = 0,
                                  correction_method = "fdr",
                                  domain_scope = "annotated",
                                  # custom_bg = all_de_genes$target,
                                  sources = "HP",
                             evcodes = TRUE,
                             highlight = TRUE)

p <- gostplot(gprofiler_results_oe, capped = FALSE, interactive = TRUE)
p

# get only driver terms 
filtered_data <- gprofiler_results_oe_tbl %>%
  filter(if_any(everything(), ~ .x == TRUE))

gprofiler_results_oe_tbl <- (gprofiler_results_oe$result)

idx <- gprofiler_results_oe_tbl[is.element(gprofiler_results_oe_tbl$highlighted, TRUE),]

c(idx$term_id[c(1:5)])
publish_gosttable(idx, highlight_terms = idx$term_id[c(1:5)],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = paste0(number, "_top5table_", source, '.tiff'))

saveWidget(ggplotly(p), file = paste0(number, "_top5manhat_", source, '.html'))

# publish_gostplot(
#   p,
#   highlight_terms = (idx$term_name[(1:5)]),
#   filename = paste0(1, "_top5manhat_", 'source', '.tiff'),
#   width = NA,
#   height = NA)


# function to run gprofiler on all clusters ----------------------------------------------------------------



gprofiler_analysis <- function(filename, source, folder) {
  #' run gprofiler analysis on cluster
  #'
  #' @param filename string. cluster csv file name
  #' @param source string. one of GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MIRNA, CO-RUM, HP, 
  #' HPA, WP. Defines what data sources gprofiler should pull from
  #' @param folder string. folder name that will be appended to beginning of filenames
  #'
  #' @returns manhattan plot as html and tiff, table with top 5 enriched terms as tiff, gprofiler results saved to rdata object in wd
  #' 
  #'
  
  
  
  string <- strsplit(filename, "\\.")[[1]]
  number <- strsplit(string, "cluster")[[1]][2]
  
  print(paste0("cluster ", number))
  
  cluster <- get_cluster_genes(filename)
  
  cluster <- gconvert(cluster, 
                       organism = "rnorvegicus",
                       target = "ENTREZGENE")
  
  print("convert gene names complete")
  
  gprofiler_results_oe <- gost(query = cluster$target, 
                               organism = "rnorvegicus",
                               ordered_query = F, 
                               exclude_iea = F, 
                               user_threshold = 0.05, 
                               # max_set_size = 0,
                               correction_method = "fdr",
                               domain_scope = "annotated",
                               # custom_bg = all_de_genes$target,
                               sources = source,
                               evcodes = TRUE,
                               highlight = TRUE)
  
  print("gprofiler complete")
  
  p_inter <- gostplot(gprofiler_results_oe, capped = FALSE, interactive = TRUE)
  p_stat <- gostplot(gprofiler_results_oe, capped = FALSE, interactive = FALSE)
  
  saveWidget(ggplotly(p_inter), file = paste0(folder, "/", source, "_", number, "_manhat_inter", '.html'))
  print("saved interactive Manhattan")
  
  gprofiler_results_oe_tbl <- (gprofiler_results_oe$result)
  
  idx <- gprofiler_results_oe_tbl #[is.element(gprofiler_results_oe_tbl$highlighted, TRUE),]
  
  publish_gosttable(idx, highlight_terms = idx$term_id[c(1:5)],
                    use_colors = TRUE, 
                    show_columns = c("source", "term_name", "term_size", "intersection_size"),
                    filename = paste0(folder, "/", source, "_", number, "_top5table", '.tiff'))
  
  publish_gostplot(p_stat, filename = paste0(folder, "/", source, "_", number, "_manhat_static", '.tiff'))
  
  save(gprofiler_results_oe, file = paste0(folder, "/", source, "_", number, '.RData')) 
}

cluster_path2 <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs_HP"

cluster_files <- list.files(cluster_path2)


# gprofiler_analysis(cluster_files[2], "CC", "gprofiler_GO_CC")

for( i in 1:length(cluster_files)) {
  gprofiler_analysis(cluster_files[i], "HP", "gprofiler_HP")
}


# all sig DE gprofiler ----

upreg <- subset(et.CABvsDMSO_toptags,logFC > 1 & FDR < 0.05)

downreg <- subset(et.CABvsDMSO_toptags,logFC < (-1) & FDR < 0.05)

sig_de <- rbind(upreg, downreg)


# write.csv(sig_de, file = "sig_de.csv")

# write.table(rownames(sig_de), file = "sig_de_genes.txt", sep = " ",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)


sig_de_genes <- gconvert(rownames(sig_de), organism = "rnorvegicus", target = "ENSG")

sum(is.na(sig_de$RefSeq))


sig_source <- NULL

gprofiler_results_oe <- gost(query = rownames(sig_de), 
                             organism = "rnorvegicus",
                             ordered_query = F, #T, 
                             exclude_iea = F, 
                             user_threshold = 0.05, 
                             # max_set_size = 0,
                             correction_method = "fdr",
                             domain_scope = "annotated",
                             # custom_bg = all_de_genes$target,
                             sources = NULL,
                             evcodes = TRUE,
                             highlight = TRUE)

gprofiler_results_oe_tbl <- gprofiler_results_oe$result

# gprofiler_results_oe_tbl_save <- apply(gprofiler_results_oe$result, 2, as.character)
# write.table(gprofiler_results_oe_tbl_save, "gprofiler_all.tsv", sep='\t')
# 
# save(gprofiler_results_oe, file = "gprofiler_all_ordered.RData")





idx <- gprofiler_results_oe_tbl #[is.element(gprofiler_results_oe_tbl$highlighted, TRUE),]

publish_gosttable(idx, highlight_terms = idx$term_id[c(1:10)],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = paste0("BP", "_", "sig_de_genes", "_top10table", '.tiff'))


p_inter <- gostplot(gprofiler_results_oe, capped = FALSE, interactive = TRUE)
saveWidget(ggplotly(p_inter), file = paste0('all', "_", "sig_de", "_manhat_inter", '.html'))


p_stat <- gostplot(gprofiler_results_oe, capped = FALSE, interactive = FALSE)
publish_gostplot(p_stat, filename = paste0('all', "_", "sig_de", "_manhat_static", '.tiff'))

# ready for export for revigo ----

## Order the results by p-adjusted value ----
gprofiler_results_oe_reordered <- gprofiler_results_oe_tbl[order(gprofiler_results_oe_tbl$p_value), ]

write.table(subset(gprofiler_results_oe_reordered, select=c("term_id")), "sig_de_GO.csv", quote=FALSE, row.names = FALSE, sep = ";") #, "p_value", "intersection")), 

## Extract only GO IDs and p-values for downstream analysis ----

GOpval_oe <- gprofiler_results_oe_reordered[ , c("term_id", "p_value")]

# write.table(GOpval_oe, "GOs_oe.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

# generate cytoscape input ----
gem <- gprofiler_results_oe$result[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
head(gem, 3)
write.table(gem, file = "gProfiler_gem_GOMF.txt", sep = "\t", quote = F, row.names = F)

# testing out cluster profiler plotting ----

library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(org.Rn.eg.db)

dotplot(gprofiler_results_oe_reordered, x="intersection_size")

barplot(gprofiler_results_oe_reordered, )