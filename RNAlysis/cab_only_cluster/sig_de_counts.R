library(gprofiler2)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(docstring)

setwd("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_sig_cluster/")
cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_sig_cluster/cluster_csvs_5"
load("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/et.CABvsDMSO_toptags.RData")

# get only sig DE counts ----

counts <- read.csv("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/cpm_cab.csv")
rownames(counts) <- counts[,1]
counts[,1] <- NULL

sig_de <- read.csv("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/sig_de.csv")


## get list of ens genes from sig de ----

sig_de_genes <- sig_de$X

## get de info from toptags table for only genes in cluster ----

sig_de_counts <- counts[c(sig_de_genes),]

# write.csv(sig_de_counts, "sig_de_cpm_cab.csv", quote=FALSE, row.names = TRUE)

# gprofiler analysis ----

cluster_files <- list.files(cluster_path)

get_cluster_genes <- function(filename){
  
  cluster <- read.csv(file.path(cluster_path, filename))
  
  # get cluster number
  
  string <- strsplit(filename, "\\.")[[1]]
  number <- strsplit(string, "cluster")[[1]][2]
  
  genes <- cluster$X
  
  # write.table(genes, file = paste0("cluster", number, "genes.txt"), sep = "",
  #             row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  
  ## get list of ens genes from cluster ----
  
  cluster_genes <- cluster$X
  return(cluster_genes)
}

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
                               sources = paste0("GO:", source), #NULL, 
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

# cluster_files
# [1] "sig_de_cpm_cab_kmeanscluster1.csv"  
# [2] "sig_de_cpm_cab_kmeanscluster10.csv" 
# [3] "sig_de_cpm_cab_kmeanscluster2.csv" 
# [4] "sig_de_cpm_cab_kmeanscluster3.csv"  
# [5] "sig_de_cpm_cab_kmeanscluster4.csv"  
# [6] "sig_de_cpm_cab_kmeanscluster5.csv" 
# [7] "sig_de_cpm_cab_kmeanscluster6.csv"  
# [8] "sig_de_cpm_cab_kmeanscluster7.csv"  
# [9] "sig_de_cpm_cab_kmeanscluster8.csv" 
# [10] "sig_de_cpm_cab_kmeanscluster9.csv" 

gprofiler_analysis(cluster_files[2], "BP", "gprofiler_BP_10")

for( i in 1:length(cluster_files)) {
  gprofiler_analysis(cluster_files[i], "CC", "gprofiler_CC_5")
}




# volcano labeled by cluster ----

et.CABvsDMSO_toptags$cluster <- NA


for(i in 1:length(cluster_files)) {
  
  filename <- cluster_files[i]
  
  # import cluster csv
  
  # cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs"
  
  cluster <- read.csv(file.path(cluster_path, filename))
  
  # get cluster number
  
  string <- strsplit(filename, "\\.")[[1]]
  number <- strsplit(string, "cluster")[[1]][2]
  
  # get list of ens genes from cluster
  
  cluster_genes <- cluster$X
  
  print(paste0("cluster number: ", number))
  print(paste0("num of genes in cluster: ", length(cluster_genes)))
  
  
  # replace the row value with cluster # if the 
  # data element at col index  is divisible  
  # by 2 looping over the rows of data frame 
  
  for (i in 1:nrow(et.CABvsDMSO_toptags)){ 
    
    # iterate over the 2nd column only of the 
    # data frame and check if divisible by 2 
    if(rownames(et.CABvsDMSO_toptags)[i] %in% cluster_genes){ 
      
      # replace the value with cluster number 
      et.CABvsDMSO_toptags[i,10]<- number
      
    } 
  } 
}

et.CABvsDMSO_toptags$cluster <- factor(et.CABvsDMSO_toptags$cluster, levels = c(1:10, NA))


sum(is.na(et.CABvsDMSO_toptags$cluster))

colorRampPalette(c("blue", "red"))(5)


p <- ggplot(data=et.CABvsDMSO_toptags, aes(x=logFC, y=-log10(FDR), col=cluster)) + 
  geom_point() +
  theme_minimal() +
  # scale_color_manual(values = colors) +
  theme(text=element_text(size=16,  family="Times")) + 
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
  scale_color_brewer(palette="Set3")
# scale_fill_discrete(breaks=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', 'NA'))
p



ggplotly(p)


