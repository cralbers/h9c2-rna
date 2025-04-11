# pathway analysis --------------------------------------------------------


library(clusterProfiler)
library(enrichplot)
library(pathview)
library(stringr)
library(topGO)
library(GOplot)

# SET THE DESIRED ORGANISM HERE
organism = "org.Rn.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
kegg_organism = "rno"

df = read.csv("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/et.CABvsDMSO_toptags.csv", header=TRUE)
len_path = "~/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/LENvsDMSO_pathways"
cab_path = "~/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/CABvsDMSO_pathways"


kegg_pathway <- function(df, folder_path, drug, id) {
  
  # we want the log2 fold change 
  original_gene_list <- df$logFC
  
  names(original_gene_list) <- df$X
  names(original_gene_list)
  
  # Convert gene IDs for gseKEGG function
  # We will lose some genes here because not all IDs will be converted
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Rn.eg.db)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 <- df[df$X %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # get ids of problematic entrezids that are duplicated or triplicated or quadrupilated hehe
  n_occur <- data.frame(table(df2$Y))
  
  df2_aggregated <- df2 %>%
    group_by(Y) %>%
    summarise(logFC = mean(logFC, na.rm = TRUE))
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2_aggregated$logFC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2_aggregated$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  if(anyDuplicated(names(kegg_gene_list))) {
    stop("There are duplicate ENTREZ IDs in the kegg_gene_list.")
  }
  
  # Produce the native KEGG plot (PNG)
  setwd(folder_path)
  dme <- pathview(gene.data=kegg_gene_list, 
                  pathway.id=id, 
                  species = kegg_organism, 
                  kegg.native = T, 
                  out.suffix=paste(drug, "vsDMSO_pathview_et", sep=''), 
                  low=list(gene="red"), 
                  high=list(gene="green"))
  
  setwd("~/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/")
  return(dme)
}

kegg_pathway(df, cab_path, 'CAB', '05414')



kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)












# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ens

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Rn.eg.db, 
             pAdjustMethod = "none")
plotGOgraph(gse_bp, firstSigNodes = 5)
