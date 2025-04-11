library(limma)
library(Glimma)
library(edgeR)
library(biomaRt)
# library(dplyr)
library(tidyr)
library(RColorBrewer)
library(DESeq2)
library(DEFormats)
library(ggplot2) 
library(Rvisdiff)
library(org.Rn.eg.db)
library(goseq)
library(GenomicFeatures)
library(coriell)
library(tidyverse)
# library(stringr)

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)

load("CABvsDMSO_et.RData")
load("LENvsDMSO_et.RData")
load("CABvsLEN_et.RData")

samples <- read.csv("sample_info.csv")

condition <- factor(c(rep("DMSO",6),rep("CAB",6),rep("LEN",6)))
contraststatement<-c("condition","CAB","DMSO")
comparison <- "CAB vs DMSO"
sample <- samples$full_file
day <- factor(c(rep('1',6)))
mark <- "Expression"
(coldata <- data.frame(row.names=samples$file_header, condition, sample, day))





# clustering ----------------------------------------------------
## https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html

res_df <- edger_to_df(et.CABvsDMSO)

candidate_genes <- res_df %>% 
  filter(PValue < 0.01) %>%    # filter table
  pull(EntrezGene) %>%             # extract the gene column as a vector
  unique()                   # retain only unique values


# trans_cts <- read_csv("/Users/coraalbers/Documents/rna_seq_test_data/counts_transformed.csv")

# Summarise counts 
trans_cts_mean <- counts_df %>% 
  # convert to long format
  pivot_longer(cols = samples$full_file, names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(coldata, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(entrez %in% candidate_genes) %>% 
  # for each gene
  group_by(entrez) %>% 
  # scale the cts column
  mutate(cts_scaled = (cts - mean(cts))/sd(cts)) %>% 
  # for each gene, strain and minute
  group_by(entrez, condition) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()


hclust_matrix <- res_df %>% 
  select(-c(GeneID,feature_id)) %>%
  select(-EntrezGene)  %>%
  as.matrix

rownames(hclust_matrix) <- res_df$EntrezGene


hclust_matrix <- hclust_matrix[rownames(hclust_matrix) %in% candidate_genes, ]

class(hclust_matrix) <- "numeric"

hclust_matrix <- subset(hclust_matrix, select = -c(Symbol,RefSeq,Length) )
# 
# hclust_matrix['logFC'] <- as.numeric(as.character(hclust_matrix['logFC']))
# sapply(hclust_matrix, class)

# 
# hclust_matrix <- as.numeric(hclust_matrix)
# 
# hclust_matrix <- hclust_matrix[!is.na(hclust_matrix)]

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scaling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


hc <- hclust(as.dist(1-cor(hclust_matrix, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

sampleTree = as.dendrogram(hc, method="average")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")

gene_dist <- dist(hclust_matrix)


gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust, labels = FALSE)
abline(h = 3, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram


gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(entrez = name, cluster = value)

head(gene_cluster, n=20)

trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "entrez")

head(trans_cts_cluster)


# Heatmap(hclust_matrix, show_row_names = FALSE)


# work backwards from clustering ------------------------------------------

cluster.1 <- gene_cluster[gene_cluster$cluster == 1,]
cluster.2 <- gene_cluster[gene_cluster$cluster == 2,]
cluster.3 <- gene_cluster[gene_cluster$cluster == 3,]
cluster.4 <- gene_cluster[gene_cluster$cluster == 4,]
cluster.5 <- gene_cluster[gene_cluster$cluster == 5,]
# cluster.6 <- gene_cluster[gene_cluster$cluster == 6,]
# cluster.7 <- gene_cluster[gene_cluster$cluster == 7,]
# cluster.8 <- gene_cluster[gene_cluster$cluster == 8,]
# cluster.9 <- gene_cluster[gene_cluster$cluster == 9,]
# cluster.10 <- gene_cluster[gene_cluster$cluster == 10,]


res.1 <- et.CABvsDMSO[et.CABvsDMSO$genes$EntrezGene %in% cluster.1$entrez, ]
res.2 <- et.CABvsDMSO[et.CABvsDMSO$genes$EntrezGene %in% cluster.2$entrez, ]
res.3 <- et.CABvsDMSO[et.CABvsDMSO$genes$EntrezGene %in% cluster.3$entrez, ]
res.4 <- et.CABvsDMSO[et.CABvsDMSO$genes$EntrezGene %in% cluster.4$entrez, ]
res.5 <- et.CABvsDMSO[et.CABvsDMSO$genes$EntrezGene %in% cluster.5$entrez, ]

go.1 <- goana(res.1, geneid="EntrezGene", species="Rn")
topGO(go.1, n=15)

go.2 <- goana(res.2, geneid="EntrezGene", species="Rn")
topGO(go.2, n=15)

go.3 <- goana(res.3, geneid="EntrezGene", species="Rn")
topGO(go.3, n=15)

go.4 <- goana(res.4, geneid="EntrezGene", species="Rn")
topGO(go.4, n=15)

go.5 <- goana(res.5, geneid="EntrezGene", species="Rn")
topGO(go.5, n=15)


