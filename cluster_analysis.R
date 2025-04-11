library(edgeR, goseq)
library(RColorBrewer)
library(plotly)

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)


load("et.CABvsDMSO_toptags.RData")

# indiv components of analyze go terms fx ----
## import cluster csv ----

cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs"

cluster <- read.csv(file.path(cluster_path, "h9c2_cab_filt35sum_kmeanscluster1.csv"))

cluster1genes <- cluster$X
write.table(cluster1genes, file = "cluster1genes.txt", sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## get list of ens genes from cluster ----

cluster_genes <- cluster$X

## get de info from toptags table for only genes in cluster ----

cluster_de <- et.CABvsDMSO_toptags[c(cluster_genes),]
  

## goseq ----

genes <- as.integer(p.adjust(cluster_de$PValue[cluster_de$logFC != 0], method = "BH") < .05)
names(genes) <- row.names(cluster_de[cluster_de$logFC != 0, ])
table(genes)

pwf <- nullp(genes, "rn4", "ensGene")
GO.wall <- goseq(pwf, "rn4", "ensGene")
head(GO.wall)


# analye go terms in clusters ----


cluster_go_analyze <- function(filename){
  
  # import cluster csv
  
  cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs"
  
  cluster <- read.csv(file.path(cluster_path, filename))
  
  # get cluster number

  string <- strsplit(filename, "\\.")[[1]]
  number <- strsplit(string, "cluster")[[1]][2]
  
  # get list of ens genes from cluster
  
  cluster_genes <- cluster$X
  
  # get de info from toptags table for only genes in cluster
  
  cluster_de <- et.CABvsDMSO_toptags[c(cluster_genes),]
  
  # goseq
  
  genes <- as.integer(p.adjust(cluster_de$PValue[cluster_de$logFC != 0], method = "BH") < .05)
  names(genes) <- row.names(cluster_de[cluster_de$logFC != 0, ])
  table(genes)
  
  pwf <- nullp(genes, "rn4", "ensGene")
  GO.wall <- goseq(pwf, "rn4", "ensGene", test.cats = c("GO:MF"))
  head(GO.wall)
  
  return(GO.wall)
}

cluster_files <- list.files(cluster_path)

cluster1go <- cluster_go_analyze(cluster_files[1])
cluster2go <- cluster_go_analyze(cluster_files[2])
cluster3go <- cluster_go_analyze(cluster_files[3])
cluster4go <- cluster_go_analyze(cluster_files[4])
cluster5go <- cluster_go_analyze(cluster_files[5])
cluster6go <- cluster_go_analyze(cluster_files[6])
cluster7go <- cluster_go_analyze(cluster_files[7])
cluster8go <- cluster_go_analyze(cluster_files[8])
cluster9go <- cluster_go_analyze(cluster_files[9])
cluster10go <- cluster_go_analyze(cluster_files[10])
cluster11go <- cluster_go_analyze(cluster_files[11])
cluster12go <- cluster_go_analyze(cluster_files[12])
cluster13go <- cluster_go_analyze(cluster_files[13])
cluster14go <- cluster_go_analyze(cluster_files[14])
cluster15go <- cluster_go_analyze(cluster_files[15])
cluster16go <- cluster_go_analyze(cluster_files[16])
cluster17go <- cluster_go_analyze(cluster_files[17])


cluster_df_list <- list(cluster1go, cluster2go, cluster3go, cluster4go, cluster5go, cluster6go, cluster7go, cluster8go, cluster9go, cluster10go, cluster11go, cluster12go, cluster13go, cluster14go, cluster15go, cluster16go, cluster17go)

# save(cluster17go, file= "cluster17go.RData")

# for (i in 1:18) {
#   paste0("clustergo", i) <- variable
#   variable <- cluster_go_analyze(cluster_files[i])
#   return(variable)
# }
# 

setwd("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_rdata_MF")
save(cluster_df_list,file='clustergo_MF.Rdata')

nums <- c(1:17)

for( i in 1:length(cluster_df_list)) {
  write.csv(cluster_df_list[i],paste0('cluster', nums[i],'.txt'))
}
  



# for(i in 1:length(cluster_files)) {
#   
#   filename <- cluster_files[i]
#   string <- strsplit(filename, "\\.")[[1]]
#   number <- strsplit(string, "cluster")[[1]][2]
#   print(i)
#   
#   save(list=cluster_df_list[i], file= paste0("cluster", number, "go.RData"))
# }

cluster1go_exp <- select(cluster1go, c('category', 'over_represented_pvalue', 'term'))
write.csv(cluster1go_exp,'cluster1go_exp.txt', row.names = FALSE, quote = FALSE, col.names = FALSE)

setwd("/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_rdata")
file_names=as.list(dir(pattern="cluster*"))
for(i in 1:length(file_names)) load(file_names[[i]]) 

setwd(workdir)


# volcano labeled by cluster ----

et.CABvsDMSO_toptags$cluster <- NA


for(i in 1:length(cluster_files)) {
  
  filename <- cluster_files[i]
  
  # import cluster csv
  
  cluster_path <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/RNAlysis/cab_only_cluster/cluster_csvs"
  
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
  
et.CABvsDMSO_toptags$cluster <- factor(et.CABvsDMSO_toptags$cluster, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, NA))


sum(is.na(et.CABvsDMSO_toptags$cluster))

colorRampPalette(c("blue", "red"))(17)

colors = c("#7e65cf",
            "#64be4b",
            "#c151b6",
            "#a8b635",
            "#6a85c9",
            "#d99836",
            "#4cbad1",
            "#d04c33",
            "#5fc48d",
            "#cb4270",
            "#4e8c39",
            "#be73a7",
            "#388864",
            "#d4786d",
            "#72742b",
            "#a26431",
            "#b9ac64")
  
p <- ggplot(data=et.CABvsDMSO_toptags, aes(x=logFC, y=-log10(PValue), col=cluster)) + 
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(text=element_text(size=16,  family="Times")) + 
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
  # scale_fill_discrete(breaks=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', 'NA'))
p

ggplotly(p)


# c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', 'NA')
# c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, NA)

rownames(et.CABvsDMSO_toptags)[1]


# interactive volcano? ----

fig <- plot_ly(data = et.CABvsDMSO_toptags, x = ~logFC, y = ~-log10(PValue), color = ~cluster, colors = ~colors, text = ~Symbol, mode='markers')


fig













