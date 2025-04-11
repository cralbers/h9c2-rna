setwd("~/Documents/smith_lab/h9c2_RNA_seq")

library(tximport)
library(tximportData)
library(ggplot2)
library(readr)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(Glimma)
library(limma)
library(edgeR)
library("tximeta")
library("GenomicFeatures")
library(txdbmaker)

### set up gencode VM23, tx2gene, orthologs
#gencode_vm23<-read.table("~/Documents/smith_lab/RENCA_RNA-seq/Processing_files/gencode_info_VM23.txt",header=T)
#tx2gene<-gencode_vm23[,grep("transcript_id|gene_name",names(gencode_vm23))]
#orthologs<-read.table("~/Documents/smith_lab/RENCA_RNA-seq/Processing_files/Mouse_Human_Orthologs_MGI.txt",header=T)
#head(tx2gene)

#txdb <- makeTxDbFromEnsembl(organism = "Rattus norvegicus", release = 113)
#saveDb(txdb, "txdb")
txdb <- loadDb("txdb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

quants_dir <- "quants"
samples <- read.csv(file.path("sample_table.csv"), header = TRUE)
files <- file.path(quants_dir, samples$file_name, "quant.sf")

### import data
txi<-tximport(files,type="salmon",tx2gene=tx2gene, ignoreTxVersion = TRUE)

### set up comparison
#condition <- factor(c(rep("RENCA_P1B9",2),rep("RENCA_P3D4",2),rep("RENCA_WT",2)))
condition <- factor(c(rep("DMSO",6),rep("CAB",6), rep("LEN",6)))
contraststatement<-c("condition","CAB","DMSO")
comparison <- "CAB vs DMSO"
mark <- "Expression"
(coldata <- data.frame(row.names=samples$file_header, condition))

ddsTxi<-DESeqDataSetFromTximport(txi,colData=coldata,design= ~ condition)
dds <- DESeq(ddsTxi)

res<-results(dds,contrast=contraststatement)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE) #merge all sample dataframes into one huge df
names(resdata)<-gsub('log2FoldChange','log2FC',names(resdata))
#resdata_RENCA<-merge(resdata,orthologs,by.x="Row.names",by.y="Mouse") #add in human ortholog info

ggplot(resdata,aes(x=log2FC,y=-log(padj)))+
  geom_point(alpha=0.5)+labs(title=paste0(comparison),subtitle=paste0(mark))+
  xlim(-10,10)+
  #ylim(0,60)+
  geom_text(aes(label=paste0(nrow(subset(resdata,log2FC > 1 & padj < 0.05))),x=Inf,y=Inf,vjust=2,hjust=2))+
  geom_text(aes(label=paste0(nrow(subset(resdata,log2FC < (-1) & padj < 0.05))),x=-10,y=Inf,vjust=2))+
  theme_classic()


glimmaMDS(dds, groups=condition) #works

glimmaMA(dds, groups=condition, main="CAB vs DMSO")

glimmaVolcano(dds, groups=condition, main="CAB vs DMSO")



## save resdata
#write.csv(resdata, "allRENCA_DESeq_results.csv")

#write.csv(res, "RENCASETD2vsWT_ensembl.csv")







