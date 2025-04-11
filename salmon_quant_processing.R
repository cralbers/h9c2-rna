library("tximeta")
library("GenomicFeatures")
library("tximportData")
library("tximport")
library(AnnotationHub)
library(DESeq2)
library(limma)
library(edgeR)
library(apeglm)
library(RMariaDB)


setwd("~/Documents/smith_lab/h9c2_RNA_seq")
dir <- "quants"
list.files(dir)

csvfile <- file.path("sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors = FALSE)
rownames(coldata) <- coldata$run
## coldata$names <- coldata$file_name
## coldata$files <- file.path("quants", coldata$names, "quant.sf")
files <- file.path("quants", coldata$file_name, "quant.sf")
names(files) <- coldata$run


## file.exists(coldata$files)

txdb <- makeTxDbFromEnsembl(organism = "Rattus norvegicus", release = 104)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

saveDb(z, "whatevs")
z <- loadDb("whatevs")


coldata$names

## summarized to gene level
txi <- tximport(files, type = "salmon",tx2gene=tx2gene, ignoreTxVersion=TRUE)
names(txi)
head(txi$counts)
rownames(coldata) <- colnames(txi$counts)

## original transcript level estimates as a list of matrices
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)

## create table that corresponds sample name with its drug condition
sampleTable <- data.frame(condition = factor(rep(c("DMSO", "CAB", "LEN"), each = 6)))
rownames(sampleTable) <- colnames(txi$counts)

ddstxi <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
dds <- DESeq(ddstxi)
res <- results(dds)
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_DMSO_vs_CAB", type="normal")

resLFC
plotMA(resLFC, ylim=c(-2,2))
