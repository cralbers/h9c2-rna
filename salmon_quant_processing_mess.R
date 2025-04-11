install.packages("BiocManager")
BiocManager::install("tximeta")
BiocManager::install("RMariaDB")

library("tximeta")
library("GenomicFeatures")
library("tximportData")
library("tximport")
library(AnnotationHub)

setwd("~/Documents/smith_lab/h9c2_RNA_seq")
dir <- "quants"
list.files(dir)

csvfile <- file.path("sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors = FALSE)
coldata$names <- coldata$file_name
coldata$files <- file.path("quants", coldata$names, "quant.sf")
file.exists(coldata$files)

txdb <- makeTxDbFromEnsembl(organism = "Rattus norvegicus", release = 104)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

ah <- AnnotationHub()
query(ah, "EnsDb")
ahDb <- query(ah, pattern = c("Rattus norvegicus", "EnsDb", 87))
ahEdb <- ahDb[[1]] ## download database


## gns <- genes(ahEdb) ## retrieve all genes
## tx <- transcripts(ahEdb, return.type="DataFrame", )



txi <- tximport(coldata$files, type = "salmon",tx2gene=tx2gene, ignoreTxVersion=TRUE)

head(txi$counts)


fastaFTP <- c("ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/cdna/Drosophila_melanogast
"ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.
gtfPath <- file.path(dir, "Drosophila_melanogaster.BDGP6.22.98.gtf.gz")
makeLinkedTxome(indexDir=indexDir, source="LocalEnsembl", organism="Drosophila melanogaster", release="98", genome="BDGP6.22", fasta=fastaFTP, gtf=gtfPath, write=FALSE)