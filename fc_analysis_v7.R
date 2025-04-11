# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#using-select

## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#useful-graphical-representations-of-differential-expression-results

library(limma)
library(Glimma)
library(edgeR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(DESeq2)
library(DEFormats)
library(ggplot2)
library(Rvisdiff)
library(org.Rn.eg.db)
library(goseq)
library(GenomicFeatures)


workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)

# counts_file <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/counts/h9c2_fc_frac_gene_rnaseq.txt"
# counts <- read.delim(counts_file, sep="\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# load in featurecounts data (from Rsubread.r file)
load("fc_basic.RData")


x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])


# get condensed sample names and reassign to column names of DGE list
samplenames <- substring(colnames(x), 1, 6)
colnames(x) <- samplenames

# add dosing info to sample info table
group <- as.factor(rep(c("DMSO", "CAB", "LEN"), c(6, 6, 6)))
x$samples$group <- group

# remove gene version info from ensembl gene ids in DGEList object
geneid_withversion <- rownames(x)

split_characters <- function(x) {
  # Use regexpr and substr to split characters at the first period
  before_period <- character(length(x))
  after_period <- character(length(x))
  
  for (i in seq_along(x)) {
    # Find the position of the first period
    split_pos <- regexpr("\\.", x[i])
    
    # Extract the part before the period
    before_period[i] <- substr(x[i], 1, split_pos - 1)
    
    # Extract the part after the period
    after_period[i] <- substr(x[i], split_pos + 1, nchar(x[i]))
  }
  
  # Return as a list of two vectors
  return(list(before_period = before_period, after_period = after_period))
}

result <- split_characters(geneid_withversion)
geneid_noversion <- result$before_period

rownames(x) <- geneid_noversion
head(rownames(x))

# idfound <- rownames(x$genes) %in% mappedRkeys(org.Rn.egENSEMBL)





# add Entrez gene ids to annotation
## https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
## step 3, example 2
annots <- select(org.Rn.eg.db, keys=rownames(x), columns=c("REFSEQ", "ENTREZID", "ENSEMBL", "SYMBOL"), keytype="ENSEMBL")
m <- match(rownames(x$genes), annots$ENSEMBL)
x$genes$EntrezGene <- annots$ENTREZID[m]
x$genes$Symbol <- annots$SYMBOL[m]
x$genes$RefSeq <- annots$REFSEQ[m]




# look at how many samples have low/0 reads and remove genes with low # of reads
table(rowSums(x$counts==0)==9)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
#rownames(x) <- x$genes$external_gene_name

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

lcpm <- cpm(x, log=TRUE)

# design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

# contrast matrix
contr.matrix <- makeContrasts(
  CABvsDMSO = CAB - DMSO, 
  LENvsDMSO = LEN - DMSO, 
  CABvsLEN = CAB - LEN, 
  levels = colnames(design))

v <- voom(x, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))



# dispersion estimate
x <- estimateDisp(x, design, robust=TRUE)
x$common.dispersion 
plotBCV(x)

et.CABvsDMSO <- exactTest(x, pair=c("DMSO", "CAB"))
et.LENvsDMSO <- exactTest(x, pair=c("DMSO", "LEN"))
et.CABvsLEN <- exactTest(x, pair=c("LEN", "CAB"))

summary(et.CABvsDMSO)

commongenes_cd <- rownames(et.CABvsDMSO$genes)
commongenes_ld <- rownames(et.LENvsDMSO$genes)
commongenes_cl <- rownames(et.CABvsLEN$genes)

commongenes <- data.frame(CABvsDMSO = commongenes_cd, LENvsDMSO = commongenes_ld, CABvsLEN = commongenes_cl)
dim(commongenes)
common_genes <- Reduce(intersect, commongenes)
length(common_genes)
# write.csv(data.frame(common_genes), 'common_genes.csv', row.names = FALSE)
# had to delete title from the csv

# write.table(et.CABvsLEN$table$logFC, "cabvslen_de.csv", sep=",",  col.names=FALSE, row.names = FALSE)


topTags(et.CABvsDMSO)


plotMD(et.CABvsDMSO)

common_genes <- Reduce(intersect, df)

# goseq -------------------------------------------------------------------

# https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf

genes <- as.integer(p.adjust(et.CABvsDMSO$table$PValue[et.CABvsDMSO$table$logFC != 0], method = "BH") < .05)
names(genes) <- row.names(et.CABvsDMSO$table[et.CABvsDMSO$table$logFC != 0, ])
table(genes)

pwf <- nullp(genes, "rn4", "ensGene")
GO.wall <- goseq(pwf, "rn4", "ensGene")
head(GO.wall)



# ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror="useast")
# genes <- AnnotationDbi::select(ensembl, keys=geneid_noversion, columns=(c("external_gene_name", "rgd_symbol", "go_id", "wikigene_id", "ensembl_gene_id")), keytype="ensembl_gene_id")

## https://www.biostars.org/p/191626/
# rn_gene_ensembl <- makeTxDbFromBiomart(dataset="rnorvegicus_gene_ensembl")
# txsByGene=transcriptsBy(rn_gene_ensembl,"gene")
# lengthData=median(width(txsByGene))
# select_genes <- as.vector(names(lengthData)%in%names(genes))
# new_lengthData <- lengthData[select_genes]
# 


# clustering --------------------------------------------------------------
## https://github.com/hamidghaedi/Unsupervised_Clustering_for_Gene_Expression_data

dds = as.DESeqDataSet(x)

vsd <- assay(vst(dds))
mads=apply(vsd,1,mad)
# check data distribution
hist(mads, breaks=nrow(vsd)*0.1)
# selecting features
mad2k=vsd[rev(order(mads))[1:2000],]
mad4k=vsd[rev(order(mads))[1:4000],]
#mad6k=vsd[rev(order(mads))[1:6000],]

hist(mad4k, breaks=nrow(vsd)*0.1)
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(mad4k,
                               maxK=10,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title= "geneExp",
                               clusterAlg="pam",
                               distance="spearman",
                               # seed=1262118388.71279,
                               plot="pdf")



