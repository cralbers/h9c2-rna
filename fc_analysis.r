# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#using-select

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

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)
#getwd()

# bamdir <- '/fs/ess/PAS2905/ILMN_2341_Yan_OSUCM_GSR2Pool_Nov2024_SmithS_BAM-transfer'
# filelist <- list.files(path=bamdir, pattern="\\.bam$", full.names=TRUE)

#genes_file <- "/fs/ess/PAS2905/rn6.ensGene.gtf"
#genes <- read.table(genes_file, header=TRUE, sep="\t")

counts_file <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/counts/h9c2_fc_frac_gene_rnaseq.txt"
counts <- read.delim(counts_file, sep="\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# load in featurecounts data (from Rsubread.r file)
load("fc_basic.RData")

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])



# get condensed sample names and reassign to column names of DGE list
samplenames <- substring(colnames(x), 1, 6)
colnames(x) <- samplenames

# add dosing info to sample info table
group <- as.factor(rep(c("DMSO", "CAB", "LEN"), c(6, 6, 6)))
x$samples$group <- group

x

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

ensembl <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl", mirror='useast')
#keytypes(ensembl)
## https://support.bioconductor.org/p/97134/
## have to specify AnnotationDbi select command because dplyr also has a select command
genes <- AnnotationDbi::select(ensembl, keys=geneid_noversion, columns=(c("external_gene_name", "rgd_symbol", "go_id", "wikigene_id", "ensembl_gene_id")), keytype="ensembl_gene_id")
#keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
nrow(genes)


#remove duplicated gene ids by keeping only first occurence of gene id
genes_nodups <- genes[!duplicated(genes$ensembl_gene_id),]
nrow(genes_nodups)

# add gene names and info to DGEList object
## am i insane
reordered_indices <- match(geneid_noversion, genes$ensembl_gene_id)
y2 <- genes[reordered_indices, ]
head(y2)
head(geneid_noversion)
## no i am not

x$genes <- genes[reordered_indices, ]
x$genes

dds_ds = as.DESeqDataSet(x)
dds_ds

# calculate counts per million and log cpm
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# look at stats for each sample's library
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

# look at how many samples have low/0 reads and remove genes with low # of reads
table(rowSums(x$counts==0)==9)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)


# # gene filtering graph
# lcpm.cutoff <- log2(10/M + 2/L)
# library(RColorBrewer)
# nsamples <- ncol(x)
# #col <- brewer.pal(nsamples, "Paired")
# par(mfrow=c(1,2))
# plot(density(lcpm[,1]),  lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
# title(main="A. Raw data", xlab="Log-cpm")
# abline(v=lcpm.cutoff, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y,  lwd=2)
# }
# #legend("topright", samplenames, bty="n", ncol=3)
# lcpm <- cpm(x, log=TRUE)
# plot(density(lcpm[,1]),  lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
# title(main="B. Filtered data", xlab="Log-cpm")
# abline(v=lcpm.cutoff, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y,  lwd=2)
# }
# #legend("topright", samplenames,  bty="n", ncol=3)


# normalize gene expression distribution 
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

plotMDS(x)



# # MDS
# lcpm <- cpm(x, log=TRUE)
# col.group <- group
# levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
# col.group <- as.character(col.group)
# plotMDS(lcpm, labels=group, col=col.group)
# title(main="Sample groups")
#plotMDS(x)

#glMDSPlot(lcpm, labels=group, 
        #  groups=x$samples[,c(2,5)], launch=FALSE)


# design matrix
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

x <- estimateDisp(x, design, robust=TRUE)
plotBCV(x)

fit <- glmQLFit(x, design, robust=TRUE)
fit$dispersion
plotQLDisp(fit)

#qlf <- glmQLFTest(fit, coef=1:3)


con <- makeContrasts(CAB - DMSO, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
# get only genes with diff exp significantly above 1.2
tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
topTags(tr)
summary(decideTests(tr))


con_all <- makeContrasts(
  LvsD = LEN - DMSO,
  CvsD = CAB - DMSO,
  CvsL = CAB - LEN, levels=design)

anov <- glmQLFTest(fit, contrast=con_all)
topTags(anov)
plotMD(anov)

# contrast matrix
contr.matrix <- makeContrasts(
  CABvsDMSO = CAB - DMSO, 
  LENvsDMSO = LEN - DMSO, 
  CABvsLEN = CAB - LEN, 
  levels = colnames(design))
contr.matrix

# mean variance
# par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


# quick look at DE
summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)







### DESeq2 analysis

dds_ds$condition <- relevel(dds_ds$condition, ref="DMSO")
dds_ds$group <- relevel(dds_ds$group, ref="DMSO")


dds <- DESeq(dds_ds)

plotDispEsts(dds)
glimmaMDS(x)

res<-results(dds, contrast=c("group","CAB","DMSO"))
res2<-results(dds, contrast = c("group", 'LEN', "DMSO"))
summary(res2)

plotMA(res, main="CAB vs DMSO")
plotMA(res2, main="LEN vs DMSO")



get_resdata <- function(res, dds) {
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE) #merge all sample dataframes into one huge df
  names(resdata)<-gsub('log2FoldChange','log2FC',names(resdata))
  return(resdata)
}

resdata <- get_resdata(res, dds)
resdata2 <- get_resdata(res2, dds)



contraststatement<-c("group","LEN","DMSO")
comparison <- "LEN vs DMSO"
mark <- "Expression"


ggplot(resdata2,aes(x=log2FC,y=-log(padj)))+
  geom_point(alpha=0.5)+labs(title=paste0(comparison),subtitle=paste0(mark))+
  xlim(-10,10)+
  #ylim(0,60)+
  geom_text(aes(label=paste0(nrow(subset(resdata,log2FC > 1 & padj < 0.05))),x=Inf,y=Inf,vjust=2,hjust=2))+
  geom_text(aes(label=paste0(nrow(subset(resdata,log2FC < (-1) & padj < 0.05))),x=-10,y=Inf,vjust=2))+
  theme_classic()


glimmaMA(dds, groups=group)

glimmaVolcano(dds, groups=condition, main="RENCA SETD2 mutant vs WT")

dge <- estimateDisp(dge, design)
gfit <- glmFit(dge, design)
glrt <- glmLRT(gfit, design, contrast = contr.matrix)

glimmaMA(efit, dge = x)


# pathway analysis
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

# we want the log2 fold change 
original_gene_list <- resdata$log2FoldChange



