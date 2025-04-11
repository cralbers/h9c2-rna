# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#using-select

## https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#useful-graphical-representations-of-differential-expression-results

## https://genomicsclass.github.io/book/pages/rnaseq_gene_level.html not helpful i don't think
## https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

library(limma)
library(Glimma)
library(edgeR)
library(biomaRt)
# library(dplyr)
# library(tidyr)
library(RColorBrewer)
library(DESeq2)
library(DEFormats)
library(ggplot2) 
library(Rvisdiff)
library(org.Rn.eg.db)
library(goseq)
library(GenomicFeatures)
library(coriell)
# library(tidyverse)
# library(stringr)
library(EnhancedVolcano)

workdir <- '/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R'
setwd(workdir)



# samples <- read.csv("")

# counts_file <- "/Users/coraalbers/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/counts/h9c2_fc_frac_gene_rnaseq.txt"
# counts <- read.delim(counts_file, sep="\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# load in featurecounts data (from Rsubread.r file) ----
load("fc_basic.RData")



geneid_withversion <- rownames(fc$counts)


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

rownames(fc$counts) <- geneid_noversion
head(rownames(fc$counts))

counts <- fc$counts

# write.csv(counts, "raw_counts.csv", quote=FALSE)

group <- as.factor(rep(c("DMSO", "CAB", "LEN"), c(6, 6, 6)))

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")], group=group)


## get condensed sample names and reassign to column names of DGE list ----
samplenames <- substring(colnames(x), 1, 6)
colnames(x) <- samplenames

## add dosing info to sample info table ----

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

# result <- split_characters(geneid_withversion)
# geneid_noversion <- result$before_period
# 
# rownames(x) <- geneid_noversion
head(rownames(x))

counts_df <- as.data.frame(counts)
counts_df$ensembl <- rownames(counts_df)
# counts_ens_version <- rownames(counts_df)
# counts_ens_split_result <- split_characters(counts_ens_version)
# counts_ens_noversion <- counts_ens_split_result$before_period
# rownames(counts_df) <- counts_ens_noversion
counts_df$gene <- mapIds(org.Rn.eg.db, rownames(counts_df), keytype="ENSEMBL", column="SYMBOL")
counts_df$entrez <- mapIds(org.Rn.eg.db, rownames(counts_df), keytype="ENSEMBL", column="ENTREZID")

# write.csv(counts_df, "raw_counts_df.csv", row.names = TRUE, quote = FALSE)

# idfound <- rownames(x$genes) %in% mappedRkeys(org.Rn.egENSEMBL)





# add Entrez gene ids to annotation ----

## https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
## step 3, example 2
annots <- select(org.Rn.eg.db, keys=rownames(x), columns=c("REFSEQ", "ENTREZID", "ENSEMBL", "SYMBOL"), keytype="ENSEMBL")
m <- match(rownames(x$genes), annots$ENSEMBL)
x$genes$EntrezGene <- annots$ENTREZID[m]
x$genes$Symbol <- annots$SYMBOL[m]
x$genes$RefSeq <- annots$REFSEQ[m]




# look at how many samples have low/0 reads and remove genes with low # of reads
x.full <- x # backup of full df

table(rowSums(x$counts==0)==9)
head(cpm(x))
apply(x$counts, 2, sum)

### oh hell no this cut it down to ~2200 de genes thats cazy
# # keep only genes with at least 100 cpm for at least 2 samples
# keep <- rowSums(cpm(x)>100) >= 2
# x1 <- x[keep,]
# dim(x1)

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
#rownames(x) <- x$genes$external_gene_name

# reset lib sizes
x$samples$lib.size <- colSums(x$counts)
x$samples

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

plotMDS(x, col=as.numeric(x$samples$group))







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

# v <- voom(x, design, plot=TRUE)
# vfit <- lmFit(v, design)
# vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
# efit <- eBayes(vfit)
# plotSA(efit, main="Final model: Mean-variance trend")
# 
# summary(decideTests(efit))



# dispersion estimate ----
## using common dispersion (assuming everything has the same common dispersion) this est is not great ----
x1 <- estimateDisp(x, design, robust=TRUE)
x1$common.dispersion 
x1 <- estimateTagwiseDisp(x1)
names(x1)
plotBCV(x1)

## estimate dispersion ----
x2 <- estimateGLMCommonDisp(x,design)
x2 <- estimateGLMTrendedDisp(x2,design, method="auto")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
x2 <- estimateGLMTagwiseDisp(x2,design)
plotBCV(x2)

fit <- glmFit(x1, design)

getCounts(x2)
getOffset(x2)
getDispersion(x2)

lcpm <- cpm(x2, log=TRUE)
# save(lcpm, file = "lcpm.RData")
cpm <- cpm(x2)
# save(cpm, file = "cpm.RData")
# write.csv(cpm, "cpm.csv", quote=FALSE, row.names = TRUE)

cpm_cab <- subset(cpm, select = -c(X12668, X12669, X12670, X12671, X12672, X12673))
# write.csv(cpm_cab, "cpm_cab.csv", quote=FALSE, row.names = TRUE)


# df_new <- df[1:(length(df)-1)]

# test <- exactTest(x2)

et.CABvsDMSO <- exactTest(x2, pair=c("DMSO", "CAB"))
et.LENvsDMSO <- exactTest(x2, pair=c("DMSO", "LEN"))
et.CABvsLEN <- exactTest(x2, pair=c("LEN", "CAB"))
summary(et.CABvsDMSO)

# et.CABvsLEN$genes$Ensembl <- rownames(et.CABvsLEN)
# write.csv(et.LENvsDMSO, "LENvsDMSO_et.csv", quote=FALSE, row.names = TRUE)

# lrtCD <- glmLRT(fit, contrast = c(1, -1, 0))
# lrtLD <- glmLRT(fit, contrast = c(0, -1, 1))
# lrtCL <- glmLRT(fit, contrast = c(1, 0, -1))

# save(et.CABvsDMSO, file = "CABvsDMSO_et.RData")
# save(et.LENvsDMSO, file = "LENvsDMSO_et.RData")
# save(et.CABvsLEN, file = "CABvsLEN_et.RData")

CABvsDMSO_df <- as.data.frame(et.CABvsDMSO, row.names = rownames(et.CABvsDMSO))
CABvsDMSO_df$padj <- p.adjust(CABvsDMSO_df$PValue,"fdr")
CABvsDMSO_df$Symbol <- as.character(CABvsDMSO_df$Symbol)
# write.csv(CABvsDMSO_df, 'CABvsDMSO_df_et_padj.csv', row.names = TRUE, quote = FALSE)


LENvsDMSO_df <- as.data.frame(et.LENvsDMSO, row.names = rownames(et.LENvsDMSO))
LENvsDMSO_df$padj <- p.adjust(LENvsDMSO_df$PValue,"fdr")

# avereps(et.CABvsDMSO)

et.CABvsDMSO_toptags <- topTags(et.CABvsDMSO, n = "Inf")$table
# write.csv(out, "et.CABvsDMSO_toptags.csv", quote=FALSE, row.names = TRUE)
# save(et.CABvsDMSO_toptags, file= "et.CABvsDMSO_toptags.RData")



glimmaVolcano(et.CABvsDMSO)

load("et.CABvsDMSO_toptags.RData")


# make volcano plot ------------------------------------------------------


make_volcano <- function(df, x, y, title){
  upreg <- nrow(subset(df,logFC > 1 & FDR < 0.05))
  downreg <- nrow(subset(df,logFC < (-1) & FDR < 0.05))
  
  EnhancedVolcano(df,
                  lab = NA, #CABvsDMSO_df$Symbol,
                  selectLab=c("Nfe2l2"),
                  x = 'logFC',
                  y = 'FDR',
                  title = title,
                  pCutoff = 0.05,
                  # FCcutoff = 1,
                  drawConnectors = TRUE,
                  subtitle = "",
                  caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                  boxedLabels = TRUE,
  ) + theme(text=element_text(size=16,  family="Times")) + 
    geom_text(x=x, y=y, label= upreg, family='Times') + 
    geom_text(x=-x, y=y, label= downreg, family="Times") +
    theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
}

# len 2, 20; cab 8, 300
# save at w 700px, h 600px
make_volcano(et.CABvsDMSO_toptags, 8, 300, "CAB vs DMSO")
# make_volcano(LENvsDMSO_df, 2, 20, "LEN vs DMSO")
  

# specific volcanoes ----

# define different cell-types that will be shaded
caspase <- c('Casp3', 'Casp7', 'Casp9')
calpain <- c('Capn1', 'Capn3', 'Capn6')
mito <- c("Tp53", 'Bcl2l1', 'Bax', 'Sod2', 'Cycs')
ion_channels <- c("Scn2a", "Kcnh2", 'Kcnq5', 'Cacna1c')
select_genes <- c(caspase, calpain, mito) #, ion_channels)

antiox <- c('Sod3', 'Gsta1', 'Gsta3', 'Nqo1', 'Gss')

## mito caspase volcano ----------------------------------------------------


keyvals.colour <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'gold',
         'black'))

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements

# keyvals <- ifelse(
#   CABvsDMSO_df$Symbol %in% select_genes, 'gold', 'black')
# keyvals[is.na(keyvals)] <- 'black'

# keyvals.shape <- ifelse(
#   CABvsDMSO_df$Symbol %in% caspase, 16,
#   ifelse(CABvsDMSO_df$Symbol %in% calpain, 2,
#          ifelse(CABvsDMSO_df$Symbol %in% mito, 5, 
#                 ifelse(CABvsDMSO_df$Symbol %in% ion_channels, 8, 1
#          ))))

keyvals.shape <- ifelse(
  CABvsDMSO_df$Symbol %in% caspase, 15,
  ifelse(CABvsDMSO_df$Symbol %in% calpain, 17,
         ifelse(CABvsDMSO_df$Symbol %in% mito, 18, 
                16
                )))

keyvals.shape[is.na(keyvals.shape)] <- 16
names(keyvals.shape)[keyvals.shape == 16] <- 'All'
names(keyvals.shape)[keyvals.shape == 15] <- 'Caspase'
names(keyvals.shape)[keyvals.shape == 17] <- 'Calpain'
names(keyvals.shape)[keyvals.shape == 18] <- 'Mitochondrial stress markers'
# names(keyvals.shape)[keyvals.shape == 8] <- 'Ion channel'

keyvals.color <- ifelse(
  CABvsDMSO_df$Symbol %in% select_genes, 'gold',
         'black')
keyvals.color[is.na(keyvals.color)] <- 'black'
names(keyvals.color)[keyvals.color == 'gold'] <- 'Gene of Interest'
# names(keyvals.colour)[keyvals.colour == 'black'] <- 'mid'
# names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'

keyvals.shape_basic <- ifelse(
  CABvsDMSO_df$Symbol %in% select_genes, 19, 17)
keyvals.shape_basic[is.na(keyvals.shape_basic)] <- 17
names(keyvals.shape_basic)[keyvals.shape_basic == 19] <- 'Gene of Interest'

EnhancedVolcano(et.CABvsDMSO_toptags,
                lab = et.CABvsDMSO_toptags$Symbol,
                selectLab=select_genes,
                x = 'logFC',
                y = 'FDR',
                title = "CAB vs DMSO",
                pCutoff = 0.05,
                # FCcutoff = 1,
                drawConnectors = TRUE,
                # widthConnectors = 0.75,
                subtitle = "",
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                boxedLabels = TRUE,
                labSize = 3,
                # shapeCustom = keyvals.shape_basic,
                # colCustom = keyvals.color,
                # colAlpha = 0.9,
                max.overlaps = Inf
                ) + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

## mito only ----

keyvals.color_basic <- ifelse(et.CABvsDMSO_toptags$Symbol %in% mito, 'red', 'lightgray')
keyvals.color_basic[is.na(keyvals.color_basic)] <- 'lightgray'
names(keyvals.color_basic)[keyvals.color_basic == 'lightgray'] <- 'All'
names(keyvals.color_basic)[keyvals.color_basic == 'red'] <- 'Gene of Interest'

keyvals.shape_basic <- ifelse(et.CABvsDMSO_toptags$Symbol %in% mito, 19, 1)
keyvals.shape_basic[is.na(keyvals.shape_basic)] <- 1
names(keyvals.shape_basic)[keyvals.shape_basic == 1] <- 'All'
names(keyvals.shape_basic)[keyvals.shape_basic == 19] <- 'Gene of Interest'

mito_volcano = EnhancedVolcano(et.CABvsDMSO_toptags,
                lab = et.CABvsDMSO_toptags$Symbol,
                selectLab=mito,
                x = 'logFC',
                y = 'FDR',
                title = "CAB vs DMSO",
                pCutoff = 0.05,
                # FCcutoff = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                subtitle = "",
                caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                boxedLabels = TRUE,
                labSize = 3,
                shapeCustom = keyvals.shape_basic,
                colCustom = keyvals.color_basic,
                colAlpha = 0.3,
                max.overlaps = Inf
) + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

ggsave(file="mito_volcano.pdf", plot=mito_volcano, width=7, height=7)
  
  # geom_text(aes(label=upreg,x=Inf,y=Inf,vjust=2,hjust=2))+
  # geom_text(aes(label=downreg,x=-10,y=Inf,vjust=2))
  


## nrf2 volcano ----
nrf2 <- c("Nfe2l2")
keyvals.color <- ifelse(et.CABvsDMSO_toptags$Symbol %in% nrf2, 'red', 'lightgray')
keyvals.color[is.na(keyvals.color)] <- 'lightgray'
names(keyvals.color)[keyvals.color == 'lightgray'] <- 'All'
names(keyvals.color)[keyvals.color == 'red'] <- 'Nrf2'

nrf2_volcano = EnhancedVolcano(et.CABvsDMSO_toptags,
                               lab = et.CABvsDMSO_toptags$Symbol,
                               selectLab=c("Nfe2l2"),
                               x = 'logFC',
                               y = 'FDR',
                               title = "CAB vs DMSO",
                               pCutoff = 0.05,
                               # FCcutoff = 1,
                               drawConnectors = TRUE,
                               widthConnectors = 0.2,
                               subtitle = "",
                               caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                               boxedLabels = TRUE,
                               labSize = 3,
                               # shapeCustom = keyvals.shape_basic,
                               colCustom = keyvals.color,
                               colAlpha = 0.3,
                               max.overlaps = Inf
) + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

ggsave(file="nrf2_volcano.pdf", plot=nrf2_volcano, width=7, height=7)

# default p-value cutoff: 10e-6 
# default log2FC cutoff: >|2|
dds = as.DESeqDataSet(x2, colData=coldata)

## antiox only ----

keyvals.color_basic <- ifelse(et.CABvsDMSO_toptags$Symbol %in% antiox, 'red', 'lightgray')
keyvals.color_basic[is.na(keyvals.color_basic)] <- 'lightgray'
names(keyvals.color_basic)[keyvals.color_basic == 'lightgray'] <- 'All'
names(keyvals.color_basic)[keyvals.color_basic == 'red'] <- 'Gene of Interest'

keyvals.shape_basic <- ifelse(et.CABvsDMSO_toptags$Symbol %in% antiox, 19, 1)
keyvals.shape_basic[is.na(keyvals.shape_basic)] <- 1
names(keyvals.shape_basic)[keyvals.shape_basic == 1] <- 'All'
names(keyvals.shape_basic)[keyvals.shape_basic == 19] <- 'Gene of Interest'

mito_volcano = EnhancedVolcano(et.CABvsDMSO_toptags,
                               lab = et.CABvsDMSO_toptags$Symbol,
                               selectLab=antiox,
                               x = 'logFC',
                               y = 'FDR',
                               title = "CAB vs DMSO",
                               pCutoff = 0.05,
                               # FCcutoff = 1,
                               drawConnectors = TRUE,
                               widthConnectors = 0.2,
                               subtitle = "",
                               caption = bquote(~Log[2]~ "fold change cutoff, 1; p-value cutoff, 0.05"),
                               boxedLabels = TRUE,
                               labSize = 3,
                               # shapeCustom = keyvals.shape_basic,
                               # colCustom = keyvals.color_basic,
                               colAlpha = 0.3,
                               max.overlaps = Inf
) + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

mito_volcano

ggsave(file="antiox_volcano_sigcolor.pdf", plot=mito_volcano, width=7, height=7)

geneCounts<- plotCounts(dds, gene="ENSRNOG00000001548.6", main="Nrf2", intgroup="condition", returnData = TRUE)
p1 <- ggplot(geneCounts, aes(x = condition, y = count, color = condition)) +
  geom_point(aes(fill=condition),size=4, pch=21, col="black") 
p1 <- p1 + theme(
  plot.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white", colour="black"),
  axis.title.x = element_blank(),
  legend.position = "none",
  plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(color = "black")) + labs(title="c-Met", y="Normalized Count")

commongenes_cd <- rownames(et.CABvsDMSO$genes)
commongenes_ld <- rownames(et.LENvsDMSO$genes)
commongenes_cl <- rownames(et.CABvsLEN$genes)

commongenes <- data.frame(CABvsDMSO = commongenes_cd, LENvsDMSO = commongenes_ld, CABvsLEN = commongenes_cl)
dim(commongenes)
common_genes <- Reduce(intersect, commongenes)
length(common_genes)
# common_genes <- Reduce(intersect, df)

# write.csv(data.frame(common_genes), 'common_genes.csv', row.names = FALSE)
# had to delete title from the csv

# write.table(et.LENvsDMSO$table$logFC, "lenvsdmso_de.csv", sep=",",  col.names=FALSE, row.names = FALSE, quote=FALSE)


cabtop20 <- topTags(et.CABvsDMSO, n=20)$table
# write.csv(cabtop20, 'cabtop20.csv', row.names = TRUE, quote = FALSE)
# save(cabtop20, file = "cabtop20.RData")

plotMD(et.CABvsDMSO)



pheatmap(lcpm)

# number of de genes with fdr < 0.05
de1 <- decideTests(et.CABvsDMSO, adjust.method="BH", p.value=0.05)
summary(de1)

de2 <- decideTests(et.LENvsDMSO, adjust.method="BH", p.value=0.05)
summary(de2)

de3 <- decideTests(et.CABvsLEN, adjust.method="BH", p.value=0.05)
summary(de3)


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

dds = as.DESeqDataSet(x2)

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


