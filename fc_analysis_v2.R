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
## https://support.bioconductor.org/p/97134/
## have to specify AnnotationDbi select command because dplyr also has a select command
genes <- AnnotationDbi::select(ensembl, keys=geneid_noversion, columns=(c("external_gene_name", "rgd_symbol", "go_id", "wikigene_id", "ensembl_gene_id")), keytype="ensembl_gene_id")


#remove duplicated gene ids by keeping only first occurence of gene id
genes_nodups <- genes[!duplicated(genes$ensembl_gene_id),]

# add gene names and info to DGEList object
reordered_indices <- match(geneid_noversion, genes$ensembl_gene_id)
x$genes <- genes[reordered_indices, ]


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


# normalize gene expression distribution 
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

lendmso <- x[,c(1:6,13:18)]
cabdmso <- x[,c(1:12)]
lencab <- x[,c(7:18)]

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

# edgeR dispersion estimate
x <- estimateDisp(x, design, robust=TRUE)
plotBCV(x)

# estimate technical and biological variation
fit <- glmQLFit(x, design, robust=TRUE)
fit$dispersion
plotQLDisp(fit)


con <- makeContrasts(CAB - DMSO, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
# get only genes with diff exp significantly above 1.2
tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
topTags(tr)
summary(decideTests(tr))


anov <- glmQLFTest(fit, contrast=contr.matrix)
topTags(anov)
# plotMD(anov)



# mean variance with limma
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

dds_ds = as.DESeqDataSet(x)

dds_ds$group <- relevel(dds_ds$group, ref="DMSO")


dds <- DESeq(dds_ds)

plotDispEsts(dds)
glimmaMDS(x)

res<-results(dds, contrast=c("group","CAB","DMSO"))
res2<-results(dds, contrast = c("group", 'LEN', "DMSO"))
summary(res2)

# write.csv(res, "CABvsDMSO.csv")


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

df = read.csv("CABvsDMSO.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
names(original_gene_list)

names(original_gene_list) <- df$X

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
  summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE))

# Create a vector of the gene unuiverse
kegg_gene_list <- df2_aggregated$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2_aggregated$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

if(anyDuplicated(names(kegg_gene_list))) {
  stop("There are duplicate ENTREZ IDs in the kegg_gene_list.")
}



kegg_organism = "rno"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


# Produce the native KEGG plot (PNG)
setwd("~/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/pathways")
dme <- pathview(gene.data=kegg_gene_list, pathway.id="rno04210", species = kegg_organism, kegg.native = T)
setwd("~/Documents/smith_lab/h9c2_RNA_seq/h9c2_R/")
 
