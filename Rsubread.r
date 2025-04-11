workdir <- '/users/PAS2905/coraalbers/h9c2_R'
setwd(workdir)
getwd()

library(Rsubread)
library(limma)
library(edgeR)

bamdir <- '/fs/ess/PAS2905/ILMN_2341_Yan_OSUCM_GSR2Pool_Nov2024_SmithS_BAM-transfer'
filelist <- list.files(path=bamdir, pattern="\\.bam$", full.names=TRUE)

genes_file <- "/fs/ess/PAS2905/rn6.ensGene.gtf"
genes <- read.table(genes_file, header=TRUE, sep="\t")

fc <- featureCounts(files=filelist,
annot.ext=genes_file,
isGTFAnnotationFile=TRUE,
GTF.featureType="exon",
useMetaFeatures=TRUE,
allowMultiOverlap=TRUE,
largestOverlap=TRUE,
countMultiMappingReads=TRUE,
isPairedEnd=TRUE,
countReadPairs=TRUE,
fraction=TRUE,
reportReadsPath="counts",
nthreads=8)

save(fc, file = "fc_frac.RData")

#load("fc.RData")

outdir <- paste(workdir,'/', 'counts', '/', sep='')
# dir.create(outdir)
fileroot <- "h9c2"

write.table(x=data.frame(fc$annotation[,c("GeneID", "Length")],
    fc$counts,stringsAsFactors=FALSE),
    paste0(outdir, fileroot, file="_fc_frac_gene_rnaseq.txt"),
    quote=FALSE,sep="\t",
    row.names=FALSE)

x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

#get only genes that had more than 10 reads per million mapped in at least 2 samples
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

targets <- readTargets()
drug <- factor(targets$treatment)
#design <- model.matrix(~0+drug)
#colnames(design) <- levels(drug)

#design2 <- modelMatrix(targets, ref="DMSO")

design <- model.matrix(~ 0+drug)
colnames(design) <- levels(drug)

y <- voom(x,design,plot=TRUE)

fit <- lmFit(y, design)

contrast.matrix <- makeContrasts(CAB-DMSO, LEN-CAB, LEN-DMSO, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
vennDiagram(results)
