### Workshop ATACseq ####

install.packages("BiocManager")
BiocManager::install("DESeq2",force = TRUE)
BiocManager::install("apeglm",force = TRUE)
BiocManager::install("vsn",force = TRUE)
BiocManager::install("hexbin",force = TRUE)
BiocManager::install("ashr",force = TRUE)
BiocManager::install("EnhancedVolcano",force = TRUE)
install.packages("magrittr")

#BiocManager::install("pheatmap",force = TRUE)

library("apeglm")
library("dplyr")
library("DESeq2")
library("ggplot2")
library("vsn")
library("ashr")
library("pheatmap")

##Set working directory al directorio donde esten los archivos del github o por el pipeline
# Modificar con el wd que tengas 
setwd("/home/scclab/Atacseq/Hands_on")
directory = getwd()

sampleFiles <- data.frame(
  files = list.files(pattern = "_htseq_counts", full.names = TRUE),
  stringsAsFactors = FALSE
)
sampleFiles
colData <- data.frame(condition=factor(c("A459","A459","HepG2","HepG2","K562","K562")))

sampleCondition <- colData
                      
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design= ~condition)

ddsHTSeq$condition

dds <- DESeq(dds)
res <- results(dds)
res
dds <- estimateSizeFactors(dds)
normalizedcounts <- counts(dds, normalized=TRUE)
res_A459_vs_HepG2 <- results(dds, contrast = c("condition", "A459", "HepG2"))
res_A459_vs_HepG2_LFC <- lfcShrink(dds, contrast=c("condition", "A459", "HepG2"), type = "ashr", parallel=TRUE)
res_K562_vs_HepG2 <- results(dds, contrast = c("condition", "K562", "HepG2"))
res_K562_vs_HepG2_LFC <- lfcShrink(dds, contrast=c("condition", "K562", "HepG2"), type = "ashr", parallel=TRUE)
head(res_A459_vs_HepG2)
head(res_K562_vs_HepG2)

resultsNames(dds)

res_K562_vs_HepG2_LFCO <- res_K562_vs_HepG2_LFC[order(res_K562_vs_HepG2_LFC$pvalue),]
res_A459_vs_HepG2_LFCO <- res_A459_vs_HepG2_LFC[order(res_A459_vs_HepG2_LFC$pvalue),]
summary(res)

plotMA(res, ylim=c(-2,2))

plotMA(res_A459_vs_HepG2_LFCO, ylim=c(-2,2))

plotCounts(dds,which.min(res$padj), intgroup = "condition")

vsd <- vst(dds, blind=FALSE)

rld <- rlog(dds, blind=FALSE)

ntd <- normTransform(dds)

meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

plotPCA(vsd, intgroup=c("condition"),ntop = 10000)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)["condition"]) #can also include [,c("condition","type)] for multi parameter comparison

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

library(EnhancedVolcano)
library(magrittr)

EnhancedVolcano(res_A459_vs_HepG2,
                lab = rownames(res_A459_vs_HepG2),
                x = 'log2FoldChange',
                y = 'pvalue')



