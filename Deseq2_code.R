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

##Cargar archivo de htseq_counts que esten en nuestro WD
sampleFiles <- data.frame(
  files = list.files(pattern = "_htseq_counts", full.names = TRUE),
  stringsAsFactors = FALSE
)
sampleFiles

##Selección de fenotipos para las muestras
colData <- data.frame(condition=factor(c("A459","A459","HepG2","HepG2","K562","K562")))

sampleCondition <- colData

##Tabla de archivos y fenotipo                      
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

##Crear objeto clase Deseq con nuestros datos de Htseq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,design= ~condition)

ddsHTSeq$condition

##DGEA
dds <- DESeq(ddsHTSeq)
##Extrae tabla de resultados del objeto Deseq
res <- results(dds)
res
dds <- estimateSizeFactors(dds)
resultsNames(dds)
#normalizedcounts <- counts(dds, normalized=TRUE)
##Resultados por comparativa y log2FC
res_A459_vs_K562 <- results(dds, contrast = c("condition", "A459", "K562"))
res_A459_vs_K562_LFC <- lfcShrink(dds, contrast=c("condition", "A459", "K562"), type = "ashr", parallel=TRUE)
res_A459_vs_HepG2 <- results(dds, contrast = c("condition", "A459", "HepG2"))
res_A459_vs_HepG2_LFC <- lfcShrink(dds, contrast=c("condition", "A459", "HepG2"), type = "ashr", parallel=TRUE)
res_K562_vs_HepG2 <- results(dds, contrast = c("condition", "K562", "HepG2"))
res_K562_vs_HepG2_LFC <- lfcShrink(dds, contrast=c("condition", "K562", "HepG2"), type = "ashr", parallel=TRUE)
head(res_A459_vs_HepG2)
head(res_K562_vs_HepG2)
##Ordenador por log2FC
res_K562_vs_HepG2_LFCO <- res_K562_vs_HepG2_LFC[order(res_K562_vs_HepG2_LFC$pvalue),]
res_A459_vs_HepG2_LFCO <- res_A459_vs_HepG2_LFC[order(res_A459_vs_HepG2_LFC$pvalue),]
summary(res)

plotMA(res_K562_vs_HepG2, ylim=c(-2,2))

plotMA(res_K562_vs_HepG2_LFC, ylim=c(-2,2))

plotCounts(dds,which.min(res$padj), intgroup = "condition")

vsd <- vst(dds, blind=FALSE)

#rld <- rlog(dds, blind=FALSE)

#ntd <- normTransform(dds)

#meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

#meanSdPlot(assay(rld))

plotPCA(vsd, intgroup=c("condition"),ntop = 10000)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

df <- as.data.frame(colData(dds)["condition"]) #can also include [,c("condition","type)] for multi parameter comparison

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE, annotation_col=df)

library(EnhancedVolcano)
library(magrittr)

EnhancedVolcano(res_A459_vs_HepG2,
                lab = rownames(res_A459_vs_HepG2),
                x = 'log2FoldChange',
                y = 'pvalue')



