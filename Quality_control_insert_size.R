BiocManager::install("Rsubread")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
BiocManager::install("soGGi")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicAlignments")
library(Rsubread)
library(Rsamtools)
library(GenomicAlignments)
library(soGGi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)
library(ggplot2)
library(magrittr)
library (dplyr)
setwd("/home/scclab/Atacseq/Hands_on")

sortedBAM = file.path("HepG2_sort_dedup.bam")

#Número de lecturas mapeadas
pmapped <- propmapped(sortedBAM)
pmapped
#Revisión de lecturas mapeadas
chromosomes_of_interest <- paste0("", c(1:22, "X", "Y"))
chromosome_lengths <- c(
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
  58617616, 64444167, 46709983, 50818468, 156040895, 57227415)

names(chromosome_lengths) <- chromosomes_of_interest

all_chromosomes <- GRanges(
  seqnames = Rle(chromosomes_of_interest),
  ranges = IRanges(start = 1, end = chromosome_lengths)
)

atacReads <- readGAlignmentPairs(
  sortedBAM, 
  param = ScanBamParam(
    mapqFilter = 1,
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE),
    what = c("qname", "mapq", "isize"),
    which = all_chromosomes
  )
)

# length(atacReads)
atacReads

#Recuperar tamaño de inserto

atacReads_read1 <- GenomicAlignments::first(atacReads)
insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
head(insertSizes)

#Gráficos de tamaño de inserto
fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                            Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                     Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
  geom_line()

fragLenPlot + theme_bw()

#Gráficos de insertos de zonas abiertas, mono y dinucleosomales

fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315,
                                        437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()



#Graficar señal de Atac-seq en los sitios de comienzo de transcripción (TSS)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), fix = "start", 1)
TSSs
seqlevelsStyle(TSSs) <- "NCBI" 

#Gráficos de la señal de atac-seq en los TSS (creando perfiles de nucleosomas)
# Nucleosome free
nucFree <- regionPlot(bamFile = sortedBAM, testRanges = TSSs, style = "point", 
                      format = "bam", paired = TRUE)
plotRegion(nucFree)
