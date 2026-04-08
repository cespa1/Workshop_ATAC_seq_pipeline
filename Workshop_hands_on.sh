#!/bin/bash
mkdir -p peak_calling
mkdir -p fastqc
## Reducción de tamaño de los fastq para poder hacer los analisis durante el workshop 
#zcat ENCFF047YJC.fastq.gz | seqtk sample -s100 - 500000 > R1_subset_A549.fastq
#zcat ENCFF048FVL.fastq.gz | seqtk sample -s100 - 500000 > R2_subset_A549.fastq
#zcat ENCFF261UXP.fastq.gz | seqtk sample -s100 - 500000 > R1_subset_K562.fastq
#zcat ENCFF123TMX.fastq.gz | seqtk sample -s100 - 500000 > R2_subset_K562.fastq
#zcat ENCFF875IIN.fastq.gz | seqtk sample -s100 - 500000 > R1_subset_HepG2.fastq
#zcat ENCFF461JXR.fastq.gz | seqtk sample -s100 - 500000 > R2_subset_HepG2.fastq

## Analisis de calidad de los fastq pre-trimming
for muestra in $(ls | grep subset | cut -d '_' -f3 | uniq )
do
fastqc "R1_subset"*$muestra -o fastqc/ --memory 10000 --threads 20
fastqc "R2_subset"*$muestra -o fastqc/ --memory 10000 --threads 20
done

## Trimming de los adaptadores
for muestra in $(ls | grep subset | cut -d '_' -f3 | uniq )
do 
fastp -i "R1_subset"*$muestra -o "R1_subset_trim_"$muestra \
    -I "R2_subset"*$muestra -O "R2_subset_trim_"$muestra \
    --detect_adapter_for_pe --thread 16 2> "fastp_log_"$muestra".txt"
done

## Analisis de calidad de los fastq post-trimming

fastqc R1_subset_trim_* -o fastqc/ --memory 10000 --threads 20
fastqc R2_subset_trim_* -o fastqc/ --memory 10000 --threads 20

##Multiqc
multiqc .

## Creación del index para el genoma al que vamos a alinear la muestra
#bowtie2-build Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz human_index

## Alineamiento de las muestras usando bowtie2
for muestra in $(ls | grep R1_subset_trim_ | grep ".fastq\b" | cut -d '_' -f4 | cut -d '.' -f1 | uniq )
do 
bowtie2 -x GRCh38_noalt_as -1 "R1_subset_trim"*$muestra*".fastq.gz" -2 "R2_subset_trim"*$muestra*".fastq.gz" -S $muestra".sam" -X 500 -p 23 --very-sensitive 2> "bowtie_log_"$muestra".txt"

##Ordenamiento del SAM
samtools sort $muestra".sam" -o $muestra"_sort.sam"
samtools view -h $muestra"_sort.sam" | grep -E "^@|$(paste -sd'|' chroms.txt)" > $muestra"_sort_filt.sam"


##Marcar los duplicados de PCR de las muestras 
picard MarkDuplicates I=$muestra"_sort_filt.sam" O=$muestra"_sort_dedup.sam" M=picard_log.txt REMOVE_DUPLICATES=false

##Transformar de SAM a BAM y hacer el index del bam 
samtools view -bhS $muestra"_sort_dedup.sam" > $muestra"_sort_dedup.bam"
samtools index $muestra"_sort_dedup.bam" $muestra"_sort_dedup.bai"

## Visualización de la muestras en IGV
igvtools count -z 5 -w 25 -e 250 $muestra"_sort_dedup.bam" $muestra"_dedup_sorted.tdf" Homo_sapiens.GRCh38.dna.primary_assembly.fa

##Peak calling 
macs2 callpeak -f BAMPE -g hs --keep-dup all --cutoff-analysis -n $muestra \
  -t $muestra"_sort_dedup.bam" --broad --outdir peak_calling/ 2> "macs2_"$muestra".log"

##Htseq count para poder hacer algunos analisis downstream como deseq2
touch peak_calling/$muestra"_peaks.bed"
touch peak_calling/merge_peaks_all.bed
cut -f1-4 peak_calling/$muestra"_peaks.broadPeak" > peak_calling/$muestra"_peaks.bed"
cat peak_calling/$muestra"_peaks.bed" >> peak_calling/merge_peaks_all.bed
done

    ##Ordenamiento de los picos 
sort -k1,1 -k2,2n peak_calling/merge_peaks_all.bed > peak_calling/merge_sorted_all.bed
    ##Merge de los picos que estan cerca entre ellos
bedtools merge -i peak_calling/merge_sorted_all.bed > peak_calling/merge_final_all.bed   
    ##Comando awk necesario para generar el archivo GTF para htseq  
awk '{print $1"\tunknown\texon\t"$2"\t"$3"\t.\t+\t0\tgene_id \42"$1":"$2"-"$3"\42;\ttranscript_id \42"$1":"$2"-"$3"\42;\tgene_name \42"$1":"$2"-"$3"\42;"}' peak_calling/merge_final_all.bed    > peak_calling/merge_final.gtf
    ##Obtener el número de counts por pico del archivo BED
for muestra in $(ls | grep R1_subset_trim_ | grep ".fastq\b" | cut -d '_' -f4 | cut -d '.' -f1 | uniq )
do 
htseq-count --format=bam --type=exon --idattr=gene_name $muestra"_sort_dedup.bam" peak_calling/merge_final.gtf  > $muestra"_htseq_counts"
done

###Analisis downstream Heatmap de accesibilidad con deeptools
#for muestra in $(ls | grep R1_subset_trim_ | grep ".fastq\b" | cut -d '_' -f4 | cut -d '.' -f1 | uniq )
#    do
#        bamCoverage -p 23 -b $muestra"_sort_dedup.bam" -o bigwig/$muestra".bw" 
#    done

#computeMatrix reference-point --referencePoint center -p 23 -S bigwig/A549.bw bigwig/HepG2.bw bigwig/K562.bw -R peak_calling/merge_final_all.bed -b 3000 -a 3000 -bs=1 -p=max --outFileName Hands_on_matrix.matrix.gz

#plotHeatmap -m Hands_on_matrix.matrix.gz --outFileName Hands_on_heatmap.pdf --yMax 10 --colorMap=Purples --zMax 10

#plotHeatmap -m Hands_on_matrix.matrix.gz  -out Hands_on_heatmap_kmeans4.pdf --kmeans 4 --outFileSortedRegions ac_kmeans10.bed --zMin 0 --zMax 10 --yMin 0 --yMax 10

##Analisis downstream HOMER y predicción de factores de transcripción

#bedtools getfasta -fi /home/scclab/Atacseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed peak_calling/merge_final_all.bed > HOMER/merge_final_all.fa

#findMotifs.pl HOMER/merge_final_all.fa human HOMER/ -fasta ~/Atacseq/scrambled.fasta

###Analisis downstream TOBIAS y footprinting

#TOBIAS ATACorrect --bam $muestra".bam" --genome Homo_sapiens.GRCh38.dna.primary_assembly.fa --peaks peak_calling/merge_final_all.bed --blacklist ../hg38-blacklist.v2.bed --outdir TOBIAS/$muestra --cores 22 
    
#TOBIAS FootprintScores --signal TOBIAS/$muestra/$muestra*"_corrected.bw" --regions peak_calling/merge_final_A431.bed \
#--output TOBIAS/$muestra/$muestra"_footprints.bw" --cores 22
