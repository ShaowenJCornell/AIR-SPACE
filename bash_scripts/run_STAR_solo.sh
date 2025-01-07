#!/bin/bash

BeadBarcode=${1}
fastq1=${2}
fastq2=${3}

cut -b 1-8 ${BeadBarcode}.txt > ${BeadBarcode}_1.txt
cut -b 9-14 ${BeadBarcode}.txt > ${BeadBarcode}_2.txt



/programs/STAR-2.7.10a/STAR --runThreadN 16 --genomeDir ../References/refdata-STARsolo-mm10-3.0.0/ --readFilesIn ${fastq2} ${fastq1} --readFilesCommand zcat --limitBAMsortRAM 50000000000 --outFileNamePrefix ./STARsolo_output/ --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --clipAdapterType CellRanger4 --soloType CB_UMI_Complex --soloBarcodeReadLength 0 --soloCBwhitelist ${BeadBarcode}_1.txt ${BeadBarcode}_2.txt --soloFeatures Gene GeneFull Velocyto --soloAdapterSequence TCTTCAGCGTTCCCGAGA --soloAdapterMismatchesNmax 2 --soloCBposition 0_0_2_-1 3_1_3_6 --soloUMIposition 3_7_3_14 --soloCellFilter EmptyDrops_CR --soloMultiMappers EM --soloCBmatchWLtype 1MM

gzip ./STARsolo_output/Solo.out/Gene/raw/*
