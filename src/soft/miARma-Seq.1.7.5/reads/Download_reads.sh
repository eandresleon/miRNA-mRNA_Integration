#!/bin/bash

checkcurl() {
    if hash curl 2>/dev/null; then
        echo "[`date`] curl: program found."
    else
        echo "[`date`] ERROR :: curl: program not found. Please install curl before continue."
		exit
    fi
}

checksra() {
    if hash fastq-dump 2>/dev/null; then
        echo "[`date`] fastq-dump: program found."
    else
        echo "[`date`] ERROR :: fastq-dump: program not found. Please install the sra-toolbox before continue."
		exit
    fi
}

checkcurl
checksra

echo "[`date`] Downloading SRA files (more info at sra_download.log)"
curl --retry 3 -K fastq_urls.txt >>sra_download.log 2>&1
echo "[`date`] All files downloaded, descompressing (more info at sra_descompress.log)"
for file in `ls *.sra`; do echo "[`date`] descompressing $file"; fastq-dump --split-3 --gzip $file >>sra_descompress.log 2>&1 ; done
echo "[`date`] Reordering files"
mkdir -p sra miRNA_reads mRNA_reads
mv *.sra sra

echo "[`date`] Renaming and sorting miRNA files"
mv SRR3624074.fastq.gz miRNA_reads/CC-CR_rep3.fastq.gz
mv SRR3624073.fastq.gz miRNA_reads/CC-CR_rep2.fastq.gz
mv SRR3624072.fastq.gz miRNA_reads/CC-CR_rep1.fastq.gz
mv SRR3624071.fastq.gz miRNA_reads/CC_rep3.fastq.gz
mv SRR3624070.fastq.gz miRNA_reads/CC_rep2.fastq.gz
mv SRR3624069.fastq.gz miRNA_reads/CC_rep1.fastq.gz

echo "[`date`] Renaming and sorting mRNA files"
cat SRR3624065_1.fastq.gz SRR3624066_1.fastq.gz SRR3624067_1.fastq.gz SRR3624068_1.fastq.gz > mRNA_reads/CC-CR_rep3_1.fastq.gz && rm SRR362406[5-8]_1.fastq.gz
cat SRR3624065_2.fastq.gz SRR3624066_2.fastq.gz SRR3624067_2.fastq.gz SRR3624068_2.fastq.gz > mRNA_reads/CC-CR_rep3_2.fastq.gz && rm SRR362406[5-8]_2.fastq.gz
cat SRR3624061_1.fastq.gz SRR3624062_1.fastq.gz SRR3624063_1.fastq.gz SRR3624064_1.fastq.gz > mRNA_reads/CC-CR_rep2_1.fastq.gz && rm SRR362406[1-4]_1.fastq.gz
cat SRR3624061_2.fastq.gz SRR3624062_2.fastq.gz SRR3624063_2.fastq.gz SRR3624064_2.fastq.gz > mRNA_reads/CC-CR_rep2_2.fastq.gz && rm SRR362406[1-4]_2.fastq.gz
cat SRR3624057_1.fastq.gz SRR3624058_1.fastq.gz SRR3624059_1.fastq.gz SRR3624060_1.fastq.gz > mRNA_reads/CC-CR_rep1_1.fastq.gz && rm SRR362405[7-9]_1.fastq.gz && rm SRR3624060_1.fastq.gz
cat SRR3624057_2.fastq.gz SRR3624058_2.fastq.gz SRR3624059_2.fastq.gz SRR3624060_2.fastq.gz > mRNA_reads/CC-CR_rep1_2.fastq.gz && rm SRR362405[7-9]_2.fastq.gz && rm SRR3624060_2.fastq.gz
cat SRR3624053_1.fastq.gz SRR3624054_1.fastq.gz SRR3624055_1.fastq.gz SRR3624056_1.fastq.gz > mRNA_reads/CC_rep3_1.fastq.gz && rm SRR362405[3-6]_1.fastq.gz
cat SRR3624053_2.fastq.gz SRR3624054_2.fastq.gz SRR3624055_2.fastq.gz SRR3624056_2.fastq.gz > mRNA_reads/CC_rep3_2.fastq.gz && rm SRR362405[3-6]_2.fastq.gz
cat SRR3624049_1.fastq.gz SRR3624050_1.fastq.gz SRR3624051_1.fastq.gz SRR3624052_1.fastq.gz > mRNA_reads/CC_rep2_1.fastq.gz && rm SRR362405[0-2]_1.fastq.gz && rm SRR3624049_1.fastq.gz
cat SRR3624049_2.fastq.gz SRR3624050_2.fastq.gz SRR3624051_2.fastq.gz SRR3624052_2.fastq.gz > mRNA_reads/CC_rep2_2.fastq.gz && rm SRR362405[0-2]_2.fastq.gz && rm SRR3624049_2.fastq.gz
cat SRR3624045_1.fastq.gz SRR3624046_1.fastq.gz SRR3624047_1.fastq.gz SRR3624048_1.fastq.gz > mRNA_reads/CC_rep1_1.fastq.gz && rm SRR362404[5-8]_1.fastq.gz
cat SRR3624045_2.fastq.gz SRR3624046_2.fastq.gz SRR3624047_2.fastq.gz SRR3624048_2.fastq.gz > mRNA_reads/CC_rep1_2.fastq.gz && rm SRR362404[5-8]_2.fastq.gz

echo "[`date`] Done!"
