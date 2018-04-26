<div class="align-justify">
# miRNA-mRNA_Integration
A detailed guide to analise and integrate small-RNASeq and RNASeq samples using [miARma-Seq](https://miarmaseq.com "miARma's Homepage") focus on differentially expresed genes and miRNAs + miRNA-mRNA interactions and a statistical correlation.


----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [miARma Installation](#installation)
   * [Indexes and other needed files](#indexes_and_other_needed_files)
       * [Genome Indexes](#genome_indexes)
       * [Organism Anottation](#organism_anottation)
   * [A complete Analysis](#analysis)
       * [Summary of experiment](#summary)
       * [Data retrieval](#sra)
       * [Analysis](#study)
         * [mRNA Analysis](#mRNA)
         * [miRNA Analysis](#miRNA)
         * [Integration](#integration)
       * [Results](#results)
   * [Conclusions](#conclusions)

----------------------------
# Overview
----------------------------

miARma-Seq, miRNA-Seq And RNA-Seq Multiprocess Analysis tool, is a comprehensive pipeline analysis suite designed for NGS transcriptomic analysis:
- Perform a quality analysis from fastq files
- Remove (and predict) adapter sequences over-represented in the samples
- Aligns reads using different aligners: HISAT2, STAR, BWA, TopHat, Bowtie1 or Bowtie2
- Identify mRNA, miRNA and circRNAs. Also De Novo miRNAs is offered
- Differential expresion analysis among any type of conditions (edgeR and NOISeq)
- Functional enrichment for RNASeq (Gene Ontology and KEGG)
- miRNA-mRNA target prediction for RNASeq and/or small RNASeq
- Data integration for miRNA and mRNA based on statistical correlation (Pearson or Spearman)

You can get more information about miARma at its [webpage](https://miarmaseq.com "miARma's Homepage") or you can read the [miARma-Seq article](https://www.nature.com/articles/srep25749).

<p align="center">
  <img src="https://github.com/eandresleon/miRNA-mRNA_Integration/blob/master/src/images/Figure1.png">
</p>


**Fig 1.** miARma-Seq is presented as a stand-alone tool that provides different well-established softwares at ease of installation process. Our suite can analyse a large number of samples due to its multithread design. Here we show that the analyses of miRNA, mRNA and circRNAs against validated datasets can be easily accessible to research community.

----------------------------
# Installation
----------------------------

miARma has been mainly developed in Perl and R. Our tool has been designed to reduce prerequisites to minimun so, in order to run a complete analysis, you just need:

- If you are using an OSX computer, please install Xcode.
- Java JDK runtime (javac)
- C compiler (gcc) [Included in Xcode]
- [R](https://cran.r-project.org)
- [Bioconductor](https://www.bioconductor.org/install/)
- miARma-Seq. You can use the version included in this [repository](https://github.com/eandresleon/miRNA-mRNA_Integration/raw/master/src/soft/miARma-Seq.1.7.5.tar.gz) or find the latest version at the [installation webpage](http://miarmaseq.idoproteins.com/installation). There you'll find docker containers, virtualbox images and source code. In this guide we will use the source code as it only implies to uncompress a file.

So, to install miARma-Seq, please download the [source code](https://github.com/eandresleon/miRNA-mRNA_Integration/raw/master/src/soft/miARma-Seq.1.7.5.tar.gz) and uncompress it:

```
mkdir ~/bin/
mv miARma-Seq.1.7.5.tar.gz ~/bin/
cd ~/bin/
tar -xzf miARma-Seq.1.7.5.tar.gz
```
In that way miARma will be installed inside the folder bin in your home directory.

----------------------------
# Indexes and other needed files
----------------------------

In order to analyse transcriptomic data, we need to download the following information.
## Genome Indexes ##
Sequenced reads must be aligned (placed) into the reference genome. As we will study human samples, miARma-Seq will need to download human genome data. Besides we will use the [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) aligner, so human HISAT2 indexes will be downloaded following these steps:

```
cd ~/bin/miARma-Seq.1.7.5
curl -L -O https://sourceforge.net/projects/miarma/files/Genomes/Index_hisat2_hg19.tar.bz2
tar -xjf Index_hisat2_hg19.tar.bz2
```
Once uncompressed, a Genome folder will be created including the whole human genome in fasta format and its HISAT indexes.



## Organism anottation ##

When sequences are localized in the reference genome, these must be quantified. The process of counting reads is called read summarization, and to perform this task miARma-Seq includes [FeatureCounts](http://bioinf.wehi.edu.au/featureCounts/). As a result of this part, we obtain the total number of sequences (proportional to the expression levels of the RNAs to be studied). This lets us to find out which genes (coding or not) have been transcribed and expressed in the form of RNA from the genomic DNA.
To perfom this step, miARma-Seq needs to know the position of each genomic feature (exon, transcript, gene) in the genome. To do that we will use the GENCODE v26 anottation file. For miRNAs we will use the microRNA annotation file from miRBase version 20.
These files can be obtained from this git-hub repository:

```
cd ~/bin/miARma-Seq.1.7.5
mkdir data
cd data

curl -L -O https://github.com/eandresleon/miRNA-mRNA_Integration/raw/master/src/data/gencode.v26_GRCh37.annotation.gtf.gz
gunzip gencode.v26_GRCh37.annotation.gtf.gz

curl -L -O https://raw.githubusercontent.com/eandresleon/miRNA-mRNA_Integration/master/src/data/miRBase_Annotation_20_for_hsa_mature_miRNA.gtf
```

## Analysis

At this moment, we have already downloaded all needed software to perform a complete anaylysis. Besides we have obtained the human genome, its HISAT indexes and all human gene and miRNA anottation.
So we are ready to perfom the complete study.

## Summary

Briefly, the work by [Lu et al](https://www.nature.com/articles/nm.4424)  is focused on colorectal cancer (CRC), as it remains the leading cause of cancer-related death worldwide. Cetuximab and panitumumab are typical CRC treatments that bind the extracellular domain of the EGF receptor enhancing its internalization and degradation. When these are combined with chemotherapy, up to 72% response rates are reported. However, de novo and acquired drug resistance frequently arises, and little is known about non-genetic resistance mechanisms.
In order to increase our knowledge about these mechanisms, CRC samples from colon adenocarcinoma cancer cell line HCA-7 were treated with cetuximab for approximately four months to induce resistance and deposited at the [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) repository with accession number: [GSE82236](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82236):
- mRNA samples. Three cetuximab-resistant and three non-resistant RNASeq colon cancer samples were downloaded from GEO and analysed using miARma-Seq. Sequenced reads were inspected using FastQC to identify sequencing errors such as bad quality reads or possible accumulation of adapter sequences. Subsequently, sequences were aligned using Hisat2 and those with a quality score above 25 were counted and summarized as expression values for each gene using the GENCODE version 26 annotation file. To obtain those genes with a different expression pattern amongst both types of samples, we use edgeR as it enables the calculation of differential expression values. Low expressed genes (those having a count per million (CPM) value smaller than 1) were removed for the analysis as recommended by Anders et al [64]. Finally, as it was done by Lu et al, genes having a |log2FC|≥1 and a FDR ≤0.01 were considered as differentially expressed. See Figure 2.
- miRNA analysis. Correspondingly, three CRC-CR and three CRC miRNASeq colon cancer samples were obtained and studied. Single-end sequenced reads were inspected using FastQC. Once this process, reads were pre-processed using Minion to predict (as it was not facilitated by the authors) and to remove the accumulated adapter sequences detected in the previous step. Afterwards, trimmed reads with a size comprised between 18 and 35 nucleotides were aligned using Bowtie1 and quantified using FeatureCounts and the microRNA annotation file from miRBase version 20. Low expressed miRNAs (having a CPM value lower than 1) were removed and differentially expressed microRNAs were acquired using edgeR. As above, miRNAs having a |log2FC|≥1 and a FDR ≤0.01 were considered as differentially expressed. See Figure 3 for more detail.


## sra

For transcripts, the *generateEvents* operation outputs one single file: An *ioi* file that shows the transcript "events" in each gene. This is a tab separated file with the following fields

## study

For transcripts, the *generateEvents* operation outputs one single file: An *ioi* file that shows the transcript "events" in each gene. This is a tab separated file with the following fields
### **mRNA**
-------

1. **seqname**: field 1 from the input GTF file of the generateEvents operation (generally the chromosome name)

2. **gene_id**: ID of the gene where the event appears taken from the GTF file. (comma separated list in case of using the --pool-genes option)

3. **event_id**: ID of the event, formatted as **gene_id**;**transcript_id**

4. **transcript_id**: ID of the transcript that defines the event, for which the relative inclusion (PSI) is calculated.

5. **Total transcripts**: IDs of the all transcripts in the gene (including the transcript in 4.) 

Note that we call a transcript "event", a transcript in a gene. 


An example of an *ioi* file is the following one:
```
seqname	gene_id	event_id	inclusion_transcripts	total_transcripts
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000556772.1    ENST00000556772.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000554818.1    ENST00000554818.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000556600.1    ENST00000556600.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000557774.1    ENST00000557774.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
```



For local AS events, the *generateEvents* operation outputs two files:

1. An *ioe* file that shows the relationship between each event and the transcripts that define that particular event.

2. A *GTF* file (with a track header) to load into the [UCSC genome browser](http://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=auto&source=genome.ucsc.edu) to visualize the different local AS events.

The name of the output is generated as follows:

**&lt;output-file&gt;**.**&lt;event-type&gt;**.ioe/gtf
    
**&lt;output-file&gt;**: -o option specified when launching the program.

**&lt;event_type&gt;**: a two letter code referring to the event type (SE, A3, A5, MX, RI, AF, AL).


The ioe file has the following fields:

### **miRNA**
-------


1. **seqname**: field 1 from the input GTF file of the generateEvents operation (generally the chromosome name)

2. **gene_id**: ID of the gene where the event appears. (comma separated list in case of using the --pool-genes option)

3. **event_id**: ID of the event.

4. **Inclusion transcripts**: IDs of the transcripts that define the form of the event for which we calculate the PSI (e.g. exon inclusion) and contribute to the numerator of the PSI formula. For more details see Figures above. 

5. **Total transcripts**: IDs of the all transcripts that define either of the two forms of the event (e.g. inclusion and skipping) and contribute to the denominator of the PSI formula. For more details see Figures above.


**Event ID:** The event_id is formatted as follows:
```
<gene_id>;<event-type>:<seqname>:<coordinates-of-the-event>:<strand>
```
where:

- &lt;gene_id&gt;: is the gene where the event takes place
- &lt;event-type&gt;: correspond to the two letter code of the event from the following list.
  - **SE**: Skipping Exon
  - **A5**: Alternative 5' Splice Site
  - **A3**: Alternative 3' Splice Site
  - **MX**: Mutually Exclusive Exon
  - **RI**: Retained Intron
  - **AF**: Alternative First Exon
  - **AL**: Alternative Last Exon  
- &lt;seqname&gt;: coordinate reference system (e.g. chr1) 
- &lt;coordinates-of-the-event&gt;: the coordinates of the event depends on the type of event (see above)
- &lt;strand&gt;: either '+' or '-'


**Inclusion transcripts**: The IDs must be the same as those provided in the expression (TPM) files. These transcripts define the inclusion form of the event. This form is chosen according to the convention described above in Figures 3 and 4: 

- **SE**: transcripts including the middle exon
- **A5/A3**: transcripts minimizing the intron length
- **MX**: transcripts containing the alternative exon with the smallest (left most) start coordinate.
- **RI**: transcripts that have the retain intron
- **AF/AL**: transcripts maximizing the intron length

**Total transcripts**: The IDs must be the same as those provided in the expression (TPM) files. Total transcripts are the transcripts that define both forms of the event and therefore contribute to the denominator of the PSI formula. 

An example of an *ioe* file is the following one:
```
seqname	gene_id	event_id	inclusion_transcripts	total_transcripts
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73763993-73789838:73789937-73822334:-	ENST00000556112.1    ENST00000556112.1,ENST00000555859.1
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73876776-73925740:73926011-73929235:-	ENST00000553415.1    ENST00000553415.1,ENST00000556700.1
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73744001-73745989:73746132-73749067:-	ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000535282.1   ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000535282.1,ENST00000554546.1,ENST00000356296.4,ENST00000555394.1,ENST00000454166.4,ENST00000560335.1,ENST00000555738.2,ENST00000554014.2
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73749213-73750789:73751082-73753818:-	ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000555394.1,ENST00000535282.1,ENST00000553997.1,ENST00000557486.1    ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000555394.1,ENST00000535282.1,ENST00000553997.1,ENST00000557486.1,ENST00000454166.4,ENST00000560335.1,ENST00000555738.2

```


**Important: --pool-genes**
-------


This option is important when creating ioe/ioi from annotations that are not loci-based, e.g.: RefSeq and UCSC genes.
Unlike Ensembl or Gencode, which annotate gene loci, i.e. a set of transcripts will be uniquely be related to a gene at a locus, other annotations, like UCSC and Refseq 
dowloaded from UCSC, do not have this unequivocal link of transcripts to a genomic locus. 

We thus re-cluster transcripts according to their overlap in genomic extent on the same strand, and according to them sharing 
at least one splice-site position in the genome. This method borrows from the (pre-Gencode) Ensembl function to create
genes from transcripts (see Curwen et al. 2004 https://www.ncbi.nlm.nih.gov/pubmed/15123590). If transcripts are not appropriately clustered, 
transcript relative abundances or event inclusion levels may not make sense.

Using the **--pool-genes** option is also advisable to use with Ensembl and Gencode. The annotation contains genes with overlapping transcripts
that share a great deal of sequence, hence their relative contribution to alternative splicing events should be taken into account. 


**GTF for local events**
-------

This output file is aimed for visualization and it is a regular *gtf* in half-open coordinate convention with a track header to be directly uploaded in UCSC. In this file, each of the two possible forms of an event is described as a separated *transcript_id*. Moreover, the *gene_id* is defined as the event_id. Note that SUPPA generally operates with closed coordinates.

The **-l** | **--exon-length** option of the *eventGenerator* operation defines the length used to display the constitutive exons in the GTF file and it is used only for visualization purposes. For instance, in a *Skipping exon* event the flanking exons will be shown with a length given by this parameter.


-------------------
**PSI calculation for Transcripts and Events**
==============
-------------------

### **Integration**



SUPPA reads the ioi or ioe file generated in the previous step and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each transcript or each local event. 

## Results

An **ioi/ioe** file and a "**transcript expression file**" are required as input. 

The transcript expression file is a tab separated file where each line provides the estimated abundance of each transcript (in TPM units). This file might contain multiple columns with the expression values in different samples. The expression file must have a header with the naming of the different expression fields, i.e., the sample name of each expression value. 

An example of a transcript expression file for one single sample:
```
sample1
transcript1 <expression>
transcript1 <expression>
transcript1 <expression>
```
A transcript expression file with multiple samples:

```
sample1 sample2 sample3 sample4
transcript1 <expression>  <expression>  <expression>  <expression>
transcript2 <expression>  <expression>  <expression>  <expression>
transcript3 <expression>  <expression>  <expression>  <expression>
```

**Note:** these files have a header with only the sample names (1 less column)

# Conclusions
</div>