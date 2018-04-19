# miRNA-mRNA_Integration
A detailed guide to analise and integrate small-RNASeq and RNASeq samples using [miARma-Seq](https://miarmaseq.com "miARma's Homepage") focus on differentially expresed genes and miRNAs + miRNA-mRNA interactions and a statistical correlation.


----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [miARma Installation](#installation)
   * [Indexes and other needed files](#files)
   * [A complete Analysis](#analisis)
       * [Summary of experiment](#summary)
       * [Data retrieval](#sra)
       * [Analysis](#study)
         * [mRNA Analysis](#mRNA)
         * [miRNA Analysis](#miRNA)
         * [Integration](#Integration)
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
  <img src="https://github.com/eandresleon/miRNA-mRNA_Integration/blob/master/src/images/Figure1.tiff">
</p>
 
**Fig 1.** miARma-Seq is presented as a stand-alone tool that provides different well-established softwares at ease of installation process. Our suite can analyse a large number of samples due to its multithread design. Here we show that the analyses of miRNA, mRNA and circRNAs against validated datasets can be easily accessible to research community.

----------------------------
# Installation
----------------------------

miARma has been mainly developed in Perl and R. Our tool has been designed to reduce prerequisites to minimun so, im order to run a complete analysis, you just need:

- If you are using an OSX computer, please install Xcode.
- Java JDK runtime (javac)
- C compiler (gcc) [Included in Xcode]
- [R](https://cran.r-project.org)
- [Bioconductor](https://www.bioconductor.org/install/)
- miARma-Seq. You can use the version included in this [repository](https://github.com/eandresleon/miRNA-mRNA_Integration/raw/master/src/soft/miARma-Seq.1.7.5.tar.gz) or find the latest version at the [intallation webpage](http://miarmaseq.idoproteins.com/installation).

If necessary, to install python3 we recommend to download from the official site https://www.python.org/downloads/ the corresponding version for your OS.

A installation using pip is available using the next command:

```
pip install SUPPA==2.2.1 
```
By default SUPPA is installed into the Python package library directory. The following command can be executed to obtain the directory location:

```
pip show SUPPA 
```

SUPPA is ready to use. Once downloaded, it can be used directly from the command line by specifying the absolute path to the SUPPA executable (suppa.py).

Another option is via bioconda (thanks to Devon Ryan)

```
conda install -c bioconda suppa
```

----------------------------
# files
----------------------------

SUPPA works with a command/subcommand structure:

```
python3.4 suppa.py subcommand options

```
where the subcommand can be one of these five:

- **generateEvents**    : Generates events from an annotation.
- **psiPerEvent**       : Quantifies event inclusion levels (PSIs) from multiple samples.
- **psiPerIsoform**       : Quantifies isoform inclusion levels (PSIs) from multiple samples.
- **diffSplice**        : Calculate differential splicing across multiple conditions with replicates.
- **clusterEvents**     : Cluster events according to PSI values across conditions.

----------------------------
**Generation of transcript events and local alternative splicing events**
==============

----------------------------

SUPPA can work with local alternative splicing events or with transcripts "events" per gene. The local alternative splicing events are
standard local splicing variations (see below), whereas a transcript event is an isoform-centric approach, where each isoform in a gene is described separately. 

SUPPA generates the AS events or transcript events from an input annotation file (GTF format). The method reads transcript and gene information solely from the "exon" lines. It then generates the events and outputs an event file: **.ioe** format for local AS events, and **.ioi** format for transcripts. 


Different local event types generated by SUPPA:

- Skipping Exon (SE)
- Alternative 5'/3' Splice Sites (A5/A3)
- Mutually Exclusive Exons (MX)
- Retained Intron (RI)
- Alternative First/Last Exons (AF/AL)


![Events_legend1.jpg](https://cloud.githubusercontent.com/assets/23315833/22699555/19788560-ed58-11e6-8340-1ce83445bdb3.png)

![Events_legend2.jpg](https://cloud.githubusercontent.com/assets/23315833/22699557/199ac7c4-ed58-11e6-8512-2d3950001a8d.png)


**Fig 3.** The figure describes the nomenclature for events in the forward (upper panel) and reverse (lower panel) strands. Each event is identified uniquely by a set of coordinates: The start (s) and end (e) coordinates for the different exonic regions involved in the event. The external coordinates of the event are only used for the RI, AF and AL events. The form of the alternative splicing event that includes event. 


## Analysis

An annotation file in GTF format is required (see e.g. http://mblab.wustl.edu/GTF22.html):

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 

```

The *generateEvents* operation uses the lines where the feature (column 3) is "exon". It then reads the different transcripts and genes. For that purpose "gene_id" and "transcript_id" tags are required in the attributes field (column 9). For local AS events, this command also generates a GTF file with the calculated events and with a track header ready to be uploaded into the UCSC browser for visualization (see below). 

## Summary

To generate the events from the GTF file one has to run the following command:

```
python3.4 suppa.py generateEvents [options]
```
List of options available:

- **-i**  | **--input-file**: a GTF format file containing at least "exon" lines

- **-o**  | **--output-file**: name of the output file without any extension

- **-f**  | **--format**: [ioe,ioi]
  	  - Required. Format of the event annotation file: ioe for local events, ioi for transcript events.

- **-p**  | **--pool-genes**: 
  	  - Optional. Redefine genes by clustering together transcripts by genomic stranded overlap and sharing at least one exon.
            It is crucial when creating ioe/ioi from annotations that are not loci-based, e.g.: RefSeq and UCSC genes.

- **-e**  | **--event-type**: (only used for local AS events) space separated list of events to generate from the following list:

  - **SE**: Skipping exon (SE)
  - **SS**: Alternative 5' (A5) or 3' (A3) splice sites (generates both)
  - **MX**: Mutually Exclusive (MX) exons
  - **RI**: Retained intron (RI)
  - **FL**: Alternative First (AF) and Last (AL) exons (generates both)

- **-b** | **--boundary**: [S,V]
        - Boundary type (only used for local AS events). Options: S -- Strict (Default) V -- Variable

- **-t** | **--threshold**: THRESHOLD
        - Variability treshold (Default: 10nt. Only used forlocal AS events). In case of strict boundaries this argument is ignored

- **-l**  | **--exon-length**: (only used for local AS events). Defines the number of nucleotides to display in the output GTF. (Default: 100 nt)

- **-h**  | **--help**: display the help message describing the different paramenters


The command line to generate local AS events will be of the form:

```
python3.4 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioe -e <list-of-events>
```

The command to generate the transcript "events" would be of the form:

```
python3.4 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioi 
```

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

### **integration**



SUPPA reads the ioi or ioe file generated in the previous step and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each transcript or each local event. 

## results

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
