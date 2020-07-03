# awesome-vdj

[Antigen] presentation and recognition is central to immunology. [HLA genes] encode the proteins that present antigens. [VDJ genes][vdj] encode the receptors: T cell receptors (TCRs) in T cells and the repertoires of antibodies/immunoglobulins in B cells.

Here, researchers can find links to tools and resources for computational analysis of HLA and VDJ data.

[Contributions are welcome!](https://github.com/slowkow/awesome-vdj/blob/master/CONTRIBUTING.md)

[Antigen]: https://en.wikipedia.org/wiki/Antigen
[vdj]: https://en.wikipedia.org/wiki/V(D)J_recombination
[HLA genes]: https://en.wikipedia.org/wiki/Human_leukocyte_antigen

**Table of Contents**

- [VDJ Databases](#vdj-databases)
- [VDJ Analysis](#vdj-analysis)
- [HLA Databases](#hla-databases)
- [HLA Analysis](#hla-analysis)

---

## Literature

- Robson, K. J., Ooi, J. D., Holdsworth, S. R., Rossjohn, J. & Kitching, A. R. [HLA and kidney disease: from associations to mechanisms.](https://pubmed.ncbi.nlm.nih.gov/30206339/) Nat. Rev. Nephrol. 14, 636–655 (2018)

- La Gruta, N. L., Gras, S., Daley, S. R., Thomas, P. G. & Rossjohn, J. [Understanding the drivers of MHC restriction of T cell receptors.](https://pubmed.ncbi.nlm.nih.gov/29636542/) Nat. Rev. Immunol. 18, 467–478 (2018)

- Nemazee, D. [Mechanisms of central tolerance for B cells.](https://www.ncbi.nlm.nih.gov/pubmed/28368006) Nat. Rev. Immunol. 17, 281–294 (2017)

---

## VDJ Databases

### McPAS-TCR: A manually curated catalogue of pathology associated T-cell receptor sequences

http://friedmanlab.weizmann.ac.il/McPAS-TCR/

> McPAS-TCR is a manually curated catalogue of T cell receptor (TCR) sequences that were found in T cells associated with various pathological conditions in humans and in mice. It is meant to link TCR sequences to their antigen target or to the pathology and organ with which they are associated.

### STCRDab: The Structural T-Cell Receptor Database

http://opig.stats.ox.ac.uk/webapps/stcrdab/

> An automated, curated set of T-Cell Receptor structural data from the PDB.

### TCR3d: T cell receptor structural repertoire database

https://tcr3d.ibbr.umd.edu/

> Welcome to the T cell receptor (TCR) structural repertoire database. Here we provide an easy-to-use interface to view all experimentally determined T cell receptor structures and their complexes. This includes complementarity determining region loops and analysis of interfaces with antigenic peptide and MHC.
> 
> We have also assembled a set of known TCR sequences from recent studies including TCR repertoire sequencing efforts.
> 
> The major goal of this site is to enable insights into the basis of TCR structure and recognition, to assist efforts in predictive modeling of this key component of the adaptive immune response, and to facilitate rational engineering of improved and novel immunotherapeutics.

### VDJDB: A curated database of T-cell receptor sequences of known antigen specificity

https://vdjdb.cdr3.net

https://github.com/antigenomics/vdjdb-db

> The primary goal of VDJdb is to facilitate access to existing information on T-cell receptor antigen specificities, i.e. the ability to recognize certain epitopes in certain MHC contexts.
>
> Our mission is to both aggregate the scarce TCR specificity information available so far and to create a curated repository to store such data.

---

## VDJ Analysis

### DeepTCR: Deep Learning Methods for Parsing T-Cell Receptor Sequencing (TCRSeq) Data

Python package

https://github.com/sidhomj/DeepTCR

> DeepTCR is a python package that has a collection of unsupervised and supervised deep learning methods to parse TCRSeq data. It has the added functionality of being able to analyze paired alpha/beta chain inputs as well as also being able to take in v/d/j gene usage and the contextual HLA information the TCR-Sequences were seen in (i.e. HLA alleles for a repertoire from a given human sample).

### dkm: Dynamic Kernel Matching

Python scripts

https://github.com/jostmey/dkm

> DKM is analogous to a convolutional network, but for sequences. Consider the problem of classifying a sequence. Because some sequences are longer than others, the number of features is irregular. Given a specific sequence, the challenge is to determine the appropriate permutation of features with weights, allowing us to run the features through the statistical classifier to generate a prediction. To find the permutation of features that exhibit the maximal response, like how max-pooling identifies the image patch that exhibit the maximal response, we use a sequence alignment algorithm.

### immunarch: An R Package for Painless Bioinformatics Analysis of T-cell and B-cell Immune Repertoire Data

R package

https://github.com/immunomind/immunarch

> immunarch is an R package designed to analyse T-cell receptor (TCR) and B-cell receptor (BCR)
> repertoires, aimed at medical scientists and bioinformaticians. The mission of immunarch is to
> make immune sequencing data analysis as effortless as possible and help you focus on research
> instead of coding. Follow us on Twitter for news and updates.

### immuneSIM: Tunable Simulation of B- And T-Cell Receptor Repertoires

R package

https://CRAN.R-project.org/package=immuneSIM 

https://github.com/GreiffLab/immuneSIM

> Simulate full B-cell and T-cell receptor repertoires using an in silico recombination process
> that includes a wide variety of tunable parameters to introduce noise and biases. Additional
> post-simulation modification functions allow the user to implant motifs or codon biases as
> well as remodeling sequence similarity architecture. The output repertoires contain records
> of all relevant repertoire dimensions and can be analyzed using provided repertoire analysis
> functions. Preprint is available at bioRxiv (Weber et al., 2019 <doi:10.1101/759795>).

### MiGMAP: mapper for full-length T- and B-cell repertoire sequencing

Groovy and Java tools

https://github.com/mikessh/migmap

> In a nutshell, this software is a smart wrapper for IgBlast V-(D)-J mapping tool designed to facilitate analysis immune receptor libraries profiled using high-throughput sequencing. This package includes additional experimental modules for contig assembly, error correction and immunoglobulin lineage tree construction.

### MiXCR: a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data

> MiXCR is a universal framework that processes big immunome data from raw sequences to quantitated clonotypes. MiXCR efficiently handles paired- and single-end reads, considers sequence quality, corrects PCR errors and identifies germline hypermutations. The software supports both partial- and full-length profiling and employs all available RNA or DNA information, including sequences upstream of V and downstream of J gene segments.

https://github.com/milaboratory/mixcr

### PRESTO: The REpertoire Sequencing TOolkit

https://presto.readthedocs.io/en/stable

https://bitbucket.org/kleinstein/presto

> pRESTO is a toolkit for processing raw reads from high-throughput sequencing of B cell and T cell repertoires.

> Dramatic improvements in high-throughput sequencing technologies now enable large-scale characterization of lymphocyte repertoires, defined as the collection of trans-membrane antigen-receptor proteins located on the surface of B cells and T cells. The REpertoire Sequencing TOolkit (pRESTO) is composed of a suite of utilities to handle all stages of sequence processing prior to germline segment assignment. pRESTO is designed to handle either single reads or paired-end reads. It includes features for quality control, primer masking, annotation of reads with sequence embedded barcodes, generation of unique molecular identifier (UMI) consensus sequences, assembly of paired-end reads and identification of duplicate sequences. Numerous options for sequence sorting, sampling and conversion operations are also included.

### scirpy: A scanpy extension to analyse single-cell TCR data.

Python package

https://github.com/icbi-lab/scirpy

> Scirpy is a scalable python-toolkit to analyse T cell receptor (TCR) repertoires from single-cell RNA sequencing (scRNA-seq) data. It seamlessly integrates with the popular scanpy library and provides various modules for data import, analysis and visualization.

### tcr-dist

Python scripts

https://github.com/phbradley/tcr-dist

> Software tools for the analysis of epitope-specific T cell receptor (TCR) repertoires

### VDJtools

Groovy and Java tools

https://github.com/mikessh/vdjtools

> A comprehensive analysis framework for T-cell and B-cell repertoire sequencing data

---

## HLA Databases

### IEDB: Immune Epitope Database and Analysis Resource

https://www.iedb.org/

> The Immune Epitope Database (IEDB) is a freely available resource funded by NIAID. It catalogs experimental data on antibody and T cell epitopes studied in humans, non-human primates, and other animal species in the context of infectious disease, allergy, autoimmunity and transplantation. The IEDB also hosts tools to assist in the prediction and analysis of epitopes.

### IMGTHLA

https://www.ebi.ac.uk/ipd/imgt/hla/

https://github.com/ANHIG/IMGTHLA

> The IPD-IMGT/HLA Database provides a specialist database for sequences of the human major histocompatibility complex (MHC) and includes the official sequences named by the WHO Nomenclature Committee For Factors of the HLA System. The IPD-IMGT/HLA Database is part of the international ImMunoGeneTics project (IMGT).
>
> The database uses the 2010 naming convention for HLA alleles in all tools herein. To aid in the adoption of the new nomenclature, all search tools can be used with both the current and pre-2010 allele designations. The pre-2010 nomenclature designations are only used where older reports or outputs have been made available for download.

---

## HLA Analysis

### arcasHLA: Fast and accurate in silico inference of HLA genotypes from RNA-seq

Python scripts

https://github.com/RabadanLab/arcasHLA

> arcasHLA performs high resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.

### HLAProfiler: Using k-mers to call HLA alleles in RNA sequencing data

Perl scripts

https://github.com/ExpressionAnalysis/HLAProfiler

> HLAProfiler uses the k-mer content of next generation sequencing reads to call HLA types in a sample. Based on the k-mer content each each read pair is assigned to an HLA gene and the aggregate k-mer profile for the gene is compared to reference k-mer profiles to determin the HLA type. Currently HLAProfiler only supports paired-end RNA-seq data.

### MATER: Minimizer RNAseq HLA typer

Python scripts and C code

https://github.com/cschin/MATER

> MATER is a minimizer-based HLA typer for RNAseq read dataset. In a typical RNAseq dataset, the reads sampled from HLA genes are less uniform and may miss regions that makes assembly or variant calling base methods for HLA typing more challenge. Here we adopt a slight different approach. We try to assign each reads to possible HLA types by using minimizers. Namely, we will generate dense minimizer for each reads and compare to those from the HLA type seqeunces.

> We annotate each each reads to possible HLA serotype or 4 digit type sequence according the minimizer matches. Some reads may be able to assign to single HLA type-sequence, some other may be more ambiguous. We derive a simple score to summarize the results from all reads that are mapped to HLA-type sequences for each HLA allele.

### OptiType: Precision HLA typing from next-generation sequencing data

Python scripts

https://github.com/FRED-2/OptiType

> OptiType is a novel HLA genotyping algorithm based on integer linear programming, capable of producing accurate 4-digit HLA genotyping predictions from NGS data by simultaneously selecting all major and minor HLA Class I alleles.

### PHLAT: Inference of High Resolution HLA Types

Python scripts

https://sites.google.com/site/phlatfortype/home

> PHLAT is a bioinformatics algorithm that offers HLA typing at four-digit resolution (or higher) using genome-wide transcriptome and exome sequencing data over a wide range of read lengths and sequencing depths.

### seq2HLA: HLA typing from RNA-Seq sequence reads

Python scripts

https://github.com/TRON-Bioinformatics/seq2HLA

> In-silico method written in Python and R to determine HLA genotypes of a sample. seq2HLA takes standard RNA-Seq sequence reads in fastq format as input, uses a bowtie index comprising all HLA alleles and outputs the most likely HLA class I and class II genotypes (in 4 digit resolution), a p-value for each call, and the expression of each class.

### SNP2HLA: Imputation of Amino Acid Polymorphisms in Human Leukocyte Antigens

C shell script

http://software.broadinstitute.org/mpg/snp2hla/

>SNP2HLA is a tool to impute amino acid polymorphisms and single nucleotide polymorphisms in human luekocyte antigenes (HLA) within the major histocompatibility complex (MHC) region in chromosome 6.
