# awesome-vdj <a href="https://www.rcsb.org/structure/6py2"><img width="40%" align="right" src="https://github.com/slowkow/awesome-vdj/assets/209714/7f322310-89b2-4398-bc2c-2aa730cd095c"></a>

[Antigen] presentation and recognition is central to immunology. [HLA genes] encode the proteins that present antigens. [VDJ genes][vdj] encode the receptors: T cell receptors (TCRs) in T cells and the repertoires of antibodies/immunoglobulins in B cells.

Here, researchers can find links to tools and resources for computational analysis of HLA and VDJ data.

[Contributions are welcome!](https://github.com/slowkow/awesome-vdj/blob/master/CONTRIBUTING.md)

[Antigen]: https://en.wikipedia.org/wiki/Antigen
[vdj]: https://en.wikipedia.org/wiki/V(D)J_recombination
[HLA genes]: https://en.wikipedia.org/wiki/Human_leukocyte_antigen

[![CI](https://github.com/slowkow/awesome-vdj/workflows/CI/badge.svg)](https://github.com/slowkow/awesome-vdj/actions)

**Table of Contents**

- [📚 Literature](#literature)
- [🗄️ VDJ Databases](#vdj-databases)
- [🔬 VDJ Analysis](#vdj-analysis)
- [🗃️ HLA Databases](#hla-databases)
- [🧬 HLA Analysis](#hla-analysis)

**Related Work**

- Ming Tang's list: [TCR-BCR-seq-analysis](https://github.com/crazyhottommy/TCR-BCR-seq-analysis)

---

## 📚 Literature
- [**Why must T cells be cross-reactive?**](https://pubmed.ncbi.nlm.nih.gov/22918468/) — This perspective article discusses the immunological necessity of T cell cross-reactivity, explaining why each T cell must be capable of recognizing multiple different peptide-MHC complexes to prov...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/22918468/) · 🪝 [499](https://www.semanticscholar.org/paper/09322c7a1bebd9d48fa41ee1c7b85660b381bec1)

- [**Mechanisms of central tolerance for B cells**](https://pubmed.ncbi.nlm.nih.gov/28368006/) — An in-depth review of the mechanisms by which developing B cells are tolerized to self-antigens in the bone marrow, including receptor editing, clonal deletion, and anergy, and how failures in thes...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28368006/) · 🪝 [418](https://www.semanticscholar.org/paper/8942bac84eab2aaf30a15c28f68c40af880e0477)

- [**Understanding the drivers of MHC restriction of T cell receptors**](https://pubmed.ncbi.nlm.nih.gov/29636542/) — A comprehensive review examining how T cell receptors (TCRs) are restricted to recognizing peptide antigens presented by major histocompatibility complex (MHC) molecules, exploring the evolutionary...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/29636542/) · 🪝 [254](https://www.semanticscholar.org/paper/b3b63f41169bef73eb8d34b911526d4f88d93752)

- [**High-Throughput and Single-Cell T Cell Receptor Sequencing Technologies**](https://doi.org/10.1038/s41592-021-01201-8) — A comprehensive review of current technologies for T cell receptor sequencing, covering both bulk and single-cell approaches, their applications in immunology research, and future directions in the...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34282327/) · 🪝 [208](https://www.semanticscholar.org/paper/471eba65c98c464da909d699c05a3998033e204f)

- [**Mining adaptive immune receptor repertoires for biological and clinical information using machine learning**](http://dx.doi.org/10.1016/j.coisb.2020.10.010) — A review of machine learning approaches for analyzing adaptive immune receptor repertoire data, discussing how these methods can extract biological insights and clinical information from large-scal...<br>[Paper](http://dx.doi.org/10.1016/j.coisb.2020.10.010) · 🪝 [72](https://www.semanticscholar.org/paper/b1606c027b9224972cbb4d25a6c36cb03f573099)

- [**HLA and kidney disease: from associations to mechanisms**](https://pubmed.ncbi.nlm.nih.gov/30206339/) — This review explores the associations between HLA genes and kidney diseases, discussing how advances in understanding HLA biology are revealing the mechanisms underlying these genetic associations ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/30206339/) · 🪝 [69](https://www.semanticscholar.org/paper/8a6ca498d44b884d163b3de4a6c3466658840c4a)

- [**SweHLA: the high confidence HLA typing bio-resource drawn from 1000 Swedish genomes**](https://www.nature.com/articles/s41431-019-0559-2) — This paper presents SweHLA, a high-confidence HLA typing resource derived from whole-genome sequencing of 1000 Swedish individuals, providing a valuable reference for HLA research and clinical appl...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31844174/) · 🪝 [16](https://www.semanticscholar.org/paper/d8bc4828bd344a35e19221e1001d01dc156f1e1d)


---

## 🗄️ VDJ Databases
### Structure Databases

- [**STCRDab: The Structural T-Cell Receptor Database**](http://opig.stats.ox.ac.uk/webapps/stcrdab/) — An automated, curated set of T-Cell Receptor structural data from the PDB.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/29087479/) · 🪝 [102](https://www.semanticscholar.org/paper/c64329d5e91e3d2fd4d2cc564a09ae1b5301f57c) · [Homepage](http://opig.stats.ox.ac.uk/webapps/stcrdab/)

- [**TCR3d: T cell receptor structural repertoire database**](https://tcr3d.ibbr.umd.edu/) — Welcome to the T cell receptor (TCR) structural repertoire database. Here we provide an easy-to-use interface to view all experimentally determined T cell receptor structures and their complexes. T...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32124321/) · 🪝 [3](https://www.semanticscholar.org/paper/8e0ccc15dd5550c762e62da5f8e3c04449f87d8f) · [Homepage](https://tcr3d.ibbr.umd.edu/)

- [**Coronavirus-Binding Antibody Sequences & Structures**](http://opig.stats.ox.ac.uk/webapps/covabdab/) — The Oxford Protein Informatics Group (Dept. of Statistics, University of Oxford) is collaborating in efforts to understand the immune response to SARS-CoV2 infection and vaccination. As part of our...<br>[Homepage](http://opig.stats.ox.ac.uk/webapps/covabdab/)

### Specificity Databases

- [**VDJDB: A curated database of T-cell receptor sequences of known antigen specificity**](https://github.com/antigenomics/vdjdb-db) — The primary goal of VDJdb is to facilitate access to existing information on T-cell receptor antigen specificities, i.e. the ability to recognize certain epitopes in certain MHC contexts. > Our mis...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28977646/) · 🪝 [491](https://www.semanticscholar.org/paper/cfd86d8ddd03ccacd18343d091ac93745e4187d6) · ⭐ [149](https://github.com/antigenomics/vdjdb-db/stargazers) · [Homepage](https://vdjdb.cdr3.net)

- [**McPAS-TCR: A manually curated catalogue of pathology associated T-cell receptor sequences**](https://friedmanlab.weizmann.ac.il/McPAS-TCR/) — McPAS-TCR is a manually curated catalogue of T cell receptor (TCR) sequences that were found in T cells associated with various pathological conditions in humans and in mice. It is meant to link TC...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28481982/) · 🪝 [449](https://www.semanticscholar.org/paper/af8ef665e0cd44a7cf69b811626519a3fcf323b4) · [Homepage](https://friedmanlab.weizmann.ac.il/McPAS-TCR/)

- [**vdjmatch**](https://github.com/antigenomics/vdjmatch) — Matching T-cell repertoire against a database of TCR antigen specificities<br>⭐ [39](https://github.com/antigenomics/vdjmatch/stargazers) · [Homepage](https://vdjdb.cdr3.net) · `Groovy`

### Sequence Repositories

- [**immuneACCESS**](https://github.com/slowkow/awesome-vdj/blob/master/download-from-immuneaccess.md) — Dive into the world’s largest collection of TCR and BCR sequences. Easily incorporate millions of sequences worth of public data into your next papers and projects using immunoSEQ Analyzer. Constru...<br>⭐ [235](https://github.com/slowkow/awesome-vdj/blob/master/download-from-immuneaccess.md/stargazers) · [Homepage](https://clients.adaptivebiotech.com/immuneaccess)

- [**iReceptor**](https://gateway.ireceptor.org/home) — iReceptor facilitates the curation, analysis and sharing of antibody/B-cell and T-cell receptor repertoires (Adaptive Immune Receptor Repertoire or AIRR-seq data) from multiple labs and institution...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/29944754/) · 🪝 [133](https://www.semanticscholar.org/paper/5d764e3cb11d09a8f2ec8bdce3b390d6e42d3f8a) · [Homepage](https://gateway.ireceptor.org/home)

- [**A Public Database of Memory and Naive B-Cell Receptor Sequences**](https://datadryad.org/stash/dataset/doi:10.5061/dryad.35ks2) — We present a public database of more than 37 million unique BCR sequences from three healthy adult donors that is many fold deeper than any existing resource, together with a set of online tools de...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/27513338/) · 🪝 [104](https://www.semanticscholar.org/paper/0296d2dce034afee35366908584f9daa81ea7319) · [Homepage](https://datadryad.org/stash/dataset/doi:10.5061/dryad.35ks2)

- [**PIRD: Pan immune repertoire database**](https://db.cngb.org/pird/) — Pan immune repertoire database (PIRD) collects raw and processed sequences of immunoglobulins (IGs) and T cell receptors (TCRs) of human and other vertebrate species with different phenotypes. You ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31373607/) · [Homepage](https://db.cngb.org/pird/)

### Standards & Resources

- [**Adaptive Immune Receptor Repertoire (AIRR) Community**](https://github.com/airr-community) — The Adaptive Immune Receptor Repertoire (AIRR) Community of The Antibody Society is a research-driven group that is organizing and coordinating stakeholders in the use of next-generation sequencing...<br>[Docs](https://docs.airr-community.org/en/stable/index.html) · [Homepage](http://airr-community.org)

- [**Human Vaccines Project (Human Immunome Program)**](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP174305) — The Human Immunome Program (HIP) is open-source effort with the goal sequencing all of the adaptive receptors on the surface of human B and T cells. Under a targeted 7-to-10-year effort, the progra...<br>[Homepage](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP174305)


---

## 🔬 VDJ Analysis
### Single-Cell

- [**TRUST4: TCR and BCR assembly from RNA-seq data**](https://github.com/liulab-dfci/TRUST4) — Tcr Receptor Utilities for Solid Tissue (TRUST) is a computational tool to analyze TCR and BCR sequences using unselected RNA sequencing data, profiled from solid tissues, including tumors. TRUST4 ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33986545/) · 🪝 [227](https://www.semanticscholar.org/paper/7564c0e07f7135c0ec2eddb4009e6a51febdc991) · ⭐ [337](https://github.com/liulab-dfci/TRUST4/stargazers) · `C` `C++` `Perl`

- [**Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data**](https://github.com/scverse/scirpy) — A scalable Python toolkit that provides simplified access to the analysis and visualization of immune repertoires from single cells and seamless integration with transcriptomic data.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32614448/) · 🪝 [212](https://www.semanticscholar.org/paper/fcd27b7bd7ba5b02c64910cf80c2b5b7fabd12e4) · ⭐ [243](https://github.com/scverse/scirpy/stargazers) · [Homepage](https://scirpy.scverse.org/en/latest/) · `Python`

- [**scirpy: A scanpy extension to analyse single-cell TCR data.**](https://github.com/icbi-lab/scirpy) — Scirpy is a scalable python-toolkit to analyse T cell receptor (TCR) repertoires from single-cell RNA sequencing (scRNA-seq) data. It seamlessly integrates with the popular scanpy library and provi...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32614448/) · 🪝 [175](https://www.semanticscholar.org/paper/fcd27b7bd7ba5b02c64910cf80c2b5b7fabd12e4) · ⭐ [243](https://github.com/icbi-lab/scirpy/stargazers) · `Python`

- [**scRepertoire: A toolkit for single-cell immune profiling**](https://github.com/BorchLab/scRepertoire) — R package for analyzing and visualizing single-cell immune receptor data. This new version introduces an array of features designed to enhance both the depth and breadth of immune receptor analysis...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/40577285/) · 🪝 [9](https://www.semanticscholar.org/paper/1a0dc99021ccfd16d1d3a19f75068de450bc25f6) · ⭐ [358](https://github.com/BorchLab/scRepertoire/stargazers) · `R`

- [**DeepTCR: Deep Learning Methods for Parsing T-Cell Receptor Sequencing (TCRSeq) Data**](https://github.com/sidhomj/DeepTCR) — DeepTCR is a python package that has a collection of unsupervised and supervised deep learning methods to parse TCRSeq data. It has the added functionality of being able to analyze paired alpha/bet...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33707415/) · 🪝 [217](https://www.semanticscholar.org/paper/a1a242bdb47b7fe9e8d519aacb41157bb78842fb) · ⭐ [123](https://github.com/sidhomj/DeepTCR/stargazers) · `Python`

- [**dandelion**](https://github.com/zktuong/dandelion) — dandelion - A single cell BCR/TCR V(D)J-seq analysis package for 10X Chromium 5' data<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/37055623/) · 🪝 [44](https://www.semanticscholar.org/paper/a3773c9b4d58198f7195df8f54a8c58ca892abce) · ⭐ [122](https://github.com/zktuong/dandelion/stargazers) · [Homepage](https://sc-dandelion.readthedocs.io/) · `Python`

- [**STARTRAC**](https://github.com/Japrin/STARTRAC) — STARTRAC(Single T-cell Analysis by Rna-seq and Tcr TRACking)<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33900375/) · 🪝 [27](https://www.semanticscholar.org/paper/11c92257f87b515bc46af2b874ff14890120fadd) · ⭐ [114](https://github.com/Japrin/STARTRAC/stargazers) · `HTML`

- [**TCRGP**](https://github.com/emmijokinen/TCRGP) — TCRGP is a novel Gaussian process method that can predict if TCRs recognize certain epitopes. This method can utilize different CDR sequences from both TCRα and TCRβ chains from single-cell data an...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33764977/) · 🪝 [109](https://www.semanticscholar.org/paper/ee044c6c9cbe0cbf79bfcc41d1e4c5c601d25ed1) · ⭐ [30](https://github.com/emmijokinen/TCRGP/stargazers) · `Python`

- [**CONGA: Clonotype Neighbor Graph Analysis**](https://github.com/phbradley/conga) — CONGA was developed to detect correlation between T cell gene expression profile and TCR sequence in single-cell datasets.<br>[Paper](https://doi.org/10.1101/2020.06.04.134536) · 🪝 [9](https://www.semanticscholar.org/paper/d0a9125325f851f69dbc486e2b2e75f9ba63d4f5) · ⭐ [93](https://github.com/phbradley/conga/stargazers) · `Python`

- [**airrflow**](https://github.com/nf-core/airrflow) — B-cell and T-cell Adaptive Immune Receptor Repertoire (AIRR) sequencing analysis pipeline using the Immcantation framework<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/38293151/) · 🪝 [10](https://www.semanticscholar.org/paper/https://www.semanticscholar.org/paper/04c2e0be97ba6d6035506595694eb22e2093037b) · ⭐ [73](https://github.com/nf-core/airrflow/stargazers) · [Homepage](https://nf-co.re/airrflow) · `Nextflow`

- [**Platypus**](https://github.com/alexyermanos/Platypus) — R package for the analysis of single-cell immune repertoires<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33884369/) · 🪝 [40](https://www.semanticscholar.org/paper/309e51a1d8ff5c00c98dc3ea8e102292a30dfba7) · ⭐ [43](https://github.com/alexyermanos/Platypus/stargazers) · `R`

- [**mvTCR**](https://github.com/SchubertLab/mvTCR) — A multi-view Variational Autoencoder (mvTCR) to jointly embed transcriptomic and TCR sequence information at a single-cell level to better capture the phenotypic behavior of T cells.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/38956082/) · 🪝 [17](https://www.semanticscholar.org/paper/62559a2f08e304d5a6149f4605e45529ac2c150e) · ⭐ [56](https://github.com/SchubertLab/mvTCR/stargazers) · [Homepage](https://zenodo.org/record/5006839) · `Python`

- [**enclone**](https://github.com/10XGenomics/enclone) — enclone is standalone software (primarily written in Rust) developed by 10x Genomics for analysis of single cell TCR and BCR sequences. enclone performs SHM-aware clonotyping, phylogenetic/lineage ...<br>⭐ [50](https://github.com/10XGenomics/enclone/stargazers) · [Homepage](https://10xgenomics.github.io/enclone/) · `Rust`

- [**covid19**](https://github.com/immunomind/covid19) — Regularly updated list of publicly available datasets with single-cell (scRNAseq) and T-cell/antibody immune repertoire (AIRR / RepSeq / immunosequencing) data of COVID-19 patients with SARS-CoV-2.<br>⭐ [46](https://github.com/immunomind/covid19/stargazers)

- [**TCRconvert**](https://github.com/seshadrilab/tcrconvert) — TCRconvert converts T cell receptor (TCR) gene names between the 10X, Adaptive, and IMGT naming conventions. It supports alpha-beta and gamma-delta TCRs for human, mouse, and rhesus macaque.<br>⭐ [15](https://github.com/seshadrilab/tcrconvert/stargazers) · `Python`

- [**TCRconvertR**](https://github.com/seshadrilab/tcrconvertr) — TCRconvertR converts T cell receptor (TCR) gene names between the 10X, Adaptive, and IMGT naming conventions. It supports alpha-beta and gamma-delta TCRs for human, mouse, and rhesus macaque.<br>⭐ [6](https://github.com/seshadrilab/tcrconvertr/stargazers) · `R`

### Repertoire Analysis

- [**VDJtools**](https://github.com/mikessh/vdjtools) — A comprehensive analysis framework for T-cell and B-cell repertoire sequencing data<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/26606115/) · 🪝 [529](https://www.semanticscholar.org/paper/c5994d9f6ed808f510cb95a3225c9f8ab0d6b460) · ⭐ [142](https://github.com/mikessh/vdjtools/stargazers) · `Java` `Groovy`

- [**immunarch: An R Package for Painless Bioinformatics Analysis of T-cell and B-cell Immune Repertoire Data**](https://github.com/immunomind/immunarch) — immunarch is an R package designed to analyse T-cell receptor (TCR) and B-cell receptor (BCR) repertoires, aimed at medical scientists and bioinformaticians. The mission of immunarch is to make imm...<br>⭐ [334](https://github.com/immunomind/immunarch/stargazers) · `R`

- [**msm: Max Snippet Model**](https://github.com/jostmey/msm) — Improved statistical classifier for immune repertoires<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/33868241/) · 🪝 [8](https://www.semanticscholar.org/paper/77b6920d6e016f551c73c7d6eb5ac385128772f5) · ⭐ [177](https://github.com/jostmey/msm/stargazers) · `Python`

- [**DeepRC**](https://github.com/ml-jku/DeepRC) — DeepRC: Immune repertoire classification with attention-based deep massive multiple instance learning<br>⭐ [124](https://github.com/ml-jku/DeepRC/stargazers) · `Python`

- [**Recon: Reconstruction of Estimated Communities from Observed Numbers**](https://github.com/ArnaoutLab/Recon) — Recon uses the distribution of species counts in a sample to estimate the distribution of species counts in the population from which the sample was drawn.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/27302887/) · 🪝 [91](https://www.semanticscholar.org/paper/42674800f3ce4fc5230fbb08afd3cbfb72a47e02) · ⭐ [14](https://github.com/ArnaoutLab/Recon/stargazers) · `Python` `R`

- [**dkm: Dynamic Kernel Matching**](https://github.com/jostmey/dkm) — DKM is analogous to a convolutional network, but for sequences. Consider the problem of classifying a sequence. Because some sequences are longer than others, the number of features is irregular. G...<br>⭐ [94](https://github.com/jostmey/dkm/stargazers) · `Python`

- [**immuneML**](https://github.com/uio-bmi/immuneML) — immuneML is a platform for machine learning analysis of adaptive immune receptor repertoire data.<br>⭐ [73](https://github.com/uio-bmi/immuneML/stargazers) · [Homepage](https://immuneml.uio.no) · `Python`

- [**abstar**](https://github.com/brineylab/abstar) — VDJ assignment and antibody sequence annotation. Scalable from a single sequence to billions of sequences.<br>⭐ [44](https://github.com/brineylab/abstar/stargazers) · `Pkl`

- [**vdjer**](https://github.com/mozack/vdjer) — V'DJer -  B Cell Receptor Repertoire Reconstruction from short read mRNA-Seq data<br>⭐ [29](https://github.com/mozack/vdjer/stargazers) · `C`

- [**CATT**](https://github.com/GuoBioinfoLab/CATT) — An ultra-sensitive and precise tool for characterizing T cell CDR3 sequences in TCR-seq and RNA-seq data.<br>⭐ [21](https://github.com/GuoBioinfoLab/CATT/stargazers) · [Homepage](http://bioinfo.life.hust.edu.cn/CATT/) · `Julia`

- [**epitopefindr**](https://github.com/brandonsie/epitopefindr) — R package to BLAST peptide sequences against each other and identify the minimal overlap of aligning regions.<br>⭐ [16](https://github.com/brandonsie/epitopefindr/stargazers) · [Homepage](https://brandonsie.github.io/epitopefindr/) · `R`

### Sequence Processing

- [**PRESTO: The REpertoire Sequencing TOolkit**](https://github.com/immcantation/presto) — pRESTO is a toolkit for processing raw reads from high-throughput sequencing of B cell and T cell repertoires. > The REpertoire Sequencing TOolkit (pRESTO) is composed of a suite of utilities to ha...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/24618469/) · 🪝 [424](https://www.semanticscholar.org/paper/ca321103129928fcc5522d2a314fdd167bf8805e) · [Docs](https://presto.readthedocs.io/en/stable) · `Python`

- [**MiXCR: a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data**](https://github.com/milaboratory/mixcr) — MiXCR is a universal framework that processes big immunome data from raw sequences to quantitated clonotypes. MiXCR efficiently handles paired- and single-end reads, considers sequence quality, cor...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/39433438/) · 🪝 [4](https://www.semanticscholar.org/paper/0eaa429866419cf0a165ccae340447d6f2a4c925) · ⭐ [380](https://github.com/milaboratory/mixcr/stargazers) · `Java`

- [**IMSEQ: IMmunogenetic SEQuence Analysis**](https://github.com/lkuchenb/imseq) — IMSEQ is a fast, PCR and sequencing error aware tool to analyze high throughput data from recombined T-cell receptor or immunoglobolin gene sequencing experiments. It derives immune repertoires fro...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/25987567/) · 🪝 [98](https://www.semanticscholar.org/paper/0b738ce3da3d79cefd8c1460687878f60bd05183) · ⭐ [15](https://github.com/lkuchenb/imseq/stargazers)

- [**vidjil**](https://github.com/vidjil/vidjil) — Vidjil -- High-throughput Analysis of V(D)J Immune Repertoire (mirror, please go to http://gitlab.vidjil.org)<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/27835690/) · 🪝 [81](https://www.semanticscholar.org/paper/03f774017c20e5297317851016bb37be66c35291) · ⭐ [31](https://github.com/vidjil/vidjil/stargazers) · [Homepage](http://gitlab.vidjil.org) · `JavaScript`

- [**stitchr**](https://github.com/JamieHeather/stitchr) — Stitchr - a Python script to stitch together coding TCR nucleotide sequences from V, J, and CDR3 info<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/35325179/) · 🪝 [22](https://www.semanticscholar.org/paper/aeaf4f3d97ca02f2b25cf040d6ad39a62db41a1e) · ⭐ [64](https://github.com/JamieHeather/stitchr/stargazers) · [Homepage](https://jamieheather.github.io/stitchr/) · `Python`

- [**pyIR: An IgBLAST wrapper and parser**](https://github.com/crowelab/PyIR) — PyIR is a minimally-dependent high-speed wrapper for the IgBLAST immunoglobulin and T-cell analyzer. This is achieved through chunking the input data set and running IgBLAST single-core in parallel...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32677886/) · 🪝 [29](https://www.semanticscholar.org/paper/4c463450d5f5f6d9e1c1f5a21b3b33cbd4ed141c) · ⭐ [50](https://github.com/crowelab/PyIR/stargazers) · `Python`

- [**vdjviz**](https://github.com/antigenomics/vdjviz) — A lightweight immune repertoire browser<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/27297497/) · 🪝 [33](https://www.semanticscholar.org/paper/f7da52e1bc9aabf2f08b78dd976bb7c39f2388c6) · ⭐ [27](https://github.com/antigenomics/vdjviz/stargazers) · [Homepage](https://vdjviz.cdr3.net) · `JavaScript`

- [**MiGMAP: mapper for full-length T- and B-cell repertoire sequencing**](https://github.com/mikessh/migmap) — In a nutshell, this software is a smart wrapper for IgBlast V-(D)-J mapping tool designed to facilitate analysis immune receptor libraries profiled using high-throughput sequencing. This package in...<br>⭐ [53](https://github.com/mikessh/migmap/stargazers) · `Java` `Groovy`

- [**BepiPred-3.0**](https://github.com/UberClifford/BepiPred-3.0) — BepiPred3.0 predicts B-cell epitopes from proteins sequences in fasta format.<br>⭐ [17](https://github.com/UberClifford/BepiPred-3.0/stargazers) · `HTML`

### Clustering & Similarity

- [**tcr-dist**](https://github.com/phbradley/tcr-dist) — Software tools for the analysis of epitope-specific T cell receptor (TCR) repertoires<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28636592/) · 🪝 [815](https://www.semanticscholar.org/paper/b3e8d6f21fbdcd58888af31e791b5a8d24a1c592) · ⭐ [86](https://github.com/phbradley/tcr-dist/stargazers) · `Python`

- [**tcrdist3**](https://github.com/kmayerb/tcrdist3) — tcrdist3 is a Python API-enabled toolkit for analyzing T-cell receptor repertoires. Some of the functionality and code is adapted from the original tcr-dist package.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34845983/) · 🪝 [126](https://www.semanticscholar.org/paper/fe7c08b0dae5d9fc6f667a9f222559ff976c5217) · ⭐ [69](https://github.com/kmayerb/tcrdist3/stargazers) · `Python`

- [**GIANA: Geometry Isometry based TCR AligNment Algorithm**](https://github.com/s175573/GIANA) — Geometric Isometry- based TCR AligNment Algorithm (GIANA), a mathematical framework to transform the CDR3 sequences, which converted the sequence alignment and clustering problem into a classic nea...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34349111/) · 🪝 [102](https://www.semanticscholar.org/paper/6488fe7fe980684b1ae8cbe4ed3623977b2b6628) · ⭐ [71](https://github.com/s175573/GIANA/stargazers) · `Python`

- [**ClusTCR: a Python interface for rapid clustering of large sets of CDR3 sequences with unknown antigen specificity**](https://github.com/svalkiers/clusTCR) — CDR3 clustering module providing a new method for fast and accurate clustering of large data sets of CDR3 amino acid sequences, and offering functionalities for downstream analysis of clustering re...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34132766/) · 🪝 [58](https://www.semanticscholar.org/paper/d8937080ab6f1d2fcbcce7c49937d29bc027ee8b) · ⭐ [54](https://github.com/svalkiers/clusTCR/stargazers) · `Python`

- [**immuneSIM: Tunable Simulation of B- And T-Cell Receptor Repertoires**](https://github.com/GreiffLab/immuneSIM) — Simulate full B-cell and T-cell receptor repertoires using an in silico recombination process that includes a wide variety of tunable parameters to introduce noise and biases. Additional post-simul...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32154832/) · 🪝 [63](https://www.semanticscholar.org/paper/https://www.semanticscholar.org/paper/530b5f57cc806f6ee93c66c3db94df8875693c73) · ⭐ [38](https://github.com/GreiffLab/immuneSIM/stargazers) · `R`

- [**ALICE: Antigen-specific Lymphocyte Identification by Clustering of Expanded sequences**](https://github.com/pogorely/ALICE) — Detecting TCR involved in immune responses from single RepSeq datasets.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31194732/) · ⭐ [27](https://github.com/pogorely/ALICE/stargazers) · `R`

- [**ImReP: Rapid and accurate profiling of the adaptive immune repertoires from regular RNA-Seq data**](https://github.com/Mangul-Lab-USC/imrep) — ImReP is a method to quantify individual immune response based on a recombination landscape of genes encoding B and T cell receptors (BCR and TCR). ImReP is able to efficiently extract TCR and BCR ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32561710/) · ⭐ [11](https://github.com/Mangul-Lab-USC/imrep/stargazers) · `Python`

### Epitope Prediction

- [**epitopepredict**](https://github.com/dmnfarrell/epitopepredict) — Python package and command line tool for epitope prediction<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/36824339/) · 🪝 [6](https://www.semanticscholar.org/paper/073554e81b4370b4f409fa7bdedaa9c36e78d83f) · ⭐ [52](https://github.com/dmnfarrell/epitopepredict/stargazers) · `Jupyter Notebook`

- [**MuPeXI**](https://github.com/ambj/MuPeXI) — MuPeXI: the mutant peptide extractor and informer, a tool for predicting neo-epitopes from tumor sequencing data.<br>⭐ [52](https://github.com/ambj/MuPeXI/stargazers) · `Python`

- [**epitopeprediction**](https://github.com/nf-core/epitopeprediction) — A bioinformatics best-practice analysis pipeline for epitope prediction and annotation<br>⭐ [49](https://github.com/nf-core/epitopeprediction/stargazers) · [Homepage](https://nf-co.re/epitopeprediction) · `Nextflow`

- [**neoantigens**](https://github.com/umccr/neoantigens) — Exploring novel tumor epitope identification<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/36604431/) · ⭐ [37](https://github.com/umccr/neoantigens/stargazers) · `Python`

- [**MixTCRpred**](https://github.com/GfellerLab/MixTCRpred) — Predictor of TCR-epitope interactions<br>⭐ [34](https://github.com/GfellerLab/MixTCRpred/stargazers) · `Python`

- [**topiary**](https://github.com/openvax/topiary) — Predict mutated T-cell epitopes from sequencing data<br>⭐ [30](https://github.com/openvax/topiary/stargazers) · `Python`

- [**AsEP-dataset**](https://github.com/biochunan/AsEP-dataset) — NeurIPS 2024 Dataset and Benchmark Submission "AsEP: Benchmarking Deep Learning Methods for Antibody-specific Epitope Prediction"<br>⭐ [30](https://github.com/biochunan/AsEP-dataset/stargazers) · `Jupyter Notebook`

- [**EpiDope**](https://github.com/rnajena/EpiDope) — Prediction of B-cell epitopes from amino acid sequences using deep neural networks.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34109374/) · 🪝 [11](https://www.semanticscholar.org/paper/73df19c2feb455fb3df7aaec1c6ebde0c85305c6) · ⭐ [18](https://github.com/rnajena/EpiDope/stargazers) · `Python`

- [**Repitope**](https://github.com/masato-ogishi/Repitope) — Epitope immunogenicity prediction through in silico TCR-peptide contact potential profiling.<br>⭐ [25](https://github.com/masato-ogishi/Repitope/stargazers) · `R`

- [**pyrepseq**](https://github.com/andim/pyrepseq) — Python library for immune repertoire analysis<br>[Paper](https://www.pnas.org/doi/10.1073/pnas.2213264120) · ⭐ [17](https://github.com/andim/pyrepseq/stargazers) · `Python`

- [**ImRex**](https://github.com/pmoris/ImRex) — Generic TCR-epitope recognition prediction using CNN approach on both known and novel epitopes<br>⭐ [17](https://github.com/pmoris/ImRex/stargazers) · `Jupyter Notebook`

### Structure & Modeling

- [**TITAN - Tcr epITope bimodal Attention Networks**](https://github.com/PaccMann/TITAN) — a bimodal neural network that explicitly encodes both TCR sequences and epitopes to enable the independent study of generalization capabilities to unseen TCRs and/or epitopes.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34252922/) · 🪝 [150](https://www.semanticscholar.org/paper/a732443cae8cd2d6a76f4f3cf785a562baf41137) · ⭐ [30](https://github.com/PaccMann/TITAN/stargazers) · `Python`

- [**TCRdock**](https://github.com/phbradley/TCRdock) — Python tools for TCR:peptide-MHC modeling and analysis: - Set up and run TCR-specialized AlphaFold simulations starting from a TSV file with TCR, peptide, and MHC information. - Parse a TCR:peptide...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/36661395/) · 🪝 [93](https://www.semanticscholar.org/paper/daa04b1951a9f7d7390754ea48144ab18ccb9b0c) · ⭐ [86](https://github.com/phbradley/TCRdock/stargazers) · `Python`

- [**Absolut: Unconstrained lattice antibody-antigen bindings generator - One tool to simulate them all!**](https://github.com/csi-greifflab/Absolut) — Absolut! is a database and C++ user interface that allows the high-throughput computation for the 3D-lattice binding of any CDRH3 sequence to any antigen, enabling the custom generation of new anti...<br>[Paper](https://doi.org/10.1101/2021.07.06.451258) · 🪝 [20](https://www.semanticscholar.org/paper/77e244e5e7df68c8afeebababc8774a07290964a) · ⭐ [111](https://github.com/csi-greifflab/Absolut/stargazers) · `C++`

- [**tcr-bert**](https://github.com/wukevin/tcr-bert) — TCR-BERT is a large language model trained on T-cell receptor sequences, built using a lightly modified BERT architecture with tweaked pre-training objectives.<br>[Paper](http://dx.doi.org/10.1101/2021.11.18.469186) · 🪝 [74](https://www.semanticscholar.org/paper/7ef95e6164999fe9fc6d30ce2b64e8f0cabaf225) · ⭐ [57](https://github.com/wukevin/tcr-bert/stargazers) · `Python`

- [**TCRmodel2: high-resolution modeling of T cell receptor recognition using deep learning**](https://github.com/piercelab/tcrmodel2) — This method, named TCRmodel2, allows users to submit sequences through an easy-to-use interface and shows similar or greater accuracy than AlphaFold and other methods to model TCR–peptide–MHC compl...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/37140040/) · 🪝 [70](https://www.semanticscholar.org/paper/ae85735476489e8e49f3dac80dd9bc27bf9d6b52) · ⭐ [45](https://github.com/piercelab/tcrmodel2/stargazers) · `Python` `R`

- [**vampire: Deep generative models for TCR sequences**](https://github.com/matsengrp/vampire/) — Fit and test variational autoencoder (VAE) models for T cell receptor sequences.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31487240/) · 🪝 [66](https://www.semanticscholar.org/paper/6597f693534cafff625af2122f929ff8a2577e80) · ⭐ [17](https://github.com/matsengrp/vampire/stargazers) · `Python`

- [**TEINet**](https://github.com/jiangdada1221/TEINet) — TEINet: a deep learning framework for prediction of TCR-epitope binding specificity<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/36907658/) · 🪝 [55](https://www.semanticscholar.org/paper/f7d928737a616310666da93d023f1477a6e709d9) · ⭐ [16](https://github.com/jiangdada1221/TEINet/stargazers) · `Python`

- [**TEIM**](https://github.com/pengxingang/TEIM) — TEIM: TCR-Epitope Interaction Modeling<br>⭐ [55](https://github.com/pengxingang/TEIM/stargazers) · `Python`

- [**TCRconv**](https://github.com/emmijokinen/TCRconv) — TCRconv is a deep learning model for predicting recognition between T cell receptors and epitopes. It uses protBERT embeddings for the TCRs and convolutional neural networks for the prediction.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/36477794/) · 🪝 [20](https://www.semanticscholar.org/paper/e7408fe04e2819eecccfe2a2425d7cdfe40145ff) · ⭐ [26](https://github.com/emmijokinen/TCRconv/stargazers) · `Python` `R`

- [**compairr**](https://github.com/uio-bmi/compairr) — Comparison of Adaptive Immune Receptor Repertoires<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/35852318/) · 🪝 [15](https://www.semanticscholar.org/paper/d01242ad2ed3bd695117f33a02773fc20f8c0c4d) · ⭐ [28](https://github.com/uio-bmi/compairr/stargazers) · `C++`


---

## 🗃️ HLA Databases
- [**Nomenclature of HLA Alleles**](https://hla.alleles.org/nomenclature/index.html) — A Nomenclature Committee composed of geneticists and immunologists, including specialists in tissue typing, has met after each of the Histocompatibility Workshops beginning with the Third Workshop ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/26760826/) · 🪝 [2605](https://www.semanticscholar.org/paper/68a5ac30187f0350760786cee36df3b2fb3f2b11) · [Homepage](https://hla.alleles.org/nomenclature/index.html)

- [**IEDB: Immune Epitope Database and Analysis Resource**](https://www.iedb.org/) — The Immune Epitope Database (IEDB) is a freely available resource funded by NIAID. It catalogs experimental data on antibody and T cell epitopes studied in humans, non-human primates, and other ani...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/30357391/) · 🪝 [1835](https://www.semanticscholar.org/paper/288b317e427c6bf4c94d455049bd1368ff2071eb) · [Homepage](https://www.iedb.org/)

- [**Allele Frequency Net Database**](http://www.allelefrequencies.net/collaborators.asp) — AFND is a public resource that collects information on allele, genotype and haplotype frequencies from different polymorphic areas in the human genome such as human leukocyte antigens (HLA), killer...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31722398/) · 🪝 [641](https://www.semanticscholar.org/paper/3ee580308b6c1e1f3fadc00690fe871587a70885) · [Homepage](http://www.allelefrequencies.net/collaborators.asp)

- [**IMGTHLA**](https://github.com/ANHIG/IMGTHLA) — The IPD-IMGT/HLA Database provides a specialist database for sequences of the human major histocompatibility complex (MHC) and includes the official sequences named by the WHO Nomenclature Committe...<br>⭐ [246](https://github.com/ANHIG/IMGTHLA/stargazers) · [Homepage](https://www.ebi.ac.uk/ipd/imgt/hla/)

- [**pHLA3D: An online database of predicted three-dimensional structures of HLA molecules**](https://www.phla3d.com.br/) — The limited number of solved HLA structures available in the literature led our research group to develop, in 2019, the pHLA3D, an online database of predicted three-dimensional structures of HLA m...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31239187/) · 🪝 [98](https://www.semanticscholar.org/paper/810def2f2bb6693affa9a449255555ea71fdc064) · [Homepage](https://www.phla3d.com.br/)


---

## 🧬 HLA Analysis
### Association Studies

- [**BIGDAWG: Case-Control Analysis of Multi-Allelic Loci**](https://github.com/IgDAWG/BIGDAWG) — Data sets and functions for chi-squared Hardy-Weinberg and case-control association tests of highly polymorphic genetic data [e.g., human leukocyte antigen (HLA) data]. Performs association tests a...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/26708359/) · 🪝 [82](https://www.semanticscholar.org/paper/8e4c6ba503b7d729681e09f60b22782fd36b3ad4) · ⭐ [3](https://github.com/IgDAWG/BIGDAWG/stargazers) · `R`

- [**HLA_analyses_tutorial**](https://github.com/immunogenomics/HLA_analyses_tutorial) — A thorough tutorial on HLA imputation and association, accompanying our manuscript "Tutorial: A statistical genetics guide to identifying HLA alleles driving complex disease"<br>⭐ [70](https://github.com/immunogenomics/HLA_analyses_tutorial/stargazers) · `Jupyter Notebook`

- [**HLA-TAPAS: HLA-Typing At Protein for Association Studies**](https://github.com/immunogenomics/HLA-TAPAS) — An HLA-focused pipeline that can handle HLA reference panel construction (MakeReference), HLA imputation (SNP2HLA), and HLA association (HLAassoc). It is an updated version of the SNP2HLA.<br>⭐ [54](https://github.com/immunogenomics/HLA-TAPAS/stargazers) · `Python` `R`

- [**HLA Electrostatic Potential**](https://pubmed.ncbi.nlm.nih.gov/30429288/) — A method for predicting humoral alloimmunity from differences in donor and recipient HLA surface electrostatic potential, enabling assessment of immunological compatibility in transplantation.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/30429288/) · 🪝 [53](https://www.semanticscholar.org/paper/1c38ce8e6f0a43a1c806a8cd65d105879459bdb4)

- [**HATK: HLA Analysis Toolkit**](https://github.com/WansonChoi/HATK) — HATK(HLA Analysis Tool-Kit) is a collection of tools and modules to perform HLA fine-mapping analysis, which is to identify which HLA allele or amino acid position of the HLA gene is driving the di...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32735319/) · 🪝 [15](https://www.semanticscholar.org/paper/55a244823d8fd8527819dbd02ffafdc4e661a795) · ⭐ [28](https://github.com/WansonChoi/HATK/stargazers) · `Python`

- [**MATER: Minimizer RNAseq HLA typer**](https://github.com/genentech/midasHLA) — MATER is a minimizer-based HLA typer for RNAseq read dataset. In a typical RNAseq dataset, the reads sampled from HLA genes are less uniform and may miss regions that makes assembly or variant call...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/34228721/) · 🪝 [24](https://www.semanticscholar.org/paper/16f5a5d9b7119ea85eab9aae4650445770a03b3b) · ⭐ [14](https://github.com/genentech/midasHLA/stargazers) · `Python` `R` `C`

- [**PyHLA**](https://github.com/felixfan/PyHLA) — Python for HLA analysis: summary, association analysis, zygosity test and interaction test<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28166716/) · ⭐ [38](https://github.com/felixfan/PyHLA/stargazers) · `Python`

- [**cdr3-QTL**](https://github.com/immunogenomics/cdr3-QTL) — Trans-association between HLA and TCR-CDR3<br>⭐ [19](https://github.com/immunogenomics/cdr3-QTL/stargazers) · `HTML`

- [**hlabud: HLA genotype analysis in R**](https://github.com/slowkow/hlabud) — hlabud provides methods to retrieve sequence alignment data from IMGTHLA and convert the data into convenient R matrices ready for downstream analysis. See the usage examples to learn how to use th...<br>⭐ [17](https://github.com/slowkow/hlabud/stargazers) · `R`

### HLA Typing

- [**OptiType: Precision HLA typing from next-generation sequencing data**](https://github.com/FRED-2/OptiType) — OptiType is a novel HLA genotyping algorithm based on integer linear programming, capable of producing accurate 4-digit HLA genotyping predictions from NGS data by simultaneously selecting all majo...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/25143287/) · 🪝 [643](https://www.semanticscholar.org/paper/1d4162253d3e32a2b3a62d0f8faab4fba3386c10) · ⭐ [205](https://github.com/FRED-2/OptiType/stargazers) · `Python`

- [**arcasHLA: Fast and accurate in silico inference of HLA genotypes from RNA-seq**](https://github.com/RabadanLab/arcasHLA) — arcasHLA performs high resolution genotyping for HLA class I and class II genes from RNA sequencing, supporting both paired and single-end samples.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31173059/) · 🪝 [157](https://www.semanticscholar.org/paper/9fccbea05592fb2c8d0cf8ff0fe330729ad81db8) · ⭐ [152](https://github.com/RabadanLab/arcasHLA/stargazers) · `Python`

- [**xHLA: Fast and accurate HLA typing from short read sequence data**](https://github.com/humanlongevity/HLA) — xHLA iteratively refines the mapping results at the amino acid level to achieve 99 to 100% 4-digit typing accuracy for both class I and II HLA genes, taking only about 3 minutes to process a 30X wh...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28674023/) · 🪝 [125](https://www.semanticscholar.org/paper/0e7b3c0eb32913f710c93bfe149610bc2d2ce8e3) · ⭐ [113](https://github.com/humanlongevity/HLA/stargazers) · [Homepage](https://pubmed.ncbi.nlm.nih.gov/28674023) · `Python` `R` `Perl` `Bash`

- [**HLA-LA: Fast HLA type inference from whole-genome data**](https://github.com/DiltheyLab/HLA-LA) — HLA typing based on a population reference graph and employs a new linear projection method to align reads to the graph.<br>⭐ [141](https://github.com/DiltheyLab/HLA-LA/stargazers) · `Perl`

- [**Kourami: Graph-guided assembly for HLA alleles**](https://github.com/Kingsford-Group/kourami) — Kourami is a graph-guided assembler for HLA haplotypes covering typing exons (exons 2 and 3 for Class I and exon 3 for Class II) using high-coverage whole genome sequencing data. Kourami constructs...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/29415772/) · 🪝 [82](https://www.semanticscholar.org/paper/42b1d2f74a6fcdd9e45b64a4acabddafb43d7426) · ⭐ [38](https://github.com/Kingsford-Group/kourami/stargazers) · `Java` `Bash`

- [**T1K: efficient and accurate inference of KIR or HLA alleles from RNA-seq, whole-genome sequencing, or whole-exome sequencing data**](https://github.com/mourisl/T1K) — T1K (The ONE genotyper for Kir and HLA) is a computational tool to infer the alleles for the polymorphic genes such as KIR and HLA. T1K calculates the allele abundances based on the RNA-seq/WES/WGS...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/37169596/) · 🪝 [25](https://www.semanticscholar.org/paper/c8ff85a07e0dc87973ba73daccc10731f225a914) · ⭐ [94](https://github.com/mourisl/T1K/stargazers) · `C` `C++` `Python` `Perl`

- [**scHLAcount**](https://github.com/10XGenomics/scHLAcount) — Count HLA alleles in single-cell RNA-seq data<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32330223/) · 🪝 [20](https://www.semanticscholar.org/paper/f38834e9b0845e743e8e41604c0ee9d327d3fe48) · ⭐ [63](https://github.com/10XGenomics/scHLAcount/stargazers) · `TeX`

- [**hlatyping**](https://github.com/nf-core/hlatyping) — Precision HLA typing from next-generation sequencing data<br>⭐ [76](https://github.com/nf-core/hlatyping/stargazers) · [Homepage](https://nf-co.re/hlatyping) · `Nextflow`

- [**HLAProfiler: Using k-mers to call HLA alleles in RNA sequencing data**](https://github.com/ExpressionAnalysis/HLAProfiler) — HLAProfiler uses the k-mer content of next generation sequencing reads to call HLA types in a sample. Based on the k-mer content each each read pair is assigned to an HLA gene and the aggregate k-m...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/28954626/) · 🪝 [44](https://www.semanticscholar.org/paper/9781ba0fad5a5841a506d25b8e3a9328996ebe52) · ⭐ [23](https://github.com/ExpressionAnalysis/HLAProfiler/stargazers) · `Perl`

- [**SpecHLA**](https://github.com/deepomicslab/SpecHLA) — SpecHLA reconstructs entire diploid sequences of HLA genes and infers LOH events. It supports HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, and -DRB1 genes. Also, it supports both short- and long-read...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/37714157/) · 🪝 [11](https://www.semanticscholar.org/paper/78feeabf537fe7503ff5fe2421b2668f635fa504) · ⭐ [53](https://github.com/deepomicslab/SpecHLA/stargazers) · `C++`

- [**seq2HLA: HLA typing from RNA-Seq sequence reads**](https://github.com/TRON-Bioinformatics/seq2HLA) — In-silico method written in Python and R to determine HLA genotypes of a sample. seq2HLA takes standard RNA-Seq sequence reads in fastq format as input, uses a bowtie index comprising all HLA allel...<br>⭐ [50](https://github.com/TRON-Bioinformatics/seq2HLA/stargazers) · `Python` `R`

- [**PHLAT: Inference of High Resolution HLA Types**](https://sites.google.com/site/phlatfortype/home) — PHLAT is a bioinformatics algorithm that offers HLA typing at four-digit resolution (or higher) using genome-wide transcriptome and exome sequencing data over a wide range of read lengths and seque...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/29858810/) · 🪝 [46](https://www.semanticscholar.org/paper/6d09c67998833e831897289ba1ea08efd0b81fee) · [Homepage](https://sites.google.com/site/phlatfortype/home) · `Python`

- [**MultiHLA: WES HLA Typing based on multiple alternative tools**](https://github.com/lkuchenb/MultiHLA) — This workflow enables the concurrent analysis of WES or WGS data using publicly available software to derive HLA haplotypes from this type of data. It includes automated Snakemake workflows for the...<br>⭐ [18](https://github.com/lkuchenb/MultiHLA/stargazers) · `Snakemake`

- [**hla3: weight of evidence of HLA allele expression based on bulk TCR beta-chain repertoires**](https://github.com/kmayerb/hla3) — This repository contains Python functions for inferring HLA-alleles from bulk TCR beta chain data using a simple weight of evidence predictor.<br>⭐ [3](https://github.com/kmayerb/hla3/stargazers) · `Python`

- [**SNP2HLA: Imputation of Amino Acid Polymorphisms in Human Leukocyte Antigens**](http://software.broadinstitute.org/mpg/snp2hla/) — SNP2HLA is a tool to impute amino acid polymorphisms and single nucleotide polymorphisms in human luekocyte antigenes (HLA) within the major histocompatibility complex (MHC) region in chromosome 6.<br>[Homepage](http://software.broadinstitute.org/mpg/snp2hla/)

### Peptide Prediction

- [**HLAMatchmaker**](http://www.epitopes.net) — A molecularly based algorithm for histocompatibility determination that identifies acceptable HLA antigens for highly alloimmunized patients based on amino acid triplets (eplets) on exposed parts o...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/11975978/) · 🪝 [267](https://www.semanticscholar.org/paper/7a00824f5126ab1433ac8fdcfba4dab4854ab3b2) · [Homepage](http://www.epitopes.net)

- [**High-Throughput Prediction of MHC Class I and II Neoantigens with MHCnuggets**](https://github.com/KarchinLab/mhcnuggets) — MHC Class I and Class II neoantigen binding prediction<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31871119/) · 🪝 [131](https://www.semanticscholar.org/paper/33d23cc483e4b077b1f637444b10e98cb1f6bab7) · ⭐ [33](https://github.com/KarchinLab/mhcnuggets/stargazers) · `Python`

- [**NeoBert**](https://github.com/CHB-learner/NeoBert) — NeoBERT is an advanced model designed specifically for predicting the binding affinity between neoantigens and HLA. It is a variant of the original BERT model, enhanced to integrate biological feat...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/41224698/) · ⭐ [155](https://github.com/CHB-learner/NeoBert/stargazers) · `Python`

- [**bigmhc**](https://github.com/KarchinLab/bigmhc) — BigMHC predicts MHC-I (neo)epitope presentation and immunogenicity<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/37829001/) · 🪝 [60](https://www.semanticscholar.org/paper/ef7763384b5f987dc546cdd9ece14b3e81b89190) · ⭐ [59](https://github.com/KarchinLab/bigmhc/stargazers) · `Jupyter Notebook`

- [**HLA-EMMA**](https://hla-emma.lumc.nl) — A user-friendly tool to analyze HLA class I and class II compatibility on the amino acid level, facilitating the assessment of donor-recipient compatibility in transplantation.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32227681/) · 🪝 [82](https://www.semanticscholar.org/paper/2739806ab693e3eab964e75041527569045cb62c) · [Homepage](https://hla-emma.lumc.nl)

- [**PIRCHE-II**](https://www.pirche.com) — An algorithm to predict indirectly recognizable HLA epitopes in solid organ transplantation, helping to evaluate immunological compatibility between donors and recipients.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/31741009/) · 🪝 [81](https://www.semanticscholar.org/paper/75defe0480b1dbd3c64501d05370a69ea3a6b260) · [Homepage](https://www.pirche.com)

- [**MHCAttnNet**](https://github.com/gopuvenkat/MHCAttnNet) — MHCAttnNet: Allele-Peptide predictions for class I & class II MHC alleles<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/32657386/) · 🪝 [39](https://www.semanticscholar.org/paper/64fd328e9f126c6277e2ab50f4a4b86be9bfda94) · ⭐ [30](https://github.com/gopuvenkat/MHCAttnNet/stargazers) · `Python`

- [**MixMHC2pred**](https://github.com/GfellerLab/MixMHC2pred) — HLA-II ligand predictor.<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/38907900/) · ⭐ [46](https://github.com/GfellerLab/MixMHC2pred/stargazers) · `C++`

- [**MixMHCpred**](https://github.com/GfellerLab/MixMHCpred) — HLA-I ligand predictor<br>⭐ [43](https://github.com/GfellerLab/MixMHCpred/stargazers) · `Python`

- [**EpVix: epitope reactivity analysis and epitope virtual crossmatching**](https://pubmed.ncbi.nlm.nih.gov/26531328/) — Performs automated epitope virtual crossmatching at the initiation of the organ donation process. EpViX is a free, web-based application developed for use over the internet on a tablet, smartphone ...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/26531328/) · 🪝 [33](https://www.semanticscholar.org/paper/f87da5d55ad244c92cacc84804775bd1e5df2fd0) · `Ruby`

- [**immunogenetr**](https://github.com/k96nb01/immunogenetr_package) — immunogenetr is a comprehensive toolkit for clinical HLA informatics. It is built on tidyverse principles and makes use of genotype list string (GL string, https://glstring.org/) for storing and us...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/41344288/) · ⭐ [6](https://github.com/k96nb01/immunogenetr_package/stargazers) · [Homepage](https://glstring.org/) · `R`

### Data & Nomenclature

- [**MHC-PRG**](https://github.com/AlexanderDilthey/MHC-PRG) — Population Reference Graphs for the HLA and MHC.<br>⭐ [35](https://github.com/AlexanderDilthey/MHC-PRG/stargazers) · `C++`

- [**py-ard**](https://github.com/nmdp-bioinformatics/py-ard) — HLA ARD Reduction in Python. Although HLA nomenclature has not always conformed to the same standard, it is now defined by The WHO Nomenclature Committee for Factors of the HLA System. py-ard is aw...<br>⭐ [19](https://github.com/nmdp-bioinformatics/py-ard/stargazers) · `Python`

- [**HLAtools: Functions and Datasets for HLA Informatics**](https://github.com/sjmack/HLAtools) — We have developed HLAtools, an R package that automates the consumption of IPD-IMGT/HLA resources, renders them computable, and makes them available alongside tools for data analysis, visualization...<br>[PubMed](https://pubmed.ncbi.nlm.nih.gov/40947766/) · 🪝 [1](https://www.semanticscholar.org/paper/a47f89a247c3149305dba16cfcbbd94b66810d49) · ⭐ [4](https://github.com/sjmack/HLAtools/stargazers) · `R`


---
