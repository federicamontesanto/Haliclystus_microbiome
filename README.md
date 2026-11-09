# Haliclystus_microbiome
This repository contains the R analysis code and sample metadata for the microbiome study of the Antarctic stalked jellyfish *Haliclystus antarcticus*.

The study investigates bacterial community variation across three anatomical compartments (gut, peduncle, and tentacles) and two Antarctic sampling sites (Byers Bay and Risopatron).

The workflow includes preprocessing of 16S rRNA gene amplicon sequencing data, taxonomic assignment, contaminant identification, alpha and beta diversity analyses, bacterial community composition, differential abundance testing, indicator taxa identification, co-abundance analyses, and figure generation.

Additional scripts support the genomic and phylogenomic characterization of the bacterial isolate FM48a.

## Workflow

### 01 — Pre-processing

Quality filtering, denoising, ASV inference using DADA2, taxonomy assignment with SILVA v138.1, and contaminant identification using the *decontam* package.

### 02 — Phyloseq and diversity preparation

Construction of phyloseq objects, phylogenetic tree inference, and preparation of bacterial community data for downstream analyses.

### 03 — Alpha diversity and bacterial composition

Calculation and statistical analysis of alpha diversity metrics, including Observed richness, Fisher's alpha, Shannon diversity, and Inverse Simpson diversity.

Analysis of bacterial community composition, ASV sharing across anatomical compartments, and differential abundance.

### 04 — Beta diversity

Bray–Curtis and UniFrac distance analyses, principal coordinates analysis (PCoA), PERMANOVA, and dispersion assessment.

### 05 — Indicator taxa

Identification of bacterial taxa associated with anatomical compartments and sampling sites using marker and indicator analyses.

### 06 — Co-abundance analysis

Exploratory analysis of bacterial co-abundance patterns among prevalent bacterial genera.

### 07 — FM48a genomic analysis

Genome-based taxonomic identification using average nucleotide identity (ANI), phylogenomic inference, and visualization of the phylogenetic relationships between FM48a and reference *Oceanisphaera* genomes.

Sequences submitted: BioProject submission PRJNA1558909.
