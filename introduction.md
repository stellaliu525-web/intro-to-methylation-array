---
title: "Introduction: loading packages and downloading data"
teaching: 15
exercises: 5
---

:::::::::::::::::::::::::::::::::::::: questions

- FIXME
- FIXME

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Load the required R packages.
- Create a sample table and define the variable of interest.
- Download the six tutorial samples from GEO.
- Import the IDAT files into a `minfi` object.

::::::::::::::::::::::::::::::::::::::::::::::::

## Load packages

Complete the [setup instructions](../learners/setup.md) first. Open RStudio and
create a project using **File > New Project > New Directory > New Project**,
or open your existing workshop project. Create an R script, and run the chunks below in order.

Installing a package downloads it to your computer. Loading it with `library()`
makes its functions available in the current session. Load the packages again
whenever you restart R.


``` r
# Data handling
library(tidyverse)
library(readxl)
library(reshape2)

# Methylation analysis and annotation
library(minfi)
library(GEOquery)
library(limma)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(AnnotationHub)
library(DMRcate)
library(methylKit)
library(missMethyl)
library(clusterProfiler)
library(org.Hs.eg.db)

# Plotting and parallel computing
library(Polychrome)
library(ComplexHeatmap)
library(ggsci)
library(gridExtra)
library(RColorBrewer)
library(parallel)
```

## Load annotation data


``` r
epic <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(epic)
```

## Download data

For this tutorial we are using the publicly available LNCaP and PREC cell lines
(500 ng DNA input) from Peters *et al.* (2024), available from GEO
([GSE240469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240469)).
We will download these data in the chunk below, so you do not have to manually
get them from the above link.

If you would like to use your own data, create a sample table and variables of
interest in the same format, with your own IDAT files, samples and variables.
The GEO download chunk is only needed for data hosted on GEO.

### Create the sample table

Each row describes one sample. `Accession` identifies the GEO sample, `Name`
matches its filename label, and `Type` identifies the cell line to compare.


``` r
# Create the sample table for IDAT downloads. Replace this table with your own sample information when using your own data.

targets <- data.frame(
  Accession = c(
    "GSM7698438",
    "GSM7698446",
    "GSM7698462",
    "GSM7698435",
    "GSM7698443",
    "GSM7698459"
  ),
  Name = c(
    "LNCAP_500_1",
    "LNCAP_500_2",
    "LNCAP_500_3",
    "PREC_500_1",
    "PREC_500_2",
    "PREC_500_3"
  ),
  Type = c("LNCaP", "LNCaP", "LNCaP", "PREC", "PREC", "PREC")
)

targets$Name2 <- paste(
  targets$Accession, targets$Name, sep = "_"
)

# Define the variable of interest and convert it to a factor.
variable.of.choice <- as.factor(targets$Type)

targets
```

### Download the IDAT files

The following chunk downloads the IDAT files for the six accessions in
`targets`. Each sample is saved in its own accession-named folder within
your current project directory.


``` r
# Download only IDAT files for the selected tutorial samples.
gsmlist <- lapply(
  targets$Accession,
  GEOquery::getGEOSuppFiles,
  makeDirectory = TRUE,
  baseDir = getwd(),
  filter_regex = "[.]idat([.]gz)?$"
)
names(gsmlist) <- targets$Accession
```

## Read in .idats

Add the sample names and file paths to `targets`, preserving its sample order.
`Basename` is the path to each file pair without the `_Grn.idat` or `_Red.idat`
suffix (or their compressed `.gz` equivalents).


``` r
# Link the sample information to the downloaded IDAT files.
targets$Sample <- targets$Name
targets$Basename <- file.path(
  getwd(), targets$Accession, targets$Name2
)

# Read in raw data from the red/green IDAT files.
rgSet <- read.metharray.exp(targets = targets)
sampleNames(rgSet) <- targets$Sample
rgSet@annotation <- c(
  array = "IlluminaHumanMethylationEPICv2",
  annotation = "20a1.hg38"
)
rgSet
```

`rgSet` contains the raw intensities for the six samples. For your own data,
set `targets$Basename` to your file paths and use the annotation matching your
array type. Quality control, filtering and normalisation follow in the next
episode.

::::::::::::::::::::::::::::::::::::: keypoints

- Raw Illumina methylation array data are stored in `.idat` files, with one red-channel and one green-channel file per sample.
- The `targets` data frame holds the sample metadata, linking each sample's name and cell-line group to its IDAT file paths.
- `read.metharray.exp(targets = targets)` reads the raw intensities into `rgSet`, an `RGChannelSet` object, ready for quality control and preprocessing.

::::::::::::::::::::::::::::::::::::::::::::::::
