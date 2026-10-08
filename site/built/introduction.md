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
```

``` output
── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.2.1     ✔ readr     2.2.0
✔ forcats   1.0.1     ✔ stringr   1.6.0
✔ ggplot2   4.0.3     ✔ tibble    3.3.1
✔ lubridate 1.9.5     ✔ tidyr     1.3.2
✔ purrr     1.2.2     
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(readxl)
library(reshape2)
```

``` output

Attaching package: 'reshape2'

The following object is masked from 'package:tidyr':

    smiths
```

``` r
# Methylation analysis and annotation
library(minfi)
```

``` output
Loading required package: BiocGenerics
Loading required package: generics

Attaching package: 'generics'

The following object is masked from 'package:lubridate':

    as.difftime

The following object is masked from 'package:dplyr':

    explain

The following objects are masked from 'package:base':

    as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    setequal, union


Attaching package: 'BiocGenerics'

The following object is masked from 'package:dplyr':

    combine

The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs

The following objects are masked from 'package:base':

    anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    unsplit, which.max, which.min

Loading required package: GenomicRanges
Loading required package: stats4
Loading required package: S4Vectors

Attaching package: 'S4Vectors'

The following objects are masked from 'package:lubridate':

    second, second<-

The following objects are masked from 'package:dplyr':

    first, rename

The following object is masked from 'package:tidyr':

    expand

The following object is masked from 'package:utils':

    findMatches

The following objects are masked from 'package:base':

    expand.grid, I, unname

Loading required package: IRanges

Attaching package: 'IRanges'

The following object is masked from 'package:lubridate':

    %within%

The following objects are masked from 'package:dplyr':

    collapse, desc, slice

The following object is masked from 'package:purrr':

    reduce

Loading required package: Seqinfo
Loading required package: SummarizedExperiment
Loading required package: MatrixGenerics
Loading required package: matrixStats

Attaching package: 'matrixStats'

The following object is masked from 'package:dplyr':

    count


Attaching package: 'MatrixGenerics'

The following objects are masked from 'package:matrixStats':

    colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    colWeightedMeans, colWeightedMedians, colWeightedSds,
    colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    rowWeightedSds, rowWeightedVars

Loading required package: Biobase
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.


Attaching package: 'Biobase'

The following object is masked from 'package:MatrixGenerics':

    rowMedians

The following objects are masked from 'package:matrixStats':

    anyMissing, rowMedians

Loading required package: Biostrings
Loading required package: XVector

Attaching package: 'XVector'

The following object is masked from 'package:purrr':

    compact


Attaching package: 'Biostrings'

The following object is masked from 'package:base':

    strsplit

Loading required package: bumphunter
Loading required package: foreach

Attaching package: 'foreach'

The following objects are masked from 'package:purrr':

    accumulate, when

Loading required package: iterators
Loading required package: parallel
Loading required package: locfit
locfit 1.5-9.12 	 2025-03-05

Attaching package: 'locfit'

The following object is masked from 'package:purrr':

    none

Setting options('download.file.method.GEOquery'='auto')
Setting options('GEOquery.inmemory.gpl'=FALSE)
```

``` r
library(GEOquery)
library(limma)
```

``` output

Attaching package: 'limma'

The following object is masked from 'package:BiocGenerics':

    plotMA
```

``` r
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(AnnotationHub)
```

``` output
Loading required package: BiocFileCache
Loading required package: dbplyr

Attaching package: 'dbplyr'

The following objects are masked from 'package:dplyr':

    ident, sql, sql_escape_ident, sql_escape_string


Attaching package: 'AnnotationHub'

The following object is masked from 'package:Biobase':

    cache
```

``` r
library(DMRcate)
```

``` output

```

``` r
library(methylKit)
```

``` output

Attaching package: 'methylKit'

The following object is masked from 'package:dplyr':

    select

The following object is masked from 'package:tidyr':

    unite
```

``` r
library(missMethyl)
```

``` output
Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19

Attaching package: 'IlluminaHumanMethylation450kanno.ilmn12.hg19'

The following objects are masked from 'package:IlluminaHumanMethylationEPICv2anno.20a1.hg38':

    Islands.UCSC, Locations, Manifest, Other, SNPs.141CommonSingle,
    SNPs.142CommonSingle, SNPs.144CommonSingle, SNPs.146CommonSingle,
    SNPs.147CommonSingle, SNPs.Illumina

Loading required package: IlluminaHumanMethylationEPICanno.ilm10b4.hg19

Attaching package: 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'

The following objects are masked from 'package:IlluminaHumanMethylation450kanno.ilmn12.hg19':

    Islands.UCSC, Locations, Manifest, Other, SNPs.132CommonSingle,
    SNPs.135CommonSingle, SNPs.137CommonSingle, SNPs.138CommonSingle,
    SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle,
    SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.Illumina

The following objects are masked from 'package:IlluminaHumanMethylationEPICv2anno.20a1.hg38':

    Islands.UCSC, Locations, Manifest, Other, SNPs.141CommonSingle,
    SNPs.142CommonSingle, SNPs.144CommonSingle, SNPs.146CommonSingle,
    SNPs.147CommonSingle, SNPs.Illumina
```

``` r
library(clusterProfiler)
```

``` output

clusterProfiler v4.18.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

Please cite:

S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
multiomics data. Nature Protocols. 2024, 19(11):3292-3320

Attaching package: 'clusterProfiler'

The following object is masked from 'package:methylKit':

    select

The following object is masked from 'package:XVector':

    slice

The following object is masked from 'package:IRanges':

    slice

The following object is masked from 'package:S4Vectors':

    rename

The following object is masked from 'package:purrr':

    simplify

The following object is masked from 'package:stats':

    filter
```

``` r
library(org.Hs.eg.db)
```

``` output
Loading required package: AnnotationDbi

Attaching package: 'AnnotationDbi'

The following object is masked from 'package:clusterProfiler':

    select

The following object is masked from 'package:methylKit':

    select

The following object is masked from 'package:dplyr':

    select
```

``` r
# Plotting and parallel computing
library(Polychrome)
library(ComplexHeatmap)
```

``` output
Loading required package: grid

Attaching package: 'grid'

The following object is masked from 'package:Biostrings':

    pattern

========================================
ComplexHeatmap version 2.26.1
Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
Github page: https://github.com/jokergoo/ComplexHeatmap
Documentation: http://jokergoo.github.io/ComplexHeatmap-reference

If you use it in published research, please cite either one:
- Gu, Z. Complex Heatmap Visualization. iMeta 2022.
- Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    genomic data. Bioinformatics 2016.


The new InteractiveComplexHeatmap package can directly export static 
complex heatmaps into an interactive Shiny app with zero effort. Have a try!

This message can be suppressed by:
  suppressPackageStartupMessages(library(ComplexHeatmap))
========================================
```

``` r
library(ggsci)
library(gridExtra)
```

``` output

Attaching package: 'gridExtra'

The following object is masked from 'package:minfi':

    combine

The following object is masked from 'package:Biobase':

    combine

The following object is masked from 'package:BiocGenerics':

    combine

The following object is masked from 'package:dplyr':

    combine
```

``` r
library(RColorBrewer)
library(parallel)
```

## Load annotation data


``` r
epic <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(epic)
```

``` output
DataFrame with 6 rows and 42 columns
                        chr       pos      strand            Name    AddressA
                <character> <integer> <character>     <character> <character>
cg25324105_BC11       chr19  37692358           + cg25324105_BC11     1754126
cg25383568_TC11       chr19  38727081           - cg25383568_TC11    79792482
cg25455143_BC11       chr19   1591515           - cg25455143_BC11    80699190
cg25459778_BC11       chr16   1614581           + cg25459778_BC11    60797262
cg25487775_BC11        chr2 161237458           + cg25487775_BC11     5799427
cg25595446_BC11       chr19  49677322           + cg25595446_BC11    65640459
                   AddressB              ProbeSeqA              ProbeSeqB
                <character>            <character>            <character>
cg25324105_BC11    99753217 ATTTATAAACTAATAACCCA.. GTTTATAAACTAATAACCCG..
cg25383568_TC11    69667133 AAACCAAAAAAATAACAAAC.. AAACCGAAAAAATAACAAAC..
cg25455143_BC11     7659147 ATAAAAAAAAATATACAACT.. ATAAAAAAAAATATACGACT..
cg25459778_BC11    65710482 AAAAATTTAAAACAAACAAC.. AAAAATTTAAAACAAACAAC..
cg25487775_BC11    89606481 AAAAACAACCTAAAAAACAA.. AAAAACAACCTAAAAAACAA..
cg25595446_BC11    39619855 AATAAAAATAACAACAACCA.. AATAAAAATAACGACAACCG..
                       Type    NextBase       Color    Probe_rs Probe_maf
                <character> <character> <character> <character> <numeric>
cg25324105_BC11           I           A         Red          NA        NA
cg25383568_TC11           I           C         Grn          NA        NA
cg25455143_BC11           I           T         Red          NA        NA
cg25459778_BC11           I           C         Grn          NA        NA
cg25487775_BC11           I           A         Red          NA        NA
cg25595446_BC11           I           C         Grn          NA        NA
                     CpG_rs   CpG_maf      SBE_rs   SBE_maf
                <character> <numeric> <character> <numeric>
cg25324105_BC11          NA        NA          NA        NA
cg25383568_TC11          NA        NA          NA        NA
cg25455143_BC11          NA        NA          NA        NA
cg25459778_BC11          NA        NA          NA        NA
cg25487775_BC11          NA        NA          NA        NA
cg25595446_BC11          NA        NA          NA        NA
                          Islands_Name Relation_to_Island         col
                           <character>        <character> <character>
cg25324105_BC11 chr19:37691892-37692..             Island           R
cg25383568_TC11 chr19:38726890-38727..             Island           G
cg25455143_BC11 chr19:1591428-159162..             Island           R
cg25459778_BC11  chr16:1610053-1615094             Island           G
cg25487775_BC11 chr2:161237949-16123..              Shore           R
cg25595446_BC11 chr19:49676753-49677..             Island           G
                 Probe_Type   Strand_FR   Strand_TB   Strand_CO Infinium_Design
                <character> <character> <character> <character>     <character>
cg25324105_BC11          cg           F           B           C               1
cg25383568_TC11          cg           R           T           C               1
cg25455143_BC11          cg           R           B           C               1
cg25459778_BC11          cg           F           B           C               1
cg25487775_BC11          cg           F           B           C               1
cg25595446_BC11          cg           F           B           C               1
                    Species     Rep_Num     UCSC_RefGene_Group
                <character> <character>            <character>
cg25324105_BC11       Human           1 TSS200;TSS200;TSS200..
cg25383568_TC11       Human           1        exon_18;exon_18
cg25455143_BC11       Human           1                       
cg25459778_BC11       Human           1                       
cg25487775_BC11       Human           1                       
cg25595446_BC11       Human           1                       
                     UCSC_RefGene_Name UCSC_RefGene_Accession
                           <character>            <character>
cg25324105_BC11 ZNF781;ZNF781;ZNF781.. NR_173332.1;NR_17333..
cg25383568_TC11            ACTN4;ACTN4 NM_001322033.2_2;NM_..
cg25455143_BC11                                              
cg25459778_BC11                                              
cg25487775_BC11                                              
cg25595446_BC11                                              
                      GencodeV41_Group        GencodeV41_Name
                           <character>            <character>
cg25324105_BC11 TSS200;TSS200;TSS200.. ENSG00000120784.17;E..
cg25383568_TC11 exon_18;exon_18;exon.. ACTN4;ACTN4;ACTN4;AC..
cg25455143_BC11                                              
cg25459778_BC11                 TSS200     ENSG00000007545.16
cg25487775_BC11                                              
cg25595446_BC11  exon_1;TSS200;TSS1500 PRMT1;ENSG0000012645..
                  GencodeV41_Accession Phantom5_Enhancers  HMM_Island
                           <character>        <character> <character>
cg25324105_BC11 ENST00000587199.5;EN..                               
cg25383568_TC11 ENST00000252699.7;EN..                               
cg25455143_BC11                                                      
cg25459778_BC11      ENST00000293925.9                               
cg25487775_BC11                                                      
cg25595446_BC11 ENST00000524771.5;EN..                               
                Regulatory_Feature_Group Methyl450_Enhancer         DMR
                             <character>        <character> <character>
cg25324105_BC11   Unclassified_Cell_ty..               TRUE         DMR
cg25383568_TC11   Gene_Associated_Cell..               TRUE            
cg25455143_BC11                                       FALSE            
cg25459778_BC11                                       FALSE            
cg25487775_BC11                                       FALSE            
cg25595446_BC11      Promoter_Associated               TRUE            
                Methyl450_Loci Methyl27_Loci EPICv1_Loci Manifest_probe_match
                   <character>   <character> <character>          <character>
cg25324105_BC11     cg25324105                cg25324105                 TRUE
cg25383568_TC11     cg25383568                cg25383568                 TRUE
cg25455143_BC11     cg25455143                cg25455143                 TRUE
cg25459778_BC11     cg25459778    cg25459778  cg25459778                 TRUE
cg25487775_BC11     cg25487775                cg25487775                 TRUE
cg25595446_BC11     cg25595446    cg25595446  cg25595446                 TRUE
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

``` output
   Accession        Name  Type                  Name2
1 GSM7698438 LNCAP_500_1 LNCaP GSM7698438_LNCAP_500_1
2 GSM7698446 LNCAP_500_2 LNCaP GSM7698446_LNCAP_500_2
3 GSM7698462 LNCAP_500_3 LNCaP GSM7698462_LNCAP_500_3
4 GSM7698435  PREC_500_1  PREC  GSM7698435_PREC_500_1
5 GSM7698443  PREC_500_2  PREC  GSM7698443_PREC_500_2
6 GSM7698459  PREC_500_3  PREC  GSM7698459_PREC_500_3
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
suffix.


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

``` output
class: RGChannelSet 
dim: 1105209 6 
metadata(0):
assays(2): Green Red
rownames(1105209): 1600157 1600179 ... 99810982 99810990
rowData names(0):
colnames(6): LNCAP_500_1 LNCAP_500_2 ... PREC_500_2 PREC_500_3
colData names(7): Accession Name ... Basename filenames
Annotation
  array: IlluminaHumanMethylationEPICv2
  annotation: 20a1.hg38
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
