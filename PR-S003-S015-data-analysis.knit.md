---
title: "PR-S003-S015_data-analysis"
author: "Jenna Hanlon"
date: "2024-01-11"
output: html_document
---

# Load Libraries


```r
library(qiime2R) 
library(ggplot2)
library(vegan) 
```

```
## Warning: package 'vegan' was built under R version 4.1.3
```

```
## Loading required package: permute
```

```
## Warning: package 'permute' was built under R version 4.1.3
```

```
## Loading required package: lattice
```

```
## Warning: package 'lattice' was built under R version 4.1.3
```

```
## This is vegan 2.6-4
```

```r
library(plyr) 
```

```
## Warning: package 'plyr' was built under R version 4.1.3
```

```r
library(dplyr) 
```

```
## Warning: package 'dplyr' was built under R version 4.1.3
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(scales)
```

```
## Warning: package 'scales' was built under R version 4.1.3
```

```r
library(grid) 
library(reshape2)
```

```
## Warning: package 'reshape2' was built under R version 4.1.3
```

```r
library(phyloseq) 
library(picante) 
```

```
## Warning: package 'picante' was built under R version 4.1.3
```

```
## Loading required package: ape
```

```
## Warning: package 'ape' was built under R version 4.1.3
```

```
## 
## Attaching package: 'ape'
```

```
## The following object is masked from 'package:dplyr':
## 
##     where
```

```
## Loading required package: nlme
```

```
## Warning: package 'nlme' was built under R version 4.1.3
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```r
library(tidyr) 
```

```
## Warning: package 'tidyr' was built under R version 4.1.3
```

```
## 
## Attaching package: 'tidyr'
```

```
## The following object is masked from 'package:reshape2':
## 
##     smiths
```

```r
library(viridis) 
```

```
## Loading required package: viridisLite
```

```
## 
## Attaching package: 'viridis'
```

```
## The following object is masked from 'package:scales':
## 
##     viridis_pal
```

```r
library(qiime2R) 
library(DESeq2) 
```

```
## Loading required package: S4Vectors
```

```
## Warning: package 'S4Vectors' was built under R version 4.1.3
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:plyr':
## 
##     rename
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following object is masked from 'package:nlme':
## 
##     collapse
```

```
## The following object is masked from 'package:phyloseq':
## 
##     distance
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:plyr':
## 
##     desc
```

```
## The following object is masked from 'package:grDevices':
## 
##     windows
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## Warning: package 'matrixStats' was built under R version 4.1.3
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## The following object is masked from 'package:plyr':
## 
##     count
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:phyloseq':
## 
##     sampleNames
```

```r
library(patchwork) 
library(RColorBrewer)
```

```
## Warning: package 'RColorBrewer' was built under R version 4.1.3
```

```r
library(microViz) 
```

```
## microViz version 0.10.10 - Copyright (C) 2023 David Barnett
## ! Website: https://david-barnett.github.io/microViz
## v Useful?  For citation details, run: `citation("microViz")`
## x Silence? `suppressPackageStartupMessages(library(microViz))`
```

```r
library(speedyseq) 
```

```
## 
## Attaching package: 'speedyseq'
```

```
## The following objects are masked from 'package:phyloseq':
## 
##     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom,
##     tip_glom, transform_sample_counts
```

```r
library(ComplexHeatmap) 
```

```
## ========================================
## ComplexHeatmap version 2.10.0
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite:
## Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##   genomic data. Bioinformatics 2016.
## 
## The new InteractiveComplexHeatmap package can directly export static 
## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(ComplexHeatmap))
## ========================================
```

```r
library(ggVennDiagram) 
```

```
## Warning: package 'ggVennDiagram' was built under R version 4.1.3
```

```r
library(SuperExactTest) 
```

```
## 
## Attaching package: 'SuperExactTest'
```

```
## The following objects are masked from 'package:GenomicRanges':
## 
##     intersect, union
```

```
## The following object is masked from 'package:GenomeInfoDb':
## 
##     intersect
```

```
## The following objects are masked from 'package:IRanges':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:dplyr':
## 
##     intersect, union
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, union
```

```r
library(nVennR)  
library(vegan) 
library(plyr) 
library(dplyr)
library(scales) 
library(grid) 
library(reshape2)
library(phyloseq) 
library(picante) 
library(tidyr) 
library(viridis) 
library(qiime2R)
library(DESeq2) 
library(patchwork) 
library(RColorBrewer)
library(speedyseq) 
library(RColorBrewer)
library(microViz) 
library(ComplexHeatmap) 
library(plotly) 
```

```
## 
## Attaching package: 'plotly'
```

```
## The following object is masked from 'package:ComplexHeatmap':
## 
##     add_heatmap
```

```
## The following object is masked from 'package:microViz':
## 
##     add_paths
```

```
## The following object is masked from 'package:IRanges':
## 
##     slice
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     rename
```

```
## The following objects are masked from 'package:plyr':
## 
##     arrange, mutate, rename, summarise
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
library(seecolor)
library(htmlwidgets)
```

```
## Warning: package 'htmlwidgets' was built under R version 4.1.3
```

```r
library(htmltools)
```

```
## Warning: package 'htmltools' was built under R version 4.1.3
```

```r
set.seed(17384)
```

# Loading Data


```r
PRphyseq <- qza_to_phyloseq(features="S003-S015-table.qza", tree="S003-S015-rooted-tree.qza",taxonomy="S003-S015-taxonomy.qza", metadata = "S003-S015-Field-Measurements-2.txt")
```


```r
rank_names(PRphyseq)
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```


```r
sample_variables(PRphyseq)
```

```
##  [1] "Project"                    "Sub_Location"              
##  [3] "Location"                   "Season"                    
##  [5] "Month"                      "Year"                      
##  [7] "Sample_Type"                "Sample_ID"                 
##  [9] "Day"                        "Station_ID"                
## [11] "Station_Name"               "Collected_by"              
## [13] "Time"                       "Samples_Collected_by"      
## [15] "Temperature_C"              "mm_Hg"                     
## [17] "DO_mg_L"                    "DO_LDO"                    
## [19] "Specific_Conductance_uS_cm" "Salinity_ppt"              
## [21] "pH"                         "pH_mv"                     
## [23] "IVCH_ug_L"                  "Turbidity_.Fluorometer"    
## [25] "Comments"
```

# Removing chloroplast sequences and any contaminant sequences


```r
## Removing chloroplast sequences and any contaminant sequences
PRphyseq <- subset_taxa(PRphyseq, Kingdom != "d__Eukaryota") 
PRphyseq <- subset_taxa(PRphyseq, Kingdom != "d__Archaea") 
PRphyseq <- subset_taxa(PRphyseq, Order != "o__Chloroplast") 
PRphyseq <- subset_taxa(PRphyseq, Order != "Chloroplast") 
PRphyseq <- subset_taxa(PRphyseq, Family != "f__Mitochondria") 
PRphyseq <- subset_taxa(PRphyseq, Family != "Mitochondria") 
PRphyseq <- subset_taxa(PRphyseq, Family != "Chloroplast")  
##Check that contaminant sequences are removed (easiest to save as data frame and search) taxtabl <- as.data.frame(tax_table(S003physeq))  #At this point, blanks and positive controls should also be removed #S003physeq = subset_samples(S003physeq, Sample != "mock")
```


```r
View((tax_table(PRphyseq)))

##This will allow you to open your data frame. You can search in the top right corner of the dataframe for the sequences you wanted to remove in above code to ensure they were removed. 
```

# Creating a subset of samples


```r
PRphyseq <- PRphyseq %>% subset_samples(Project %in% c("S003","S015"))  ##This will create a subset (based on project column) of the samples we would like to use in the analysis. If you would like to keep all samples in the metadata file, you can ignore this step. 


PR_Normal <- PRphyseq %>% subset_samples(Sample_Type %in% c("N"))
```

# Examine the number of reads


```r
 readsumsdf = data.frame(nreads = sort(taxa_sums(PR_Normal), TRUE),                                           sorted = 1:ntaxa(PR_Normal), type = "ASVs") 
 readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(PR_Normal),                          TRUE), sorted = 1:nsamples(PR_Normal), type = "Samples"))
 title = "Total number of reads" 
 p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity") 
 p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
```

```
## Warning: Transformation introduced infinite values in continuous y-axis
```

```
## Warning: Removed 1908 rows containing missing values (`geom_bar()`).
```

<img src="PR-S003-S015-data-analysis_files/figure-html/readsumsdf-1.png" width="672" />

Plot a rarefaction curve using vegan:


```r
taxa_are_rows(PR_Normal)
```

```
## [1] TRUE
```


```r
mat <- t(otu_table(PR_Normal))
class(mat) <- "matrix"
```

```
## Warning in class(mat) <- "matrix": Setting class(x) to "matrix" sets attribute
## to NULL; result will no longer be an S4 object
```

```r
class(mat)
```

```
## [1] "matrix" "array"
```

```r
mat <- as(t(otu_table(PR_Normal)), "matrix")
class(mat)
```

```
## [1] "matrix" "array"
```

```r
raremax <- min(rowSums(mat))

system.time(rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE))
```

<img src="PR-S003-S015-data-analysis_files/figure-html/mat-1.png" width="672" />

```
##    user  system elapsed 
##    4.28    0.17    6.32
```

# Creating table for number of reads per sample


```r
#Sum reads and sort by samples with least to greatest number of reads
samples_reads <- sort(sample_sums(PR_Normal))
#Create table
knitr::kable(samples_reads) %>% 
kableExtra::kable_styling("striped", latex_options="scale_down") %>% 
kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PR109C </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr242a </td>
   <td style="text-align:right;"> 679 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR325 </td>
   <td style="text-align:right;"> 699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR164C </td>
   <td style="text-align:right;"> 700 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR151C </td>
   <td style="text-align:right;"> 855 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR566A </td>
   <td style="text-align:right;"> 1040 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR321 </td>
   <td style="text-align:right;"> 1088 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR147C </td>
   <td style="text-align:right;"> 1128 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR564A </td>
   <td style="text-align:right;"> 1219 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lsj3-1 </td>
   <td style="text-align:right;"> 1251 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR323 </td>
   <td style="text-align:right;"> 1323 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR559A </td>
   <td style="text-align:right;"> 1335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR145C </td>
   <td style="text-align:right;"> 1351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR422 </td>
   <td style="text-align:right;"> 1363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR425 </td>
   <td style="text-align:right;"> 1372 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR162C </td>
   <td style="text-align:right;"> 1450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR141C </td>
   <td style="text-align:right;"> 1461 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR421 </td>
   <td style="text-align:right;"> 1479 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR317 </td>
   <td style="text-align:right;"> 1481 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR399 </td>
   <td style="text-align:right;"> 1579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR561A </td>
   <td style="text-align:right;"> 1624 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR311 </td>
   <td style="text-align:right;"> 1634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR420 </td>
   <td style="text-align:right;"> 1700 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR139C </td>
   <td style="text-align:right;"> 1730 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR333 </td>
   <td style="text-align:right;"> 1785 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR558A </td>
   <td style="text-align:right;"> 1793 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR554A </td>
   <td style="text-align:right;"> 1972 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR419 </td>
   <td style="text-align:right;"> 2037 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR417 </td>
   <td style="text-align:right;"> 2075 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR160C </td>
   <td style="text-align:right;"> 2191 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR107C </td>
   <td style="text-align:right;"> 2201 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR556A </td>
   <td style="text-align:right;"> 2311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR334 </td>
   <td style="text-align:right;"> 2388 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR557A </td>
   <td style="text-align:right;"> 2532 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR392 </td>
   <td style="text-align:right;"> 2952 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR395 </td>
   <td style="text-align:right;"> 3001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR149C </td>
   <td style="text-align:right;"> 3054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR394 </td>
   <td style="text-align:right;"> 3055 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR332 </td>
   <td style="text-align:right;"> 3152 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR393 </td>
   <td style="text-align:right;"> 3167 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR555A </td>
   <td style="text-align:right;"> 3238 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR416 </td>
   <td style="text-align:right;"> 3539 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR562A </td>
   <td style="text-align:right;"> 3614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR415 </td>
   <td style="text-align:right;"> 3728 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR328 </td>
   <td style="text-align:right;"> 3792 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR400 </td>
   <td style="text-align:right;"> 3815 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR316 </td>
   <td style="text-align:right;"> 4073 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR565A </td>
   <td style="text-align:right;"> 4122 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR331 </td>
   <td style="text-align:right;"> 4337 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR402 </td>
   <td style="text-align:right;"> 4429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR560A </td>
   <td style="text-align:right;"> 4453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR401 </td>
   <td style="text-align:right;"> 4461 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR563A </td>
   <td style="text-align:right;"> 4999 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR315 </td>
   <td style="text-align:right;"> 5085 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 239 </td>
   <td style="text-align:right;"> 5356 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:right;"> 5464 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr248a </td>
   <td style="text-align:right;"> 5585 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR318 </td>
   <td style="text-align:right;"> 5696 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR102C </td>
   <td style="text-align:right;"> 6357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PR322 </td>
   <td style="text-align:right;"> 6779 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr250a </td>
   <td style="text-align:right;"> 7781 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr262a </td>
   <td style="text-align:right;"> 8718 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr271a </td>
   <td style="text-align:right;"> 9101 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 101 </td>
   <td style="text-align:right;"> 9859 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr267a </td>
   <td style="text-align:right;"> 10277 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 220 </td>
   <td style="text-align:right;"> 10514 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr247a </td>
   <td style="text-align:right;"> 10711 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 44 </td>
   <td style="text-align:right;"> 10755 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 104 </td>
   <td style="text-align:right;"> 10896 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs1-1 </td>
   <td style="text-align:right;"> 11179 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpw-3 </td>
   <td style="text-align:right;"> 11765 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr295a </td>
   <td style="text-align:right;"> 12285 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lsj3-2 </td>
   <td style="text-align:right;"> 12316 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpe-3 </td>
   <td style="text-align:right;"> 12325 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr249a </td>
   <td style="text-align:right;"> 12448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmp2-2 </td>
   <td style="text-align:right;"> 12550 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 142 </td>
   <td style="text-align:right;"> 12638 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 198 </td>
   <td style="text-align:right;"> 12984 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lp1-2 </td>
   <td style="text-align:right;"> 13656 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 42 </td>
   <td style="text-align:right;"> 13896 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmp2-1 </td>
   <td style="text-align:right;"> 13943 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpe-2 </td>
   <td style="text-align:right;"> 14167 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmp2-3 </td>
   <td style="text-align:right;"> 14339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lp1-3 </td>
   <td style="text-align:right;"> 14775 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:right;"> 14776 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 218 </td>
   <td style="text-align:right;"> 15212 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 181 </td>
   <td style="text-align:right;"> 15472 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpe-1 </td>
   <td style="text-align:right;"> 15558 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpw-1 </td>
   <td style="text-align:right;"> 15630 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 183 </td>
   <td style="text-align:right;"> 16414 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 63 </td>
   <td style="text-align:right;"> 16771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lsj3-3 </td>
   <td style="text-align:right;"> 16802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs1-2 </td>
   <td style="text-align:right;"> 19707 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 95 </td>
   <td style="text-align:right;"> 21854 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lp1-1 </td>
   <td style="text-align:right;"> 23460 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 106 </td>
   <td style="text-align:right;"> 25062 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 90 </td>
   <td style="text-align:right;"> 25753 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cs1-3 </td>
   <td style="text-align:right;"> 29344 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cmpw-2 </td>
   <td style="text-align:right;"> 51074 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pr293a </td>
   <td style="text-align:right;"> 73188 </td>
  </tr>
</tbody>
</table></div>

# Remove samples with \<1,000 reads


```r
#Remove samples with less than 1,000 reads
PR_Normal = subset_samples(PR_Normal, 
                            Sample_ID != "PR105C")
PR_Normal = subset_samples(PR_Normal, 
                            Sample_ID != "PR109C")
PR_Normal = subset_samples(PR_Normal, 
                           Sample_ID != "pr242a")
PR_Normal = subset_samples(PR_Normal, 
                            Sample_ID != "PR325")
PR_Normal = subset_samples(PR_Normal, 
                            Sample_ID != "PR164C")
PR_Normal = subset_samples(PR_Normal, 
                            Sample_ID != "PR151C")
```

# Examining abundance of individual ASVs


```r
#For each ASV, find the number of samples in which each ASV is 0, then divide by the total number of samples 
test_filter_method <- as.data.frame(rowSums(otu_table(PR_Normal) == 0)/ncol(otu_table(PR_Normal)))    
#Change the name of the column in the test_filter_method4 dataframe; this column contains the proportion of samples in which each ASV is 0  
colnames(test_filter_method) <- c('samplew0')   
#Create a bar plot 
ggplot(data = test_filter_method) + geom_bar(mapping = aes(x= samplew0)) +labs(x = "Proportion of samples in which ASV = 0", y = "# of ASV's") +theme(text = element_text(size = 18), axis.title = element_text(size = 15),panel.spacing = unit(1, "lines"),panel.border = element_rect(colour = "black", fill = NA, size = 0.5), panel.background = element_blank()) 
```

```
## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
## i Please use the `linewidth` argument instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

<img src="PR-S003-S015-data-analysis_files/figure-html/filter method-1.png" width="672" />

## Removing singletons and doubletons


```r
##How many singletons and doubletons are present 
#Create a dataframe that includes only the total number of times that each ASV occurs
readsumsdf_ASVs <- readsumsdf %>% 
    filter(type == "ASVs")    
#How many singletons are present?
length(which(readsumsdf_ASVs$nreads <= 0))
```

```
## [1] 1908
```

```r
length(which(readsumsdf_ASVs$nreads == 1))
```

```
## [1] 7
```

```r
#How many doubletons are present?
length(which(readsumsdf_ASVs$nreads == 2))
```

```
## [1] 158
```

```r
#Filter out low abundance ASVs

PR_Normal = prune_taxa(taxa_sums(PR_Normal) > 2, PR_Normal)
```

Examine the abundance of individual ASVs again:


```r
#For each ASV, find the number of samples in which each ASV is 0, then divide by the total number of samples 
test_filter_method_filtered <- as.data.frame(rowSums(otu_table(PR_Normal) == 0)/ncol(otu_table(PR_Normal)))   

#Change the name of the column in the test_filter_method4 dataframe; this column contains the proportion of samples in which each ASV is 0 
colnames(test_filter_method_filtered) <- c('samplew0')  

#Create a bar plot 
ggplot(data = test_filter_method_filtered) + geom_bar(mapping = aes(x= samplew0)) + 
  labs(x = "Proportion of samples in which ASV = 0", y = "# of ASV's") +   
  theme(text = element_text(size = 18), axis.title = element_text(size = 15), 
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank()) 
```

<img src="PR-S003-S015-data-analysis_files/figure-html/filt abundance check-1.png" width="672" />

# Label variables


```r
sample_data(PR_Normal)$Sub_Location <- factor(sample_data(PR_Normal)$Sub_Location,              levels = c("BSJ1","BSJ1-ENR","BSJ2","BSJ2-ENR","BSJ3","BSJ3-ENR","BSJ4-ENR","BSJ5-ENR","BSJ6-ENR","BSJ7-ENR","CMP1","CMP1-ENR","CMP2","CMP2-ENR","CMPE","CMPW","CS1","CS1-ENR","CS2","CS2-ENR","CSA1-ENR","LLC1","LLC1-ENR","LP1","LP1-ENR","LP1-ENR","LP1L-ENR","LP2-ENR","LP3-ENR","LP4-ENR","LP5-ENR","LSJ1","LSJ1-ENR","LSJ2","LSJ2-ENR","LSJ3","LT1","LT1-ENR","LT2","LT2-ENR","LT3-ENR","QB1-ENR","QSA1-ENR","RPN1-ENR"),                         
labels = c("BSJ1","BSJ1-ENR","BSJ2","BSJ2-ENR","BSJ3","BSJ3-ENR","BSJ4-ENR","BSJ5-ENR","BSJ6-ENR","BSJ7-ENR","CMP1","CMP1-ENR","CMP2","CMP2-ENR","CMPE","CMPW","CS1","CS1-ENR","CS2","CS2-ENR","CSA1-ENR","LLC1","LLC1-ENR","LP1","LP1-ENR","LP1-ENR","LP1L-ENR","LP2-ENR","LP3-ENR","LP4-ENR","LP5-ENR","LSJ1","LSJ1-ENR","LSJ2","LSJ2-ENR","LSJ3","LT1","LT1-ENR","LT2","LT2-ENR","LT3-ENR","QB1-ENR","QSA1-ENR","RPN1-ENR"))

sample_data(PR_Normal)$Location <- factor(sample_data(PR_Normal)$Location,              levels = c(
  "BSJ","BSJ-ENR","CMP","CMP-ENR","CMPE","CMPW", "CS","CS-ENR", "CSA","CSA-ENR","LLC","LLC-ENR","LP","LP-ENR","LSJ","LSJ-ENR","LT","LT-ENR","QB","QB-ENR","QSA","QSA-ENR","RPN","RPN-ENR"),
labels = c("BSJ","BSJ-ENR","CMP","CMP-ENR","CMPE","CMPW", "CS","CS-ENR", "CSA","CSA-ENR","LLC","LLC-ENR","LP","LP-ENR","LSJ","LSJ-ENR","LT","LT-ENR","QB","QB-ENR","QSA","QSA-ENR","RPN","RPN-ENR"))
sample_data(PR_Normal)$Season <- factor(sample_data(PR_Normal)$Season,                                             levels = c("W", "D"), 
labels = c("W", "D"))

sample_data(PR_Normal)$Project <- factor(sample_data(PR_Normal)$Project,
                               levels = c("S003","S015"), labels = c("S003","S015"))

sample_data(PR_Normal)$Month <- factor(sample_data(PR_Normal)$Month,
                               levels = c("April","August","December","January","July","June","March","May","November","October","September"), labels = c("April","August","December","January","July","June","March","May","November","October","September"))
  
sample_data(PR_Normal)$Year <- factor(sample_data(PR_Normal)$Year,
                                     levels = c("2021","2022","2023"), labels = c("2021","2022","2023"))
  
sample_data(PR_Normal)$Sample_Type <- factor(sample_data(PR_Normal)$Sample_Type,
                                            levels = c("N","E"), labels = c("N","E"))
sample_data(PR_Normal)$Sample_ID <- factor(sample_data(PR_Normal)$Sample_ID,
                                            levels = c("PR419","PR417","PR422","PR421","PR420","PR416","PR415","90","44","42","63","5","18","pr293a", "pr271a","pr262a","pr285a","pr295a","pr267a","pr290a","PR425","PR325","PR323","PR321","PR322","PR332","PR333","PR331","PR315","PR318","PR328","PR317","PR316","PR334",	"PR311","PR563A","PR562A","PR560A","PR561A","PR565A","cmp2-3","cmpe-3","cmpw-3","cs1-3","PR556A","PR559A","PR564A","lp1-3","PR558A","PR557A","PR566A","lsj3-3","PR555A","PR554A","cmp2-2","cmpe-2","cmpw-2","cs1-2","lp1-2","lsj3-2","PR395","PR394","PR392",	"PR393","PR400","PR401","PR399","PR402","PR464A","PR465A","PR466A","PR467A","PR468A","PR469A","PR470A","PR473A","cmp2-1","PR472A","cmpe-1","cmpw-1","cs1-1","PR476A","PR477A","PR474A","PR478A","lp1-1","PR488A","PR489A","PR490A","PR491A","PR492A","PR493A","PR479A","PR480A","lsj3-1","PR481A","PR482A","PR483A","PR484A","PR485A","PR471A","Mock","Undeter","B53","B84","B59","B43","B36","B44","B107","B90","B126","B141","B128","B152","B134","B132","B140","B89","B100","B102","B60","B76","B56","blank","blank2","Mock-com","B130","B136","B138","B33","B25","B41","B81","B68","B57","B87","B91","B117","B143","B139","B154","B69","B52","B75","B99","B118","B124","B26","B34","B42","pr249a","pr248a","pr250a","pr242a","pr247a","183","218","181","198","239","220","PR160C","PR162C","PR164C","PR102C","101","PR109C","108","PR107C","106","PR143C","142","PR145C","PR151C","PR96C","95","PR147C","PR149C","PR105C","104","PR139C","PR141C"
),
labels = c("PR419","PR417","PR422","PR421","PR420","PR416","PR415","90","44","42","63","5","18","pr293a", "pr271a","pr262a","pr285a","pr295a","pr267a","pr290a","PR425","PR325","PR323","PR321","PR322","PR332","PR333","PR331","PR315","PR318","PR328","PR317","PR316","PR334",	"PR311","PR563A","PR562A","PR560A","PR561A","PR565A","cmp2-3","cmpe-3","cmpw-3","cs1-3","PR556A","PR559A","PR564A","lp1-3","PR558A","PR557A","PR566A","lsj3-3","PR555A","PR554A","cmp2-2","cmpe-2","cmpw-2","cs1-2","lp1-2","lsj3-2","PR395","PR394","PR392",	"PR393","PR400","PR401","PR399","PR402","PR464A","PR465A","PR466A","PR467A","PR468A","PR469A","PR470A","PR473A","cmp2-1","PR472A","cmpe-1","cmpw-1","cs1-1","PR476A","PR477A","PR474A","PR478A","lp1-1","PR488A","PR489A","PR490A","PR491A","PR492A","PR493A","PR479A","PR480A","lsj3-1","PR481A","PR482A","PR483A","PR484A","PR485A","PR471A","Mock","Undeter","B53","B84","B59","B43","B36","B44","B107","B90","B126","B141","B128","B152","B134","B132","B140","B89","B100","B102","B60","B76","B56","blank","blank2","Mock-com","B130","B136","B138","B33","B25","B41","B81","B68","B57","B87","B91","B117","B143","B139","B154","B69","B52","B75","B99","B118","B124","B26","B34","B42","pr249a","pr248a","pr250a","pr242a","pr247a","183","218","181","198","239","220","PR160C","PR162C","PR164C","PR102C","101","PR109C","108","PR107C","106","PR143C","142","PR145C","PR151C","PR96C","95","PR147C","PR149C","PR105C","104","PR139C","PR141C"))
```

# Tax Fix


```r
#View the tax_table
knitr::kable(head(tax_table(PR_Normal))) %>%   
  kableExtra::kable_styling("striped") %>%    
  kableExtra::scroll_box(width = "100%")  
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Kingdom </th>
   <th style="text-align:left;"> Phylum </th>
   <th style="text-align:left;"> Class </th>
   <th style="text-align:left;"> Order </th>
   <th style="text-align:left;"> Family </th>
   <th style="text-align:left;"> Genus </th>
   <th style="text-align:left;"> Species </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1779c60c03d948a73687f440be973f68 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Flammeovirgaceae </td>
   <td style="text-align:left;"> Algivirga </td>
   <td style="text-align:left;"> Algivirga_pacifica </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25feebceed837c1dfa810cd8c0f26b41 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Cyclobacteriaceae </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0872b535e6e15b9cd031dae6a67d6257 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Amoebophilaceae </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured_bacterium </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6590e42d5aeb18c338a836bc26028592 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Chitinophagales </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured_Bacteroidetes </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0f22f07d74171ba4bbe6b3977217cbdc </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> uncultured_organism </td>
  </tr>
  <tr>
   <td style="text-align:left;"> d817bb1b166e13e13db6db8d56b43983 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table></div>


```r
#Now fix the labels
PR_Normal <- tax_fix(PR_Normal)

#View the relabeled tax_table
knitr::kable(head(tax_table(PR_Normal))) %>%   
  kableExtra::kable_styling("striped") %>%    
  kableExtra::scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Kingdom </th>
   <th style="text-align:left;"> Phylum </th>
   <th style="text-align:left;"> Class </th>
   <th style="text-align:left;"> Order </th>
   <th style="text-align:left;"> Family </th>
   <th style="text-align:left;"> Genus </th>
   <th style="text-align:left;"> Species </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1779c60c03d948a73687f440be973f68 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Flammeovirgaceae </td>
   <td style="text-align:left;"> Algivirga </td>
   <td style="text-align:left;"> Algivirga_pacifica </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25feebceed837c1dfa810cd8c0f26b41 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Cyclobacteriaceae </td>
   <td style="text-align:left;"> Cyclobacteriaceae Family </td>
   <td style="text-align:left;"> Cyclobacteriaceae Family </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0872b535e6e15b9cd031dae6a67d6257 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Cytophagales </td>
   <td style="text-align:left;"> Amoebophilaceae </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured_bacterium </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6590e42d5aeb18c338a836bc26028592 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Chitinophagales </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured </td>
   <td style="text-align:left;"> uncultured_Bacteroidetes </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 0f22f07d74171ba4bbe6b3977217cbdc </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> uncultured_organism </td>
  </tr>
  <tr>
   <td style="text-align:left;"> d817bb1b166e13e13db6db8d56b43983 </td>
   <td style="text-align:left;"> d__Bacteria </td>
   <td style="text-align:left;"> Bacteroidota </td>
   <td style="text-align:left;"> Bacteroidia </td>
   <td style="text-align:left;"> Bacteroidales </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 </td>
   <td style="text-align:left;"> Bacteroidetes_BD2-2 Genus </td>
  </tr>
</tbody>
</table></div>

# Bacterial community composition bar plots


```r
library(seecolor) 
bcc_hex <- print_color(c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF","#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF", "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF","#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF","#DAE319FF", "#E4E419FF","#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3","#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299","#300018"), type = "r")
```

```
## 
##  ------ c black #440154FF #450659FF #460B5EFF #472D7AFF #3B518BFF #3A548CFF #38598CFF #365C8DFF #34608DFF #33638DFF #31678EFF #306A8EFF #2E6E8EFF #2D718EFF #2B748EFF #2A778EFF #297B8EFF #1F958BFF #1F988BFF #1E9C89FF #1F9F88FF #1FA287FF #21A585FF #23A983FF #25AC82FF #29AF7FFF #2DB27DFF #32B67AFF #37B878FF #3CBC74FF #57C766FF #5EC962FF #A2DA37FF #DAE319FF #E4E419FF #ECE51BFF #F5E61FFF darkorchid1 darkorchid2 darkorchid3 #A21C9AFF #A62098FF #AB2394FF #AE2892FF #B22B8FFF #B6308BFF #BA3388FF #BE3885FF #C13B82FF #C53F7EFF #C8437BFF #CC4678FF #CE4B75FF #D14E72FF #D5536FFF #D7566CFF #DA5B69FF #DD5E66FF #E06363FF #FF6699 #F68D45FF #FCA537FF #F6E726FF #F4ED27FF #00489C #CCCCCC #999999 #A1C299 #300018 ------
## black               
## #440154FF           
## #450659FF           
## #460B5EFF           
## #472D7AFF           
## #3B518BFF           
## #3A548CFF           
## #38598CFF           
## #365C8DFF           
## #34608DFF           
## #33638DFF           
## #31678EFF           
## #306A8EFF           
## #2E6E8EFF           
## #2D718EFF           
## #2B748EFF           
## #2A778EFF           
## #297B8EFF           
## #1F958BFF           
## #1F988BFF           
## #1E9C89FF           
## #1F9F88FF           
## #1FA287FF           
## #21A585FF           
## #23A983FF           
## #25AC82FF           
## #29AF7FFF           
## #2DB27DFF           
## #32B67AFF           
## #37B878FF           
## #3CBC74FF           
## #57C766FF           
## #5EC962FF           
## #A2DA37FF           
## #DAE319FF           
## #E4E419FF           
## #ECE51BFF           
## #F5E61FFF           
## darkorchid1         
## darkorchid2         
## darkorchid3         
## #A21C9AFF           
## #A62098FF           
## #AB2394FF           
## #AE2892FF           
## #B22B8FFF           
## #B6308BFF           
## #BA3388FF           
## #BE3885FF           
## #C13B82FF           
## #C53F7EFF           
## #C8437BFF           
## #CC4678FF           
## #CE4B75FF           
## #D14E72FF           
## #D5536FFF           
## #D7566CFF           
## #DA5B69FF           
## #DD5E66FF           
## #E06363FF           
## #FF6699             
## #F68D45FF           
## #FCA537FF           
## #F6E726FF           
## #F4ED27FF           
## #00489C             
## #CCCCCC             
## #999999             
## #A1C299             
## #300018
```

## Creating Family Barplot


```r
##Subsetting non-enriched samples to view composition in bar plot

PR_N <- PR_Normal %>% subset_samples(Sample_Type %in% c("N"))
```


```r
PRSeqR_family <- PR_N %>%  
  tax_glom(taxrank = "Family") %>%   
  transform_sample_counts(function(x) {x/sum(x)}) %>%   
  psmelt() %>%   
  group_by(Location,             
           Kingdom, Phylum, Class, Order, Family) %>%  
  filter(Abundance > 0.01) %>%    
  arrange(Class)
PRSeqR_family$Class_family <- paste(PRSeqR_family$Class, PRSeqR_family$Family, sep="_") 
PRSeqR_family_bar <- ggplot(PRSeqR_family, aes(x = Location, y = Abundance, fill = Class_family)) +   
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") + scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (Family > 1%) \n") + xlab("Location") +   
  theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(),                         
panel.grid.major = element_blank(), axis.text.x  =element_text(angle = 90, vjust = 0.5, hjust=1, size=8, colour="black"),                          
axis.text.y = element_text(size=8, colour="black"), plot.title = element_text(hjust = 0.5),     
axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),                        
legend.text = element_text(size = 7))
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## i Please use `linewidth` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```
## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
## i Please use the `linewidth` argument instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```


```r
print(PRSeqR_family_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/family bar plot1-1.png" width="672" />


```r
p <- ggplotly(PRSeqR_family_bar)
```


```r
htmlwidgets::saveWidget(p, "PRSeqR_family_bar.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "PRSeqR_family_bar.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0")
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/PRSeqR_family_bar.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

## Creating Species Barplot


```r
PRSeqR_Species <- PR_N %>%  
  tax_glom(taxrank = "Species") %>%   
  transform_sample_counts(function(x) {x/sum(x)}) %>%   
  psmelt() %>%   
  group_by(Location,             
           Kingdom, Phylum, Class, Order, Family,Genus,Species) %>%  
  filter(Abundance > 0.05) %>%    
  arrange(Species)
PRSeqR_Species_bar <- ggplot(PRSeqR_Species, aes(x = Location, y = Abundance, fill = Species)) +   
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") + scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +  
  ylab("Relative Abundance (Species > 5%) \n") + xlab("Sample ID") +   
  theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(),                         
panel.grid.major = element_blank(), axis.text.x  =element_text(angle = 90, vjust = 0.5, hjust=1, size=8, colour="black"),                          
axis.text.y = element_text(size=8, colour="black"), plot.title = element_text(hjust = 0.5),     
axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),                        
legend.text = element_text(size = 7))
```


```r
p2 <- ggplotly(PRSeqR_Species_bar)
```


```r
print(PRSeqR_Species_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/species bar plot1-1.png" width="672" />


```r
htmlwidgets::saveWidget(p, "PRSeqR_Species_bar.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "PRSeqR_Species_bar.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0")
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/PRSeqR_Species_bar.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

## Creating Season vs Location Barplot with non-enriched samples


```r
PRSeqR_family_Season <- PR_N %>%    
  tax_glom(taxrank = "Family") %>%     
  transform_sample_counts(function(x) {x/sum(x)}) %>%   
  psmelt() %>%    
  group_by(Season,             
           Kingdom, Phylum, Class, Order, Family)%>%   
  dplyr::summarize(Mean =                      
                     mean(Abundance, na.rm=TRUE)) %>%    
  filter(Mean > 0.01) %>%                                
  arrange(Class)  
```

```
## `summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument.
```

```r
## `summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument.
```


```r
## Creaitng the plot
PRSeqR_family_Season$Class_family <- paste(PRSeqR_family_Season$Class, PRSeqR_family_Season$Family, sep="_")  
PRSeqR_Season_family_bar <- ggplot(PRSeqR_family_Season, aes(x = Season, y = Mean, fill = Class_family)) +   
  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black", "#440154FF", "#450659FF","#460B5EFF","#472D7AFF", "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF","#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF","#2E6E8EFF", "#2D718EFF", "#2B748EFF","#2A778EFF", "#297B8EFF","#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF", "#1FA287FF", "#21A585FF", "#23A983FF","#25AC82FF", "#29AF7FFF", "#2DB27DFF","#32B67AFF", "#37B878FF", "#3CBC74FF","#57C766FF", "#5EC962FF","#A2DA37FF", "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF","darkorchid1", "darkorchid2", "darkorchid3", "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF","#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#FF6699","#F68D45FF", "#FCA537FF","#F6E726FF", "#F4ED27FF", "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999","#999999", "#A1C299", "#300018","#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample ID") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(),                          panel.grid.major = element_blank(), axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8, colour="black"),                          axis.text.y = element_text(size=8, colour="black"), plot.title = element_text(hjust = 0.5),                          axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),                         legend.text = element_text(size = 7))    
print(PRSeqR_Season_family_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/creating the plot-1.png" width="672" />


```r
p3 <- ggplotly(PRSeqR_Season_family_bar)

p
```

```{=html}
<div class="plotly html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-443b8eb1c663449b913e" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-443b8eb1c663449b913e">{"x":{"data":[{"orientation":"v","width":[0.899999999999999,0.899999999999999],"base":[0.99824563711163,0.999371270662033],"x":[8,8],"y":[0.00112563355040352,0.000628729337966871],"text":["Location: LSJ<br />Abundance: 0.0011256336<br />Class_family: ABY1_Candidatus_Uhrbacteria","Location: LSJ<br />Abundance: 0.0006287293<br />Class_family: ABY1_Candidatus_Uhrbacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"ABY1_Candidatus_Uhrbacteria","legendgroup":"ABY1_Candidatus_Uhrbacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.998240857565648,"x":[7],"y":[0.00175914243435238],"text":"Location: LP<br />Abundance: 0.0017591424<br />Class_family: Acidimicrobiia_Actinomarinaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(68,1,84,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Acidimicrobiia_Actinomarinaceae","legendgroup":"Acidimicrobiia_Actinomarinaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9],"base":[0.981859864856417,0.98670289563782,0.989787377929104,0.978376639759518,0.983218764634239,0.986525042893824,0.98977178398411,0.987974873399688,0.992184078654571,0.992249498432893,0.993761189266698,0.991338725093522,0.998319942368078,0.990271718050602,0.995206730838505,0.992292926745932,0.99429163847261,0.996403157385922,0.994997404744755,0.996202047946317,0.997720380568005,0.997520444146993,0.995613293020271,0.997893777165094,0.998845413626852,0.998974052531538,0.997014963076727],"x":[3,8,8,7,7,5,3,7,7,8,8,6,2,5,8,5,5,8,3,5,3,8,7,5,4,5,7],"y":[0.00791191912769262,0.00308448229128333,0.00246212050378891,0.00484212487472069,0.0047561087654493,0.0037466751567784,0.00522562076064514,0.00420920525488311,0.00342921436569943,0.00151169083380576,0.00144554157180654,0.00866127490647806,0.00168005763192169,0.0020212086953294,0.00119642654741714,0.00199871172667854,0.00191040947370713,0.00111728676107103,0.00272297582325065,0.00169172921877714,0.00227961943199462,0.000725192964636467,0.00140167005645686,0.00108027536644373,0.00115458637314791,0.00102594746846185,0.00122589448892019],"text":["Location: CMPE<br />Abundance: 0.0079119191<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0030844823<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0024621205<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0048421249<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0047561088<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0037466752<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CMPE<br />Abundance: 0.0052256208<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0042092053<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0034292144<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0015116908<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0014455416<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LLC<br />Abundance: 0.0086612749<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CMP<br />Abundance: 0.0016800576<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0020212087<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0011964265<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0019987117<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0019104095<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0011172868<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CMPE<br />Abundance: 0.0027229758<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0016917292<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CMPE<br />Abundance: 0.0022796194<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LSJ<br />Abundance: 0.0007251930<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0014016701<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0010802754<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CMPW<br />Abundance: 0.0011545864<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: CS<br />Abundance: 0.0010259475<br />Class_family: Acidimicrobiia_Ilumatobacteraceae","Location: LP<br />Abundance: 0.0012258945<br />Class_family: Acidimicrobiia_Ilumatobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(69,6,89,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Acidimicrobiia_Ilumatobacteraceae","legendgroup":"Acidimicrobiia_Ilumatobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999],"base":[0.944985853118349,0.967027860223369,0.97213251523338,0.958426085828719,0.961264175721351,0.974276398415957,0.95877474980568,0.97901221739384,0.968017061748783,0.976404170849428,0.984147240740396,0.967529772846724,0.988727519174256,0.983049764480657,0.96594330051754,0.969995955245991,0.985623852780994,0.988152584845649,0.975458981110386,0.971642506592005,0.97476328705902,0.990651121222103,0.980251180126586,0.984558292696944,0.977829977861487,0.992956527401484,0.992845357039413,0.978960686066314,0.988610507294832,0.980437094808233,0.974032102592817,0.98031801410371,0.994613211688185,0.981509758300314,0.995456621544787,0.992075878975267,0.994772147808899,0.997450188530823,0.996147436234393,0.982480761377955,0.976483891570658,0.982622222240834,0.98594023152979,0.98405466478126,0.997373833039571,0.98008475136429,0.997419774433295,0.983535777014467,0.984206149071109,0.985435586355665,0.984861783048617,0.9854883636104,0.986114790334579],"x":[7,8,8,3,5,2,7,2,3,8,4,5,4,2,7,7,2,2,3,5,5,2,9,9,5,4,2,8,9,5,7,8,2,8,4,9,9,9,2,5,7,8,6,5,4,3,2,8,8,5,8,8,8],"y":[0.0137888966873307,0.00510465501001023,0.00427165561604803,0.00959097592006464,0.00626559712537345,0.00473581897788322,0.00716855071186029,0.00403754708681658,0.00744191936160232,0.0025565152168866,0.00458027843385944,0.00411273374528065,0.00422900822722783,0.00257408830033701,0.00405265472845107,0.00403614734682634,0.00252873206465476,0.00249853637645447,0.00462577025390454,0.00312078046701503,0.00306669080246735,0.0021942358173096,0.00430711257035798,0.00405221459788807,0.00260711694674554,0.00250009414330354,0.00176785464877205,0.00135732803739597,0.00346537168043493,0.00204366656972188,0.00245178897784037,0.00119174419660395,0.00153422454620789,0.00111246394051945,0.00191721149478352,0.00269626883363216,0.00267804072192412,0.00254981146917688,0.00127233819890249,0.00157390340330477,0.0018927481888602,0.000913554773633085,0.00539849356373223,0.00138092157440506,0.00147158058728147,0.00177511349212689,0.000900167934782981,0.000670372056642554,0.000655633977507786,0.00108945653815917,0.000626580561783219,0.000626426724179008,0.000588105303241271],"text":["Location: LP<br />Abundance: 0.0137888967<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0051046550<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0042716556<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPE<br />Abundance: 0.0095909759<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0062655971<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0047358190<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LP<br />Abundance: 0.0071685507<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0040375471<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPE<br />Abundance: 0.0074419194<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0025565152<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPW<br />Abundance: 0.0045802784<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0041127337<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPW<br />Abundance: 0.0042290082<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0025740883<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LP<br />Abundance: 0.0040526547<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LP<br />Abundance: 0.0040361473<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0025287321<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0024985364<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPE<br />Abundance: 0.0046257703<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0031207805<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0030666908<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0021942358<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0043071126<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0040522146<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0026071169<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPW<br />Abundance: 0.0025000941<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0017678546<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0013573280<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0034653717<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0020436666<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LP<br />Abundance: 0.0024517890<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0011917442<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0015342245<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0011124639<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPW<br />Abundance: 0.0019172115<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0026962688<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0026780407<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LT<br />Abundance: 0.0025498115<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0012723382<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0015739034<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LP<br />Abundance: 0.0018927482<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0009135548<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LLC<br />Abundance: 0.0053984936<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0013809216<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPW<br />Abundance: 0.0014715806<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMPE<br />Abundance: 0.0017751135<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CMP<br />Abundance: 0.0009001679<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0006703721<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0006556340<br />Class_family: Actinobacteria_Microbacteriaceae","Location: CS<br />Abundance: 0.0010894565<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0006265806<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0006264267<br />Class_family: Actinobacteria_Microbacteriaceae","Location: LSJ<br />Abundance: 0.0005881053<br />Class_family: Actinobacteria_Microbacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(70,11,94,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Microbacteriaceae","legendgroup":"Actinobacteria_Microbacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.956874230308496,"x":[3],"y":[0.00155185552022274],"text":"Location: CMPE<br />Abundance: 0.0015518555<br />Class_family: Actinobacteria_Nitriliruptoraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(71,45,122,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Nitriliruptoraceae","legendgroup":"Actinobacteria_Nitriliruptoraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9],"base":[0.951872740305619,0.920064226388251,0.928563273014011,0.947005404521682,0.959474577232501,0.94454849085076,0.952325509384671,0.96220201244575,0.935472276843272,0.955882775618675,0.950615703971498,0.963989671263464,0.96515382878885,0.938946026460474,0.957978843895313,0.979378659630943,0.941242267045628,0.959823808566093,0.943324458071523,0.966271628349058,0.953541333230934,0.973412368776683,0.955248600057811],"x":[8,7,7,5,8,3,5,8,7,5,3,8,8,7,5,6,7,5,7,8,3,2,3],"y":[0.00760183692688288,0.00849904662576006,0.00690900382926085,0.00532010486298917,0.00272743521324881,0.00606721312073788,0.00355726623400376,0.00178765881771337,0.00347374961720215,0.00209606827663789,0.00292562925943618,0.00116415752538634,0.00111779956020841,0.00229624058515376,0.0018449646707801,0.0065615718988471,0.00208219102589591,0.0014403671552583,0.00166139504682561,0.000756231874310953,0.0017072668268775,0.000864029639274122,0.00162563025068463],"text":["Location: LSJ<br />Abundance: 0.0076018369<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0084990466<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0069090038<br />Class_family: Actinobacteria_PeM15","Location: CS<br />Abundance: 0.0053201049<br />Class_family: Actinobacteria_PeM15","Location: LSJ<br />Abundance: 0.0027274352<br />Class_family: Actinobacteria_PeM15","Location: CMPE<br />Abundance: 0.0060672131<br />Class_family: Actinobacteria_PeM15","Location: CS<br />Abundance: 0.0035572662<br />Class_family: Actinobacteria_PeM15","Location: LSJ<br />Abundance: 0.0017876588<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0034737496<br />Class_family: Actinobacteria_PeM15","Location: CS<br />Abundance: 0.0020960683<br />Class_family: Actinobacteria_PeM15","Location: CMPE<br />Abundance: 0.0029256293<br />Class_family: Actinobacteria_PeM15","Location: LSJ<br />Abundance: 0.0011641575<br />Class_family: Actinobacteria_PeM15","Location: LSJ<br />Abundance: 0.0011177996<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0022962406<br />Class_family: Actinobacteria_PeM15","Location: CS<br />Abundance: 0.0018449647<br />Class_family: Actinobacteria_PeM15","Location: LLC<br />Abundance: 0.0065615719<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0020821910<br />Class_family: Actinobacteria_PeM15","Location: CS<br />Abundance: 0.0014403672<br />Class_family: Actinobacteria_PeM15","Location: LP<br />Abundance: 0.0016613950<br />Class_family: Actinobacteria_PeM15","Location: LSJ<br />Abundance: 0.0007562319<br />Class_family: Actinobacteria_PeM15","Location: CMPE<br />Abundance: 0.0017072668<br />Class_family: Actinobacteria_PeM15","Location: CMP<br />Abundance: 0.0008640296<br />Class_family: Actinobacteria_PeM15","Location: CMPE<br />Abundance: 0.0016256303<br />Class_family: Actinobacteria_PeM15"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(59,81,139,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_PeM15","legendgroup":"Actinobacteria_PeM15","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.899999999999999],"base":[0.948275830126438,0.94225886795207,0.950315347350724,0.951206804774509],"x":[8,3,8,8],"y":[0.00203951722428586,0.00228962289868928,0.000891457423785491,0.000665935531109496],"text":["Location: LSJ<br />Abundance: 0.0020395172<br />Class_family: Actinobacteria_Sporichthyaceae","Location: CMPE<br />Abundance: 0.0022896229<br />Class_family: Actinobacteria_Sporichthyaceae","Location: LSJ<br />Abundance: 0.0008914574<br />Class_family: Actinobacteria_Sporichthyaceae","Location: LSJ<br />Abundance: 0.0006659355<br />Class_family: Actinobacteria_Sporichthyaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(58,84,140,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Actinobacteria_Sporichthyaceae","legendgroup":"Actinobacteria_Sporichthyaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9],"base":[0.891557635239227,0.899375964871347,0.973313348484188,0.978061943134645,0.98280443663277,0.906724281369274,0.963349827282501,0.9597467950688,0.911250661020672,0.987289098590552,0.965062905755604,0.969402213409731,0.966116071842374,0.96812133494279,0.990502104131744,0.97368833101448,0.993071301629599,0.970072981382768,0.914967985967992,0.971825939992608,0.94355883462761,0.977258067071328,0.917722383819543,0.995583235008276,0.944772610907075,0.94586413434989,0.946789461736586,0.997584970062226,0.98291801050292,0.998855728304806,0.947663640373263,0.94600382962867],"x":[7,7,1,1,1,7,2,9,7,1,9,9,2,2,1,9,1,2,7,2,8,9,7,1,8,8,8,1,4,1,8,5],"y":[0.00781832963212026,0.00734831649792744,0.00474859465045718,0.00474249349812517,0.00448466195778185,0.0045263796513979,0.00276624455987351,0.0053161106868046,0.00371732494732013,0.00321300554119219,0.00433930765412671,0.00428611760474862,0.00200526310041615,0.00195164643997725,0.00256919749785456,0.00356973605684774,0.00251193337867683,0.00175295860984059,0.00275439785155041,0.00158642878407489,0.0012137762794644,0.00299311305525829,0.00234184256870784,0.00200173505395063,0.0010915234428156,0.000925327386695374,0.000874178636677048,0.00127075824257972,0.0012292302374759,0.00114427169519393,0.000612189753174897,0.001001574893012],"text":["Location: LP<br />Abundance: 0.0078183296<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LP<br />Abundance: 0.0073483165<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0047485947<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0047424935<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0044846620<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LP<br />Abundance: 0.0045263797<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMP<br />Abundance: 0.0027662446<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LT<br />Abundance: 0.0053161107<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LP<br />Abundance: 0.0037173249<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0032130055<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LT<br />Abundance: 0.0043393077<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LT<br />Abundance: 0.0042861176<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMP<br />Abundance: 0.0020052631<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMP<br />Abundance: 0.0019516464<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0025691975<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LT<br />Abundance: 0.0035697361<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0025119334<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMP<br />Abundance: 0.0017529586<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LP<br />Abundance: 0.0027543979<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMP<br />Abundance: 0.0015864288<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LSJ<br />Abundance: 0.0012137763<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LT<br />Abundance: 0.0029931131<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LP<br />Abundance: 0.0023418426<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0020017351<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LSJ<br />Abundance: 0.0010915234<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LSJ<br />Abundance: 0.0009253274<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LSJ<br />Abundance: 0.0008741786<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0012707582<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CMPW<br />Abundance: 0.0012292302<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: BSJ<br />Abundance: 0.0011442717<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: LSJ<br />Abundance: 0.0006121898<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group","Location: CS<br />Abundance: 0.0010015749<br />Class_family: Alphaproteobacteria_AEGEAN-169_marine_group"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(56,89,140,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_AEGEAN-169_marine_group","legendgroup":"Alphaproteobacteria_AEGEAN-169_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.941426612074948,"x":[8],"y":[0.00213222255266232],"text":"Location: LSJ<br />Abundance: 0.0021322226<br />Class_family: Alphaproteobacteria_Candidatus_Hepatincola","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(54,92,141,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Candidatus_Hepatincola","legendgroup":"Alphaproteobacteria_Candidatus_Hepatincola","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.962234578617884,"x":[2],"y":[0.00111524866461721],"text":"Location: CMP<br />Abundance: 0.0011152487<br />Class_family: Alphaproteobacteria_Caulobacteraceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(52,96,141,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Caulobacteraceae","legendgroup":"Alphaproteobacteria_Caulobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9],"base":[0.681183239599809,0.728427926921468,0.772267631105151,0.813643467236753,0.854298540017431,0.88253404325821,0.834317659026866,0.906286349658444,0.929144943404404,0.867392929077353,0.817462835582447,0.838572325589171,0.949562640954699,0.857903738683683,0.894203962604928,0.912300784615951,0.872489033550782,0.92023549327365,0.937343878442655,0.964388153784798,0.945069525599912,0.930247824039815,0.926265220597441,0.942678729480388,0.930126176242765,0.930052446202672,0.933699655964787,0.952266020604494,0.935862881072854,0.936989708895259,0.951344819537685,0.971386777468897,0.954218672537212,0.937801336313141,0.956944893848053,0.940678792991312,0.975339596271761,0.977671100893266,0.885331396024341,0.939400047433051,0.979926571947787,0.942886810369088,0.887847242843352,0.959124247795423,0.944543232597635,0.889838045250908,0.940596246687365,0.960365625354435,0.961365099009997,0.981730693937685],"x":[1,1,1,1,1,1,9,1,1,9,7,7,1,7,9,9,7,8,2,1,2,9,8,9,8,5,8,9,5,8,2,4,2,3,2,5,4,4,7,8,4,5,7,2,5,7,8,2,2,4],"y":[0.0472446873216595,0.0438397041836831,0.041375836131602,0.0406550727806776,0.0282355032407788,0.0237523064002345,0.0330752700504865,0.0228585937459599,0.020417697550295,0.0268110335275755,0.0211094900067244,0.0193314130945118,0.0148255128300985,0.0145852948670985,0.0180968220110229,0.017947039423864,0.0128423624735594,0.00602972732379126,0.00772564715725621,0.00892519469939035,0.00627529393777304,0.0124309054405731,0.00386095564532363,0.0095872911241055,0.00357347972202149,0.00581043487018273,0.00329005293047235,0.00748077446430595,0.00481591191845754,0.00241033853779227,0.00287385299952692,0.00395281880286413,0.00272622131084188,0.00445753163892948,0.00217935394736934,0.00220801737777632,0.00233150462150467,0.00225547105452151,0.00251584681901063,0.00119619925431336,0.00180412198989766,0.00165642222854701,0.00199080240755622,0.00124137755901232,0.00146059703103452,0.00171958998831889,0.00083036538758352,0.000999473655561745,0.000869479607886636,0.00118731656523552],"text":["Location: BSJ<br />Abundance: 0.0472446873<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0438397042<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0413758361<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0406550728<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0282355032<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0237523064<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0330752701<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0228585937<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0204176976<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0268110335<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0211094900<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0193314131<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0148255128<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0145852949<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0180968220<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0179470394<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0128423625<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0060297273<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0077256472<br />Class_family: Alphaproteobacteria_Clade_I","Location: BSJ<br />Abundance: 0.0089251947<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0062752939<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0124309054<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0038609556<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0095872911<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0035734797<br />Class_family: Alphaproteobacteria_Clade_I","Location: CS<br />Abundance: 0.0058104349<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0032900529<br />Class_family: Alphaproteobacteria_Clade_I","Location: LT<br />Abundance: 0.0074807745<br />Class_family: Alphaproteobacteria_Clade_I","Location: CS<br />Abundance: 0.0048159119<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0024103385<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0028738530<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPW<br />Abundance: 0.0039528188<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0027262213<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPE<br />Abundance: 0.0044575316<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0021793539<br />Class_family: Alphaproteobacteria_Clade_I","Location: CS<br />Abundance: 0.0022080174<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPW<br />Abundance: 0.0023315046<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPW<br />Abundance: 0.0022554711<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0025158468<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0011961993<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPW<br />Abundance: 0.0018041220<br />Class_family: Alphaproteobacteria_Clade_I","Location: CS<br />Abundance: 0.0016564222<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0019908024<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0012413776<br />Class_family: Alphaproteobacteria_Clade_I","Location: CS<br />Abundance: 0.0014605970<br />Class_family: Alphaproteobacteria_Clade_I","Location: LP<br />Abundance: 0.0017195900<br />Class_family: Alphaproteobacteria_Clade_I","Location: LSJ<br />Abundance: 0.0008303654<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0009994737<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMP<br />Abundance: 0.0008694796<br />Class_family: Alphaproteobacteria_Clade_I","Location: CMPW<br />Abundance: 0.0011873166<br />Class_family: Alphaproteobacteria_Clade_I"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(51,99,141,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_I","legendgroup":"Alphaproteobacteria_Clade_I","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.865207040343424,0.890836115445565,0.875504561682106,0.902521155864039,0.88220659789663,0.888735296680537,0.895170407574376,0.901208288687652,0.906500093948874,0.626921052027027,0.635818087806634,0.644318920194283,0.652061766721834,0.913645080516823,0.911548389871896,0.811465918490365,0.955627456564556,0.818852631809187,0.659301711189626,0.804429192116198,0.664346083525581,0.914433955681266,0.810019396957629,0.66911690993523,0.927609400582945,0.825967426562955,0.672994248616113,0.676259762008392,0.966749710309116,0.917041400289539,0.918659548224317,0.919029450604316,0.921642049991983,0.814545776609027,0.830863064583775,0.932509113134076,0.932495634037589,0.924194251119896,0.926262563897647,0.934184289122525,0.972851221744083,0.679241803395867,0.928279077915009,0.969675496157901,0.935523024314745,0.935652984174013,0.936518181871896],"x":[8,5,8,5,8,8,8,8,8,1,1,1,1,5,8,9,6,9,1,7,1,8,7,1,3,9,1,1,4,8,8,5,5,7,9,2,3,5,5,2,6,1,5,4,3,2,2],"y":[0.0102975213386826,0.0116850404184737,0.00670203621452392,0.0111239246527841,0.00652869878390627,0.00643511089383952,0.00603788111327574,0.00529180526122175,0.00504829592302203,0.00889703577960754,0.00850083238764854,0.00774284652755119,0.00723994446779186,0.00538437008749348,0.00288556580937005,0.00738671331882179,0.0172237651795266,0.00711479475376864,0.00504437233595523,0.00559020484143169,0.00477082640964888,0.00260744460827289,0.0045263796513979,0.0038773386808828,0.00488623345464423,0.00489563802081971,0.00326551339227987,0.00298204138747482,0.002925785848785,0.00161814793477799,0.00157594504933367,0.00261259938766656,0.00255220112791266,0.00291705897341976,0.00345459444309126,0.00167517598844935,0.00302739027715571,0.00206831277775144,0.00201651401736158,0.00146869505148806,0.00652743788685994,0.00194143620394138,0.00177336828766306,0.00171128131099618,0.0022783119983959,0.000865197697882403,0.000825696570759593],"text":["Location: LSJ<br />Abundance: 0.0102975213<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0116850404<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0067020362<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0111239247<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0065286988<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0064351109<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0060378811<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0052918053<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0050482959<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0088970358<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0085008324<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0077428465<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0072399445<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0053843701<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0028855658<br />Class_family: Alphaproteobacteria_Clade_II","Location: LT<br />Abundance: 0.0073867133<br />Class_family: Alphaproteobacteria_Clade_II","Location: LLC<br />Abundance: 0.0172237652<br />Class_family: Alphaproteobacteria_Clade_II","Location: LT<br />Abundance: 0.0071147948<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0050443723<br />Class_family: Alphaproteobacteria_Clade_II","Location: LP<br />Abundance: 0.0055902048<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0047708264<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0026074446<br />Class_family: Alphaproteobacteria_Clade_II","Location: LP<br />Abundance: 0.0045263797<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0038773387<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMPE<br />Abundance: 0.0048862335<br />Class_family: Alphaproteobacteria_Clade_II","Location: LT<br />Abundance: 0.0048956380<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0032655134<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0029820414<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMPW<br />Abundance: 0.0029257858<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0016181479<br />Class_family: Alphaproteobacteria_Clade_II","Location: LSJ<br />Abundance: 0.0015759450<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0026125994<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0025522011<br />Class_family: Alphaproteobacteria_Clade_II","Location: LP<br />Abundance: 0.0029170590<br />Class_family: Alphaproteobacteria_Clade_II","Location: LT<br />Abundance: 0.0034545944<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMP<br />Abundance: 0.0016751760<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMPE<br />Abundance: 0.0030273903<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0020683128<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0020165140<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMP<br />Abundance: 0.0014686951<br />Class_family: Alphaproteobacteria_Clade_II","Location: LLC<br />Abundance: 0.0065274379<br />Class_family: Alphaproteobacteria_Clade_II","Location: BSJ<br />Abundance: 0.0019414362<br />Class_family: Alphaproteobacteria_Clade_II","Location: CS<br />Abundance: 0.0017733683<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMPW<br />Abundance: 0.0017112813<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMPE<br />Abundance: 0.0022783120<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMP<br />Abundance: 0.0008651977<br />Class_family: Alphaproteobacteria_Clade_II","Location: CMP<br />Abundance: 0.0008256966<br />Class_family: Alphaproteobacteria_Clade_II"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(49,103,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_II","legendgroup":"Alphaproteobacteria_Clade_II","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.824089296931328,0.716061265372769,0.834465126417073,0.73498874726137,0.752873176294845,0.76775520810244,0.75737393765381,0.858843152066175,0.843231195523751,0.867414432123213,0.781680112168893,0.772512140257815,0.847630419564009,0.850869065417442,0.854052567275963,0.781123096278327,0.960454278806537,0.914433580911874,0.788262568392023,0.794581948470848,0.800869695121071,0.789495946288049,0.79400350517211,0.874352550076753,0.877803576775075,0.920543708405826,0.806705425649443,0.856990685732687,0.858813465548784,0.613377683727353,0.616472749883455,0.798324302350881,0.860596621990068,0.801494803909991,0.881124513182479,0.88347091843455,0.885779300238589,0.928480170054694,0.619521147554104,0.930142260369047,0.621844365844346,0.862126259908283,0.863200339193976,0.86425892873981,0.925301035984214,0.623961574893504,0.888020920237236,0.950208451526408,0.889525138380784,0.965407168014758,0.625655322571405,0.931753087199711],"x":[8,7,8,7,7,7,9,5,8,5,7,9,8,8,8,9,4,3,9,9,9,7,7,5,5,3,9,8,8,1,1,7,8,7,5,5,5,2,1,2,1,8,8,8,3,1,5,6,5,4,1,2],"y":[0.0103758294857447,0.0189274818886012,0.00876606910667799,0.0178844290334749,0.0148820318075953,0.0139249040664532,0.0151382026040052,0.00857128005703822,0.00439922404025772,0.00693811795354049,0.00781583411915576,0.00861095602051221,0.0032386458534337,0.00318350185852101,0.00293811845672343,0.00713947211369548,0.00495288920822046,0.00611012749395268,0.00631938007882538,0.0062877466502228,0.00583573052837194,0.00450755888406051,0.00432079717877121,0.00345102669832176,0.00332093640740438,0.00475732757838743,0.00476049284092195,0.0018227798160968,0.00178315644128491,0.00309506615610144,0.00304839767064913,0.00317050155911058,0.00152963791821425,0.00293438820620617,0.00234640525207097,0.00230838180403892,0.00224161999864625,0.00166209031435327,0.00232321829024185,0.00161082683066438,0.00211720904915802,0.00107407928569359,0.00105858954583404,0.000948111603613411,0.0023083645987314,0.00169374767790176,0.0015042181435484,0.00541900503814796,0.001310977064781,0.00134254229435815,0.0012657294556212,0.000756025934364635],"text":["Location: LSJ<br />Abundance: 0.0103758295<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0189274819<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0087660691<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0178844290<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0148820318<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0139249041<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0151382026<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0085712801<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0043992240<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0069381180<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0078158341<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0086109560<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0032386459<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0031835019<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0029381185<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0071394721<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMPW<br />Abundance: 0.0049528892<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMPE<br />Abundance: 0.0061101275<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0063193801<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0062877467<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0058357305<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0045075589<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0043207972<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0034510267<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0033209364<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMPE<br />Abundance: 0.0047573276<br />Class_family: Alphaproteobacteria_Clade_III","Location: LT<br />Abundance: 0.0047604928<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0018227798<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0017831564<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0030950662<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0030483977<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0031705016<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0015296379<br />Class_family: Alphaproteobacteria_Clade_III","Location: LP<br />Abundance: 0.0029343882<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0023464053<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0023083818<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0022416200<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMP<br />Abundance: 0.0016620903<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0023232183<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMP<br />Abundance: 0.0016108268<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0021172090<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0010740793<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0010585895<br />Class_family: Alphaproteobacteria_Clade_III","Location: LSJ<br />Abundance: 0.0009481116<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMPE<br />Abundance: 0.0023083646<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0016937477<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0015042181<br />Class_family: Alphaproteobacteria_Clade_III","Location: LLC<br />Abundance: 0.0054190050<br />Class_family: Alphaproteobacteria_Clade_III","Location: CS<br />Abundance: 0.0013109771<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMPW<br />Abundance: 0.0013425423<br />Class_family: Alphaproteobacteria_Clade_III","Location: BSJ<br />Abundance: 0.0012657295<br />Class_family: Alphaproteobacteria_Clade_III","Location: CMP<br />Abundance: 0.0007560259<br />Class_family: Alphaproteobacteria_Clade_III"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,106,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Clade_III","legendgroup":"Alphaproteobacteria_Clade_III","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.714595328566746,"x":[7],"y":[0.00146593680602247],"text":"Location: LP<br />Abundance: 0.0014659368<br />Class_family: Alphaproteobacteria_Fokiniaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(46,110,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Fokiniaceae","legendgroup":"Alphaproteobacteria_Fokiniaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.823507161241316,"x":[8],"y":[0.000582135690012531],"text":"Location: LSJ<br />Abundance: 0.0005821357<br />Class_family: Alphaproteobacteria_Holosporaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,113,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Holosporaceae","legendgroup":"Alphaproteobacteria_Holosporaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.71332569989242,"x":[7],"y":[0.00126962867432678],"text":"Location: LP<br />Abundance: 0.0012696287<br />Class_family: Alphaproteobacteria_Hyphomicrobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(43,116,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Hyphomicrobiaceae","legendgroup":"Alphaproteobacteria_Hyphomicrobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9],"base":[0.940696820009356,0.819309811630791,0.853526328888502,0.85557717271705,0.820878876579839,0.821887463226161,0.912258828549768,0.857590458923679,0.822789224691131,0.612207166389728,0.755735136615021,0.712116158152511],"x":[6,8,5,5,8,8,3,5,8,1,9,7],"y":[0.00951163151705203,0.00156906494904741,0.00205084382854848,0.00201328620662811,0.00100858664632197,0.000901761464969786,0.00217475236210529,0.00125269314249599,0.00071793655018515,0.00117051733762508,0.00163880103878944,0.00120954173990895],"text":["Location: LLC<br />Abundance: 0.0095116315<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LSJ<br />Abundance: 0.0015690649<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: CS<br />Abundance: 0.0020508438<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: CS<br />Abundance: 0.0020132862<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LSJ<br />Abundance: 0.0010085866<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LSJ<br />Abundance: 0.0009017615<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: CMPE<br />Abundance: 0.0021747524<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: CS<br />Abundance: 0.0012526931<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LSJ<br />Abundance: 0.0007179366<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: BSJ<br />Abundance: 0.0011705173<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LT<br />Abundance: 0.0016388010<br />Class_family: Alphaproteobacteria_Hyphomonadaceae","Location: LP<br />Abundance: 0.0012095417<br />Class_family: Alphaproteobacteria_Hyphomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(42,119,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Hyphomonadaceae","legendgroup":"Alphaproteobacteria_Hyphomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.899999999999999],"base":[0.850357904464474,0.852074889995352,0.754016394062144],"x":[5,5,9],"y":[0.00171698553087773,0.00145143889315036,0.00171874255287663],"text":["Location: CS<br />Abundance: 0.0017169855<br />Class_family: Alphaproteobacteria_Kiloniellaceae","Location: CS<br />Abundance: 0.0014514389<br />Class_family: Alphaproteobacteria_Kiloniellaceae","Location: LT<br />Abundance: 0.0017187426<br />Class_family: Alphaproteobacteria_Kiloniellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,123,142,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Kiloniellaceae","legendgroup":"Alphaproteobacteria_Kiloniellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.709844860325879,"x":[7],"y":[0.00227129782663205],"text":"Location: LP<br />Abundance: 0.0022712978<br />Class_family: Alphaproteobacteria_Micavibrionaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,149,139,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Micavibrionaceae","legendgroup":"Alphaproteobacteria_Micavibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.909439028161433,"x":[3],"y":[0.00281980038833518],"text":"Location: CMPE<br />Abundance: 0.0028198004<br />Class_family: Alphaproteobacteria_Rhizobiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,152,139,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhizobiaceae","legendgroup":"Alphaproteobacteria_Rhizobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.707209729798846,0.70863616738719],"x":[7,7],"y":[0.00142643758834449,0.00120869293868808],"text":["Location: LP<br />Abundance: 0.0014264376<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis","Location: LP<br />Abundance: 0.0012086929<br />Class_family: Alphaproteobacteria_Rhizobiales_Incertae_Sedis"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(30,156,137,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","legendgroup":"Alphaproteobacteria_Rhizobiales_Incertae_Sedis","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999],"base":[0.807227245989656,0.769440327177658,0.848186517501639,0.789910536393574,0.807526975729911,0.798924005180634,0.824446492089783,0.874621944662386,0.651689666846922,0.8400170262142,0.766712784835454,0.830291100464972,0.751303547861372,0.511231295270065,0.853206724682178,0.528739519945609,0.895847538336186,0.854181300454588,0.865191120623489,0.875807886323868,0.783551705241729,0.545425934532534,0.760898549348352,0.556938001284181,0.79468699140094,0.693930143525139,0.767140284072141,0.884795715943335,0.567648019754821,0.77288061020131,0.7783041173078,0.911505489287543,0.891943183786097,0.898536737479237,0.897883273676774,0.905025720213356,0.911319250232927,0.577670574458235,0.921122408204268,0.917490951243004,0.7836706028274,0.929410878732639,0.788140203132033,0.586187705638541,0.804472498777342,0.811583357891598,0.818578815662937,0.825562333339214,0.79250304832162,0.796599507435682,0.923294665559367,0.937343771803868,0.708917375266766,0.944438633561892,0.718948714341709,0.874996787349025,0.594010167361839,0.832430275462725,0.728861462088533,0.737708594198888,0.674690259911026,0.884071048185417,0.600655927740954,0.800632377009387,0.681487884489314,0.606579355305453,0.803841516338593,0.806763763760281,0.838395661231581,0.892128006531434,0.842812376528783,0.898564390902194,0.951467547628086,0.687941831683433,0.956053015270405,0.693030649892336,0.809539973106075,0.927703523838342,0.811796562811681,0.904926502810395,0.746411646005716,0.697718060165064,0.813742659723054,0.815250996850499,0.750389351897632,0.700891156024851,0.703138925794085,0.81674180208462,0.847205987224813,0.705339249235737,0.848845518951453,0.817809050849643,0.818702218358759],"x":[4,2,4,2,2,3,2,4,7,2,5,3,8,1,2,1,4,3,2,2,5,1,8,1,5,9,8,2,1,8,8,4,2,2,6,2,2,1,4,2,8,4,8,1,5,5,5,5,8,8,2,4,9,4,9,3,1,5,9,9,7,3,1,8,7,1,8,8,5,3,5,3,4,7,4,7,8,6,8,3,9,7,8,8,9,7,7,8,5,7,5,8,8],"y":[0.0409592715119835,0.0204702092159161,0.0264354271607466,0.0176164393363375,0.016919516359872,0.0313670952843382,0.015570534124417,0.0212255936738007,0.0230005930641041,0.013189698467978,0.0168389204062755,0.0238901999896165,0.00959500148698034,0.0175082246755447,0.0119843959413108,0.0166864145869252,0.0156579509513568,0.0208154868944362,0.0106167657003792,0.00898782961946709,0.0111352861592109,0.0115120667516465,0.00624173472378908,0.0107100184706403,0.00978550737640149,0.0149872317416267,0.00574032612916864,0.00714746784276177,0.0100225547034137,0.00542350710648987,0.00536648551960073,0.00961691891672445,0.00659355369313996,0.00648898273411902,0.0298202501615685,0.00629353001957134,0.00617170101007614,0.00851713118030639,0.00828847052837112,0.00580371431636395,0.00446960030463261,0.00793289307122957,0.00436284518958718,0.00782246172329748,0.00711085911425668,0.00699545777133848,0.00698351767627747,0.00686794212351105,0.00409645911406231,0.00403286957370452,0.00518550449532629,0.00709486175802398,0.0100313390749434,0.00702891406619399,0.00991274774682349,0.00907426083639218,0.00664576037911557,0.00596538576885597,0.00884713211035548,0.00870305180682773,0.00679762457828814,0.00805695834601727,0.00592342756449815,0.0032091393292063,0.00645394719411829,0.00562781108427546,0.0029222474216879,0.00277620934579359,0.00441671529720133,0.00643638437075955,0.00439361069603028,0.00636211190820091,0.0045854676423186,0.00508881820890383,0.00440126353613213,0.00468741027272745,0.00225658970560672,0.0129932961710139,0.00194609691137226,0.00451252535103852,0.00397770589191626,0.00317309585978698,0.00150833712744558,0.00149080523412026,0.00362704216451193,0.00224776976923458,0.00220032344165166,0.00106724876502295,0.00163953172663978,0.0018704805631089,0.00151238551302113,0.000893167509116721,0.000607593272032192],"text":["Location: CMPW<br />Abundance: 0.0409592715<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0204702092<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0264354272<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0176164393<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0169195164<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0313670953<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0155705341<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0212255937<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0230005931<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0131896985<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0168389204<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0238902000<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0095950015<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0175082247<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0119843959<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0166864146<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0156579510<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0208154869<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0106167657<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0089878296<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0111352862<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0115120668<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0062417347<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0107100185<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0097855074<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0149872317<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0057403261<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0071474678<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0100225547<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0054235071<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0053664855<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0096169189<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0065935537<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0064889827<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LLC<br />Abundance: 0.0298202502<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0062935300<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0061717010<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0085171312<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0082884705<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0058037143<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0044696003<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0079328931<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0043628452<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0078224617<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0071108591<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0069954578<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0069835177<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0068679421<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0040964591<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0040328696<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMP<br />Abundance: 0.0051855045<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0070948618<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0100313391<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0070289141<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0099127477<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0090742608<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0066457604<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0059653858<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0088471321<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0087030518<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0067976246<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0080569583<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0059234276<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0032091393<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0064539472<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: BSJ<br />Abundance: 0.0056278111<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0029222474<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0027762093<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0044167153<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0064363844<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0043936107<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0063621119<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0045854676<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0050888182<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPW<br />Abundance: 0.0044012635<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0046874103<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0022565897<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LLC<br />Abundance: 0.0129932962<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0019460969<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CMPE<br />Abundance: 0.0045125254<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0039777059<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0031730959<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0015083371<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0014908052<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LT<br />Abundance: 0.0036270422<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0022477698<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0022003234<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0010672488<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0016395317<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LP<br />Abundance: 0.0018704806<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: CS<br />Abundance: 0.0015123855<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0008931675<br />Class_family: Alphaproteobacteria_Rhodobacteraceae","Location: LSJ<br />Abundance: 0.0006075933<br />Class_family: Alphaproteobacteria_Rhodobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,159,136,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rhodobacteraceae","legendgroup":"Alphaproteobacteria_Rhodobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.746807339435106,"x":[8],"y":[0.00449620842626608],"text":"Location: LSJ<br />Abundance: 0.0044962084<br />Class_family: Alphaproteobacteria_Rickettsiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,162,135,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Rickettsiaceae","legendgroup":"Alphaproteobacteria_Rickettsiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9],"base":[0.473026398609364,0.481487074173827,0.48748787698601,0.493193554231862,0.675510258319621,0.498662917713192,0.682917740674678,0.503733809992252,0.689279249954029,0.507992016956289,0.509673474401608],"x":[1,1,1,1,9,1,9,1,9,1,1],"y":[0.00846067556446306,0.00600080281218268,0.0057056772458518,0.00546936348132998,0.00740748235505695,0.00507089227906071,0.00636150927935086,0.00425820696403678,0.00465089357111004,0.00168145744531845,0.00155782086845679],"text":["Location: BSJ<br />Abundance: 0.0084606756<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0060008028<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0057056772<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0054693635<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: LT<br />Abundance: 0.0074074824<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0050708923<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: LT<br />Abundance: 0.0063615093<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0042582070<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: LT<br />Abundance: 0.0046508936<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0016814574<br />Class_family: Alphaproteobacteria_SAR116_clade","Location: BSJ<br />Abundance: 0.0015578209<br />Class_family: Alphaproteobacteria_SAR116_clade"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(33,165,133,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_SAR116_clade","legendgroup":"Alphaproteobacteria_SAR116_clade","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.76427018538112,0.789832595163351,0.795130821740097,0.763926321746054,0.471790939206856,0.765582743974601],"x":[2,3,3,5,1,5],"y":[0.00517014179653774,0.00529822657674606,0.00379318344053625,0.00165642222854701,0.00123545940250824,0.00113004086085278],"text":["Location: CMP<br />Abundance: 0.0051701418<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Location: CMPE<br />Abundance: 0.0052982266<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Location: CMPE<br />Abundance: 0.0037931834<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Location: CS<br />Abundance: 0.0016564222<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Location: BSJ<br />Abundance: 0.0012354594<br />Class_family: Alphaproteobacteria_Sphingomonadaceae","Location: CS<br />Abundance: 0.0011300409<br />Class_family: Alphaproteobacteria_Sphingomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(35,169,131,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Sphingomonadaceae","legendgroup":"Alphaproteobacteria_Sphingomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.649617318665878,"x":[7],"y":[0.00207234818104396],"text":"Location: LP<br />Abundance: 0.0020723482<br />Class_family: Alphaproteobacteria_Terasakiellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(37,172,130,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Terasakiellaceae","legendgroup":"Alphaproteobacteria_Terasakiellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.647424016310843,"x":[7],"y":[0.00219330235503479],"text":"Location: LP<br />Abundance: 0.0021933024<br />Class_family: Alphaproteobacteria_Thalassobaculaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(41,175,127,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Alphaproteobacteria_Thalassobaculaceae","legendgroup":"Alphaproteobacteria_Thalassobaculaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.637698059911068,0.644270327765183],"x":[7,7],"y":[0.00657226785411558,0.00315368854565989],"text":["Location: LP<br />Abundance: 0.0065722679<br />Class_family: Anaerolineae_Anaerolineaceae","Location: LP<br />Abundance: 0.0031536885<br />Class_family: Anaerolineae_Anaerolineaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(45,178,125,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Anaerolineae_Anaerolineaceae","legendgroup":"Anaerolineae_Anaerolineaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9],"base":[0.757272761920848,0.739864856841606,0.742026308599841,0.743966117932553,0.785114472605909,0.761124546740885,0.745833933487151,0.636036664864242,0.787929155752924,0.762819745441991],"x":[5,8,8,8,3,5,8,7,3,5],"y":[0.00385178482003712,0.00216145175823501,0.00193980933271198,0.00186781555459836,0.00281468314701427,0.0016951987011059,0.000973405947954586,0.00166139504682561,0.00190343941042759,0.00110657630406308],"text":["Location: CS<br />Abundance: 0.0038517848<br />Class_family: Anaerolineae_Caldilineaceae","Location: LSJ<br />Abundance: 0.0021614518<br />Class_family: Anaerolineae_Caldilineaceae","Location: LSJ<br />Abundance: 0.0019398093<br />Class_family: Anaerolineae_Caldilineaceae","Location: LSJ<br />Abundance: 0.0018678156<br />Class_family: Anaerolineae_Caldilineaceae","Location: CMPE<br />Abundance: 0.0028146831<br />Class_family: Anaerolineae_Caldilineaceae","Location: CS<br />Abundance: 0.0016951987<br />Class_family: Anaerolineae_Caldilineaceae","Location: LSJ<br />Abundance: 0.0009734059<br />Class_family: Anaerolineae_Caldilineaceae","Location: LP<br />Abundance: 0.0016613950<br />Class_family: Anaerolineae_Caldilineaceae","Location: CMPE<br />Abundance: 0.0019034394<br />Class_family: Anaerolineae_Caldilineaceae","Location: CS<br />Abundance: 0.0011065763<br />Class_family: Anaerolineae_Caldilineaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(50,182,122,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Anaerolineae_Caldilineaceae","legendgroup":"Anaerolineae_Caldilineaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.763350646448518,"x":[2],"y":[0.000919538932601549],"text":"Location: CMP<br />Abundance: 0.0009195389<br />Class_family: Anaerolineae_SBR1031 Order","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(55,184,120,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Anaerolineae_SBR1031 Order","legendgroup":"Anaerolineae_SBR1031 Order","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.633824231053749,"x":[7],"y":[0.00221243381049363],"text":"Location: LP<br />Abundance: 0.0022124338<br />Class_family: Anaerolineae_uncultured","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(60,188,116,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Anaerolineae_uncultured","legendgroup":"Anaerolineae_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.806124707585215,"x":[4],"y":[0.00110253840444019],"text":"Location: CMPW<br />Abundance: 0.0011025384<br />Class_family: Bacilli_Erysipelotrichaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(87,199,102,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacilli_Erysipelotrichaceae","legendgroup":"Bacilli_Erysipelotrichaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.75588546584446,"x":[5],"y":[0.00138729607638755],"text":"Location: CS<br />Abundance: 0.0013872961<br />Class_family: Bacilli_Mycoplasmataceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(94,201,98,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacilli_Mycoplasmataceae","legendgroup":"Bacilli_Mycoplasmataceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0.782776346305262,0.805015990116119,0.762567362301366],"x":[3,4,2],"y":[0.00233812630064789,0.00110871746909602,0.000783284147152119],"text":["Location: CMPE<br />Abundance: 0.0023381263<br />Class_family: Bacteroidia_Bacteroidaceae","Location: CMPW<br />Abundance: 0.0011087175<br />Class_family: Bacteroidia_Bacteroidaceae","Location: CMP<br />Abundance: 0.0007832841<br />Class_family: Bacteroidia_Bacteroidaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(162,218,55,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Bacteroidaceae","legendgroup":"Bacteroidia_Bacteroidaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.622741698886065,0.631127006148213],"x":[7,7],"y":[0.00838530726214759,0.00269722490553592],"text":["Location: LP<br />Abundance: 0.0083853073<br />Class_family: Bacteroidia_Bacteroidetes_BD2-2","Location: LP<br />Abundance: 0.0026972249<br />Class_family: Bacteroidia_Bacteroidetes_BD2-2"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(218,227,25,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Bacteroidetes_BD2-2","legendgroup":"Bacteroidia_Bacteroidetes_BD2-2","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.899999999999999],"base":[0.801218630151981,0.803725263338379,0.893683867661512,0.739246510104525],"x":[4,4,6,8],"y":[0.00250663318639877,0.00129072677773989,0.00419940601526214,0.000618346737080722],"text":["Location: CMPW<br />Abundance: 0.0025066332<br />Class_family: Bacteroidia_Bacteroidetes_VC2.1_Bac22","Location: CMPW<br />Abundance: 0.0012907268<br />Class_family: Bacteroidia_Bacteroidetes_VC2.1_Bac22","Location: LLC<br />Abundance: 0.0041994060<br />Class_family: Bacteroidia_Bacteroidetes_VC2.1_Bac22","Location: LSJ<br />Abundance: 0.0006183467<br />Class_family: Bacteroidia_Bacteroidetes_VC2.1_Bac22"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(228,228,25,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Bacteroidetes_VC2.1_Bac22","legendgroup":"Bacteroidia_Bacteroidetes_VC2.1_Bac22","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9],"base":[0.469428588797526,0.672295218704235,0.673998990751075,0.470720440249417],"x":[1,9,9,1],"y":[0.00129185145189092,0.00170377204683925,0.00151126756854658,0.00107049895743944],"text":["Location: BSJ<br />Abundance: 0.0012918515<br />Class_family: Bacteroidia_Crocinitomicaceae","Location: LT<br />Abundance: 0.0017037720<br />Class_family: Bacteroidia_Crocinitomicaceae","Location: LT<br />Abundance: 0.0015112676<br />Class_family: Bacteroidia_Crocinitomicaceae","Location: BSJ<br />Abundance: 0.0010704990<br />Class_family: Bacteroidia_Crocinitomicaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(236,229,27,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Crocinitomicaceae","legendgroup":"Bacteroidia_Crocinitomicaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.614438369537265,0.677690193928774,0.684343664663483,0.70106386996514,0.825840209270226,0.689841612846934,0.695085134704875,0.700303981165836,0.705237394430757,0.710200014267599,0.709907120920091,0.636402416055166,0.714441663658684,0.718239121608885,0.858488241774702,0.721982032258473,0.647510232683984,0.7178310610715,0.723674520302344,0.745012859139227,0.729488522097543,0.760029329724656,0.608931445758987,0.725636521535406,0.783464817589077,0.734779250966124,0.739237915587176,0.788641517813503,0.728616255071827,0.749466321553416,0.656535619742917,0.73103938966621,0.662611174813551,0.743559017052951,0.767750411313516,0.881060049106735,0.752585234742437,0.755156689777088,0.747279613264322,0.733400808722444,0.772965663471641,0.735173319453632,0.615375250123824,0.461051154873825,0.736830096107905,0.668629978258516,0.750376823081202,0.618537252980508,0.77710804892327,0.757685421841743,0.793277205102617,0.795290794439656,0.797291306394764,0.799291694413357,0.463689191319579,0.738267028697743,0.752757357766813,0.465497145114366,0.780480764558312,0.75439115105344,0.759272474990735,0.760142504835838,0.621449824174028,0.760993870285227,0.467249693591379,0.761791300404106,0.46835204592095],"x":[9,8,8,5,6,8,8,8,8,5,8,9,8,8,6,8,9,5,5,2,5,3,7,8,4,5,5,4,8,2,9,8,9,5,3,6,2,2,5,8,3,8,7,1,8,9,5,7,3,2,4,4,4,4,1,8,5,1,3,5,2,2,7,2,1,2,1],"y":[0.0219640465179012,0.00665347073470857,0.00549794818345117,0.00913614430245835,0.0326480325044759,0.00524352185794053,0.00521884646096193,0.00493341326492003,0.00466972648933495,0.00763104680390114,0.00453454273859244,0.0111078166288181,0.00379745795020148,0.00374291064958743,0.0225718073320337,0.0036544892769329,0.00902538705893252,0.00584345923084439,0.00581400179519909,0.00445346241418942,0.00529072886858073,0.00772108158885998,0.00644380436483727,0.00297973353642178,0.00517670022442673,0.00445866462105216,0.00432110146577469,0.00463568728911312,0.00242313459438259,0.00311891318902069,0.00607555507063384,0.00236141905623388,0.00601880344496586,0.00372059621137144,0.00521525215812535,0.0126238185547763,0.00257145503465139,0.00252873206465476,0.00309720981688011,0.00177251073118767,0.00414238545162915,0.00165677665427366,0.003162002856684,0.00263803644575417,0.00143693258983768,0.00366524044571881,0.00238053468561017,0.00291257119352051,0.00337271563504116,0.00158705314899243,0.00201358933703921,0.00200051195510798,0.00200038801859359,0.00192693573862335,0.0018079537947866,0.000979481406781946,0.00163379328662727,0.00175254847701384,0.0022955817469501,0.00149431479102047,0.000870029845102316,0.000851365449388908,0.00129187471203662,0.000797430118879672,0.00110235232957029,0.000776061897259761,0.00107654287657588],"text":["Location: LT<br />Abundance: 0.0219640465<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0066534707<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0054979482<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0091361443<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LLC<br />Abundance: 0.0326480325<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0052435219<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0052188465<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0049334133<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0046697265<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0076310468<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0045345427<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LT<br />Abundance: 0.0111078166<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0037974580<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0037429106<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LLC<br />Abundance: 0.0225718073<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0036544893<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LT<br />Abundance: 0.0090253871<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0058434592<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0058140018<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0044534624<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0052907289<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPE<br />Abundance: 0.0077210816<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LP<br />Abundance: 0.0064438044<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0029797335<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0051767002<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0044586646<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0043211015<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0046356873<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0024231346<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0031189132<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LT<br />Abundance: 0.0060755551<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0023614191<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LT<br />Abundance: 0.0060188034<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0037205962<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPE<br />Abundance: 0.0052152522<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LLC<br />Abundance: 0.0126238186<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0025714550<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0025287321<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0030972098<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0017725107<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPE<br />Abundance: 0.0041423855<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0016567767<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LP<br />Abundance: 0.0031620029<br />Class_family: Bacteroidia_Cryomorphaceae","Location: BSJ<br />Abundance: 0.0026380364<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0014369326<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LT<br />Abundance: 0.0036652404<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0023805347<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LP<br />Abundance: 0.0029125712<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPE<br />Abundance: 0.0033727156<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0015870531<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0020135893<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0020005120<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0020003880<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPW<br />Abundance: 0.0019269357<br />Class_family: Bacteroidia_Cryomorphaceae","Location: BSJ<br />Abundance: 0.0018079538<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LSJ<br />Abundance: 0.0009794814<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0016337933<br />Class_family: Bacteroidia_Cryomorphaceae","Location: BSJ<br />Abundance: 0.0017525485<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMPE<br />Abundance: 0.0022955817<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CS<br />Abundance: 0.0014943148<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0008700298<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0008513654<br />Class_family: Bacteroidia_Cryomorphaceae","Location: LP<br />Abundance: 0.0012918747<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0007974301<br />Class_family: Bacteroidia_Cryomorphaceae","Location: BSJ<br />Abundance: 0.0011023523<br />Class_family: Bacteroidia_Cryomorphaceae","Location: CMP<br />Abundance: 0.0007760619<br />Class_family: Bacteroidia_Cryomorphaceae","Location: BSJ<br />Abundance: 0.0010765429<br />Class_family: Bacteroidia_Cryomorphaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(245,230,31,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cryomorphaceae","legendgroup":"Bacteroidia_Cryomorphaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.9,0.899999999999999],"base":[0.674879838645244,0.756153372604177,0.605683083486263,0.743495619900434,0.677058547829332],"x":[8,3,7,2,8],"y":[0.00217870918408802,0.0038759571204785,0.00324836227272385,0.00151723923879288,0.000631646099441885],"text":["Location: LSJ<br />Abundance: 0.0021787092<br />Class_family: Bacteroidia_Cyclobacteriaceae","Location: CMPE<br />Abundance: 0.0038759571<br />Class_family: Bacteroidia_Cyclobacteriaceae","Location: LP<br />Abundance: 0.0032483623<br />Class_family: Bacteroidia_Cyclobacteriaceae","Location: CMP<br />Abundance: 0.0015172392<br />Class_family: Bacteroidia_Cyclobacteriaceae","Location: LSJ<br />Abundance: 0.0006316461<br />Class_family: Bacteroidia_Cyclobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(191,62,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Cyclobacteriaceae","legendgroup":"Bacteroidia_Cyclobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.509779746183748,0.545869689285711,0.721894121757551,0.56647953138514,0.682106826799172,0.734938677518094,0.690998948649216,0.69946823633215,0.403801433298289,0.568577328923666,0.746806788139818,0.707768095368097,0.584392161823532,0.59776687980517,0.757635115214766,0.715339450000101,0.415520648202082,0.42407961089275,0.721408406955273,0.432064421374288,0.766143602406385,0.4389156147019,0.653536860073649,0.445482526249012,0.73990505518094,0.451188798269934,0.581496370845364,0.656833155292048,0.681360256812569,0.659458779180564,0.726629186991384,0.772868214881989,0.587162118995464,0.661927126848195,0.729762936411475,0.664003040317154,0.591788427085865,0.732385026358438,0.747507006795212,0.68560933992058,0.688515159277059,0.734743068136966,0.61032138156177,0.666003346056141,0.691409922217215,0.667560663764344,0.669050344010099,0.59558255699772,0.736808541413477,0.670487276599937,0.456488569085176,0.738573153320289,0.77707037226333,0.458811787375418,0.671772453525903,0.693919032873779,0.672936724905928,0.598527853710269,0.740209709250539,0.751780969539432,0.779313699981105,0.695810800250671,0.60080810680936,0.697476200537109,0.82054436343749,0.67405104668566,0.602723122623972,0.741612240255527,0.754368898156901,0.742571763335788,0.781117853733825,0.604384226613447,0.699003335995803,0.782341090728788,0.700057769960931],"x":[9,9,4,9,2,4,2,2,1,7,4,2,9,9,4,2,1,1,2,1,4,1,8,1,3,1,7,8,5,8,2,4,7,8,2,8,7,2,3,5,5,2,9,8,5,8,8,7,2,8,1,2,4,1,8,5,8,7,2,3,4,5,7,5,6,8,7,2,3,2,4,7,5,4,5],"y":[0.0360899431019631,0.0206098420994288,0.013044555760543,0.0179126304383916,0.00889212185004418,0.0118681106217242,0.0084692876829342,0.0082998590359471,0.0117192149037921,0.0129190419216981,0.0108283270749476,0.00757135463200431,0.013374717981638,0.0125545017565998,0.00850848719161856,0.00606895695517151,0.00855896269066891,0.00798481048153726,0.00522078003611082,0.00685119332761225,0.0067246124756043,0.00656691154711248,0.00329629521839847,0.00570627202092183,0.00760195161427191,0.00529977081524213,0.00566574815009979,0.0026256238885165,0.00424908310801175,0.00246834766763093,0.00313374942009126,0.00420215738134067,0.00462630809040143,0.00207591346895875,0.00262208994696345,0.00200030573898691,0.00379412991185524,0.0023580417785275,0.00427396274421976,0.00290581935647882,0.00289476294015578,0.00206547327651141,0.00411698797549531,0.00155731770820344,0.00250911065656434,0.00148968024575502,0.00143693258983768,0.00294529671254884,0.00176461190681132,0.00128517692596586,0.0023232182902419,0.00163655593024981,0.00224332771777525,0.00223936749840659,0.00116427138002495,0.00189176737689201,0.00111432177973192,0.00228025309909097,0.00140253100498844,0.00258792861746937,0.00180415375272047,0.00166540028643802,0.0019150158146114,0.00152713545869332,0.00529584583273557,0.00082879195958474,0.00166110398947494,0.00095952308026126,0.00178447444727592,0.000923856564645664,0.00122323699496241,0.00129885687281628,0.00105443396512817,0.00112372686028905,0.00100610000420953],"text":["Location: LT<br />Abundance: 0.0360899431<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LT<br />Abundance: 0.0206098421<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0130445558<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LT<br />Abundance: 0.0179126304<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0088921219<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0118681106<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0084692877<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0082998590<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0117192149<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0129190419<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0108283271<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0075713546<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LT<br />Abundance: 0.0133747180<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LT<br />Abundance: 0.0125545018<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0085084872<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0060689570<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0085589627<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0079848105<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0052207800<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0068511933<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0067246125<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0065669115<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0032962952<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0057062720<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPE<br />Abundance: 0.0076019516<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0052997708<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0056657482<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0026256239<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0042490831<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0024683477<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0031337494<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0042021574<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0046263081<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0020759135<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0026220899<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0020003057<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0037941299<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0023580418<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPE<br />Abundance: 0.0042739627<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0029058194<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0028947629<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0020654733<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LT<br />Abundance: 0.0041169880<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0015573177<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0025091107<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0014896802<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0014369326<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0029452967<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0017646119<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0012851769<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0023232183<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0016365559<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0022433277<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: BSJ<br />Abundance: 0.0022393675<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0011642714<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0018917674<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0011143218<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0022802531<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0014025310<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPE<br />Abundance: 0.0025879286<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0018041538<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0016654003<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0019150158<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0015271355<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LLC<br />Abundance: 0.0052958458<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LSJ<br />Abundance: 0.0008287920<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0016611040<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0009595231<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPE<br />Abundance: 0.0017844744<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMP<br />Abundance: 0.0009238566<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0012232370<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: LP<br />Abundance: 0.0012988569<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0010544340<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CMPW<br />Abundance: 0.0011237269<br />Class_family: Bacteroidia_Flavobacteriaceae","Location: CS<br />Abundance: 0.0010061000<br />Class_family: Bacteroidia_Flavobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(178,58,238,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Flavobacteriaceae","legendgroup":"Bacteroidia_Flavobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9],"base":[0.644045580598674,0.564217335854084,0.645810965003179,0.647098418528271,0.677092716984184,0.648327433965449,0.50692395733827,0.649453067515852,0.735804902420152,0.650527146801546,0.651401325438223,0.679004044244113,0.652185465842544,0.816688296606253,0.680310987877032,0.73839389332742,0.652929393239635,0.720806242954654,0.567344667706932],"x":[8,7,8,8,5,8,9,8,3,8,8,5,8,6,5,3,8,4,7],"y":[0.00176538440450469,0.00312733185284808,0.00128745352509174,0.00122901543717802,0.0019113272599286,0.00112563355040352,0.00285578884547821,0.00107407928569359,0.00258899090726816,0.000874178636676937,0.000784140404321509,0.00130694363291906,0.000743927397091171,0.00385606683123729,0.00104926893553647,0.00151116185351996,0.000607466834013803,0.00108787880289707,0.00123266121673338],"text":["Location: LSJ<br />Abundance: 0.0017653844<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LP<br />Abundance: 0.0031273319<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0012874535<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0012290154<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CS<br />Abundance: 0.0019113273<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0011256336<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LT<br />Abundance: 0.0028557888<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0010740793<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CMPE<br />Abundance: 0.0025889909<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0008741786<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0007841404<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CS<br />Abundance: 0.0013069436<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0007439274<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LLC<br />Abundance: 0.0038560668<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CS<br />Abundance: 0.0010492689<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CMPE<br />Abundance: 0.0015111619<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LSJ<br />Abundance: 0.0006074668<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: CMPW<br />Abundance: 0.0010878788<br />Class_family: Bacteroidia_NS11-12_marine_group","Location: LP<br />Abundance: 0.0012326612<br />Class_family: Bacteroidia_NS11-12_marine_group"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(154,50,205,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS11-12_marine_group","legendgroup":"Bacteroidia_NS11-12_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999],"base":[0.55936305053155,0.633990072709779,0.502616844767912,0.67309368348812,0.680405679773859,0.636097031962054,0.804214088899509,0.637196621598744,0.810775660798356,0.638270283910024,0.399437000422726,0.675616039990643,0.639219648397575,0.401132927083604,0.640069544294344,0.402509581846398,0.640808845604305,0.641508139480536,0.642176164402126,0.642814713737627,0.643431598662091],"x":[7,8,9,5,2,8,6,8,6,8,1,5,8,1,8,1,8,8,8,8,8],"y":[0.00485428532253429,0.00210695925227533,0.0043071125703581,0.00252235650252275,0.0017011470253131,0.00109958963669021,0.0065615718988471,0.00107366231128014,0.00591263580789714,0.000949364487550342,0.00169592666087748,0.00147667699354148,0.000849895896769426,0.00137665476279486,0.000739301309960472,0.00129185145189098,0.000699293876231244,0.000668024921589905,0.000638549335501226,0.000616884924463657,0.000613981936583685],"text":["Location: LP<br />Abundance: 0.0048542853<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0021069593<br />Class_family: Bacteroidia_NS9_marine_group","Location: LT<br />Abundance: 0.0043071126<br />Class_family: Bacteroidia_NS9_marine_group","Location: CS<br />Abundance: 0.0025223565<br />Class_family: Bacteroidia_NS9_marine_group","Location: CMP<br />Abundance: 0.0017011470<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0010995896<br />Class_family: Bacteroidia_NS9_marine_group","Location: LLC<br />Abundance: 0.0065615719<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0010736623<br />Class_family: Bacteroidia_NS9_marine_group","Location: LLC<br />Abundance: 0.0059126358<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0009493645<br />Class_family: Bacteroidia_NS9_marine_group","Location: BSJ<br />Abundance: 0.0016959267<br />Class_family: Bacteroidia_NS9_marine_group","Location: CS<br />Abundance: 0.0014766770<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0008498959<br />Class_family: Bacteroidia_NS9_marine_group","Location: BSJ<br />Abundance: 0.0013766548<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0007393013<br />Class_family: Bacteroidia_NS9_marine_group","Location: BSJ<br />Abundance: 0.0012918515<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0006992939<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0006680249<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0006385493<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0006168849<br />Class_family: Bacteroidia_NS9_marine_group","Location: LSJ<br />Abundance: 0.0006139819<br />Class_family: Bacteroidia_NS9_marine_group"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(162,28,154,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_NS9_marine_group","legendgroup":"Bacteroidia_NS9_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.718750489528935,0.733073575271106],"x":[4,3],"y":[0.00205575342571884,0.00273132714904678],"text":["Location: CMPW<br />Abundance: 0.0020557534<br />Class_family: Bacteroidia_Paludibacteraceae","Location: CMPE<br />Abundance: 0.0027313271<br />Class_family: Bacteroidia_Paludibacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(166,32,152,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Paludibacteraceae","legendgroup":"Bacteroidia_Paludibacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.716774080127503,"x":[4],"y":[0.00197640940143218],"text":"Location: CMPW<br />Abundance: 0.0019764094<br />Class_family: Bacteroidia_Prevotellaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(171,35,148,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Prevotellaceae","legendgroup":"Bacteroidia_Prevotellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.715624371039342,0.731623740363035],"x":[4,3],"y":[0.00114970908816137,0.00144983490807016],"text":["Location: CMPW<br />Abundance: 0.0011497091<br />Class_family: Bacteroidia_Prolixibacteraceae","Location: CMPE<br />Abundance: 0.0014498349<br />Class_family: Bacteroidia_Prolixibacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(174,40,146,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Prolixibacteraceae","legendgroup":"Bacteroidia_Prolixibacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9],"base":[0.579650261262176,0.587299439417176,0.594187987740353,0.465211039356177,0.600106367255479,0.535630696271288,0.779542578559844,0.637229826032276,0.643765611649431,0.605802554180782,0.609365574872113,0.612882106710118,0.480399278420072,0.671775935729509,0.649819267255485,0.489056282333741,0.654762235102617,0.616256796810938,0.49617107708751,0.618882420699454,0.621253682285884,0.544328587987022,0.623571315495299,0.659349420526446,0.549061763684912,0.662734283341303,0.625648131388929,0.665933300644538,0.627529542365975,0.553152453563453,0.709361495876268,0.629291232036789,0.630845327390882,0.675668505277375,0.72636226148273,0.395434119913218,0.669069917101364,0.632380282232341,0.397694859954027,0.677553560089209,0.556515805109977,0.670850319559793,0.678589524187329,0.729893803061804,0.679529200468938,0.712307552758647,0.558113190733644,0.713488989977682,0.633401967406537,0.672110450689534,0.714579532461934],"x":[8,8,8,9,8,7,6,5,5,8,8,8,9,2,5,9,5,8,9,8,8,7,8,5,7,5,8,5,8,7,4,8,8,2,3,1,5,8,1,2,7,5,2,3,2,4,7,4,8,5,4],"y":[0.00764917815499944,0.00688854832317665,0.00591837951512686,0.0151882390638942,0.00569618692530227,0.00869789171573376,0.0246715103396647,0.00653578561715429,0.00605365560605464,0.00356302069133108,0.00351653183800527,0.00337469010081992,0.00865700391366991,0.00389256954786665,0.00494296784713188,0.00711479475376858,0.00458718542382852,0.0026256238885165,0.00644576768040184,0.0023712615864292,0.00231763320941547,0.00473317569789022,0.00207681589362962,0.00338486281485684,0.00409068987854011,0.00319901730323524,0.00188141097704619,0.00313661645682606,0.00176168967081447,0.00336335154652478,0.00294605688237959,0.00155409535409234,0.00153495484145949,0.00188505481183354,0.00353154157907354,0.00226074004080923,0.00178040245842948,0.00102168517419621,0.00174214046869903,0.00103596409812012,0.00159738562366685,0.0012601311297411,0.000939676281609181,0.00172993730123183,0.000876479304920408,0.00118143721903496,0.00124985979790593,0.00109054248425222,0.00058810530324116,0.000983232798585809,0.00104483857740734],"text":["Location: LSJ<br />Abundance: 0.0076491782<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0068885483<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0059183795<br />Class_family: Bacteroidia_Saprospiraceae","Location: LT<br />Abundance: 0.0151882391<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0056961869<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0086978917<br />Class_family: Bacteroidia_Saprospiraceae","Location: LLC<br />Abundance: 0.0246715103<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0065357856<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0060536556<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0035630207<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0035165318<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0033746901<br />Class_family: Bacteroidia_Saprospiraceae","Location: LT<br />Abundance: 0.0086570039<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMP<br />Abundance: 0.0038925695<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0049429678<br />Class_family: Bacteroidia_Saprospiraceae","Location: LT<br />Abundance: 0.0071147948<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0045871854<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0026256239<br />Class_family: Bacteroidia_Saprospiraceae","Location: LT<br />Abundance: 0.0064457677<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0023712616<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0023176332<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0047331757<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0020768159<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0033848628<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0040906899<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0031990173<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0018814110<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0031366165<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0017616897<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0033633515<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPW<br />Abundance: 0.0029460569<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0015540954<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0015349548<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMP<br />Abundance: 0.0018850548<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPE<br />Abundance: 0.0035315416<br />Class_family: Bacteroidia_Saprospiraceae","Location: BSJ<br />Abundance: 0.0022607400<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0017804025<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0010216852<br />Class_family: Bacteroidia_Saprospiraceae","Location: BSJ<br />Abundance: 0.0017421405<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMP<br />Abundance: 0.0010359641<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0015973856<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0012601311<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMP<br />Abundance: 0.0009396763<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPE<br />Abundance: 0.0017299373<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMP<br />Abundance: 0.0008764793<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPW<br />Abundance: 0.0011814372<br />Class_family: Bacteroidia_Saprospiraceae","Location: LP<br />Abundance: 0.0012498598<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPW<br />Abundance: 0.0010905425<br />Class_family: Bacteroidia_Saprospiraceae","Location: LSJ<br />Abundance: 0.0005881053<br />Class_family: Bacteroidia_Saprospiraceae","Location: CS<br />Abundance: 0.0009832328<br />Class_family: Bacteroidia_Saprospiraceae","Location: CMPW<br />Abundance: 0.0010448386<br />Class_family: Bacteroidia_Saprospiraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(178,43,143,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Saprospiraceae","legendgroup":"Bacteroidia_Saprospiraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.533062223776577,"x":[7],"y":[0.00256847249471182],"text":"Location: LP<br />Abundance: 0.0025684725<br />Class_family: Bacteroidia_SB-5","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,48,139,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_SB-5","legendgroup":"Bacteroidia_SB-5","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999],"base":[0.630812409302593,0.57784102451294,0.634845952858846,0.771931155157181,0.531458673344189,0.463473081652349],"x":[5,8,5,6,7,9],"y":[0.00403354355625318,0.00180923674923639,0.00238387317343014,0.00761142340266241,0.0016035504323878,0.00173795770382862],"text":["Location: CS<br />Abundance: 0.0040335436<br />Class_family: Bacteroidia_Sphingobacteriaceae","Location: LSJ<br />Abundance: 0.0018092367<br />Class_family: Bacteroidia_Sphingobacteriaceae","Location: CS<br />Abundance: 0.0023838732<br />Class_family: Bacteroidia_Sphingobacteriaceae","Location: LLC<br />Abundance: 0.0076114234<br />Class_family: Bacteroidia_Sphingobacteriaceae","Location: LP<br />Abundance: 0.0016035504<br />Class_family: Bacteroidia_Sphingobacteriaceae","Location: LT<br />Abundance: 0.0017379577<br />Class_family: Bacteroidia_Sphingobacteriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(186,51,136,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Sphingobacteriaceae","legendgroup":"Bacteroidia_Sphingobacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.577205303407577,"x":[8],"y":[0.000635721105363252],"text":"Location: LSJ<br />Abundance: 0.0006357211<br />Class_family: Bacteroidia_uncultured","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(190,56,133,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_uncultured","legendgroup":"Bacteroidia_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.723098634247243,0.724806855574573],"x":[3,3],"y":[0.00170822132733017,0.0015554059081565],"text":["Location: CMPE<br />Abundance: 0.0017082213<br />Class_family: Bacteroidia_Weeksellaceae","Location: CMPE<br />Abundance: 0.0015554059<br />Class_family: Bacteroidia_Weeksellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(193,59,130,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bacteroidia_Weeksellaceae","legendgroup":"Bacteroidia_Weeksellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.576527693062391,"x":[8],"y":[0.00067761034518532],"text":"Location: LSJ<br />Abundance: 0.0006776103<br />Class_family: BD2-11_terrestrial_group_BD2-11_terrestrial_group","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(197,63,126,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"BD2-11_terrestrial_group_BD2-11_terrestrial_group","legendgroup":"BD2-11_terrestrial_group_BD2-11_terrestrial_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.461619396829227,"x":[9],"y":[0.00185368482312215],"text":"Location: LT<br />Abundance: 0.0018536848<br />Class_family: BD7-11_BD7-11","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(200,67,123,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"BD7-11_BD7-11","legendgroup":"BD7-11_BD7-11","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9],"base":[0.459381034434783,0.629782355894551],"x":[9,5],"y":[0.00223836239444403,0.00103005340804219],"text":["Location: LT<br />Abundance: 0.0022383624<br />Class_family: Bdellovibrionia_Bacteriovoracaceae","Location: CS<br />Abundance: 0.0010300534<br />Class_family: Bdellovibrionia_Bacteriovoracaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,70,120,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bdellovibrionia_Bacteriovoracaceae","legendgroup":"Bdellovibrionia_Bacteriovoracaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.707345872526064,0.670926647828156],"x":[4,2],"y":[0.00201562335020411,0.000849287901352747],"text":["Location: CMPW<br />Abundance: 0.0020156234<br />Class_family: Bdellovibrionia_Bdellovibrionaceae","Location: CMP<br />Abundance: 0.0008492879<br />Class_family: Bdellovibrionia_Bdellovibrionaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(206,75,117,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Bdellovibrionia_Bdellovibrionaceae","legendgroup":"Bdellovibrionia_Bdellovibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999],"base":[0.404169288842576,0.407069161908266,0.571797446361944,0.472854005187687,0.456193152047791,0.531828559433359,0.485249700314542,0.617437771746672,0.575136496539164,0.497061002425445,0.516905905411614,0.603456997108465,0.53568425162258,0.554193743039947,0.629435483555818,0.572399124773587,0.589479248123014,0.604568905345561,0.619641428911992,0.657576405502228,0.633724742376621,0.682503039540037,0.653957422195682,0.670532370664666,0.646823292521853,0.656709035032635,0.687075602374599,0.607650561430769,0.520238161540656,0.616406015399128,0.705906509528342,0.698643976679715,0.665444654892352,0.567710216290927,0.569954461754181,0.713707189896469,0.572150054886281,0.621913160429636,0.574103523813749,0.704202101778357,0.719049643326743,0.452841739530266,0.669128781270923,0.625124822070038,0.39358749440119,0.457078639776892,0.627246395328559,0.628651976216792,0.530075612895151,0.575926303629846],"x":[2,4,6,4,2,4,7,3,4,2,2,4,2,2,4,2,2,2,2,3,2,3,4,4,2,2,4,5,7,5,3,4,2,8,8,3,8,5,8,4,3,9,2,5,1,9,5,5,7,8],"y":[0.0520238632052147,0.0657848432794209,0.200133708795237,0.0589745542456727,0.0408678503776539,0.0433079371058047,0.0349884612261138,0.0401386337555558,0.0283205005693014,0.0198449029861686,0.0187783462109663,0.0259784864473522,0.0185094914173676,0.0182053817336398,0.0245219386398645,0.0170801233494272,0.0150896572225463,0.0150725235664314,0.0140833134646285,0.0249266340378093,0.013098550145232,0.0234034699883053,0.0165749484689837,0.0165432317099332,0.00988574251078278,0.00873561985971638,0.0115683743051158,0.00875545396835919,0.00983745135449499,0.00550714503050798,0.00780068036812653,0.00555812509864229,0.00368412637857118,0.00224424546325375,0.00219559313210027,0.00534245343027473,0.00195346892746795,0.00321166164040232,0.0018227798160968,0.00314377074770644,0.00404899092049993,0.00423690024662615,0.00179786655723291,0.00212157325852058,0.00184662551202791,0.00230239465789095,0.00140558088823384,0.00113037967775842,0.00138306044903813,0.000601389432545507],"text":["Location: CMP<br />Abundance: 0.0520238632<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0657848433<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LLC<br />Abundance: 0.2001337088<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0589745542<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0408678504<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0433079371<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LP<br />Abundance: 0.0349884612<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0401386338<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0283205006<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0198449030<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0187783462<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0259784864<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0185094914<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0182053817<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0245219386<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0170801233<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0150896572<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0150725236<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0140833135<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0249266340<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0130985501<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0234034700<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0165749485<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0165432317<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0098857425<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0087356199<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0115683743<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0087554540<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LP<br />Abundance: 0.0098374514<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0055071450<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0078006804<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0055581251<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0036841264<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LSJ<br />Abundance: 0.0022442455<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LSJ<br />Abundance: 0.0021955931<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0053424534<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LSJ<br />Abundance: 0.0019534689<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0032116616<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LSJ<br />Abundance: 0.0018227798<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPW<br />Abundance: 0.0031437707<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMPE<br />Abundance: 0.0040489909<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LT<br />Abundance: 0.0042369002<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CMP<br />Abundance: 0.0017978666<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0021215733<br />Class_family: Campylobacteria_Arcobacteraceae","Location: BSJ<br />Abundance: 0.0018466255<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LT<br />Abundance: 0.0023023947<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0014055809<br />Class_family: Campylobacteria_Arcobacteraceae","Location: CS<br />Abundance: 0.0011303797<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LP<br />Abundance: 0.0013830604<br />Class_family: Campylobacteria_Arcobacteraceae","Location: LSJ<br />Abundance: 0.0006013894<br />Class_family: Campylobacteria_Arcobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(209,78,114,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Campylobacteria_Arcobacteraceae","legendgroup":"Campylobacteria_Arcobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.361637164431561,0.374877250739355,0.597003719339897,0.386525204012642,0.599448731513404,0.382123803028723,0.388603776767718,0.394439218825654,0.392180368448142,0.60952449688263,0.395217667734919,0.397391055705984,0.398949603207455,0.401550924132281,0.399511127810479,0.479274064325305,0.449819204393172,0.481823373466982,0.403973269194553,0.606000124198895,0.401207553988593,0.402368296851866,0.613681999535243,0.61563151368215,0.405825977560767,0.483958178020967,0.40341926311404],"x":[2,2,5,2,3,4,4,4,2,3,2,2,4,4,2,7,9,7,4,5,2,2,3,3,4,7,2],"y":[0.0132400863077947,0.0116479532732861,0.00899640485899778,0.00565516443550068,0.0100757653692257,0.00647997373899523,0.00583544205793601,0.00451038438180112,0.0030372992867771,0.00415750265261317,0.00217338797106437,0.00212007210449522,0.00260132092482607,0.00242234506227124,0.0016964261781145,0.00254930914167673,0.00302253513709339,0.00213480455398501,0.00185270836621415,0.0016504372318733,0.00116074286327283,0.0010509662621736,0.00194951414690747,0.00180625806452139,0.00124318434749904,0.0012915222935751,0.000750025728536441],"text":["Location: CMP<br />Abundance: 0.0132400863<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0116479533<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CS<br />Abundance: 0.0089964049<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0056551644<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPE<br />Abundance: 0.0100757654<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0064799737<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0058354421<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0045103844<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0030372993<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPE<br />Abundance: 0.0041575027<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0021733880<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0021200721<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0026013209<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0024223451<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0016964262<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: LP<br />Abundance: 0.0025493091<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: LT<br />Abundance: 0.0030225351<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: LP<br />Abundance: 0.0021348046<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0018527084<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CS<br />Abundance: 0.0016504372<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0011607429<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0010509663<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPE<br />Abundance: 0.0019495141<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPE<br />Abundance: 0.0018062581<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMPW<br />Abundance: 0.0012431843<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: LP<br />Abundance: 0.0012915223<br />Class_family: Campylobacteria_Sulfurimonadaceae","Location: CMP<br />Abundance: 0.0007500257<br />Class_family: Campylobacteria_Sulfurimonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(213,83,111,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Campylobacteria_Sulfurimonadaceae","legendgroup":"Campylobacteria_Sulfurimonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.360387287495224,"x":[2],"y":[0.00124987693633699],"text":"Location: CMP<br />Abundance: 0.0012498769<br />Class_family: Campylobacteria_Sulfurospirillaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(215,86,108,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Campylobacteria_Sulfurospirillaceae","legendgroup":"Campylobacteria_Sulfurospirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.27940148088272,0.507562049735103,0.297600461192022,0.541311932144141,0.337194874314979,0.314309309010296,0.571041368561512,0.326035180792338,0.354989352260573,0.334983699354728,0.366333834647899,0.584364764747651,0.342849137860051,0.570619309214417,0.441557815233039,0.580744530477509,0.591959006334852,0.588246720657392,0.376167676547707,0.348459330309504,0.473926388349041,0.351820703181063,0.354826126763181,0.594869103606172,0.357338292235731,0.359113479768499,0.380916942808649],"x":[2,3,2,3,4,2,5,2,4,2,4,5,2,3,9,3,5,3,4,2,7,2,2,3,2,2,4],"y":[0.0181989803093021,0.0337498824090375,0.0167088478182738,0.0293073770702759,0.0177944779455947,0.0117258717820414,0.0133233961861385,0.00894851856239093,0.0113444823873255,0.00786543850532267,0.0098338418998084,0.00759424158720123,0.00561019244945271,0.010125221263092,0.00826138916013341,0.00750219017988385,0.00504471300504539,0.00662238294877948,0.00474926626094196,0.00336137287155941,0.00534767597626445,0.00300542358211819,0.00251216547254923,0.0045796279072321,0.00177518753276845,0.00127380772672464,0.00120686022007371],"text":["Location: CMP<br />Abundance: 0.0181989803<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0337498824<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0167088478<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0293073771<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPW<br />Abundance: 0.0177944779<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0117258718<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CS<br />Abundance: 0.0133233962<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0089485186<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPW<br />Abundance: 0.0113444824<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0078654385<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPW<br />Abundance: 0.0098338419<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CS<br />Abundance: 0.0075942416<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0056101924<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0101252213<br />Class_family: Campylobacteria_Sulfurovaceae","Location: LT<br />Abundance: 0.0082613892<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0075021902<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CS<br />Abundance: 0.0050447130<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0066223829<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPW<br />Abundance: 0.0047492663<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0033613729<br />Class_family: Campylobacteria_Sulfurovaceae","Location: LP<br />Abundance: 0.0053476760<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0030054236<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0025121655<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPE<br />Abundance: 0.0045796279<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0017751875<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMP<br />Abundance: 0.0012738077<br />Class_family: Campylobacteria_Sulfurovaceae","Location: CMPW<br />Abundance: 0.0012068602<br />Class_family: Campylobacteria_Sulfurovaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(218,91,105,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Campylobacteria_Sulfurovaceae","legendgroup":"Campylobacteria_Sulfurovaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.4712877020982,"x":[7],"y":[0.00263868625084057],"text":"Location: LP<br />Abundance: 0.0026386863<br />Class_family: Chlamydiae_Simkaniaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(221,94,102,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Chlamydiae_Simkaniaceae","legendgroup":"Chlamydiae_Simkaniaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999],"base":[0.310793597543152,0.493661037850516,0.506336268816788,0.518275422928729,0.527920468476761,0.53748569087471,0.545196227749485,0.329235733439882,0.552756773369353,0.42860312491181,0.502519599613786,0.557955424004511,0.434119089383833,0.561355754211612,0.564375651282299,0.566727981237544,0.564382514361026,0.569022928433343,0.438804018846328,0.565678698449369,0.566418575301126,0.567123290060119],"x":[4,5,5,5,5,5,5,4,5,9,3,5,9,5,5,5,8,5,9,8,8,8],"y":[0.0184421358967302,0.0126752309662724,0.0119391541119409,0.00964504554803203,0.00956522239794855,0.00771053687477519,0.00756054561986819,0.00795914087509669,0.00519865063515768,0.00551596447202274,0.00504245012131699,0.00340033020710107,0.00468492946249471,0.00301989707068706,0.00235232995524481,0.00229494719579926,0.00129618408834342,0.00201844012816943,0.00275379638671114,0.000739876851757137,0.000704714758992697,0.00058692623080836],"text":["Location: CMPW<br />Abundance: 0.0184421359<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0126752310<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0119391541<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0096450455<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0095652224<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0077105369<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0075605456<br />Class_family: Chlorobia_Chlorobiaceae","Location: CMPW<br />Abundance: 0.0079591409<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0051986506<br />Class_family: Chlorobia_Chlorobiaceae","Location: LT<br />Abundance: 0.0055159645<br />Class_family: Chlorobia_Chlorobiaceae","Location: CMPE<br />Abundance: 0.0050424501<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0034003302<br />Class_family: Chlorobia_Chlorobiaceae","Location: LT<br />Abundance: 0.0046849295<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0030198971<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0023523300<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0022949472<br />Class_family: Chlorobia_Chlorobiaceae","Location: LSJ<br />Abundance: 0.0012961841<br />Class_family: Chlorobia_Chlorobiaceae","Location: CS<br />Abundance: 0.0020184401<br />Class_family: Chlorobia_Chlorobiaceae","Location: LT<br />Abundance: 0.0027537964<br />Class_family: Chlorobia_Chlorobiaceae","Location: LSJ<br />Abundance: 0.0007398769<br />Class_family: Chlorobia_Chlorobiaceae","Location: LSJ<br />Abundance: 0.0007047148<br />Class_family: Chlorobia_Chlorobiaceae","Location: LSJ<br />Abundance: 0.0005869262<br />Class_family: Chlorobia_Chlorobiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(224,99,99,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Chlorobia_Chlorobiaceae","legendgroup":"Chlorobia_Chlorobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.468141402841067,"x":[7],"y":[0.00314629925713261],"text":"Location: LP<br />Abundance: 0.0031462993<br />Class_family: Cloacimonadia_LK-44f","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,204,204,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cloacimonadia_LK-44f","legendgroup":"Cloacimonadia_LK-44f","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.308899346665072,"x":[4],"y":[0.00189425087807948],"text":"Location: CMPW<br />Abundance: 0.0018942509<br />Class_family: Clostridia_Clostridiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,102,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Clostridia_Clostridiaceae","legendgroup":"Clostridia_Clostridiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0.304735847114637,0.30728447556878,0.501054507597071],"x":[4,4,3],"y":[0.00254862845414328,0.00161487109629205,0.00146509201671496],"text":["Location: CMPW<br />Abundance: 0.0025486285<br />Class_family: Clostridia_Lachnospiraceae","Location: CMPW<br />Abundance: 0.0016148711<br />Class_family: Clostridia_Lachnospiraceae","Location: CMPE<br />Abundance: 0.0014650920<br />Class_family: Clostridia_Lachnospiraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,102,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Clostridia_Lachnospiraceae","legendgroup":"Clostridia_Lachnospiraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.303128603945357,"x":[4],"y":[0.00160724316927946],"text":"Location: CMPW<br />Abundance: 0.0016072432<br />Class_family: Clostridia_Peptostreptococcaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(246,141,69,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Clostridia_Peptostreptococcaceae","legendgroup":"Clostridia_Peptostreptococcaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.300798101349902,0.499448927304781],"x":[4,3],"y":[0.00233050259545536,0.00160558029229046],"text":["Location: CMPW<br />Abundance: 0.0023305026<br />Class_family: Clostridia_Ruminococcaceae","Location: CMPE<br />Abundance: 0.0016055803<br />Class_family: Clostridia_Ruminococcaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(252,165,55,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Clostridia_Ruminococcaceae","legendgroup":"Clostridia_Ruminococcaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9],"base":[0.0703396307332749,0.136236748709325,0.179868362500331,0.341505462684631,0.370275596782205,0.192933682250535,0.258912518044817,0.395851934774149,0.286171635761218,0.313341127757359,0.235930558069938,0.352609543355078,0.267134732693288,0.32642597282522,0.281650993142034,0.417032070776301,0.314280450072233,0.432393984818297,0.451838996478918,0.297025420606802,0.317787694819137,0.345827039185894,0.365675074979405,0.324656687772248,0.374336944351296,0.447678753067858,0.401849987128425,0.46104834696154,0.341443322512622,0.473445276369756,0.399064574637323,0.362184852501836,0.381548633766735,0.400828514576389,0.427557659598788,0.48510225745028,0.454412505362104,0.495712091200602,0.236105286263444,0.349191832706426,0.505706444845459,0.361474999815571,0.419825706728341,0.237922744239711,0.435817105312844,0.515496045723211,0.451507937465941,0.465579411059105,0.479469633708865,0.524256699254802,0.532299748498953,0.386553347502929,0.539923091721802,0.405979135426139,0.428627766118436,0.367218415085746,0.380825646278087,0.254915044581291,0.266011847326646,0.547475023833095,0.479360080603189,0.444074211099803,0.55256432445992,0.277020694140382,0.54396346593873,0.248855723648564,0.557301967357184,0.453849698804622,0.254635749335065,0.48759977851143,0.25994249943202,0.285385874003652,0.264079021011487,0.561433554138779,0.462316715498376,0.267877225236982,0.271351677044794,0.274492632023136,0.290846107407563,0.293838145151956,0.276638806246565,0.425359176471346,0.296781306326519,0.298951084605046,0.278379066258602],"x":[1,1,9,8,8,1,9,8,5,9,1,6,1,3,7,8,7,8,6,1,5,7,3,1,7,8,7,8,5,8,3,5,5,5,3,8,3,8,2,1,8,9,5,4,5,8,5,5,3,8,8,9,8,9,7,1,1,4,4,8,5,7,8,4,6,2,8,7,2,5,2,4,2,8,7,2,2,2,4,4,2,9,4,4,2],"y":[0.0658971179760504,0.0566969335412098,0.0790441555444862,0.0287701340975742,0.0255763379919436,0.0429968758194026,0.0544286097125421,0.021180136002152,0.0316160590579181,0.0481338720582117,0.03120417462335,0.09922945312384,0.0298906879135143,0.0392491021541855,0.0326294569301984,0.015361914041996,0.0315465891136612,0.015284768249561,0.092124469459812,0.027631267165446,0.0236556276934851,0.028509905165402,0.0333894996579176,0.024535144934178,0.0275130427771285,0.013369593893682,0.0267777789900115,0.0123969294082168,0.0207415299892145,0.0116569810805232,0.0284930849614652,0.0193637812648993,0.0192798808096531,0.0189971921519522,0.0268548457633159,0.0106098337503222,0.0250571283467605,0.00999435364485707,0.0127504373851198,0.0180265823793199,0.00978960087775205,0.0250783476873582,0.0159913985845028,0.0169923003415808,0.015690832153097,0.00876065353159072,0.014071473593164,0.0137806695440847,0.0199792935959161,0.00804304924415111,0.0076233432228493,0.0194257879232106,0.00755193211129324,0.019380041045207,0.0154464449813675,0.013607231192341,0.0127618481231031,0.0110968027453543,0.011008846813736,0.00508930062682456,0.00823969790824097,0.00977548770481868,0.00473764289726453,0.00836517986327018,0.0278339804232144,0.00578002568650116,0.00413158678159486,0.00846701669375349,0.00530675009695464,0.00606125933908541,0.00413652157946781,0.00546023340391111,0.00379820422549421,0.00294896022224633,0.00582468734269176,0.00347445180781242,0.00314095497834238,0.00214617422342878,0.00299203774439344,0.00294316117456278,0.00174026001203692,0.00324394844046366,0.00216977827852738,0.0018470167448556,0.00102241462411812],"text":["Location: BSJ<br />Abundance: 0.0658971180<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0566969335<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0790441555<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0287701341<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0255763380<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0429968758<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0544286097<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0211801360<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0316160591<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0481338721<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0312041746<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LLC<br />Abundance: 0.0992294531<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0298906879<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0392491022<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0326294569<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0153619140<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0315465891<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0152847682<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LLC<br />Abundance: 0.0921244695<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0276312672<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0236556277<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0285099052<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0333894997<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0245351449<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0275130428<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0133695939<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0267777790<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0123969294<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0207415300<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0116569811<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0284930850<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0193637813<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0192798808<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0189971922<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0268548458<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0106098338<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0250571283<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0099943536<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0127504374<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0180265824<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0097896009<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0250783477<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0159913986<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0169923003<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0156908322<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0087606535<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0140714736<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0137806695<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPE<br />Abundance: 0.0199792936<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0080430492<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0076233432<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0194257879<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0075519321<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0193800410<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0154464450<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0136072312<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: BSJ<br />Abundance: 0.0127618481<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0110968027<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0110088468<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0050893006<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0082396979<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0097754877<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0047376429<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0083651799<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LLC<br />Abundance: 0.0278339804<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0057800257<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0041315868<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0084670167<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0053067501<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CS<br />Abundance: 0.0060612593<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0041365216<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0054602334<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0037982042<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LSJ<br />Abundance: 0.0029489602<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LP<br />Abundance: 0.0058246873<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0034744518<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0031409550<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0021461742<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0029920377<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0029431612<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0017402600<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: LT<br />Abundance: 0.0032439484<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0021697783<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMPW<br />Abundance: 0.0018470167<br />Class_family: Cyanobacteriia_Cyanobiaceae","Location: CMP<br />Abundance: 0.0010224146<br />Class_family: Cyanobacteriia_Cyanobiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(246,231,38,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Cyanobiaceae","legendgroup":"Cyanobacteriia_Cyanobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.280110766732878,"x":[7],"y":[0.00154022640915619],"text":"Location: LP<br />Abundance: 0.0015402264<br />Class_family: Cyanobacteriia_Limnotrichaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(244,237,39,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Limnotrichaceae","legendgroup":"Cyanobacteriia_Limnotrichaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999],"base":[0.318832626716056,0.310579017964225,0.267900964590836,0.323230985262817,0.274038698746974,0.326706764618762,0.329157249836082,0.331588085089479,0.16891989934631,0.333413288172212,0.333976483574877,0.336091310245332,0.278191327195721,0.281375128638761,0.174511589349933,0.319677485859377,0.323094953856972,0.338003695314981,0.284126515593095,0.339356207560368,0.346266844276337,0.340544817458749,0.234292002024023,0.235282837908935,0.178348782026123],"x":[8,6,5,8,5,8,8,8,9,6,8,8,5,5,9,3,3,8,5,8,6,8,2,2,9],"y":[0.00439835854676096,0.0228342702079876,0.00613773415613866,0.00347577935594484,0.00415262844874686,0.00245048521732022,0.00243083525339671,0.00238839848539779,0.0055916900036227,0.0128535561041244,0.00211482667045521,0.00191238506964947,0.0031838014430397,0.00275138695433397,0.00383719267618979,0.00341746799759401,0.00333101896824833,0.00135251224538657,0.00204512016812353,0.00118860989838082,0.00634269907874135,0.00096064522588224,0.000990835884911528,0.000822448354509536,0.00151958047420803],"text":["Location: LSJ<br />Abundance: 0.0043983585<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LLC<br />Abundance: 0.0228342702<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CS<br />Abundance: 0.0061377342<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0034757794<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CS<br />Abundance: 0.0041526284<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0024504852<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0024308353<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0023883985<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LT<br />Abundance: 0.0055916900<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LLC<br />Abundance: 0.0128535561<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0021148267<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0019123851<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CS<br />Abundance: 0.0031838014<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CS<br />Abundance: 0.0027513870<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LT<br />Abundance: 0.0038371927<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CMPE<br />Abundance: 0.0034174680<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CMPE<br />Abundance: 0.0033310190<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0013525122<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CS<br />Abundance: 0.0020451202<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0011886099<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LLC<br />Abundance: 0.0063426991<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LSJ<br />Abundance: 0.0009606452<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CMP<br />Abundance: 0.0009908359<br />Class_family: Cyanobacteriia_Microcystaceae","Location: CMP<br />Abundance: 0.0008224484<br />Class_family: Cyanobacteriia_Microcystaceae","Location: LT<br />Abundance: 0.0015195805<br />Class_family: Cyanobacteriia_Microcystaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,72,156,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Microcystaceae","legendgroup":"Cyanobacteriia_Microcystaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9],"base":[0.292467482035366,0.30237362578071,0.245312676059244,0.311673990099913,0.25682857561723,0.31587830341667,0.261794978061876,0.164083843126961,0.265536473540618,0.302442668809654,0.278128948041013,0.317575771510364],"x":[8,8,5,8,5,8,5,9,5,6,7,3],"y":[0.00990614374534426,0.00930036431920289,0.0115158995579858,0.0042043133167568,0.00496640244464608,0.00295432329938605,0.00374149547874209,0.00483605621934938,0.00236449105021741,0.00813634915457023,0.00198181869186531,0.00210171434901379],"text":["Location: LSJ<br />Abundance: 0.0099061437<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LSJ<br />Abundance: 0.0093003643<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: CS<br />Abundance: 0.0115158996<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LSJ<br />Abundance: 0.0042043133<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: CS<br />Abundance: 0.0049664024<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LSJ<br />Abundance: 0.0029543233<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: CS<br />Abundance: 0.0037414955<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LT<br />Abundance: 0.0048360562<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: CS<br />Abundance: 0.0023644911<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LLC<br />Abundance: 0.0081363492<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: LP<br />Abundance: 0.0019818187<br />Class_family: Cyanobacteriia_Nodosilineaceae","Location: CMPE<br />Abundance: 0.0021017143<br />Class_family: Cyanobacteriia_Nodosilineaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,204,204,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Nodosilineaceae","legendgroup":"Cyanobacteriia_Nodosilineaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9],"base":[0.232587188964028,0.250160792682236,0.234825313738147,0.265223619149902,0.231303500493067,0.285511242314276,0.287705657389086,0.272651032300416,0.240977030899582,0.289835518372364,0.291193834982393,0.243793662327399,0.161840009034085,0.276391993426634,0.315983860550167],"x":[7,7,5,7,2,8,8,7,5,8,8,5,9,7,3],"y":[0.0175736037182082,0.015062826467666,0.00615171716143501,0.0074274131505141,0.0029885015309557,0.00219441507480933,0.00212986098327828,0.00374096112621763,0.00281663142781707,0.00135831661002916,0.00127364705297289,0.00151901373184563,0.00224383409287573,0.00173695461437878,0.00159191096019656],"text":["Location: LP<br />Abundance: 0.0175736037<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LP<br />Abundance: 0.0150628265<br />Class_family: Cyanobacteriia_Nostocaceae","Location: CS<br />Abundance: 0.0061517172<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LP<br />Abundance: 0.0074274132<br />Class_family: Cyanobacteriia_Nostocaceae","Location: CMP<br />Abundance: 0.0029885015<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LSJ<br />Abundance: 0.0021944151<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LSJ<br />Abundance: 0.0021298610<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LP<br />Abundance: 0.0037409611<br />Class_family: Cyanobacteriia_Nostocaceae","Location: CS<br />Abundance: 0.0028166314<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LSJ<br />Abundance: 0.0013583166<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LSJ<br />Abundance: 0.0012736471<br />Class_family: Cyanobacteriia_Nostocaceae","Location: CS<br />Abundance: 0.0015190137<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LT<br />Abundance: 0.0022438341<br />Class_family: Cyanobacteriia_Nostocaceae","Location: LP<br />Abundance: 0.0017369546<br />Class_family: Cyanobacteriia_Nostocaceae","Location: CMPE<br />Abundance: 0.0015919110<br />Class_family: Cyanobacteriia_Nostocaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(153,153,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Nostocaceae","legendgroup":"Cyanobacteriia_Nostocaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0.289388527602046,0.221326929649655,0.314223346733225],"x":[3,7,3],"y":[0.0248348191311785,0.0112602593143732,0.00176051381694237],"text":["Location: CMPE<br />Abundance: 0.0248348191<br />Class_family: Cyanobacteriia_Oscillatoriaceae","Location: LP<br />Abundance: 0.0112602593<br />Class_family: Cyanobacteriia_Oscillatoriaceae","Location: CMPE<br />Abundance: 0.0017605138<br />Class_family: Cyanobacteriia_Oscillatoriaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(161,194,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Oscillatoriaceae","legendgroup":"Cyanobacteriia_Oscillatoriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9],"base":[0.26918140466258,0.274328605494587,0.278595392888793,0.225122150809571,0.281937155884971,0.230603233352347,0.280784990868323,0.285979533508442,0.28454133764792,0.236687935011866,0.219946336853075],"x":[8,8,8,5,8,5,3,3,8,4,7],"y":[0.00514720083200754,0.00426678739420644,0.00334176299617733,0.00548108254277618,0.00260418176294958,0.00422208038579922,0.00519454264011882,0.00340899409360396,0.000969904666356047,0.00123480922784489,0.00138059279658032],"text":["Location: LSJ<br />Abundance: 0.0051472008<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: LSJ<br />Abundance: 0.0042667874<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: LSJ<br />Abundance: 0.0033417630<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: CS<br />Abundance: 0.0054810825<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: LSJ<br />Abundance: 0.0026041818<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: CS<br />Abundance: 0.0042220804<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: CMPE<br />Abundance: 0.0051945426<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: CMPE<br />Abundance: 0.0034089941<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: LSJ<br />Abundance: 0.0009699047<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: CMPW<br />Abundance: 0.0012348092<br />Class_family: Cyanobacteriia_Prochlorotrichaceae","Location: LP<br />Abundance: 0.0013805928<br />Class_family: Cyanobacteriia_Prochlorotrichaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Cyanobacteriia_Prochlorotrichaceae","legendgroup":"Cyanobacteriia_Prochlorotrichaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.216276453132359,0.218316122466395],"x":[7,7],"y":[0.0020396693340359,0.00163021438667951],"text":["Location: LP<br />Abundance: 0.0020396693<br />Class_family: Desulfobacteria_Desulfatiglandaceae","Location: LP<br />Abundance: 0.0016302144<br />Class_family: Desulfobacteria_Desulfatiglandaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,72,156,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfobacteria_Desulfatiglandaceae","legendgroup":"Desulfobacteria_Desulfatiglandaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0.203309681612429,0.207842280132509,0.21217180619447],"x":[7,7,7],"y":[0.00453259852007981,0.00432952606196055,0.00410464693788948],"text":["Location: LP<br />Abundance: 0.0045325985<br />Class_family: Desulfobacteria_Desulfosarcinaceae","Location: LP<br />Abundance: 0.0043295261<br />Class_family: Desulfobacteria_Desulfosarcinaceae","Location: LP<br />Abundance: 0.0041046469<br />Class_family: Desulfobacteria_Desulfosarcinaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,204,204,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfobacteria_Desulfosarcinaceae","legendgroup":"Desulfobacteria_Desulfosarcinaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.201270012278393,"x":[7],"y":[0.0020396693340359],"text":"Location: LP<br />Abundance: 0.0020396693<br />Class_family: Desulfobulbia_Desulfurivibrionaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(153,153,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfobulbia_Desulfurivibrionaceae","legendgroup":"Desulfobulbia_Desulfurivibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9],"base":[0.268343439591776,0.279216151081463],"x":[8,3],"y":[0.000837965070803137,0.0015688397868604],"text":["Location: LSJ<br />Abundance: 0.0008379651<br />Class_family: Desulfuromonadia_Bradymonadales","Location: CMPE<br />Abundance: 0.0015688398<br />Class_family: Desulfuromonadia_Bradymonadales"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(161,194,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfuromonadia_Bradymonadales","legendgroup":"Desulfuromonadia_Bradymonadales","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.267695347547605,"x":[8],"y":[0.000648092044171766],"text":"Location: LSJ<br />Abundance: 0.0006480920<br />Class_family: Desulfuromonadia_PB19","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfuromonadia_PB19","legendgroup":"Desulfuromonadia_PB19","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.198625996475014,"x":[7],"y":[0.00264401580337989],"text":"Location: LP<br />Abundance: 0.0026440158<br />Class_family: Desulfuromonadia_Sva1033","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,72,156,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Desulfuromonadia_Sva1033","legendgroup":"Desulfuromonadia_Sva1033","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.211802288734006,0.216845533332563,0.265957789752287,0.283229832765312,0.221073035331391,0.224471495399707,0.223131388613458,0.227653503794538,0.275470852984108,0.230496784652,0.225250745157681,0.233114294956255,0.227067277613353,0.235444644951412,0.228021310094064,0.228902131597924,0.229761761381471,0.230547474558703],"x":[4,2,3,6,4,4,2,4,3,4,2,4,2,4,2,2,2,2],"y":[0.00927074659738461,0.00628585528089454,0.00951306323182083,0.0192128360443428,0.00339846006831615,0.00318200839483107,0.0021193565442236,0.00284328085746197,0.00374529809735502,0.00261751030425519,0.00181653245567112,0.00233034999515669,0.00095403248071102,0.00124329006045373,0.000880821503860596,0.00085962978354634,0.000785713177232,0.000756025934364746],"text":["Location: CMPW<br />Abundance: 0.0092707466<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0062858553<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPE<br />Abundance: 0.0095130632<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: LLC<br />Abundance: 0.0192128360<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0033984601<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0031820084<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0021193565<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0028432809<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPE<br />Abundance: 0.0037452981<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0026175103<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0018165325<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0023303500<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0009540325<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMPW<br />Abundance: 0.0012432901<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0008808215<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0008596298<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0007857132<br />Class_family: Gammaproteobacteria_Aeromonadaceae","Location: CMP<br />Abundance: 0.0007560259<br />Class_family: Gammaproteobacteria_Aeromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,204,204,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aeromonadaceae","legendgroup":"Gammaproteobacteria_Aeromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.255029020999017,0.212224219405113,0.216438074256195,0.258047961555655,0.146302185869767,0.260183021599168,0.220348985743194,0.151135467424058,0.211326017091198,0.261976879945681,0.155252455399553,0.262747441027275,0.263545944894728,0.213349002742922,0.223292429723484,0.159086212647374,0.264832930109836,0.265881532554412,0.192287572943969,0.214882439231476,0.066452725053671,0.266889816748048,0.194269391635834,0.195831861726743,0.0680355469211235,0.210614972168771,0.216016431287368,0.197337235602046,0.0692453666910534],"x":[8,5,5,8,9,8,5,9,2,8,9,3,8,2,5,9,8,8,7,2,1,8,7,7,1,4,2,7,1],"y":[0.00301894055663787,0.00421385485108197,0.00391091148699937,0.00213506004351272,0.00483328155429086,0.00179385834651286,0.00294344398028995,0.00411698797549526,0.00202298565172382,0.00156906494904735,0.0038337572478207,0.00321034872501258,0.00128698521510789,0.00153343648855356,0.00182972108608681,0.00275379638671114,0.0010486024445765,0.00100828419363574,0.00198181869186528,0.00113399205589276,0.00158282186745248,0.00080553079955642,0.00156247009090918,0.00150537387530275,0.00120981976992988,0.00118731656523546,0.000829102045194846,0.00128876087296745,0.00109426404222154],"text":["Location: LSJ<br />Abundance: 0.0030189406<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CS<br />Abundance: 0.0042138549<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CS<br />Abundance: 0.0039109115<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0021350600<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LT<br />Abundance: 0.0048332816<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0017938583<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CS<br />Abundance: 0.0029434440<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LT<br />Abundance: 0.0041169880<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMP<br />Abundance: 0.0020229857<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0015690649<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LT<br />Abundance: 0.0038337572<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMPE<br />Abundance: 0.0032103487<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0012869852<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMP<br />Abundance: 0.0015334365<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CS<br />Abundance: 0.0018297211<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LT<br />Abundance: 0.0027537964<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0010486024<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0010082842<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LP<br />Abundance: 0.0019818187<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMP<br />Abundance: 0.0011339921<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: BSJ<br />Abundance: 0.0015828219<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LSJ<br />Abundance: 0.0008055308<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LP<br />Abundance: 0.0015624701<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LP<br />Abundance: 0.0015053739<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: BSJ<br />Abundance: 0.0012098198<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMPW<br />Abundance: 0.0011873166<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: CMP<br />Abundance: 0.0008291020<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: LP<br />Abundance: 0.0012887609<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1","Location: BSJ<br />Abundance: 0.0010942640<br />Class_family: Gammaproteobacteria_Alcanivoracaceae1"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(153,153,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alcanivoracaceae1","legendgroup":"Gammaproteobacteria_Alcanivoracaceae1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9],"base":[0.057027934754537,0.207023846452963,0.190198107539873,0.2095813611559,0.210456970698199],"x":[1,2,7,2,2],"y":[0.00942479029913399,0.00255751470293719,0.00208946540409577,0.000875609542299316,0.000869046392999062],"text":["Location: BSJ<br />Abundance: 0.0094247903<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Location: CMP<br />Abundance: 0.0025575147<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Location: LP<br />Abundance: 0.0020894654<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Location: CMP<br />Abundance: 0.0008756095<br />Class_family: Gammaproteobacteria_Alteromonadaceae","Location: CMP<br />Abundance: 0.0008690464<br />Class_family: Gammaproteobacteria_Alteromonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(161,194,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Alteromonadaceae","legendgroup":"Gammaproteobacteria_Alteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9],"base":[0.206948823012385,0.209069847735873,0.206005577061665,0.261212104872772],"x":[4,4,2,3],"y":[0.0021210247234881,0.00154512443289742,0.00101826939129757,0.00153533615450274],"text":["Location: CMPW<br />Abundance: 0.0021210247<br />Class_family: Gammaproteobacteria_Aquaspirillaceae","Location: CMPW<br />Abundance: 0.0015451244<br />Class_family: Gammaproteobacteria_Aquaspirillaceae","Location: CMP<br />Abundance: 0.0010182694<br />Class_family: Gammaproteobacteria_Aquaspirillaceae","Location: CMPE<br />Abundance: 0.0015353362<br />Class_family: Gammaproteobacteria_Aquaspirillaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Aquaspirillaceae","legendgroup":"Gammaproteobacteria_Aquaspirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.259147623482523,"x":[3],"y":[0.00206448139024934],"text":"Location: CMPE<br />Abundance: 0.0020644814<br />Class_family: Gammaproteobacteria_Burkholderiaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,72,156,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Burkholderiaceae","legendgroup":"Gammaproteobacteria_Burkholderiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9],"base":[0.148247566361669,0.164355745934631,0.161962527169062,0.177131272176564,0.195352196228185,0.16339388595503,0.187092936000701,0.170066088057609,0.175214264658283,0.180301277614896,0.185252201256203,0.190057311328333,0.194700982937972,0.195674377204369,0.204050446504274,0.199286803239962,0.132931871117716,0.137419539303467,0.173694146175664,0.177093595065724,0.201999806258172,0.208030198308074,0.180451060409719,0.183122359140502,0.201765980507114,0.210468476260512,0.203793083866549,0.141776351821224,0.203939670570266,0.185460459844388,0.257163502575611,0.144261974652229,0.187320360175831,0.205671387667557,0.205136097453779,0.188923452969035],"x":[4,4,7,4,5,2,4,2,2,2,2,2,2,4,5,2,9,9,7,7,2,5,7,7,4,5,4,9,2,7,3,9,7,4,2,7],"y":[0.0161081795729614,0.0127755262419337,0.0117316190066028,0.00996166382413671,0.00869825027608936,0.0066722021025781,0.00858144120366749,0.00514817660067424,0.0050870129566136,0.00495092364130634,0.00480511007213077,0.00464367160963877,0.00458582030198931,0.00609160330274511,0.00397975180379961,0.00271300301821015,0.00448766818575147,0.00435681251775713,0.00339944889005986,0.0033574653439947,0.00193986431209447,0.00243827795243867,0.00267129873078356,0.00233810070388601,0.00202710335943543,0.00175574314460081,0.00187830380100804,0.00248562283100465,0.00119642688351237,0.00185990033144293,0.00198412090691164,0.002040211217538,0.00160309279320342,0.00127743534482805,0.000869479607886636,0.00127465457083836],"text":["Location: CMPW<br />Abundance: 0.0161081796<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0127755262<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0117316190<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0099616638<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CS<br />Abundance: 0.0086982503<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0066722021<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0085814412<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0051481766<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0050870130<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0049509236<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0048051101<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0046436716<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0045858203<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0060916033<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CS<br />Abundance: 0.0039797518<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0027130030<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LT<br />Abundance: 0.0044876682<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LT<br />Abundance: 0.0043568125<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0033994489<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0033574653<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0019398643<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CS<br />Abundance: 0.0024382780<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0026712987<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0023381007<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0020271034<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CS<br />Abundance: 0.0017557431<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0018783038<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LT<br />Abundance: 0.0024856228<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0011964269<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0018599003<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPE<br />Abundance: 0.0019841209<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LT<br />Abundance: 0.0020402112<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0016030928<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMPW<br />Abundance: 0.0012774353<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: CMP<br />Abundance: 0.0008694796<br />Class_family: Gammaproteobacteria_Chromatiaceae","Location: LP<br />Abundance: 0.0012746546<br />Class_family: Gammaproteobacteria_Chromatiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,204,204,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Chromatiaceae","legendgroup":"Gammaproteobacteria_Chromatiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9],"base":[0.217325540799778,0.238308468244649,0.144955635607181,0.149645917864123,0.153471497650122,0.134590653734887,0.157060778300659,0.252595410255047,0.159531059194564,0.139111094226428,0.141027751897969,0.161278183166507,0.277933986932576,0.142939477474796,0.144494433732744,0.160365175238517,0.0556379454161074,0.145919213611027,0.162511855698272,0.194217135802049,0.131295334318343,0.147194628790667],"x":[3,3,2,2,2,4,2,3,2,4,4,2,6,4,4,7,1,4,2,5,9,4],"y":[0.0209829274448708,0.0142869420103981,0.00469028225694207,0.00382557978599818,0.00358928065053721,0.00452044049154096,0.00247028089390519,0.00456809232056404,0.00174712397194329,0.00191665767154089,0.0019117255768267,0.00123367253176437,0.00529584583273546,0.00155495625794891,0.00142477987828257,0.00159735193054425,0.00138998933842966,0.00127541517964014,0.000882030256758842,0.00113506042613523,0.00163653679937303,0.00105293757100214],"text":["Location: CMPE<br />Abundance: 0.0209829274<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPE<br />Abundance: 0.0142869420<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0046902823<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0038255798<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0035892807<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0045204405<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0024702809<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPE<br />Abundance: 0.0045680923<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0017471240<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0019166577<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0019117256<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0012336725<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: LLC<br />Abundance: 0.0052958458<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0015549563<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0014247799<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: LP<br />Abundance: 0.0015973519<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: BSJ<br />Abundance: 0.0013899893<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0012754152<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMP<br />Abundance: 0.0008820303<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CS<br />Abundance: 0.0011350604<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: LT<br />Abundance: 0.0016365368<br />Class_family: Gammaproteobacteria_Comamonadaceae","Location: CMPW<br />Abundance: 0.0010529376<br />Class_family: Gammaproteobacteria_Comamonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(153,153,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Comamonadaceae","legendgroup":"Gammaproteobacteria_Comamonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9],"base":[0.238081114887652,0.182524574896466,0.248147694175725,0.149458663957431,0.12811115294335,0.187427877472172,0.142717558833228,0.154134865365203,0.209989891612153,0.213831370031065,0.190418263212153,0.250757393272918,0.131346679872753,0.251981649934099,0.192693190207438,0.252930844892276,0.157392564451358,0.253804765201933,0.159094696673298,0.254431068216846,0.133535951500943],"x":[8,5,8,7,4,5,2,7,3,3,5,8,4,8,5,8,7,8,7,8,4],"y":[0.0100665792880727,0.00490330257570654,0.00260969909719319,0.00467620140777206,0.0032355269294029,0.00299038573998081,0.0022380767739528,0.0032576990861547,0.00384147841891183,0.00349417076871306,0.00227492699528492,0.0012242566611807,0.00218927162819024,0.000949194958177579,0.00152394559461141,0.000873920309656717,0.00170213222194013,0.000626303014913399,0.00127047856521956,0.000597952782170952,0.00105470223394391],"text":["Location: LSJ<br />Abundance: 0.0100665793<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CS<br />Abundance: 0.0049033026<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0026096991<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LP<br />Abundance: 0.0046762014<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMPW<br />Abundance: 0.0032355269<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CS<br />Abundance: 0.0029903857<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMP<br />Abundance: 0.0022380768<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LP<br />Abundance: 0.0032576991<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMPE<br />Abundance: 0.0038414784<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMPE<br />Abundance: 0.0034941708<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CS<br />Abundance: 0.0022749270<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0012242567<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMPW<br />Abundance: 0.0021892716<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0009491950<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CS<br />Abundance: 0.0015239456<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0008739203<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LP<br />Abundance: 0.0017021322<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0006263030<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LP<br />Abundance: 0.0012704786<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: LSJ<br />Abundance: 0.0005979528<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae","Location: CMPW<br />Abundance: 0.0010547022<br />Class_family: Gammaproteobacteria_Ectothiorhodospiraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(161,194,153,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Ectothiorhodospiraceae","legendgroup":"Gammaproteobacteria_Ectothiorhodospiraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.899999999999999],"base":[0.0539013582184834,0.147975107094814,0.129572664003213],"x":[1,7,9],"y":[0.00173658719762394,0.0014835568626172,0.00172267031512952],"text":["Location: BSJ<br />Abundance: 0.0017365872<br />Class_family: Gammaproteobacteria_Halieaceae","Location: LP<br />Abundance: 0.0014835569<br />Class_family: Gammaproteobacteria_Halieaceae","Location: LT<br />Abundance: 0.0017226703<br />Class_family: Gammaproteobacteria_Halieaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(48,0,24,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Halieaceae","legendgroup":"Gammaproteobacteria_Halieaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.180590768244532,"x":[5],"y":[0.0019338066519341],"text":"Location: CS<br />Abundance: 0.0019338067<br />Class_family: Gammaproteobacteria_Hydrogenophilaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(25,25,112,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Hydrogenophilaceae","legendgroup":"Gammaproteobacteria_Hydrogenophilaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.146464240921454,"x":[7],"y":[0.00151086617335996],"text":"Location: LP<br />Abundance: 0.0015108662<br />Class_family: Gammaproteobacteria_KI89A_clade","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(53,57,53,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_KI89A_clade","legendgroup":"Gammaproteobacteria_KI89A_clade","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9],"base":[0.23654314441782,0.176143432353914,0.178592056517853],"x":[8,5,5],"y":[0.00153797046983159,0.00244862416393898,0.00199871172667856],"text":["Location: LSJ<br />Abundance: 0.0015379705<br />Class_family: Gammaproteobacteria_Legionellaceae","Location: CS<br />Abundance: 0.0024486242<br />Class_family: Gammaproteobacteria_Legionellaceae","Location: CS<br />Abundance: 0.0019987117<br />Class_family: Gammaproteobacteria_Legionellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(240,255,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Legionellaceae","legendgroup":"Gammaproteobacteria_Legionellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999],"base":[0.127858922585587,0.133975519122359,0.138049963458745,0.118166621255398,0.12208463690427,0.139480118280575,0.140150383967951,0.142419444879746,0.0521011173748286,0.144651201513422,0.141731012281152,0.125881858823416,0.125690178808765,0.127776796170058,0.235269606967745,0.126978332920578,0.235930954664645],"x":[2,2,2,9,9,7,2,7,1,7,2,9,4,9,8,4,8],"y":[0.00611659653677196,0.00407444433638607,0.00210042050920609,0.00391801564887173,0.00379722191914612,0.00293932659917098,0.00158062831320094,0.00223175663367536,0.00180024084365481,0.0018130394080319,0.000986546552076389,0.00189493734664245,0.00128815411181288,0.00179586783315494,0.000661347696900849,0.00113282002277207,0.00061218975317498],"text":["Location: CMP<br />Abundance: 0.0061165965<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMP<br />Abundance: 0.0040744443<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMP<br />Abundance: 0.0021004205<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LT<br />Abundance: 0.0039180156<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LT<br />Abundance: 0.0037972219<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LP<br />Abundance: 0.0029393266<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMP<br />Abundance: 0.0015806283<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LP<br />Abundance: 0.0022317566<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: BSJ<br />Abundance: 0.0018002408<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LP<br />Abundance: 0.0018130394<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMP<br />Abundance: 0.0009865466<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LT<br />Abundance: 0.0018949373<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMPW<br />Abundance: 0.0012881541<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LT<br />Abundance: 0.0017958678<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LSJ<br />Abundance: 0.0006613477<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: CMPW<br />Abundance: 0.0011328200<br />Class_family: Gammaproteobacteria_Litoricolaceae","Location: LSJ<br />Abundance: 0.0006121898<br />Class_family: Gammaproteobacteria_Litoricolaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(137,207,240,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Litoricolaceae","legendgroup":"Gammaproteobacteria_Litoricolaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.202587205140359,0.124618801440688],"x":[3,4],"y":[0.00740268647179362,0.00107137736807741],"text":["Location: CMPE<br />Abundance: 0.0074026865<br />Class_family: Gammaproteobacteria_Methylomonadaceae","Location: CMPW<br />Abundance: 0.0010713774<br />Class_family: Gammaproteobacteria_Methylomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylomonadaceae","legendgroup":"Gammaproteobacteria_Methylomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.126971450541845,"x":[2],"y":[0.000887472043742177],"text":"Location: CMP<br />Abundance: 0.0008874720<br />Class_family: Gammaproteobacteria_Methylophagaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(115,147,179,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophagaceae","legendgroup":"Gammaproteobacteria_Methylophagaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9],"base":[0.118036230742765,0.112720178712564,0.120194671315378,0.0503361753712455,0.121577981154702,0.122790279567396,0.137587370091715,0.123998399690394,0.122048671949799,0.116185550392999,0.125176303003249,0.234602009524277,0.126113079049421,0.123447373267809],"x":[2,9,2,1,2,2,7,2,4,9,2,8,2,4],"y":[0.00215844057261218,0.00346537168043493,0.00138330983932454,0.00176494200358315,0.00121229841269355,0.00120812012299819,0.00189274818886012,0.00117790331285504,0.00139870131801041,0.00198107086239895,0.000936776046172288,0.000667597443468004,0.000858371492423626,0.00117142817287841],"text":["Location: CMP<br />Abundance: 0.0021584406<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: LT<br />Abundance: 0.0034653717<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0013833098<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: BSJ<br />Abundance: 0.0017649420<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0012122984<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0012081201<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: LP<br />Abundance: 0.0018927482<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0011779033<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMPW<br />Abundance: 0.0013987013<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: LT<br />Abundance: 0.0019810709<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0009367760<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: LSJ<br />Abundance: 0.0006675974<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMP<br />Abundance: 0.0008583715<br />Class_family: Gammaproteobacteria_Methylophilaceae","Location: CMPW<br />Abundance: 0.0011714282<br />Class_family: Gammaproteobacteria_Methylophilaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(8,143,143,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Methylophilaceae","legendgroup":"Gammaproteobacteria_Methylophilaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.136287080045197,"x":[7],"y":[0.00130029004651819],"text":"Location: LP<br />Abundance: 0.0013002900<br />Class_family: Gammaproteobacteria_Milano-WF1B-44","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,150,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Milano-WF1B-44","legendgroup":"Gammaproteobacteria_Milano-WF1B-44","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9],"base":[0.136873107317603,0.0496708196843308,0.0678130168011289,0.174554069802295,0.0829899844424684,0.079049242308303,0.0854356479732979,0.0911732043511861,0.0966513922724164,0.0932743734193047,0.0996245108248479,0.195025321110033,0.101818051649972,0.105499030466696,0.108896069848909,0.10557403591645,0.110031276333719,0.114124373452986,0.112236166659551,0.114718921777576,0.267958091294167,0.233038606509605,0.11775307477267,0.174314569475352,0.120265082977775,0.11693650685333],"x":[3,4,4,3,4,2,2,2,2,4,4,3,2,2,2,4,4,4,2,2,6,8,4,5,4,2],"y":[0.0376809624846914,0.0181421971167981,0.0151769676413394,0.0204712513077379,0.0102843889768364,0.00638640566499486,0.00573755637788818,0.00547818792123028,0.00516665937755603,0.0063501374055432,0.00594952509160239,0.00756188403032682,0.00368097881672402,0.00339703938221242,0.00334009681064208,0.00445724041726871,0.00409309711926714,0.00362870131968336,0.00248275511802472,0.00221758507575436,0.00997589563840878,0.00156340301467109,0.00251200820510514,0.00182886287856238,0.00178358897202407,0.00109972388943545],"text":["Location: CMPE<br />Abundance: 0.0376809625<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0181421971<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0151769676<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPE<br />Abundance: 0.0204712513<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0102843890<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0063864057<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0057375564<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0054781879<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0051666594<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0063501374<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0059495251<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPE<br />Abundance: 0.0075618840<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0036809788<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0033970394<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0033400968<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0044572404<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0040930971<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0036287013<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0024827551<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0022175851<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: LLC<br />Abundance: 0.0099758956<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: LSJ<br />Abundance: 0.0015634030<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0025120082<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CS<br />Abundance: 0.0018288629<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMPW<br />Abundance: 0.0017835890<br />Class_family: Gammaproteobacteria_Moraxellaceae","Location: CMP<br />Abundance: 0.0010997239<br />Class_family: Gammaproteobacteria_Moraxellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(95,158,160,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Moraxellaceae","legendgroup":"Gammaproteobacteria_Moraxellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999],"base":[0.14990550357695,0.154625538193029,0.100327368617069,0.158595687719444,0.217626319242909,0.219545116535937,0.1618834461378,0.0466909777459608,0.221299628825381,0.16470480167607,0.222944655290618,0.224539302820356,0.259217673143363,0.22607183058165,0.167444662351961,0.22737593006856,0.129080812608882,0.228568524276086,0.229721718927416,0.169495506180509,0.171253927259644,0.131955561359458,0.105630976472957,0.134496776463977,0.108275694717913,0.230778878912756,0.173003592410571,0.110713910897576,0.231672046421873,0.232431538011913],"x":[5,5,9,5,8,8,5,4,8,5,8,8,6,8,5,8,3,8,8,5,5,3,9,3,9,8,5,9,8,8],"y":[0.00472003461607934,0.00397014952641478,0.00530360785588811,0.00328775841835574,0.00191879729302732,0.00175451228944487,0.00282135553827001,0.00297984193837007,0.00164502646523623,0.00273986067589099,0.00159464752973834,0.00153252776129445,0.00874041815080451,0.00130409948690988,0.00205084382854845,0.00119259420752613,0.00287474875057639,0.00115319465132965,0.00105715998534003,0.00175842107913426,0.00174966515092698,0.00254121510451857,0.0026447182449567,0.00237633085362676,0.00243821617966224,0.000893167509116694,0.00131097706478106,0.00200626781498865,0.000759491590040295,0.000607068497692376],"text":["Location: CS<br />Abundance: 0.0047200346<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0039701495<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LT<br />Abundance: 0.0053036079<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0032877584<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0019187973<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0017545123<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0028213555<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CMPW<br />Abundance: 0.0029798419<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0016450265<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0027398607<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0015946475<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0015325278<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LLC<br />Abundance: 0.0087404182<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0013040995<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0020508438<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0011925942<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CMPE<br />Abundance: 0.0028747488<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0011531947<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0010571600<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0017584211<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0017496652<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CMPE<br />Abundance: 0.0025412151<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LT<br />Abundance: 0.0026447182<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CMPE<br />Abundance: 0.0023763309<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LT<br />Abundance: 0.0024382162<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0008931675<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: CS<br />Abundance: 0.0013109771<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LT<br />Abundance: 0.0020062678<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0007594916<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group","Location: LSJ<br />Abundance: 0.0006070685<br />Class_family: Gammaproteobacteria_MWH-UniP1_aquatic_group"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,71,171,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_MWH-UniP1_aquatic_group","legendgroup":"Gammaproteobacteria_MWH-UniP1_aquatic_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.0307355964181198,0.0137693456995667,0.118627467035291,0.0449642514804836,0.0250742669314585,0.0364868109733021,0.0557201886869856,0.0611193863931864,0.0647714314126695,0.0323815913592661,0.0458337361840427,0.0365145023659848,0.0679442889808196,0.0702145932809548,0.0724650131179122,0.0745339757162661,0.0397536471714173,0.076503070547901,0.0424564516506645,0.0778616324352238,0.0443013336758539,0.148801494888062,0.0491351343388752,0.0455498167957564],"x":[2,4,7,2,1,1,2,2,2,4,1,4,2,2,2,2,4,2,4,2,4,5,1,4],"y":[0.0142286550623638,0.0186122456596994,0.0176596130099061,0.010755937206502,0.0114125440418435,0.00934692521074062,0.0053991977062008,0.00365204501948305,0.0031728575681501,0.00413291100671876,0.0033013981548325,0.00323914480543244,0.00227030430013517,0.00225041983695745,0.00206896259835392,0.00196909483163485,0.00270280447924722,0.00135856188732281,0.0018448820251894,0.00118760987307927,0.00124848311990255,0.0011040086888883,0.0012010410323703,0.00114116095020433],"text":["Location: CMP<br />Abundance: 0.0142286551<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0186122457<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: LP<br />Abundance: 0.0176596130<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0107559372<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: BSJ<br />Abundance: 0.0114125440<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: BSJ<br />Abundance: 0.0093469252<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0053991977<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0036520450<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0031728576<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0041329110<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: BSJ<br />Abundance: 0.0033013982<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0032391448<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0022703043<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0022504198<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0020689626<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0019690948<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0027028045<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0013585619<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0018448820<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMP<br />Abundance: 0.0011876099<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0012484831<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CS<br />Abundance: 0.0011040087<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: BSJ<br />Abundance: 0.0012010410<br />Class_family: Gammaproteobacteria_Nitrincolaceae","Location: CMPW<br />Abundance: 0.0011411610<br />Class_family: Gammaproteobacteria_Nitrincolaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(100,149,237,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Nitrincolaceae","legendgroup":"Gammaproteobacteria_Nitrincolaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.0985034887852607,"x":[9],"y":[0.00182387983180785],"text":"Location: LT<br />Abundance: 0.0018238798<br />Class_family: Gammaproteobacteria_Nitrosomonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,255,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Nitrosomonadaceae","legendgroup":"Gammaproteobacteria_Nitrosomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.115779662835886,"x":[7],"y":[0.0028478041994047],"text":"Location: LP<br />Abundance: 0.0028478042<br />Class_family: Gammaproteobacteria_Oceanospirillaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,139,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Oceanospirillaceae","legendgroup":"Gammaproteobacteria_Oceanospirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.028210251711446,"x":[2],"y":[0.00252534470667386],"text":"Location: CMP<br />Abundance: 0.0025253447<br />Class_family: Gammaproteobacteria_Pseudoalteromonadaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(111,143,175,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudoalteromonadaceae","legendgroup":"Gammaproteobacteria_Pseudoalteromonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999],"base":[0.211388901346711,0.145087261509533,0.212721940034058,0.113472376713989,0.213863228809719,0.14703681531466,0.214920388795059,0.0960253018485548,0.215916976765344,0.216807676660798],"x":[8,5,8,7,8,5,8,9,8,8],"y":[0.00133303868734658,0.0019495538051271,0.00114128877566172,0.00230728612189728,0.00105715998534001,0.00176467957340215,0.000996587970284901,0.0024781869367059,0.000890699895453179,0.000818642582111673],"text":["Location: LSJ<br />Abundance: 0.0013330387<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: CS<br />Abundance: 0.0019495538<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LSJ<br />Abundance: 0.0011412888<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LP<br />Abundance: 0.0023072861<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LSJ<br />Abundance: 0.0010571600<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: CS<br />Abundance: 0.0017646796<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LSJ<br />Abundance: 0.0009965880<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LT<br />Abundance: 0.0024781869<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LSJ<br />Abundance: 0.0008906999<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae","Location: LSJ<br />Abundance: 0.0008186426<br />Class_family: Gammaproteobacteria_Pseudohongiellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(20,52,164,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudohongiellaceae","legendgroup":"Gammaproteobacteria_Pseudohongiellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9],"base":[0.0253808122337961,0.124211426370748,0.126732651431407,0.254722362145808,0.0125285639579975,0.0273867617847597],"x":[2,3,3,6,4,2],"y":[0.00200594955096358,0.00252122506065848,0.0023481611774748,0.00449531099755457,0.00124078174156915,0.000823489926686256],"text":["Location: CMP<br />Abundance: 0.0020059496<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Location: CMPE<br />Abundance: 0.0025212251<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Location: CMPE<br />Abundance: 0.0023481612<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Location: LLC<br />Abundance: 0.0044953110<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Location: CMPW<br />Abundance: 0.0012407817<br />Class_family: Gammaproteobacteria_Pseudomonadaceae","Location: CMP<br />Abundance: 0.0008234899<br />Class_family: Gammaproteobacteria_Pseudomonadaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(125,249,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Pseudomonadaceae","legendgroup":"Gammaproteobacteria_Pseudomonadaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999],"base":[0.101283580885704,0.113460367138812,0.0210872852662746,0.00782847564947748,0.120866106237002,0.0237505976241017,0.0106425891734077,0.210772016422248],"x":[3,3,2,4,3,2,4,8],"y":[0.0121767862531079,0.00740573909818966,0.00266331235782707,0.0028141135239302,0.0033453201337464,0.00163021460969443,0.00188597478458986,0.000616884924463573],"text":["Location: CMPE<br />Abundance: 0.0121767863<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMPE<br />Abundance: 0.0074057391<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMP<br />Abundance: 0.0026633124<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMPW<br />Abundance: 0.0028141135<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMPE<br />Abundance: 0.0033453201<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMP<br />Abundance: 0.0016302146<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: CMPW<br />Abundance: 0.0018859748<br />Class_family: Gammaproteobacteria_Rhodocyclaceae","Location: LSJ<br />Abundance: 0.0006168849<br />Class_family: Gammaproteobacteria_Rhodocyclaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(96,130,182,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Rhodocyclaceae","legendgroup":"Gammaproteobacteria_Rhodocyclaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.111659337305957,"x":[7],"y":[0.00181303940803192],"text":"Location: LP<br />Abundance: 0.0018130394<br />Class_family: Gammaproteobacteria_Run-SP154","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,163,108,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Run-SP154","legendgroup":"Gammaproteobacteria_Run-SP154","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0201284674307807,"x":[2],"y":[0.000958817835493984],"text":"Location: CMP<br />Abundance: 0.0009588178<br />Class_family: Gammaproteobacteria_Saccharospirillaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(63,0,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Saccharospirillaceae","legendgroup":"Gammaproteobacteria_Saccharospirillaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.9],"base":[0.209416795731877,0.020178775138848,0.0218535027164645,0.0235042017938807],"x":[8,1,1,1],"y":[0.00135522069037058,0.00167472757761652,0.00165069907741625,0.00157006513757781],"text":["Location: LSJ<br />Abundance: 0.0013552207<br />Class_family: Gammaproteobacteria_SAR86_clade","Location: BSJ<br />Abundance: 0.0016747276<br />Class_family: Gammaproteobacteria_SAR86_clade","Location: BSJ<br />Abundance: 0.0016506991<br />Class_family: Gammaproteobacteria_SAR86_clade","Location: BSJ<br />Abundance: 0.0015700651<br />Class_family: Gammaproteobacteria_SAR86_clade"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(93,63,211,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_SAR86_clade","legendgroup":"Gammaproteobacteria_SAR86_clade","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.104851516913471,0.109459658742219],"x":[7,7],"y":[0.00460814182874779,0.00219967856373804],"text":["Location: LP<br />Abundance: 0.0046081418<br />Class_family: Gammaproteobacteria_Sedimenticolaceae","Location: LP<br />Abundance: 0.0021996786<br />Class_family: Gammaproteobacteria_Sedimenticolaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(173,216,230,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Sedimenticolaceae","legendgroup":"Gammaproteobacteria_Sedimenticolaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9],"base":[0.0776151979195749,0.0164737086930737,0.0830919939878902,0.0879252755421811,0.0915970040577959,0.0189203473077824,0.00636510454862721,0.103311290504315,0.142935405216331,0.0944748985649382,0.144105920452743],"x":[9,2,9,9,9,2,4,7,5,9,5],"y":[0.00547679606831533,0.00244663861470878,0.00483328155429084,0.00367172851561483,0.00287789450714233,0.00120812012299821,0.00146337110085027,0.00154022640915621,0.00117051523641165,0.00155040328361655,0.00098134105678957],"text":["Location: LT<br />Abundance: 0.0054767961<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: CMP<br />Abundance: 0.0024466386<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: LT<br />Abundance: 0.0048332816<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: LT<br />Abundance: 0.0036717285<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: LT<br />Abundance: 0.0028778945<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: CMP<br />Abundance: 0.0012081201<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: CMPW<br />Abundance: 0.0014633711<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: LP<br />Abundance: 0.0015402264<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: CS<br />Abundance: 0.0011705152<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: LT<br />Abundance: 0.0015504033<br />Class_family: Gammaproteobacteria_Thioglobaceae","Location: CS<br />Abundance: 0.0009813411<br />Class_family: Gammaproteobacteria_Thioglobaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(25,25,112,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Thioglobaceae","legendgroup":"Gammaproteobacteria_Thioglobaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9,0.9],"base":[0.0145128110684514,0.101383173780072,0.141788608860374,0.0997326922557409,0.0157239379687204],"x":[2,7,5,3,2],"y":[0.00121112690026896,0.00192811672424284,0.00114679635595713,0.00155088862996357,0.000749770724353272],"text":["Location: CMP<br />Abundance: 0.0012111269<br />Class_family: Gammaproteobacteria_Thiomicrospiraceae","Location: LP<br />Abundance: 0.0019281167<br />Class_family: Gammaproteobacteria_Thiomicrospiraceae","Location: CS<br />Abundance: 0.0011467964<br />Class_family: Gammaproteobacteria_Thiomicrospiraceae","Location: CMPE<br />Abundance: 0.0015508886<br />Class_family: Gammaproteobacteria_Thiomicrospiraceae","Location: CMP<br />Abundance: 0.0007497707<br />Class_family: Gammaproteobacteria_Thiomicrospiraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,0,128,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Thiomicrospiraceae","legendgroup":"Gammaproteobacteria_Thiomicrospiraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9,0.9],"base":[0.0833063349932585,0.0963753273928219,0.0990438330853032,0.097344728877743],"x":[7,7,7,3],"y":[0.0130689923995634,0.00266850569248135,0.002339340694769,0.00238796337799793],"text":["Location: LP<br />Abundance: 0.0130689924<br />Class_family: Gammaproteobacteria_Thiotrichaceae","Location: LP<br />Abundance: 0.0026685057<br />Class_family: Gammaproteobacteria_Thiotrichaceae","Location: LP<br />Abundance: 0.0023393407<br />Class_family: Gammaproteobacteria_Thiotrichaceae","Location: CMPE<br />Abundance: 0.0023879634<br />Class_family: Gammaproteobacteria_Thiotrichaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(31,81,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Thiotrichaceae","legendgroup":"Gammaproteobacteria_Thiotrichaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0811060115516067,"x":[7],"y":[0.00220032344165176],"text":"Location: LP<br />Abundance: 0.0022003234<br />Class_family: Gammaproteobacteria_UBA10353_marine_group","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(167,199,231,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_UBA10353_marine_group","legendgroup":"Gammaproteobacteria_UBA10353_marine_group","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999],"base":[0.203925619129087,0.205779725696619,0.0748863633618756,0.0772328743372989,0.248295584093746,0.140023933554513,0.207013495545546,0.0794806441065335,0.207934468450422,0.208717501855646],"x":[8,8,7,7,6,5,8,7,8,8],"y":[0.00185410656753243,0.00123376984892717,0.00234651097542334,0.00224776976923453,0.00642677805206218,0.00176467530586077,0.000920972904875611,0.00162536744507326,0.000783033405223871,0.000699293876231216],"text":["Location: LSJ<br />Abundance: 0.0018541066<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LSJ<br />Abundance: 0.0012337698<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LP<br />Abundance: 0.0023465110<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LP<br />Abundance: 0.0022477698<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LLC<br />Abundance: 0.0064267781<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: CS<br />Abundance: 0.0017646753<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LSJ<br />Abundance: 0.0009209729<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LP<br />Abundance: 0.0016253674<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LSJ<br />Abundance: 0.0007830334<br />Class_family: Gammaproteobacteria_Unknown_Family","Location: LSJ<br />Abundance: 0.0006992939<br />Class_family: Gammaproteobacteria_Unknown_Family"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,204,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Unknown_Family","legendgroup":"Gammaproteobacteria_Unknown_Family","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9],"base":[0.157483429013704,0.0393618838424186,0.180402585495789,0.19020452636507,0.129047262935411,0.00809920834946079,0.016667223125973,0.198904220843589,0.200747743999356,0.01226627299606,0.0939071240506518,0.202562837570933,0.138258284002748,0.0731841347894269],"x":[6,9,8,8,5,1,1,8,8,2,3,8,5,7],"y":[0.0908121550800426,0.0382533140771563,0.00980194086928074,0.00869969447851984,0.00921102106733701,0.00856801477651224,0.00351155201287493,0.00184352315576691,0.00181509357157661,0.00224653807239139,0.00343760482709114,0.0013627815581537,0.00176564955176589,0.0017022285724487],"text":["Location: LLC<br />Abundance: 0.0908121551<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LT<br />Abundance: 0.0382533141<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LSJ<br />Abundance: 0.0098019409<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LSJ<br />Abundance: 0.0086996945<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: CS<br />Abundance: 0.0092110211<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: BSJ<br />Abundance: 0.0085680148<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: BSJ<br />Abundance: 0.0035115520<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LSJ<br />Abundance: 0.0018435232<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LSJ<br />Abundance: 0.0018150936<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: CMP<br />Abundance: 0.0022465381<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: CMPE<br />Abundance: 0.0034376048<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LSJ<br />Abundance: 0.0013627816<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: CS<br />Abundance: 0.0017656496<br />Class_family: Gammaproteobacteria_Vibrionaceae","Location: LP<br />Abundance: 0.0017022286<br />Class_family: Gammaproteobacteria_Vibrionaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(182,208,226,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Vibrionaceae","legendgroup":"Gammaproteobacteria_Vibrionaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0.0664173431225227,0.0690955524720676,0.071522181998731],"x":[7,7,7],"y":[0.00267820934954491,0.00242662952666337,0.00166195279069593],"text":["Location: LP<br />Abundance: 0.0026782093<br />Class_family: Gammaproteobacteria_Woeseiaceae","Location: LP<br />Abundance: 0.0024266295<br />Class_family: Gammaproteobacteria_Woeseiaceae","Location: LP<br />Abundance: 0.0016619528<br />Class_family: Gammaproteobacteria_Woeseiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(150,222,209,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gammaproteobacteria_Woeseiaceae","legendgroup":"Gammaproteobacteria_Woeseiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0647160406989248,"x":[7],"y":[0.00170130242359788],"text":"Location: LP<br />Abundance: 0.0017013024<br />Class_family: Gracilibacteria_Absconditabacteriales_(SR1)","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(65,105,225,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gracilibacteria_Absconditabacteriales_(SR1)","legendgroup":"Gracilibacteria_Absconditabacteriales_(SR1)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0113927110100884,"x":[2],"y":[0.000873561985971655],"text":"Location: CMP<br />Abundance: 0.0008735620<br />Class_family: Gracilibacteria_Gracilibacteria","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(15,82,186,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gracilibacteria_Gracilibacteria","legendgroup":"Gracilibacteria_Gracilibacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0105813735720441,"x":[2],"y":[0.000811337438044256],"text":"Location: CMP<br />Abundance: 0.0008113374<br />Class_family: Gracilibacteria_JGI_0000069-P22","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(159,226,191,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Gracilibacteria_JGI_0000069-P22","legendgroup":"Gracilibacteria_JGI_0000069-P22","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.060107898870177,"x":[7],"y":[0.0046081418287478],"text":"Location: LP<br />Abundance: 0.0046081418<br />Class_family: Ignavibacteria_Ignavibacteriaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(135,206,235,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Ignavibacteria_Ignavibacteriaceae","legendgroup":"Ignavibacteria_Ignavibacteriaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9],"base":[0.164792905512527,0.167418529401043,0.117816829530345,0.169070018037515,0.170613646001613,0.090498129957048,0.172053044158694,0.173379271699463,0.174594458243528,0.120580373493896,0.17572992637593,0.176798366931869,0.12246377492865,0.177852177772384,0.152234171494626,0.178784114640413,0.124145742917388,0.179629772351205,0.12548117610807,0.12675245884073,0.127976770922699],"x":[8,8,5,8,8,3,8,8,8,5,8,8,5,8,6,8,5,8,5,5,5],"y":[0.00262562388851648,0.0016514886364721,0.00276354396355023,0.0015436279640974,0.00143939815708169,0.00340899409360386,0.00132622754076892,0.00121518654406449,0.00113546813240223,0.00188340143475478,0.0010684405559386,0.00105381084051562,0.00168196798873711,0.000931936868028893,0.00524925751907759,0.000845657710791259,0.0013354331906827,0.000772813144584383,0.00127128273265922,0.00122431208196949,0.00107049201271153],"text":["Location: LSJ<br />Abundance: 0.0026256239<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0016514886<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0027635440<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0015436280<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0014393982<br />Class_family: Kapabacteria_Kapabacteriales","Location: CMPE<br />Abundance: 0.0034089941<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0013262275<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0012151865<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0011354681<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0018834014<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0010684406<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0010538108<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0016819680<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0009319369<br />Class_family: Kapabacteria_Kapabacteriales","Location: LLC<br />Abundance: 0.0052492575<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0008456577<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0013354332<br />Class_family: Kapabacteria_Kapabacteriales","Location: LSJ<br />Abundance: 0.0007728131<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0012712827<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0012243121<br />Class_family: Kapabacteria_Kapabacteriales","Location: CS<br />Abundance: 0.0010704920<br />Class_family: Kapabacteria_Kapabacteriales"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(70,130,180,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Kapabacteria_Kapabacteriales","legendgroup":"Kapabacteria_Kapabacteriales","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9],"base":[0.0150335449372424,0.157965711746064,0.0269448305362482,0.106349697361745,0.0350757662376701,0.161464532348174,0.10951826212143,0.111615542865256,0.113598151697081,0.162903857066575,0.163980238423715,0.115424155465883,0.116727372992186],"x":[9,8,9,5,9,8,5,5,5,8,8,5,5],"y":[0.0119112855990057,0.00349882060211024,0.00813093570142199,0.00316856475968422,0.00428611760474847,0.00143932471840019,0.00209728074382638,0.00198260883182454,0.00182600376880217,0.00107638135714067,0.000812667088811586,0.00130321752630341,0.00108945653815928],"text":["Location: LT<br />Abundance: 0.0119112856<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LSJ<br />Abundance: 0.0034988206<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LT<br />Abundance: 0.0081309357<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0031685648<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LT<br />Abundance: 0.0042861176<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LSJ<br />Abundance: 0.0014393247<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0020972807<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0019826088<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0018260038<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LSJ<br />Abundance: 0.0010763814<br />Class_family: Kiritimatiellae_WCHB1-41","Location: LSJ<br />Abundance: 0.0008126671<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0013032175<br />Class_family: Kiritimatiellae_WCHB1-41","Location: CS<br />Abundance: 0.0010894565<br />Class_family: Kiritimatiellae_WCHB1-41"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,128,128,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Kiritimatiellae_WCHB1-41","legendgroup":"Kiritimatiellae_WCHB1-41","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.00570992748100833,0.00697565693662945],"x":[1,1],"y":[0.00126572945562113,0.00112355141283134],"text":["Location: BSJ<br />Abundance: 0.0012657295<br />Class_family: Marinimicrobia_(SAR406_clade)_Marinimicrobia_(SAR406_clade)","Location: BSJ<br />Abundance: 0.0011235514<br />Class_family: Marinimicrobia_(SAR406_clade)_Marinimicrobia_(SAR406_clade)"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(64,224,208,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Marinimicrobia_(SAR406_clade)_Marinimicrobia_(SAR406_clade)","legendgroup":"Marinimicrobia_(SAR406_clade)_Marinimicrobia_(SAR406_clade)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.899999999999999,"base":0.156513238956672,"x":[8],"y":[0.00145247278939209],"text":"Location: LSJ<br />Abundance: 0.0014524728<br />Class_family: Oligoflexia_0319-6G20","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(4,55,242,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Oligoflexia_0319-6G20","legendgroup":"Oligoflexia_0319-6G20","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9],"base":[0.153834868416025,0.104782409008604,0.154859930328318,0.155760474019533,0.00966131167891199],"x":[8,5,8,8,2],"y":[0.00102506191229304,0.00156728835314142,0.000900543691215588,0.00075276493713855,0.000920061893132132],"text":["Location: LSJ<br />Abundance: 0.0010250619<br />Class_family: Oligoflexia_Oligoflexaceae","Location: CS<br />Abundance: 0.0015672884<br />Class_family: Oligoflexia_Oligoflexaceae","Location: LSJ<br />Abundance: 0.0009005437<br />Class_family: Oligoflexia_Oligoflexaceae","Location: LSJ<br />Abundance: 0.0007527649<br />Class_family: Oligoflexia_Oligoflexaceae","Location: CMP<br />Abundance: 0.0009200619<br />Class_family: Oligoflexia_Oligoflexaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(64,181,173,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Oligoflexia_Oligoflexaceae","legendgroup":"Oligoflexia_Oligoflexaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9],"base":[0.134220067828135,0.14097890821128,0.145708359298834,0.0996382717809685,0.148727338369244,0.1505700683194,0.152377904021563,0.102874647307629,0.146835677930894],"x":[8,8,8,5,8,8,8,5,6],"y":[0.0067588403831447,0.00472945108755415,0.00301897907041015,0.00323637552666027,0.00184272995015561,0.00180783570216325,0.00145696439446172,0.0019077617009753,0.00539849356373223],"text":["Location: LSJ<br />Abundance: 0.0067588404<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LSJ<br />Abundance: 0.0047294511<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LSJ<br />Abundance: 0.0030189791<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: CS<br />Abundance: 0.0032363755<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LSJ<br />Abundance: 0.0018427300<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LSJ<br />Abundance: 0.0018078357<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LSJ<br />Abundance: 0.0014569644<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: CS<br />Abundance: 0.0019077617<br />Class_family: Oligoflexia_Silvanigrellaceae","Location: LLC<br />Abundance: 0.0053984936<br />Class_family: Oligoflexia_Silvanigrellaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(8,24,168,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Oligoflexia_Silvanigrellaceae","legendgroup":"Oligoflexia_Silvanigrellaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9],"base":[0.132382388571993,0.00442881704648001],"x":[8,1],"y":[0.0018376792561425,0.00128111043452832],"text":["Location: LSJ<br />Abundance: 0.0018376793<br />Class_family: OM190_OM190","Location: BSJ<br />Abundance: 0.0012811104<br />Class_family: OM190_OM190"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(170,255,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"OM190_OM190","legendgroup":"OM190_OM190","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999],"base":[0.141951326611327,0.0134746854125403],"x":[6,9],"y":[0.00488435131956724,0.00155885952470209],"text":["Location: LLC<br />Abundance: 0.0048843513<br />Class_family: Parcubacteria_Candidatus_Campbellbacteria","Location: LT<br />Abundance: 0.0015588595<br />Class_family: Parcubacteria_Candidatus_Campbellbacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(95,158,160,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Parcubacteria_Candidatus_Campbellbacteria","legendgroup":"Parcubacteria_Candidatus_Campbellbacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0569033560961762,"x":[7],"y":[0.00320454277400081],"text":"Location: LP<br />Abundance: 0.0032045428<br />Class_family: Parcubacteria_Candidatus_Kaiserbacteria","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(9,121,105,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Parcubacteria_Candidatus_Kaiserbacteria","legendgroup":"Parcubacteria_Candidatus_Kaiserbacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.0975790328962366,0.0986300147722877],"x":[5,5],"y":[0.00105098187605113,0.00100825700868076],"text":["Location: CS<br />Abundance: 0.0010509819<br />Class_family: Parcubacteria_Candidatus_Moranbacteria","Location: CS<br />Abundance: 0.0010082570<br />Class_family: Parcubacteria_Candidatus_Moranbacteria"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(175,225,175,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Parcubacteria_Candidatus_Moranbacteria","legendgroup":"Parcubacteria_Candidatus_Moranbacteria","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9],"base":[0.128117489201738,0.0829211981564429,0.0555779725755292],"x":[8,3,7],"y":[0.00426489937025432,0.00757693180060505,0.00132538352064699],"text":["Location: LSJ<br />Abundance: 0.0042648994<br />Class_family: Parcubacteria_GWA2-38-13b","Location: CMPE<br />Abundance: 0.0075769318<br />Class_family: Parcubacteria_GWA2-38-13b","Location: LP<br />Abundance: 0.0013253835<br />Class_family: Parcubacteria_GWA2-38-13b"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(223,255,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Parcubacteria_GWA2-38-13b","legendgroup":"Parcubacteria_GWA2-38-13b","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999],"base":[0.122715984783122,0.0899530941371516,0.123935683404455,0.0535326276362592,0.0919593309162193,0.12506673692438,0.0935386468936168,0.12599522217438,0.0950921836648226,0.0963913462871636,0.126853037434792,0.08142022314508,0.127521329436294],"x":[8,5,8,7,5,8,5,8,5,5,8,3,8],"y":[0.00121969862133355,0.00200623677906772,0.00113105351992514,0.00204534493927006,0.00157931597739747,0.000928485249999167,0.00155353677120577,0.000857815260412231,0.00129916262234107,0.00118768660907295,0.00066829200150223,0.0015009750113629,0.000596159765444498],"text":["Location: LSJ<br />Abundance: 0.0012196986<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CS<br />Abundance: 0.0020062368<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LSJ<br />Abundance: 0.0011310535<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LP<br />Abundance: 0.0020453449<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CS<br />Abundance: 0.0015793160<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LSJ<br />Abundance: 0.0009284852<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CS<br />Abundance: 0.0015535368<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LSJ<br />Abundance: 0.0008578153<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CS<br />Abundance: 0.0012991626<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CS<br />Abundance: 0.0011876866<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LSJ<br />Abundance: 0.0006682920<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: CMPE<br />Abundance: 0.0015009750<br />Class_family: Phycisphaerae_Phycisphaeraceae","Location: LSJ<br />Abundance: 0.0005961598<br />Class_family: Phycisphaerae_Phycisphaeraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(228,208,10,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Phycisphaerae_Phycisphaeraceae","legendgroup":"Phycisphaerae_Phycisphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0,0.00186600765273142,0.00332646471690964],"x":[1,1,1],"y":[0.00186600765273142,0.00146045706417822,0.00110235232957037],"text":["Location: BSJ<br />Abundance: 0.0018660077<br />Class_family: Pla3_lineage_Pla3_lineage","Location: BSJ<br />Abundance: 0.0014604571<br />Class_family: Pla3_lineage_Pla3_lineage","Location: BSJ<br />Abundance: 0.0011023523<br />Class_family: Pla3_lineage_Pla3_lineage"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(0,255,255,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Pla3_lineage_Pla3_lineage","legendgroup":"Pla3_lineage_Pla3_lineage","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9],"base":[0.0675143895167727,0.12594109117814,0.109665783629819,0.111832622031532,0.113807300165637,0.115661406733169,0.117287429634382,0.00777625686707842,0.118829641945541,0.120089482526394,0.0782308429713651,0.0767275649964118,0.0800950870968121,0.0818399498223644,0.0834850433204768,0.0487652824039223,0.0793461431700685,0.085116982193002,0.0864133126327344,0.121331155366667,0.0876994903264622,0.0507294084292903,0.0522020567855647,0.122070456676628,0.0889186293026815],"x":[3,6,8,8,8,8,8,2,8,8,5,3,5,5,5,7,3,5,5,8,5,7,7,8,5],"y":[0.00921317547963918,0.0160102354331867,0.00216683840171313,0.00197467813410478,0.00185410656753243,0.00162602290121247,0.00154221231115895,0.00188505481183357,0.00125984058085388,0.0012416728402728,0.00186424412544695,0.00261857817365667,0.00174486272555234,0.00164509349811234,0.00163193887252525,0.00196412602536791,0.00207407997501149,0.00129633043973239,0.00128617769372778,0.000739301309960541,0.00121913897621934,0.00147264835627444,0.00133057085069447,0.000645528106493848,0.0010344648344701],"text":["Location: CMPE<br />Abundance: 0.0092131755<br />Class_family: Planctomycetes_Pirellulaceae","Location: LLC<br />Abundance: 0.0160102354<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0021668384<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0019746781<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0018541066<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0016260229<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0015422123<br />Class_family: Planctomycetes_Pirellulaceae","Location: CMP<br />Abundance: 0.0018850548<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0012598406<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0012416728<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0018642441<br />Class_family: Planctomycetes_Pirellulaceae","Location: CMPE<br />Abundance: 0.0026185782<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0017448627<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0016450935<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0016319389<br />Class_family: Planctomycetes_Pirellulaceae","Location: LP<br />Abundance: 0.0019641260<br />Class_family: Planctomycetes_Pirellulaceae","Location: CMPE<br />Abundance: 0.0020740800<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0012963304<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0012861777<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0007393013<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0012191390<br />Class_family: Planctomycetes_Pirellulaceae","Location: LP<br />Abundance: 0.0014726484<br />Class_family: Planctomycetes_Pirellulaceae","Location: LP<br />Abundance: 0.0013305709<br />Class_family: Planctomycetes_Pirellulaceae","Location: LSJ<br />Abundance: 0.0006455281<br />Class_family: Planctomycetes_Pirellulaceae","Location: CS<br />Abundance: 0.0010344648<br />Class_family: Planctomycetes_Pirellulaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(2,48,32,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Pirellulaceae","legendgroup":"Planctomycetes_Pirellulaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.9,0.9],"base":[0.106008940471595,0.00984764324802832,0.108337483231645,0.0763811148695151,0.12121675941097,0.0473864048204262],"x":[8,9,8,5,6,7],"y":[0.00232854276004993,0.00362704216451203,0.00132830039817355,0.00184972810185001,0.00472433176716984,0.00137887758349619],"text":["Location: LSJ<br />Abundance: 0.0023285428<br />Class_family: Planctomycetes_Rubinisphaeraceae","Location: LT<br />Abundance: 0.0036270422<br />Class_family: Planctomycetes_Rubinisphaeraceae","Location: LSJ<br />Abundance: 0.0013283004<br />Class_family: Planctomycetes_Rubinisphaeraceae","Location: CS<br />Abundance: 0.0018497281<br />Class_family: Planctomycetes_Rubinisphaeraceae","Location: LLC<br />Abundance: 0.0047243318<br />Class_family: Planctomycetes_Rubinisphaeraceae","Location: LP<br />Abundance: 0.0013788776<br />Class_family: Planctomycetes_Rubinisphaeraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(80,200,120,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Planctomycetes_Rubinisphaeraceae","legendgroup":"Planctomycetes_Rubinisphaeraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9],"base":[0.038062170391322,0.10279207103991,0.115304123603073,0.0442909275130458,0.103944741492849,0.0738643384436174,0.0751769786957601,0.104768088917387,0.105399440154987,0.0461616854041049],"x":[7,8,6,7,8,5,5,8,8,7],"y":[0.00622875712172373,0.00115267045293939,0.00591263580789721,0.00187075789105916,0.000823347424537602,0.00131264025214275,0.00120413617375496,0.000631351237600081,0.000609500316608697,0.00122471941632125],"text":["Location: LP<br />Abundance: 0.0062287571<br />Class_family: Rhodothermia_Balneolaceae","Location: LSJ<br />Abundance: 0.0011526705<br />Class_family: Rhodothermia_Balneolaceae","Location: LLC<br />Abundance: 0.0059126358<br />Class_family: Rhodothermia_Balneolaceae","Location: LP<br />Abundance: 0.0018707579<br />Class_family: Rhodothermia_Balneolaceae","Location: LSJ<br />Abundance: 0.0008233474<br />Class_family: Rhodothermia_Balneolaceae","Location: CS<br />Abundance: 0.0013126403<br />Class_family: Rhodothermia_Balneolaceae","Location: CS<br />Abundance: 0.0012041362<br />Class_family: Rhodothermia_Balneolaceae","Location: LSJ<br />Abundance: 0.0006313512<br />Class_family: Rhodothermia_Balneolaceae","Location: LSJ<br />Abundance: 0.0006095003<br />Class_family: Rhodothermia_Balneolaceae","Location: LP<br />Abundance: 0.0012247194<br />Class_family: Rhodothermia_Balneolaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(95,133,117,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Rhodothermia_Balneolaceae","legendgroup":"Rhodothermia_Balneolaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999],"base":[0.101165211888536,0.0983686623971382,0.0704780614243862,0.100619929497945,0.101811895534508],"x":[6,8,5,8,8],"y":[0.0141389117145368,0.0022512671008071,0.00338627701923117,0.00119196603656234,0.000980175505401915],"text":["Location: LLC<br />Abundance: 0.0141389117<br />Class_family: Saccharimonadia_Saccharimonadales","Location: LSJ<br />Abundance: 0.0022512671<br />Class_family: Saccharimonadia_Saccharimonadales","Location: CS<br />Abundance: 0.0033862770<br />Class_family: Saccharimonadia_Saccharimonadales","Location: LSJ<br />Abundance: 0.0011919660<br />Class_family: Saccharimonadia_Saccharimonadales","Location: LSJ<br />Abundance: 0.0009801755<br />Class_family: Saccharimonadia_Saccharimonadales"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(79,121,66,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Saccharimonadia_Saccharimonadales","legendgroup":"Saccharimonadia_Saccharimonadales","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.9,0.9],"base":[0.0946841791154435,0.0288067478774587,0.0964134407163563,0.0649449238193548,0.0311003133298421,0.0331565342876873,0.0974744227489715,0.0351111166957173,0.0366655788495626,0.00524344068973961,0.0695048868673758],"x":[8,7,8,3,7,7,8,7,7,4,5],"y":[0.00172926160091287,0.00229356545238343,0.00106098203261512,0.00256946569741784,0.00205622095784519,0.00195458240803004,0.000894239648166775,0.00155446215384524,0.00139659154175946,0.0011216638588876,0.000973174557010459],"text":["Location: LSJ<br />Abundance: 0.0017292616<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LP<br />Abundance: 0.0022935655<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LSJ<br />Abundance: 0.0010609820<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: CMPE<br />Abundance: 0.0025694657<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LP<br />Abundance: 0.0020562210<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LP<br />Abundance: 0.0019545824<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LSJ<br />Abundance: 0.0008942396<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LP<br />Abundance: 0.0015544622<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: LP<br />Abundance: 0.0013965915<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: CMPW<br />Abundance: 0.0011216639<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","Location: CS<br />Abundance: 0.0009731746<br />Class_family: SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(34,139,34,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","legendgroup":"SAR324_clade(Marine_group_B)_SAR324_clade(Marine_group_B)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999],"base":[0,0.0932055503494572,0.0940435154202603],"x":[4,8,8],"y":[0.00524344068973961,0.000837965070803137,0.000640663695183147],"text":["Location: CMPW<br />Abundance: 0.0052434407<br />Class_family: Sericytochromatia_Sericytochromatia","Location: LSJ<br />Abundance: 0.0008379651<br />Class_family: Sericytochromatia_Sericytochromatia","Location: LSJ<br />Abundance: 0.0006406637<br />Class_family: Sericytochromatia_Sericytochromatia"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(124,252,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Sericytochromatia_Sericytochromatia","legendgroup":"Sericytochromatia_Sericytochromatia","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0275482281746056,"x":[7],"y":[0.00125851970285307],"text":"Location: LP<br />Abundance: 0.0012585197<br />Class_family: Spirochaetia_Spirochaetaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(242,210,189,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Spirochaetia_Spirochaetaceae","legendgroup":"Spirochaetia_Spirochaetaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0257142369895912,"x":[7],"y":[0.00183399118501444],"text":"Location: LP<br />Abundance: 0.0018339912<br />Class_family: Synergistia_Synergistaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(255,172,28,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Synergistia_Synergistaceae","legendgroup":"Synergistia_Synergistaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9],"base":[0.0222911165911461,0.024481872542518],"x":[7,7],"y":[0.0021907559513719,0.0012323644470732],"text":["Location: LP<br />Abundance: 0.0021907560<br />Class_family: Syntrophobacteria_uncultured","Location: LP<br />Abundance: 0.0012323644<br />Class_family: Syntrophobacteria_uncultured"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(205,127,50,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Syntrophobacteria_uncultured","legendgroup":"Syntrophobacteria_uncultured","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9],"base":[0.0607949473071616,0.0863012077187847,0.0879166307817064,0.0584269747840707,0.0894504497014098,0.090596609246277,0.00655131907376054,0.0645454123307173,0.0612858417766541,0.0659857794859755,0.0917027394504927,0.0924639077224442,0.0673310662250861,0.0633719426399042,0.0684839912748691],"x":[5,8,8,3,8,8,2,5,3,5,8,8,5,3,5],"y":[0.00375046502355571,0.00161542306292173,0.00153381891970336,0.00285886699258341,0.00114615954486719,0.00110613020421577,0.00122493779331788,0.00144036715525822,0.00208610086325014,0.00134528673911057,0.000761168271951482,0.00074164262701297,0.00115292504978301,0.00157298117945058,0.00102089559250666],"text":["Location: CS<br />Abundance: 0.0037504650<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0016154231<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0015338189<br />Class_family: Thermoleophilia_67-14","Location: CMPE<br />Abundance: 0.0028588670<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0011461595<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0011061302<br />Class_family: Thermoleophilia_67-14","Location: CMP<br />Abundance: 0.0012249378<br />Class_family: Thermoleophilia_67-14","Location: CS<br />Abundance: 0.0014403672<br />Class_family: Thermoleophilia_67-14","Location: CMPE<br />Abundance: 0.0020861009<br />Class_family: Thermoleophilia_67-14","Location: CS<br />Abundance: 0.0013452867<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0007611683<br />Class_family: Thermoleophilia_67-14","Location: LSJ<br />Abundance: 0.0007416426<br />Class_family: Thermoleophilia_67-14","Location: CS<br />Abundance: 0.0011529250<br />Class_family: Thermoleophilia_67-14","Location: CMPE<br />Abundance: 0.0015729812<br />Class_family: Thermoleophilia_67-14","Location: CS<br />Abundance: 0.0010208956<br />Class_family: Thermoleophilia_67-14"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(218,160,109,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Thermoleophilia_67-14","legendgroup":"Thermoleophilia_67-14","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":0.9,"base":0.0202500557910061,"x":[7],"y":[0.0020410608001399],"text":"Location: LP<br />Abundance: 0.0020410608<br />Class_family: Thermotogae_Petrotogaceae","type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(204,85,0,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Thermotogae_Petrotogaceae","legendgroup":"Thermotogae_Petrotogaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.899999999999999,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9,0.9,0.9,0.899999999999999,0.9,0.9,0.899999999999999,0.9,0.899999999999999,0.9,0.9],"base":[0.0135738878470094,0,0.0262287377621523,0.0382793834065188,0.0133754384939945,0.00977025997406159,0.0249846595393397,0.0496494192691809,0.040093824657002,0.0545327003235769,0.0585400978099777,0.0732652697935088,0.0358002203237834,0.0623771970160267,0,0.0417826746473597,0.0466724673208112,0.0657145888375851,0.0512449757488052,0.0683917966924771,0.0708686361630621,0.0139430729881299,0.0732847687838777,0.0752958849538053,0.0772363372538469,0.07916234011444,0.0505495230851148,0.0810860855186913,0.00193103175846366,0.0824869471846401,0.0556441973358387,0.083745452468008,0.0181373936769264,0.0575446740046566,0.0591974548555639,0.00365319889176226,0.084815346332779,0.0963620028774503,0.0551033286280632,0.00803412216577231,0.00486728931984781,0.0856501622780473,0.00574931957660669,0.0569158129305508],"x":[8,6,8,8,3,5,5,8,3,8,8,6,5,8,9,5,5,8,5,8,8,7,8,8,8,8,3,8,2,8,5,8,7,5,5,2,8,6,3,9,2,8,2,3],"y":[0.0126548499151429,0.0732652697935088,0.0120506456443665,0.0113700358626621,0.0267183861630076,0.0152143995652781,0.0108155607844437,0.00488328105439598,0.0104556984281128,0.00400739748640082,0.00383709920604897,0.0230967330839414,0.00598245432357632,0.00333739182155837,0.00803412216577231,0.00488979267345157,0.00457250842799392,0.00267720785489205,0.00439922158703352,0.00247683947058494,0.00241613262081568,0.00419432068879646,0.00201111616992752,0.00194045230004161,0.00192600286059308,0.00192374540425136,0.00455380554294844,0.00140086166594877,0.0017221671332986,0.00125850528336795,0.00190047666881792,0.00106989386477097,0.00211266211407976,0.00165278085090727,0.00159749245159773,0.00121409042808555,0.000834815945268283,0.00480320901108569,0.00181248430248752,0.00181352108225601,0.000882030256758874,0.000651045440737408,0.000801999497153852,0.00151116185351992],"text":["Location: LSJ<br />Abundance: 0.0126548499<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LLC<br />Abundance: 0.0732652698<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0120506456<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0113700359<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMPE<br />Abundance: 0.0267183862<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0152143996<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0108155608<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0048832811<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMPE<br />Abundance: 0.0104556984<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0040073975<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0038370992<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LLC<br />Abundance: 0.0230967331<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0059824543<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0033373918<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LT<br />Abundance: 0.0080341222<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0048897927<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0045725084<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0026772079<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0043992216<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0024768395<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0024161326<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LP<br />Abundance: 0.0041943207<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0020111162<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0019404523<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0019260029<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0019237454<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMPE<br />Abundance: 0.0045538055<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0014008617<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMP<br />Abundance: 0.0017221671<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0012585053<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0019004767<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0010698939<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LP<br />Abundance: 0.0021126621<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0016527809<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CS<br />Abundance: 0.0015974925<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMP<br />Abundance: 0.0012140904<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0008348159<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LLC<br />Abundance: 0.0048032090<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMPE<br />Abundance: 0.0018124843<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LT<br />Abundance: 0.0018135211<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMP<br />Abundance: 0.0008820303<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: LSJ<br />Abundance: 0.0006510454<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMP<br />Abundance: 0.0008019995<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae","Location: CMPE<br />Abundance: 0.0015111619<br />Class_family: Verrucomicrobiae_Chthoniobacteraceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(233,116,81,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Chthoniobacteraceae","legendgroup":"Verrucomicrobiae_Chthoniobacteraceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999],"base":[0.0112545769657835,0.0128835642410621],"x":[7,8],"y":[0.00268849602234647,0.000690323605947341],"text":["Location: LP<br />Abundance: 0.0026884960<br />Class_family: Verrucomicrobiae_Puniceicoccaceae","Location: LSJ<br />Abundance: 0.0006903236<br />Class_family: Verrucomicrobiae_Puniceicoccaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(227,150,62,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Puniceicoccaceae","legendgroup":"Verrucomicrobiae_Puniceicoccaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.899999999999999,0.9,0.9,0.899999999999999,0.899999999999999,0.9,0.9,0.9,0.899999999999999,0.899999999999999,0.9],"base":[0.0032851181633813,0.00152641526579278,0,0.00782463776019727,0.00930324038011835,0.00494648232335415,0.00978787104021816,0.00740296584650991,0.0107701780427565,0.0120955689498608,0.00989101367293676],"x":[8,3,2,8,8,5,3,5,8,8,7],"y":[0.00453951959681598,0.00826145577442538,0.00193103175846366,0.00147860261992108,0.0014669376626381,0.00245648352315576,0.00358756745377632,0.00236729412755167,0.00132539090710433,0.000787995291201282,0.0013635632928467],"text":["Location: LSJ<br />Abundance: 0.0045395196<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: CMPE<br />Abundance: 0.0082614558<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: CMP<br />Abundance: 0.0019310318<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: LSJ<br />Abundance: 0.0014786026<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: LSJ<br />Abundance: 0.0014669377<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: CS<br />Abundance: 0.0024564835<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: CMPE<br />Abundance: 0.0035875675<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: CS<br />Abundance: 0.0023672941<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: LSJ<br />Abundance: 0.0013253909<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: LSJ<br />Abundance: 0.0007879953<br />Class_family: Verrucomicrobiae_Rubritaleaceae","Location: LP<br />Abundance: 0.0013635633<br />Class_family: Verrucomicrobiae_Rubritaleaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(242,140,40,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Rubritaleaceae","legendgroup":"Verrucomicrobiae_Rubritaleaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.899999999999999,0.899999999999999,0.9,0.899999999999999,0.9,0.9],"base":[0,0,0.00131012299281951,0.0026071169467455,0.00259529991878531,0.00385487822647586,0],"x":[5,8,8,5,8,5,3],"y":[0.0026071169467455,0.00131012299281951,0.0012851769259658,0.00124776127973036,0.000689818244595989,0.00109160409687829,0.00152641526579278],"text":["Location: CS<br />Abundance: 0.0026071169<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: LSJ<br />Abundance: 0.0013101230<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: LSJ<br />Abundance: 0.0012851769<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: CS<br />Abundance: 0.0012477613<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: LSJ<br />Abundance: 0.0006898182<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: CS<br />Abundance: 0.0010916041<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae","Location: CMPE<br />Abundance: 0.0015264153<br />Class_family: Verrucomicrobiae_Terrimicrobiaceae"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(210,125,45,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"Verrucomicrobiae_Terrimicrobiaceae","legendgroup":"Verrucomicrobiae_Terrimicrobiaceae","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[0.9,0.9,0.9],"base":[0,0.00674330930770362,0.00862502665183207],"x":[7,7,7],"y":[0.00674330930770362,0.00188171734412845,0.00126598702110469],"text":["Location: LP<br />Abundance: 0.0067433093<br />Class_family: WWE3_WWE3","Location: LP<br />Abundance: 0.0018817173<br />Class_family: WWE3_WWE3","Location: LP<br />Abundance: 0.0012659870<br />Class_family: WWE3_WWE3"],"type":"bar","textposition":"none","marker":{"autocolorscale":false,"color":"rgba(184,115,51,1)","line":{"width":1.13385826771654,"color":"rgba(0,0,0,1)"}},"name":"WWE3_WWE3","legendgroup":"WWE3_WWE3","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":49.7467828974678,"l":46.8244084682441},"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.4,9.6],"tickmode":"array","ticktext":["BSJ","CMP","CMPE","CMPW","CS","LLC","LP","LSJ","LT"],"tickvals":[1,2,3,4,5,6,7,8,9],"categoryorder":"array","categoryarray":["BSJ","CMP","CMPE","CMPW","CS","LLC","LP","LSJ","LT"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.65296803652968,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":10.6268161062682},"tickangle":-90,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"Location","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.05,1.05],"tickmode":"array","ticktext":["0.00","0.25","0.50","0.75","1.00"],"tickvals":[0,0.25,0.5,0.75,1],"categoryorder":"array","categoryarray":["0.00","0.25","0.50","0.75","1.00"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.65296803652968,"tickwidth":0.132835201328352,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":10.6268161062682},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"Relative Abundance (Family > 1%) <br />","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"transparent","line":{"color":"rgba(0,0,0,1)","width":0.66417600664176,"linetype":"solid"},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":null,"bordercolor":null,"borderwidth":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":9.29846409298464},"title":{"text":"Class_family","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"623c112c73ca":{"x":{},"y":{},"fill":{},"type":"bar"}},"cur_data":"623c112c73ca","visdat":{"623c112c73ca":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```


```r
htmlwidgets::saveWidget(p, "Season_family.html")

htmltools::tags$iframe(
  src=file.path(getwd(), "Season_family.html"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0")
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/Season_family.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

## Creating Season vs. Location barplot with only normal samples


```r
PRSeqR_family_Location <- PR_Normal %>%   
  tax_glom(taxrank = "Family") %>%                    
  #Agglomerate at family level   
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  #Transform to rel.abundance  
  psmelt() %>%                                    
  # Melt to long format  
  group_by(Season, Location, Kingdom, Phylum, Class, Order, Family) %>%   
  dplyr::summarize(Mean =  mean(Abundance, na.rm=TRUE)) %>%                              #Calculate average  
  filter(Mean > 0.01) %>%                         
  #Filter arrange
  arrange(Class) 
```

```
## `summarise()` has grouped output by 'Season', 'Location', 'Kingdom', 'Phylum',
## 'Class', 'Order'. You can override using the `.groups` argument.
```


```r
#Creating a new column that combines both the Class and Family level information #so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed 
PRSeqR_family_Location$Class_family <- paste(PRSeqR_family_Location$Class, PRSeqR_family_Location$Family, sep="_")                                                 #Creating plot   
PRSeqR_Location_family_bar <- ggplot(PRSeqR_family_Location, aes(x = Location, y = Mean, fill = Class_family)) +  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Family > 1%) \n") + xlab("Location") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"),    panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8,                                colour="black"),    axis.text.y = element_text(size=8, colour="black"),    plot.title = element_text(hjust = 0.5),    axis.ticks.x = element_line(colour="#000000", size=0.1),    axis.ticks.y = element_line(colour="#000000", size=0.1),   legend.text = element_text(size = 7)) +   facet_grid(. ~ Location, margins = TRUE, scale="free")

print(PRSeqR_Location_family_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-4-1.png" width="672" />


```r
p4 <- ggplotly(PRSeqR_Location_family_bar)
```


```r
htmlwidgets::saveWidget(p,"PRSeqR_Location_family_bar")
htmltools::tags$iframe(
  src=file.path(getwd(),"PRSeqR_Location_family_bar"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0")
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/PRSeqR_Location_family_bar" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

## Creating location vs. year barplot


```r
PR_Normal <- PRphyseq %>% subset_samples(Sample_Type %in% c("N"))

PR_2022 <- PRphyseq %>% subset_samples(Year %in% c("2022"))

##Creating a subset of samples for only normal samples.This will exclude the enriched samples
```


```r
PRSeqR_family_Year <- PRphyseq %>%    
  tax_glom(taxrank = "Family") %>%     
  transform_sample_counts(function(x) {x/sum(x)}) %>%   
  psmelt() %>%    
  group_by(Year,Season, Kingdom, Phylum, Class, Order, Family)%>%   
  dplyr::summarize(Mean =                      
                     mean(Abundance, na.rm=TRUE)) %>%    
  filter(Mean > 0.01) %>%                                
  arrange(Class)  
```

```
## `summarise()` has grouped output by 'Year', 'Season', 'Kingdom', 'Phylum',
## 'Class', 'Order'. You can override using the `.groups` argument.
```

```r
## `summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument.
```


```r
#Creating a new column that combines both the Class and Family level information #so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed 
PRSeqR_family_Year$Class_family <- paste(PRSeqR_family_Year$Class, PRSeqR_family_Year$Family, sep="_")                                                 #Creating plot   
PRSeqR_Year_family_bar <- ggplot(PRSeqR_family_Year, aes(x = Year, y = Mean, fill = Class_family)) +  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Family > 1%) \n") + xlab("Year") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"),    panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8,                                colour="black"),    axis.text.y = element_text(size=8, colour="black"),    plot.title = element_text(hjust = 0.5),    axis.ticks.x = element_line(colour="#000000", size=0.1),    axis.ticks.y = element_line(colour="#000000", size=0.1),   legend.text = element_text(size = 7)) +   facet_grid(. ~ Season, margins = TRUE, scale="free")

print(PRSeqR_Year_family_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-6-1.png" width="672" />


```r
p5 <- ggplotly(PRSeqR_Year_family_bar)
```

Note that 2023 is only enriched samples so this year cannot be accurately depicted and 2022 is missing from the wet season.

\*Now that we know there is an error in the 2022 data when attempting to generate the 2022 NMDS plot, we will be creating a family plot for the 2022 samples to determine what is too similar - this is what we think is preventing the NMDS plot from being generated.


```r
PRSeqR_2022_Sample <- PR_2022 %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>%
  group_by(Location,Season,Kingdom,Phylum,Class,Order,Family,Genus,Species)%>%   
  dplyr::summarize(Mean =          
 mean(Abundance, na.rm=TRUE)) %>%    
  filter(Mean > 0.00) %>%                                
  arrange(Class)  
```

```
## `summarise()` has grouped output by 'Location', 'Season', 'Kingdom', 'Phylum',
## 'Class', 'Order', 'Family', 'Genus'. You can override using the `.groups`
## argument.
```

```r
##`summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument
```


```r
#Creating a new column that combines both the Class and Family level information #so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed 
PRSeqR_2022_Sample$Family <- paste(PRSeqR_2022_Sample$Class, PRSeqR_2022_Sample$Family, sep="_")                                                 #Creating plot   
PRSeqR_2022_Sample_bar <- ggplot(PRSeqR_2022_Sample, aes(x = Location, y = Mean, fill = Family)) +  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Family > 1%) \n") + xlab("Location") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"),    panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8,                                colour="black"),    axis.text.y = element_text(size=8, colour="black"),    plot.title = element_text(hjust = 0.5),    axis.ticks.x = element_line(colour="#000000", size=0.1),    axis.ticks.y = element_line(colour="#000000", size=0.1),   legend.text = element_text(size = 7))

print(PRSeqR_2022_Sample_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-9-1.png" width="672" />


```r
 p6 <- ggplotly(PRSeqR_2022_Sample_bar)
```

This is showing the samples averaged together, but there are multiple samples from each site. We are going to plot based on sample name to give us individual bars for each sample to better understand the similarities between each other.


```r
PRSeqR_2022_Sample_ID <- PR_2022 %>%    
  tax_glom(taxrank = "Family") %>%     
  transform_sample_counts(function(x) {x/sum(x)}) %>%   
  psmelt() %>%    
  group_by(Sample_ID,Location,Season,Kingdom,Phylum,Class,Order,Family,Genus,Species)%>%   
  dplyr::summarize(Mean =                      
                     mean(Abundance, na.rm=TRUE)) %>%    
  filter(Mean > 0.00) %>%                                
  arrange(Class)  
```

```
## `summarise()` has grouped output by 'Sample_ID', 'Location', 'Season',
## 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'. You can override
## using the `.groups` argument.
```

```r
## `summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument
```


```r
#Creating a new column that combines both the Class and Family level information #so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed 
PRSeqR_2022_Sample_ID$Family <- paste(PRSeqR_2022_Sample_ID$Class, PRSeqR_2022_Sample_ID$Family, sep="_")                                                 #Creating plot   
PRSeqR_2022_Sample_ID_bar <- ggplot(PRSeqR_2022_Sample_ID, aes(x = Sample_ID, y = Mean, fill = Family)) +  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Family > 1%) \n") + xlab("Sample_ID") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"),    panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8,                                colour="black"),    axis.text.y = element_text(size=8, colour="black"),    plot.title = element_text(hjust = 0.5),    axis.ticks.x = element_line(colour="#000000", size=0.1),    axis.ticks.y = element_line(colour="#000000", size=0.1),   legend.text = element_text(size = 7))

print(PRSeqR_2022_Sample_ID_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-12-1.png" width="672" />


```r
p7 <- ggplotly(PRSeqR_2022_Sample_ID_bar)
```


```r
htmlwidgets::saveWidget(p,"PRSeqR_2022_Sample_ID_bar")
htmltools::tags$iframe(
  src=file.path(getwd(),"PRSeqR_2022_Sample_ID_bar"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0")
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/PRSeqR_2022_Sample_ID_bar" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

## Creating barplots for enriched samples


```r
#PR_ENR <- PRphyseq %>% subset_samples(Sample_Type %in% c("E"))

##Creating a subset of samples for only enriched samples.This will exclude the normal samples
```


```r
#PRSeqR_family_Location <- PR_ENR %>%    
#  tax_glom(taxrank = "Species") %>%     
#  transform_sample_counts(function(x) {x/sum(x)}) %>%   
#  psmelt() %>%    
#  group_by(Location,Season,Kingdom,Phylum,Class,Order,Family,Genus,Species)%>%   
#  dplyr::summarize(Mean =                      
#                     mean(Abundance, na.rm=TRUE)) %>%    
#  filter(Mean > 0.00) %>%                                
#  arrange(Class)  

## `summarise()` has grouped output by 'Season', 'Kingdom', 'Phylum', 'Class',
## 'Order'. You can override using the `.groups` argument
```


```r
#Creating a new column that combines both the Class and Family level information #so that Family level identifiers that are not unique to one Class (such as Ambiguous_taxa) aren't merged into a single large category when graphed 
PRSeqR_family_Location$Genus_species <- paste(PRSeqR_family_Location$Class, PRSeqR_family_Location$Species, sep="_")                                                 #Creating plot   
```

```
## Warning: Unknown or uninitialised column: `Species`.
```

```r
PRSeqR_Location_family_bar <- ggplot(PRSeqR_family_Location, aes(x = Location, y = Mean, fill = Genus_species)) +  geom_bar(stat="identity", colour = "black", size=0.3, position = "fill") +   scale_fill_manual(values = c("black",                                                                 "#440154FF", "#450659FF","#460B5EFF",                                                                "#472D7AFF",                                                                 "#3B518BFF", "#3A548CFF", "#38598CFF", "#365C8DFF",                                 "#34608DFF", "#33638DFF", "#31678EFF", "#306A8EFF",                                "#2E6E8EFF", "#2D718EFF", "#2B748EFF",                                "#2A778EFF", "#297B8EFF",                                                                "#1F958BFF", "#1F988BFF", "#1E9C89FF", "#1F9F88FF",                                "#1FA287FF", "#21A585FF", "#23A983FF",                                "#25AC82FF", "#29AF7FFF", "#2DB27DFF",                                "#32B67AFF", "#37B878FF", "#3CBC74FF",                                                                "#57C766FF", "#5EC962FF",                                                                "#A2DA37FF",                                                                "#DAE319FF", "#E4E419FF",  "#ECE51BFF", "#F5E61FFF",                                                                "darkorchid1", "darkorchid2", "darkorchid3",                                                                 "#A21C9AFF", "#A62098FF", "#AB2394FF", "#AE2892FF", "#B22B8FFF", "#B6308BFF", "#BA3388FF", "#BE3885FF", "#C13B82FF", "#C53F7EFF",                                "#C8437BFF", "#CC4678FF", "#CE4B75FF", "#D14E72FF", "#D5536FFF", "#D7566CFF", "#DA5B69FF", "#DD5E66FF", "#E06363FF","#00CCCC", "#006600",                                                                "#FF6699",                                                                "#F68D45FF", "#FCA537FF",                                                                "#F6E726FF", "#F4ED27FF",                                                                                                 "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018",                                "#00489C", "#CCCCCC", "#999999", "#A1C299", "#300018","#191970","#353935","#F0FFFF","#89CFF0","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58")) +      guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +   ylab("Relative Abundance (Species > 0) \n") + xlab("Location") +   theme_minimal()+theme(panel.border = element_rect(fill=NA, colour = "black"),    panel.grid.minor = element_blank(),   panel.grid.major = element_blank(),    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust=1, size=8,                                colour="black"),    axis.text.y = element_text(size=8, colour="black"),    plot.title = element_text(hjust = 0.5),    axis.ticks.x = element_line(colour="#000000", size=0.1),    axis.ticks.y = element_line(colour="#000000", size=0.1),   legend.text = element_text(size = 7))

print(PRSeqR_Location_family_bar)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-14-1.png" width="672" />


```r
p8 <- ggplotly(PRSeqR_Location_family_bar)
```


```r
#(2/2) This step is needed to create an interactive plot on gtihub. After it has been pushed to github you can download the html doc and it will open as the interactive plot via your web browser.
htmlwidgets::saveWidget(p, "PRSeqR_Location_family_bar")

htmltools::tags$iframe(
  src=file.path(getwd(), "PRSeqR_Location_family_bar"),
  width="100%",
  height="600",
  scrolling="no",
  seamless="seamless",
  frameBorder="0"
)
```

```{=html}
<iframe src="C:/Users/jhanlon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/PR-S003-S015/PRSeqR_Location_family_bar" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0"></iframe>
```

# Transforming Phyloseq object to relative abundance


```r
#Transform to relative abundance  
PR_RA_Normal = transform_sample_counts(PR_Normal, function(x){x / sum(x)})    

#View the otu_table to see the transformed relative abundances  
View(otu_table(PR_RA_Normal)) 
```

# NMDS Plots

## Choose colors and shapes


```r
#Choose colors
plot.colors <- c("steelblue2","purple4","darkorange","firebrick","springgreen4", "gold", "darkblue", "yellowgreen","turquoise4", "orange","indianred","darkslategrey", "lightblue","darkgreen","mediumaquamarine","gray48","mediumorchid1", "#5F7FC7","#DA5724", "#508578", "#CBD588","#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D","#8A7C64", "#599861","dodgerblue","darkmagenta", "forestgreen","steelblue1", "cyan","mediumorchid3", "cadetblue3", "yellow", "#00FF7F","#40E0D0","#E97451","#0000FF","#7393B3","#088F8F","#0096FF","#5F9EA0","#0047AB","#6495ED","#00FFFF","#00008B","#6F8FAF","#1434A4","#7DF9FF","#6082B6","#00A36C","#3F00FF","#5D3FD3","#ADD8E6","#191970","#000080","#1F51FF","#A7C7E7","#CCCCFF","#B6D0E2","#96DED1","#4169E1","#0F52BA","#9FE2BF","#87CEEB","#4682B4","#008080","#40E0D0","#0437F2","#40B5AD","#0818A8","#AAFF00","#5F9EA0","#097969","#AFE1AF","#DFFF00","#E4D00A","#00FFFF","#023020","#50C878","#5F8575","#4F7942","#228B22","#7CFC00","#F2D2BD","#FFAC1C","#CD7F32","#DAA06D","#CC5500","#E97451","#E3963E","#F28C28","#D27D2D","#B87333","#FF7F50","#F88379","#DA70D6","#C3B1E1","#CCCCFF","#673147","#A95C68","#800080","#51414F","#953553","#D8BFD8","#630330","#7F00FF","#722F37","#BDB5D5","#FCF55F","#FFFFF0","#F8DE7E","#F0E68C","#FAFA33","#FBEC5D","#F4BB44","#FFDB58", "#C4B454")
plot.shape <- c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21)
```

## Subset the samples

Subsetting from PR_physeq_RA_Normal, which includes only the non-enriched samples. This dataframe has been filtered so that ASVs that were 0 in all samples are excluded.


```r
#PR_S003 <- PR_physeq_RA_Normal %>% subset_samples(Project %in% c("S003"))

#PR_S015 <- PR_physeq_RA_Normal %>% subset_samples(Project %in% c("S015"))

#PR_2021 <- PR_Normal %>% subset_samples(Year %in% c("2021"))

#PR_2022 <- PR_Normal %>% subset_samples(Year %in% c("2022"))
```

## NMDS of both runs


```r
all.nmds.source.ord <- ordinate( 
  physeq = PR_RA_Normal, 
  method = "NMDS", 
  distance = "bray")
```

```
## Run 0 stress 0.1850681 
## Run 1 stress 0.1850731 
## ... Procrustes: rmse 0.0007928043  max resid 0.00578824 
## ... Similar to previous best
## Run 2 stress 0.2139555 
## Run 3 stress 0.2207744 
## Run 4 stress 0.2022847 
## Run 5 stress 0.2123585 
## Run 6 stress 0.1850805 
## ... Procrustes: rmse 0.001712935  max resid 0.01392611 
## Run 7 stress 0.1871644 
## Run 8 stress 0.1871644 
## Run 9 stress 0.2024291 
## Run 10 stress 0.2035301 
## Run 11 stress 0.2080844 
## Run 12 stress 0.1850833 
## ... Procrustes: rmse 0.003884105  max resid 0.03295702 
## Run 13 stress 0.1850062 
## ... New best solution
## ... Procrustes: rmse 0.005042124  max resid 0.03300655 
## Run 14 stress 0.2179661 
## Run 15 stress 0.2297205 
## Run 16 stress 0.1849911 
## ... New best solution
## ... Procrustes: rmse 0.002372687  max resid 0.01513598 
## Run 17 stress 0.1849877 
## ... New best solution
## ... Procrustes: rmse 0.001707955  max resid 0.015439 
## Run 18 stress 0.1850491 
## ... Procrustes: rmse 0.005735748  max resid 0.05125967 
## Run 19 stress 0.2237442 
## Run 20 stress 0.1870003 
## *** Best solution was not repeated -- monoMDS stopping criteria:
##      1: no. of iterations >= maxit
##     18: stress ratio > sratmax
##      1: scale factor of the gradient < sfgrmin
```

## NMDs plot of location for combined runs


```r
#Plot, color coding by particle type
all.nmds.Location <- plot_ordination(
  physeq = PR_RA_Normal,
  ordination = all.nmds.source.ord) + 
  scale_fill_manual(values = plot.colors, "Location")+
  scale_shape_manual(values = plot.shape, name = "Location") +
  geom_point(mapping = aes(fill = factor(Location), shape = Location, size = 5)) +
  guides(size=FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  stat_ellipse(aes(group=Location))
```

```
## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
## of ggplot2 3.3.4.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

```r
print(all.nmds.Location)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/combined NMDs plot-1.png" width="672" />

## NMDS plot of season


```r
#Plot, color coding by seaosn
all.nmds.Season <- plot_ordination(
  physeq = PR_RA_Normal,
  ordination = all.nmds.source.ord) + 
  scale_fill_manual(values = plot.colors, "Season") +
  scale_shape_manual(values = c("D" = 21, "W" = 21), name = "Season") +
  geom_point(mapping = aes(fill = factor(Season), shape = Season, size = 5)) +
  guides(size=FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 18), 
        axis.title = element_text(size = 15),
        panel.spacing = unit(1, "lines"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        panel.background = element_blank(), 
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(.99, .99),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  stat_ellipse(aes(group=Season))

print(all.nmds.Season)
```

<img src="PR-S003-S015-data-analysis_files/figure-html/unnamed-chunk-16-1.png" width="672" />

# PERMANOVA Analysis


```r
#filt_complete_bray <- phyloseq::distance(PR_RA_Normal, method = "bray")
#filt_complete_bray_df <- data.frame(sample_data(PR_RA_Normal))
#filt_complete <- adonis2(filt_complete_bray ~ mm_Hg*DO_mg_L*DO_LDO*Specific_Conductance_uS_cm, data = filt_complete_bray_df)
#filt_complete
```


```r
#Run PERMANOVA with Bray-Curtis dissimilarity
#filt_complete_bray <- phyloseq::distance(PR_RA_Normal, method = "bray")
#filt_complete_bray_df <- data.frame(sample_data(PR_RA_Normal))
#filt_complete <- adonis2(filt_complete_bray ~ Temperature_C*Location*pH*Salinity_ppt, data = filt_complete_bray_df)
#filt_complete
```


```r
#Run PERMANOVA with Bray-Curtis dissimilarity
#filt_complete_bray <- phyloseq::distance(PR_RA_Normal, method = "bray")
#filt_complete_bray_df <- data.frame(sample_data(PR_RA_Normal))
#filt_complete <- adonis2(filt_complete_bray ~ Temperature_C*Location*DO_LDO*pH*Salinity_ppt, data = filt_complete_bray_df)
#filt_complete
```


```r
#Convert output to dataframe
#filt.complete.pvalue <- as.data.frame(filt_complete)

#Create table of p-values
#knitr::kable(filt.complete.pvalue, row.names = NA) %>% 
 # kableExtra::kable_styling("striped", 
                           # latex_options="scale_down") %>% 
  #kableExtra::scroll_box(width = "100%")
```

## PERMANOVA analysis of location vs pH


```r
#Subset phyloseq object
#filt.mp.particle.RA <- subset_samples(MPfiltRA, sample_type == "Particle")
#sample_data(filt.mp.particle.RA)
```


```r
#Run PERMANOVA with Bray-Curtis dissimilarity
#filt_particle_bray <- phyloseq::distance(filt.mp.particle.RA, method = "bray")
#filt_particle_bray_df <- data.frame(sample_data(filt.mp.particle.RA))
#filt_all <- adonis2(filt_particle_bray ~ particle_type*effluent*week*polymer_type, data = filt_particle_bray_df)
#filt_all
```


```r
#Convert output to dataframe
#filt.all.pvalue <- as.data.frame(filt_all)

#Create table of p-values
#knitr::kable(filt.all.pvalue, row.names = NA) %>% 
  #kableExtra::kable_styling("striped", 
                          # latex_options="scale_down") %>% 
  #kableExtra::scroll_box(width = "100%")
```


```r
#Run PERMANOVA with Bray-Curtis dissimilarity, only looking at particle type (glass vs. plastic)
#filt_particle_only_bray <- phyloseq::distance(filt.mp.particle.RA, method = "bray")
#filt_particle_only_bray_df <- data.frame(sample_data(filt.mp.particle.RA))
#filt_only_particle <- adonis2(filt_particle_only_bray ~ particle_type, data = filt_particle_only_bray_df)
#filt_only_particle
```

# NMDS Plot extras in R: Envfit


```r
##Load libraries and read in your data
library(vegan)
library(ggplot2)
```


```r
##Subset your data to only environmental data and abundance data. 
##I already had a dataframe from the PR_RA_Normal for abundance data, so I just needed to do the same for the sample data (aka environmental data).
#df_env <- as.data.frame(sample_data(PR_RA_Normal))
#df_com <- as.data.frame(otu_table(PR_RA_Normal))
## Individual data frames created, now I am selecting the columns that will be used for this plot and renaming the data frame to com and env. 
#com = df_com[,1:95]
#env = df_env[,15:17, 20, 21,24]
```


```r
#convert com to a matrix and perform the NMDs ordination
#m_df_com = as.matrix(com)

#nmds code
#set.seed(123)
#nmds = metaMDS(m_df_com, distance = "bray")
#nmds
```


```r
##Running envfit function with env dataframe
#en = envfit(nmds,env, permutations = 999, na.rm = TRUE)
```

# RDA - Redundancy Analysis


```r
library(phyloseq)
library(ggplot2)
library(microViz)

# install.packages("remotes")
#remotes::install_github("statdivlab/corncob")
#library(corncob)
#> microViz version 0.12.0 - Copyright (C) 2023 David Barnett
#> ! Website: https://david-barnett.github.io/microViz
#> ✔ Useful?  For citation details, run: `citation("microViz")`
#> ✖ Silence? `suppressPackageStartupMessages(library(microViz))`
#knitr::opts_chunk$set(fig.width = 7, fig.height = 6)
```


```r
#PR_RA_Normal <- corncob::PR_phylo
#PR_phylo
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 36349 taxa and 91 samples ]
#> sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
#> tax_table()   Taxonomy Table:    [ 36349 taxa by 7 taxonomic ranks ]
```




```r
##Redundancy analysis is a constrained ordination method. It displays the microbial variation that can also be explained by selected constraint variables. Behind the scenes, a linear regression model is created for each microbial abundance variable (using the constraints as the explanatory variables) and a PCA is performed using the fitted values of the microbial abundances.

#PR_RA_Normal %>%
 # ps_mutate(
  #  IBD = as.numeric(ibd == "ibd"),
   # Female = as.numeric(gender == "female"),
   #  ) %>%
 # tax_transform("clr", rank = "Genus") %>%
 # ord_calc(
    #constraints = c("IBD", "Female", "Abx."),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
  #  scale_cc = FALSE # doesn't make a difference
 # ) %>%
 # ord_plot(
 #   colour = "DiseaseState", size = 2, alpha = 0.5, shape = "active",
 #   plot_taxa = 1:8
 # )
```
