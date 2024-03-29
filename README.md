---
title: "A brief tutorial for estimating differentiation potency of single cells using SCENT and CCAT"
author:
- name: "Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
  - UCL Cancer Institute, University College London
date: "2020-10-26"
package: SCENT
output:
  BiocStyle::html_document:
    toc_float: true
---
  


# Summary

The main purpose of the `SCENT` package is to provide a means of estimating the differentiation potency of single cells without the need to assume prior biological knowledge such as marker expression or timepoint. This may be particularly important in scenarios where a high dropout rate may preclude the use of a marker gene, in snapshot scRNA-Seq datasets of complex tissues where differentiation hierarchies are not well-established, or in cancer tissue where one may want to identify putative cancer stem-cell phenotypes.

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/SCENT")
```


# References

Teschendorff AE, Enver T.  Single-cell entropy for accurate estimation of differentiation potency from a cell's transcriptome. Nat Commun. 2017 Jun 1;8:15599. doi: 10.1038/ncomms15599


 
