# SCITO-seq: single-cell combinatorial indexed cytometry sequencing

This repository contains SCITO-seq package for demultiplexing and resolving single cell data from SCITO-seq experiments.
It also contains manuscript-related CyTOF data, SCITO-seq data and code and notebooks used for analysis. 

To get the data you need to install [Git LFS](https://git-lfs.github.com) and `git clone` this repository. Downloading ZIP will only download pointers to the data.

# News
We are currently working on integrating python SCITO-seq package into the unified workflow.

### Repository structure:
`paper/code` - R code and jupyter notebooks demonstrating SCITO-seq data analysis  
`paper/CyTOF` - CyTOF data  
`paper/SCITOseq_counts` - count matrices (outputs of `cellranger count`)  
`Tutorial.md` - detailed SCITO-seq data analysis workflow as presented in the manuscript

### Addendum  
Package `nero` featured in some Jupyter notebooks is an unreleased package for handling single cell data. For your analysis you can exclude it and filter demuxlet outputs with `pandas.DataFrame` methods. 
