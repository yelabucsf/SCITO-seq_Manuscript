# SCITO-seq: single-cell combinatorial indexed cytometry sequencing

This repository contains SCITO-seq package for demultiplexing and resolving single cell data from SCITO-seq experiments.
It also contains manuscript-related CyTOF data, SCITO-seq data and code and notebooks used for analysis. 

# News
We are currently working on integrating python SCITO-seq package into the unified workflow.

## Important branches  
`master` - package (stable, alpha) + manuscript data  
`solid` - only package (stable, alpha)  
`develop` - developer version of package only   

## Paper-related data/code  
We share the code (R and Jupyter notebooks), CyTOF and SCITO-seq data present in the original manuscript in the `master` 
branch.

## Package installation  
We recommend to create a conda environment or virtualenv with pre-installed scanpy
```bash
git clone https://github.com/yelabucsf/SCITO-seq
cd SCITO-seq
pip install .
```

after that follow specific tutorials in the `SCITO-seq/examples`
