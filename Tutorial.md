# Tutorial for data analysis
## Step 1. Alignment  
For alignment, we currently used the Feature barcoding analysis pipeline from 10x Genomics. In order to run this,
you will need two specific files called library.csv and features.csv. You can analyze multi-modal
data easily by specifying your library types in the library.csv file. Detailed information on how to specify 
these files can be found [HERE](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).
Alternatively, you can align **separately (Gene expression, Antibody expression.. etc) and combine the aligned matrix 
in the later step.

### Creating count matrix  
```bash
/path_to_cell_ranger/cellranger count \
--id run_id \
--transcriptome /path_to_reference/refdata-cellranger-GRCh38-3.0.0/ \
--library library.csv \
--feature-ref features.csv \
--chemistry SC3Pv3
```

## Step 2. Preprocessing and deconvolution of multiplets
Next step is to process the output matrix from the cell ranger run and generate a deconvoluted matrix along with meta
information using R (sctioseq.pbmc_final.R). The main function used here is HTODemux from the R package ‘Seurat’
but modified to output multiplet classification information. In the future release, this functionality will all be
integrated into one python package ‘SCITO-seq’.
Example data can be found in the `paper/SCITOseq_counts/200K_PBMC_28ab_10batch/` directory.

```R
## Specify your h5 input file and filter in only cell containing droplets
library(Matrix)
library(Seurat)

experiment <- "200k_pbmc_28ab"
raw <- Read10X_h5(paste("/dir/", experiment, "filtered_feature_bc_matrix.h5", sep=""))
tx <- raw$`Gene Expression`[,which(Matrix::colSums(raw$`Gene Expression`)>0)]
ab <- raw$`Antibody Capture`[,which(Matrix::colSums(raw$`Gene Expression`)>0)]
ccd_indices <- intersect(which(Matrix::colSums(tx)>0), which(Matrix::colSums(ab)>0))
ab <- ab[,ccd_indices]
tx <- tx[,ccd_indices]
## Read meta information file (donor assignment, droplet barcode sequences..etc)
Donor <- read.table("/dir/donor_ids.tsv", header=TRUE)
IND_filter <- donor[ccd_indices,]
IND <- as.data.frame(IND_filter$donor_id)
DBC_real <- as.data.frame(IND_filter$cell)

colnames(IND) <- "Demuxlet"
rownames(IND) <- colnames(ab)

colnames(DBC_real) <- "DBC_real"
rownames(DBC_real) <-colnames(ab)

aux <- cbind(IND, DBC_real)
```

Now, detect which batch (pool) each droplet contains. Initial manuscript version uses normalization to mean CD45
expression values (we also recommend using the sum of all batch specific antibody counts). ‘Ab’ object below contains
28 unique antibody * 10 batches = 280 antibodies (28 pool specific antibodies are divided into 10 chunks) and index 1,
29 …. 253 below contain the expression of CD45 antibody.


```R
## Normalizing to the CD45 expression values
ab_hashtag <- rbind(round(colSums(t(t(ab[1:28,])/mean(ab[1,]))),0),
                    round(colSums(t(t(ab[29:56,])/mean(ab[29,]))),0),
... ## insert lines here
... ## insert lines here
                    round(colSums(t(t(ab[253:280,])/mean(ab[253,]))),0))
## Sum of antibody counts
ab_hashtag <- rbind(Matrix::colSums(ab[1:28,]),
                    Matrix::colSums(ab[29:56,]),
... ## insert lines here
... ## insert lines here
                    Matrix::colSums(ab[253:280,]))
rownames(ab_hashtag) = paste("Batch", seq(1:10), sep="")
```

Let’s setup a Seurat object and store classification result (*Note, make sure to define the HTODemux function before running HTODemiux since it's been modified from the original version)

```R
ab_hashtag_obj <- CreateSeuratObject(counts=ab, meta.data=aux)
ab_hashtag_obj[["HTO"]] <- CreateAssayObject(counts = ab_hashtag)
ab_hashtag_obj <- NormalizeData(ab_hashtag_obj, assay="HTO", normalization.method = "CLR")

Rst <- HTODemux(ab_hashtag_obj, assay="HTO", positive.quantile = 0.99)
ab_hashtag_obj <- rst[[1]];
discrete <- rst[[2]];

## Check if classification is in the expected range
table(ab_hashtag_obj$HTO_classification.global)
Idents(ab_hashtag_obj) <- "HTO_maxID"
```

**(Optional)** You can rerun the HTOdemux in case you have many ‘Negative’ classifications. 
Some of the droplet barcodes could be mis-classified as ‘negative’ because of the cutoff value so we recommend 
running the UMAP projection to see if this is the case)  

Then, we are ready to resolve all the multiplets into singlets. Here we used discrete object (this is a binary matrix 
which stores if particular droplet barcode is positive for specific batch, for example, here we have P * D matrix where
p = number of batches, D = number of droplet barcodes) to assign new ‘single-cell’ like barcodes.


```R
ab_resolved <- cbind(ab[1:28, which(discrete[1,]==1)],
                     ab[29:56, which(discrete[2,]==1)],
...
...
                     ab[253:280, which(discrete[10,]==1)])

## Make deconvoluted barcodes unique here
colnames(ab_resolved) <- make.names(colnames(ab_resolved), unique=TRUE)

Classification <- as.character(ab_hashtag_obj$HTO_classification.global)
DBC_real <- as.character(ab_hashtag_obj@meta.data$DBC_real)
Demuxlet <- as.character(ab_hashtag_obj@meta.data@Demuxlet)

## Add meta information to the deconvoluted matrix (optional)
Meta <- data.frame(Batch=c(rep("Batch1", sum(discrete[1,])),
                   Batch=c(rep("Batch2", sum(discrete[2,])),
... ## insert lines here
... ## insert lines here

Type=c(classification[which(discrete[1,]==1)],
       classification[which(discrete[2,]==1)],
... ## insert lines here
... ## insert lines here
Demuxlet=c(Demuxlet[which(discrete[1,]==1)],
           Demuxlet[which(discrete[2,]==1)],
...
...
DBC=c(DBC_real[which(discrete[1,]==1)],
      DBC_real[which(discrete[2,]==1)]))

rownames(meta) <- colnames(ab_resolved)
ab_resolved_obj <- CreateSeuratObject(counts = ab_resolved, meta.data=meta)
ab_resolved_obj[["ADT"]] <- CreatAssayObject(counts = ab_resolved)
DefaultASsay(ab_resoled_obj) <- "ADT"

## Output intermediate files for downstream analysis
Selected_ab_cell_matrix <- ab_resolved_obj[["ADT"]]@counts
write.table(selected_ab_cell_matrix, "200K_selected_ab_cell_matrix.txt", sep=",")

Meta_ab_cell_matrix <- ab_resolved_obj@meta.data
write.table(meta_ab_cell_matrix, "200K_meta_ab_cell_matrix.txt", sep=",")
```

<<<<<<< HEAD
This is the end of processing files using R (this **step 2** will be replaced by the python package in the
future - will be posted as a stand-alone GitHub repository).
Users can will be able to choose either python Scanpy or continue to analyze the data using the Seurat package.

## Step 3. Downstream analysis and visualization
This tutorial covers processing of the intermediate files to generate the final UMAP plot using scanpy package. 

```python
import scanpy
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns

%config Completer.use_jedi =False
sc.settings.verbosity = 4
sc.settingset_figure_params(dpi=80)
sc.set_figure_params(dpi_save=1000, figsize=[5,5])
## read input files from R processed files above
## Note, this input matrix contains ‘Raw’ UMI counts
meta_df = pd.read_table(“200K_meta_ab_cell_matrix.txt”, sep=”,”)
adt_df = pd.read_table(“200K_selected_ab_cell_matrix.txt”, sep=”,”)

## convert dataframe to adata object for downstream processing
adt_adata = sc.AnnData(adt_df.T)
adt_adata.obs = meta_df
```

Now, let’s normalize and generate UMAP figures!

```python
## total counts normalization followed by ‘batch’ normalization
sc.pp.normalize_per_cell(adt_adata, counts_per_cell_after=1e4)

import re
def atoi(text):
   return int(text) if text.isdigit() else text
def natural_keys(text):
   return [atoi(c) for c in re.split('(\d+)', text)]

def normalize_byBatch(adt, batch_var='Batch', max_val=10):
   all_batches = list(set(adt.obs[batch_var]))
   all_batches.sort(key=natural_keys)
   ann_dats=[]
   For b in all_batches:
      Batch_adat = adt[adt.obs[batch_var]==b]
      sc.pp.log1p(batch_adata)
      sc.pp.scale(batch_adat, max_value=max_val)
      ann_dats.append(batch_adat)
   norm_concat = ann_dats[0].concatenate(ann_dats[1:len(ann_dats)+1])
   return norm_concat
ptn_names = pd.DataFrame(['CD45','CD33','CD3'.....'CD235', 'CD61'])
ptn_names.columns=['Antibody']
Ptn_names = ptn_names.set_index('Antibody')
adt_adata.var = ptn_names
```

Run PCA to find greatest variance across all cells. Focusing on top PCs help you to capture meaningful
biological variation and reduce the noise signals.

```python
sc.tl.pca(adt_adata, svd_solver='arpack', use_highly_variable=None)
sc.pl.pca_variance_ratio(adt_adata, log=True)
```

Here, we chose 15 top PCs for computing neighbors and run leiden clustering and UMAP for final figures.

```python
sc.pp.neighbors(adt_adata, n_neighbors=10, n_pcs=15)
sc.tl.umap(adt_adata)
sc.tl.leiden(adt_adata)
sc.pl.umap(adt_adata, color='leiden')
```
**(Optional)** if you want to compare with CyTOF, we arcsinh transform the data and z-score normalized, run PCA and cluster to generate UMAP. 
Ingest function is used to compare SCITO-seq to CyTOF.

