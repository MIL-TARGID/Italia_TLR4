# Italia_TLR4
scRNAseq reanalysis from Ke2024 (EGAS0000000382) and Adam2020 (GSE136831)

```
.
├── Adams2022_IPF
│   ├── analysis
│   │  ├── notebooks
│   │  └── plots
│   └── data


└── Ke2024_IBD
    ├── analysis
    │   ├── notebooks
    │   │   └── nb01_gene_expression_signature.ipynb
    │   └── plots
    ├── data
    │   ├── cellranger
    │   ├── fastq
    │   └── out_ambientRNA
    │       ├── anndata
    │       └── metrics
    └── preprocessing
        ├── adata
        ├── notebooks
        │   ├── nb00_ambientRNA_contamination_plots.ipynb
        │   ├── nb01_barcode_rank_plot.ipynb
        │   ├── nb02_first_qc.ipynb
        │   ├── nb03_dedoublet.ipynb
        │   ├── nb04_clustering_after_dedoublet.ipynb
        │   ├── nb05_final_qc_celltypist.ipynb
        │   └── nb06_mesenchymal.ipynb
        ├── plots
        ├── rank_data
        └── scripts_logs
            ├── s01_cellranger_run.sh
            ├── s02_ambient_removal.py
            ├── s03_add_metadata.py
            └── s04_concatenate_anndata.py
```

<br>

# Ke2024_IBD

## 1. preprocessing/scripts_logs

### s01_cellranger

`./cellranger_run.sh`
- in: `data/fastq` + reference genome
- out: `data/cellranger`

### s02_ambientRNA

cellranger output (out/ directories) screened by SoupX for ambient RNA<br>
based on cellranger's clustering as it is not sensitive to the clustering type/granularity<br>

conda env rpy2_ambient_doublet<br>

`python ambient_removal.py ../../data/cellranger`

- in: `data/cellranger`
- out: `data/out_ambientRNA`

### s03_add_metadata

add the metadata to every individual h5ad which is ambientRNA output<br>
also adds sample column to obs and prepends obs_names (barcodes) with sample name<br>

env is norm_scRNAseq<br>

`python add_metadata.py ../../data/out_ambientRNA/anndata metadata.csv`<br>

- in: `data/out_ambientRNA/anndata` directory + `metadata.csv` file
- out: updated .h5ad files in `data/out_ambientRNA/anndata`

### s04_concatenate

concatenating all .h5ad files in a given directory into one

env is norm_scRNAseq<br>

`python concatenate_anndata.py ../../data/out_ambientRNA/anndata`

- in: `data/out_ambientRNA/anndata`
- out: `data/out_ambientRNA/anndata/anndata.h5ad`


## 2. preprocessing/notebooks

### nb00_ambientRNA_contamination_plots.ipynb

- plots ambient RNA contamination per sample

### nb01_barcode_rank_plot.ipynb

- plots barcode ranks, coloured by ambient RNA contamination estimation

### nb02_first_qc.ipynb

- MAD-based lenient round of QC

### nb03_dedoublet.ipynb

- using scDblFinder for dedoubleting

### nb04_clustering_after_dedoublet.ipynb

- exploratory clustering after removing the doublets

### nb05_final_qc_celltypist.ipynb

- additional QC with mast cell & neutrophil mask

### nb06_mesenchymal.ipynb

- subsetting the mesenchymal compartment, clustering and fine annotation (JCI Ke2024-based)

## 3. analysis/notebooks

### nb01_gene_expression_signature.ipynb
<br>

# Adams2020_IPF

## 1. preprocessing



<br><br>

---

# Acknowledgements

**`MIL-TARGID/Italia_TLR4`** created by [Kamila Kwiecien](https://github.com/kamila-kwiecien)