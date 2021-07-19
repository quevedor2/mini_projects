# PDAC AhR Repo
## Description
"Question- Given the wide availability of single cell sequencing 
datasets of human pancreatic cancer, the authors should address the 
cell type specific expression of AhR in human tumors."

"dig up some public scRNA seq data (preferably from high impact papers 
if possible). And look for AhR expression data for all cells and provide 
a cell type ranked assessment of relative AhR expression"

In this repo, I try and answer the question of cell-type specific expression of
AhR across PDAC tumor singlecell datasets, with a focus on immune. For the sake
of comparisons, I am using 10X datasets that are not immune depleted.

## Dataset Details
* **Marina Pasca di Magliano:** Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer
  * **Data Availability**: GSE155698
    * [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698)
    * [https://github.com/PascaDiMagliano-Lab/MultimodalMappingPDA-scRNASeq](https://github.com/PascaDiMagliano-Lab/MultimodalMappingPDA-scRNASeq)
    * **scRNAseq with clinical data:** phs002071
  * **Platform**: 10X Chromium, HiSeq4000 or NovaSeq6000, Cell Ranger v3.0.0, Seurat v3.0, filtered to include (cells>200 genes, genes>3cells)
  * **scRNA Samples**: 1
    * 16 treatment-naive PDA samples (6 surgical, 10 fine-needle biopsy)
      * 46,244 cells from PDA (55,873 from blood)
    * Controls: three non-malignant pancreas samples (1196, 1258, 19732; duodenal adenoma, ampullary carcinoma, tissue adjacent to PDA)
      * 8,541 cells from adjacent/normal samples (14,240 from blood (n=4))
  * **Sample Processing**: 
  * **Notes**: 
  * **Use**: TRUE

## Workflow
Reproducibility: The conda environment used for this analysis is provided in `envs/environment.yaml`. Packages that had to be installed outside of conda are listed in the `envs/r-env.txt` file.

The dataset directory is set up as followed, with the corresponding `GSE155698.h5seurat` used as the input:
```
datasets/
├── GSE155698
│   ├── AdjNorm_TISSUE_1
│   │   └── filtered_feature_bc_matrix
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   ├── AdjNorm_TISSUE_2
...
│   └── PDAC_TISSUE_9
│       └── filtered_feature_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
```

The cell type markers used for annotation were as followed:
```
list("NatureCancer_auto"=list(
    "Epithelial"=c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1"),
    "T-cell"=c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1"),
    "Myeloid"=c("CD14","ITGAM","MNDA","MPEG1","ITGAX"),
    "B-cell"=c("CD79A","MS4A1","CD19"),
    "Fibroblasts"=c("CDH11","PDGFRA","PDGFRB","ACTA2"),
    "RBC"=c("HBA1","HBB","HBA2"),
    "NK"=c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1"),
    "Endothelial"=c("CDH5","PECAM1"),
    "Acinar"=c("TRY4","SPINK1","AMY2A")
  ))
```
Other cell type lists were provided based on the NatureCancer paper code, as well as a related CellResearch paper.

The main steps of the analysis were executed using the `code/pdac-ahr_main.R` code. This code relied on the `AutomatedClusterMarkerTable()` function which is found at `code/AutomatedClusterMarkerTable.R`. The single-cell processing workflow and the cell-type annotation workflow were adapted directly from the original authors source code, found here: https://github.com/PascaDiMagliano-Lab/MultimodalMappingPDA-scRNASeq/blob/master/R%20Code/MultimodalMappingPDA-scRNASeq.R . 

Copy-number variation (CNV) inference was done using the inferCNV R package. The code used is found in `code/tn-matched_infercnv.R`, and each sample was individually submitted to the cluster to speed up the processing time. **NOTE:** This step was done between Step 5 and Step 6 of the main `pdac-ahr_main.R` code and the results were saved as an rds file to be loaded in Step 6.

Finally, the previous clustering from steps 1-5 and the CNV inference were combined together in Step 6 to estimate the percent genome altered (PGA) for individual cells, as well as create the visualizations needed to investigate the relation between AHR, the genes of interest, and PGA.


 
## Related details
  * https://satijalab.org/seurat/articles/integration_introduction.html
  * https://satijalab.org/seurat/articles/merge_vignette.html
  * https://mojaveazure.github.io/seurat-disk/articles/h5Seurat-load.html
  * https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
  * https://satijalab.org/seurat/articles/integration_large_datasets.html
  * https://satijalab.org/seurat/archive/v3.0/conversion_vignette.html
  * https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/seurat-part-3-data-normalization/
  * https://satijalab.org/seurat/articles/get_started.html

## Other Datasets
* **Wenming Wu:** Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma
  * **Data Availability**:  deposited in GSA:PRJCA001063, accession number GSA:CRA001160
    * [https://zenodo.org/record/3969339](https://zenodo.org/record/3969339)
    * `wget '[https://zenodo.org/record/3969339/files/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad?download=1](https://zenodo.org/record/3969339/files/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad?download=1)'`
  * **Platform**: 10X Chromium, HiSeqXTen (150nt PE), Cell Ranger 2.1.0, Seurat v2.3.0, removed low quality cells (<200genes/cell, <3cells/gene and  >10% mitochondrial genes)
  * **scRNA Samples**: 
    * 24 PDAC tumor samples
    * Control pancreas without visible inflammation: 11 (3 non-pancreatic tumor [bile duct tumors, duodenal tumors], and 8 non-malignant pancreatic tumors [pancreatic cyst]).
  * **Sample Processing**: 
  * **Notes**: 
  * **Use**: TRUE
