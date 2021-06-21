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

## Details
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

## Related details
  * https://satijalab.org/seurat/articles/integration_introduction.html
  * https://satijalab.org/seurat/articles/merge_vignette.html
  * https://mojaveazure.github.io/seurat-disk/articles/h5Seurat-load.html
  * https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
  * https://satijalab.org/seurat/articles/integration_large_datasets.html
  * https://satijalab.org/seurat/archive/v3.0/conversion_vignette.html
  * https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/seurat-part-3-data-normalization/
  * https://satijalab.org/seurat/articles/get_started.html
