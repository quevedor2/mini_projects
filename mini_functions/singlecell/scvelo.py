# Workaround for ImportError: cannot import name 'Iterable' from 'collections' (/cluster/home/quever/miniconda3/envs/scvelo/lib/python3.10/collections/__init__.py)
# https://stackoverflow.com/questions/72032032/importerror-cannot-import-name-iterable-from-collections-in-python
import collections.abc
#hyper needs the four following aliases to be done manually.
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
collections.MutableSet = collections.abc.MutableSet
collections.MutableMapping = collections.abc.MutableMapping

import os
import numpy as np
import pandas as pd
import scvelo as scv
import scanpy as sc
import cellrank as cr
from pathlib import Path
import matplotlib.backends.backend_pdf as backend_pdf
from matplotlib import pyplot as plt
import matplotlib as mp



####################
#### Parameters ####
# Setup:
outdir = Path("xfer/")
indir = Path("/cluster/projects/mcgahalab/data/mcgahalab/teresa_ly6/scrna/results/scVelo/")
dmsodata = scv.read(indir / "seu_velo.DMSO.h5ad")
dabdata = scv.read(indir / "seu_velo.DAB.h5ad")
mergeddata = scv.read(indir / "seu_velo.Merged.h5ad")
datadict = {
    'dmso' : dmsodata,
    'dab' : dabdata,
    'merged' : mergeddata
}
    
#adata = scv.read("seu_velo.h5ad")
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')
cr.settings.verbosity = 2
label_id = 'functional_cluster'
basis_id = 'umap_orig' # umap umap_orig

s_genes_list = \
    ['Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
     'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Mlf1ip', 'Hells', 'Rfc2',
     'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
     'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
     'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8']

g2m_genes_list = \
    ['Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
     'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
     'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
     'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
     'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
     'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
     'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa']

retain_genes = \
    ['Ly6c2', 'Ly6g', 'Ly6c1', 'Itgam', 'H2-Ab1', 'H2-Aa', 'H2-Eb1', 'H2-Eb2',
    'H2-Ea', 'Csf1r', 'Il4ra', 'Ccr2', 'Arg1', 'Ahr']
all_retained = s_genes_list + g2m_genes_list + retain_genes

########################################
#### Preprocessing and loading data ####
# Preprocess the data
term = 'dab'
processed_f = "seu_velo." + term + ".processed.h5ad"
processed_f_path = indir / processed_f
adata = datadict[term]
if os.path.exists(processed_f_path):
    print("loading...")
    adata = scv.read(processed_f_path)
else:
    print("processing...")
if True:
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000,
        retain_genes=all_retained)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.recover_dynamics(adata) # learn full transcriptional dynamics of splicing kinetics
    scv.tl.velocity(adata, mode='dynamical')
    # var_names = ['Ly6c2', 'Ly6g', 'Itgam', 'H2-Ab1']
    #scv.tl.differential_kinetic_test(adata, groupby='seurat_clusters') # var_names=var_names,
    #scv.tl.velocity(adata, diff_kinetics=True)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)
    scv.tl.latent_time(adata)
    
    # Saving error _index fix: https://github.com/theislab/scvelo/issues/255
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    # del adata.raw # dotplot error: f"Could not find keys '{not_found}' in columns of `adata.{dim}` or in"
    # del(adata.var['_index']) #if still error while saving
    adata.write(filename=processed_f_path, compression='gzip')

############
#### QC ####
# Kinetic Rate parameters
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
scv.get_df(adata, 'fit*', dropna=True).head()

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

outf = outdir / (term + ".kinetics.pdf")
plt.savefig(outf, format="pdf", bbox_inches="tight")
plt.show()

# Proportions of spliced to unspliced
scv.pl.proportions(adata, groupby=label_id)
outf = outdir / (term + ".proportions.pdf")
plt.savefig(outf, format="pdf", bbox_inches="tight")
plt.show()

###########################
#### CellRank Analysis ####
# Fix to not use OpenMPI: https://github.com/theislab/cellrank/issues/661
cr.tl.terminal_states(adata, cluster_key="seurat_clusters", weight_connectivities=0.2, fit_kwargs={"method": "brandts"})

# Computer fate maps: how likely the cell is to develop towards terminal states
# Fix to not use OpenMPI: https://github.com/theislab/cellrank/issues/661
cr.tl.lineages(adata, use_petsc=False)

cr.pl.terminal_states(adata)
outf = outdir / (term + ".crank.terminal.pdf")
plt.savefig(outf, format="png", bbox_inches="tight")
plt.show()
print("xfer " + term + ".crank.terminal.pdf")

###########################
#### Velocity Analysis ####
#---------- Overall velocity and latent-time analysis
# Velocity_1: Streamline velocity projection
scv.pl.velocity_embedding_stream(adata, basis=basis_id, color=label_id)
outf = outdir / (term + ".scvelo.vel_stream.png")
plt.savefig(outf, format="png", bbox_inches="tight")
plt.show()
print("xfer " + term + ".scvelo.topheatmap.png")

# Velocity_2: Cellular and Gridlines level velocity projection
outf = outdir / (term + ".scvelo.vel_projection.pdf")
pdf = backend_pdf.PdfPages(outf)
# ---- Gridline velocity
fig = scv.pl.velocity_embedding_grid(adata, basis=basis_id, color=label_id,
    arrow_length=3, arrow_size=2, dpi=120)
pdf.savefig( fig )
# ---- Cellular velocity
fig = scv.pl.velocity_embedding(adata, basis=basis_id, color=label_id,
    arrow_length=3, arrow_size=2, dpi=120)
pdf.savefig( fig )
# ---- Velocity confidence
keys = 'velocity_length', 'velocity_confidence'
fig = scv.pl.scatter(adata, basis=basis_id, c=keys, cmap='coolwarm', perc=[5, 95])
pdf.savefig( fig )
# ---- Splicing efficiency: fraction of unspliced counts [https://github.com/theislab/scvelo/issues/173]
counts_s = scv.utils.sum_var(adata.layers['spliced'])
counts_u = scv.utils.sum_var(adata.layers['unspliced'])
fractions_u = counts_u / (counts_s + counts_u)
fig = scv.pl.scatter(adata, color=fractions_u, basis=basis_id, smooth=True)
pdf.savefig( fig )
# ---- Pseudotime projections
fig = scv.pl.scatter(adata,  basis=basis_id , color="latent_time", color_map="gnuplot")
pdf.savefig( fig )
pdf.close()


#---------- Top and selected velocity driver genes
# ---- Top 100 genes for velocity drivers
all_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index
top_genes = all_genes[:100]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color=label_id,
    n_convolve=100, font_scale=0.5, figsize=(8,15))
outf = outdir / (term + ".scvelo.topheatmap.png")
plt.savefig(outf, format="png", bbox_inches="tight")
plt.show()

outf = outdir / (term + ".scvelo.gene_velocities.pdf")
pdf = backend_pdf.PdfPages(outf)
# ---- Scatterplot of top genes velocities
fig = scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, color=label_id, frameon=False)
pdf.savefig( fig )
# ---- Scatterplot of selected genes gene_velocities
fig = scv.pl.scatter(adata, ['Ahr', 'Ly6c2', 'Ly6g'], ncols=2, color=label_id)
pdf.savefig( fig )
#fig = scv.pl.scatter(adata, x='latent_time', y=['Ahr', 'Ly6c2', 'Ly6g'], color=label_id)
#pdf.savefig( fig )
# ---- Cluster-specific top-likelihood genes
scv.tl.rank_dynamical_genes(adata, groupby=label_id)
df = scv.get_df(adata, 'rank_dynamical_genes/names')
#df.head(5)
for cluster in ['Neutrophil_precursor', 'Neutrophil_mature', 'M0_M2_macrophage', 'Mature_myeloid', 'LSK_HSC_MPP']:
    fig = scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, color=label_id, frameon=False)
    pdf.savefig( fig )

pdf.close()
