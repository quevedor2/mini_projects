module load apptainer
docker_imgs='/cluster/projects/mcgahalab/envs/containers/'
pdir='/cluster/projects/mcgahalab/data/mcgahalab/melissa_tcells/scrna/results/immcantation'
cd $pdir

apptainer shell ${docker_imgs}/immcantation/immcantation_suite-4.4.0.sif
R

library(airr); library(alakazam); library(dowser)
library(dplyr); library(ggplot2); library(scoper); library(Seurat); library(shazam)
sample='D8_WT' # sample='CD19_HPB238'
dir=file.path("results", sample)
outdir=file.path(".", "results", sample)
annocol="immgen.main.cluster" # "monaco.main.cell"
meta_f=file.path("data", "annotated_meta.csv")
# meta_f='/cluster/projects/mcgahalab/data/mcgahalab/ipmn_pdac/results/seurat_obj/3_annotated_meta.csv'
dir.create(outdir, showWarnings=F)


# --- Load and clean BCR data ---
bcr_data <- airr::read_rearrangement(file=file.path(dir, "BCR_data_sequences_igblast_db-pass.tsv"))

# drop cells with inconsistent locus calls, non-productive seqs, multiple heavy chains, or no heavy chain at all
bcr_data <- bcr_data %>%
  filter((grepl("^IGHV", v_call) & grepl("^IGHJ", j_call) & grepl("^IGH[MGADE]", c_call)) |
         (grepl("^IGKV", v_call) & grepl("^IGKJ", j_call) & grepl("^IGKC", c_call)) |
         (grepl("^IGLV", v_call) & grepl("^IGLJ", j_call) & grepl("^IGLC", c_call))) %>%
  filter(productive)

multi_heavy <- table(filter(bcr_data, locus == "IGH")$cell_id)
bcr_data <- filter(bcr_data, !cell_id %in% names(multi_heavy)[multi_heavy > 1])

heavy_cells <- filter(bcr_data, locus == "IGH")$cell_id
bcr_data <- filter(bcr_data, cell_id %in% heavy_cells)

# --- Merge in GEX cell-type annotations ---
## GEX cell-type anotations typically produced using SCT normalized in
## Seurat, followed by Azimuth + GPTcelltype  + SingleR annotations
## Followed by subclustering and manual-review/annotations
gex_meta <- read.csv(meta_f)
cell_id_unique <- paste0(sample, "_", gsub("-1", "", bcr_data$cell_id))
bcr_data$cell_id_unique <- cell_id_unique
bcr_data$subject_id <- sample
match.index <- match(bcr_data$cell_id_unique, rownames(gex_meta))
bcr_data$gex_annotation <- as.character(gex_meta[,annocol])[match.index]
bcr_data <- filter(bcr_data, !is.na(gex_annotation))

# --- Find heavy-chain clonal clustering threshold ---
dist_nearest <- shazam::distToNearest(dplyr::filter(bcr_data, locus == "IGH", subject_id == sample),
                               cellIdColumn = "cell_id")
threshold <- shazam::findThreshold(dist_nearest$dist_nearest,
                                    method = "gmm", model = "gamma-norm",
                                    cutoff = "user", spc = 0.995)@threshold
# ---repertoire_profiler: Distance_hist()
ggplot2::ggplot(
    data = subset(db, !is.na(dist_nearest)),
    ggplot2::aes(x = dist_nearest)
) +
    ggplot2::geom_histogram(
        binwidth = 0.02
    ) +
    ggplot2::geom_vline(
        xintercept = clone_distance,
        linetype = 2
    )


# --- SCOPer Define clonal groups based on V, J, CDR3 ---
results <- hierarchicalClones(bcr_data, cell_id = "cell_id_unique",
                               threshold = threshold, only_heavy = TRUE, summarize_clones=FALSE,
                               split_light = FALSE, fields = "subject_id")
results.bkup = results
results <- results.bkup %>%
  dplyr::mutate(v_germline_length = as.integer(v_germline_end - v_germline_start + 1),
    d_germline_length = as.integer(d_germline_end - d_germline_start + 1),
    j_germline_length = as.integer(j_germline_end - j_germline_start + 1))
    
# --- Reconstruct germline sequences ---
references <- readIMGT(dir = "/usr/local/share/germlines/imgt/mouse/vdj")
results <- createGermlines(results, references, fields = "subject_id")

# --- Alakazam and ShazaM: Repertoir-level somatic hypermutation load per cell ---
data_mut <- shazam::observedMutations(filter(results, locus == "IGH"),
              sequenceColumn = "sequence_alignment",
              germlineColumn = "germline_alignment_d_mask",
              regionDefinition = IMGT_V, frequency = TRUE, combine = FALSE)
abund <- estimateAbundance(filter(results, locus == "IGH"), group = "subject_id")  # alakazam
div   <- alphaDiversity(filter(results, locus == "IGH"), group = "subject_id")    # alakazam
clone_sizes <- countClones(filter(results, locus == "IGH"), groups = "subject_id") # alakazam

saveRDS(results, file=file.path(outdir, "results.rds"))
saveRDS(bcr_data, file=file.path(outdir, "bcr_data.rds"))
outf <- file.path(outdir, paste0(sample, "_immcantation_shazam.pdf"))
pdf(outf)
plot(div)
plot(abund)
plot(clone_sizes)
dev.off()
file.copy(outf, to="~/xfer", overwrite=T)
message("xfer ", basename(outf))

outf <- file.path(outdir, paste0(sample, "_immcantation_mufreq.pdf"))
pdf(outf)
ggplot(data_mut, aes(y = mu_freq_cdr_s, x = gex_annotation, fill = gex_annotation)) +
  geom_boxplot() + labs(y = "Mutation frequency CDR S", x = "Cell type")
  ggplot(data_mut, aes(y = mu_freq_cdr_r, x = gex_annotation, fill = gex_annotation)) +
    geom_boxplot() + labs(y = "Mutation frequency CDR R", x = "Cell type")
dev.off()
file.copy(outf, to="~/xfer", overwrite=T)
message("xfer ", basename(outf))
