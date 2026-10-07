module load apptainer
docker_imgs='/cluster/projects/mcgahalab/envs/containers/'
pdir='/cluster/projects/mcgahalab/data/mcgahalab/melissa_tcells/scrna/results/immcantation'
sample='D8_WT'

cd $pdir

####################
####  1. # Change-O tool: Annotate the raw V(D)J FASTA with IgBLAST#
# Preprocessing / V(D)J annotation
apptainer exec ${docker_imgs}/immcantation/immcantation_suite-4.4.0.sif \
AssignGenes.py igblast \
   -s data/${sample}/filtered_contig.fasta \
   -b /usr/local/share/igblast \
   --organism mouse --loci ig --format blast \
   --outdir results/${sample} --outname BCR_data_sequences

# convert IgBLAST output into AIRR-standard TSV
apptainer exec ${docker_imgs}/immcantation/immcantation_suite-4.4.0.sif \
MakeDb.py igblast \
   -i results/${sample}/BCR_data_sequences_igblast.fmt7 \
   -s data/${sample}/filtered_contig.fasta \
   -r /usr/local/share/germlines/imgt/mouse/vdj/ \
   --10x data/${sample}/filtered_contig_annotations.csv --extended

Rscript initializeBCR.R

threshold=$(awk -F',' '$1=="D8_WT" {print $2}' results/${sample}/threshold.csv)
# --act default 'set'; repertoire_profiler uses 'single'
# --dist default 0; repertoire_profiler uses '0.15'
apptainer exec ${docker_imgs}/immcantation/immcantation_suite-4.4.0.sif \
DefineClones.py \
  -d results/${sample}/BCR_data_sequences_igblast_db-pass.tsv \
  --act set \
  --model ham \
  --norm len \
  --dist ${threshold} \
  --outdir results/${sample} \
  --fail
