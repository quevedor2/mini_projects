#!/bin/bash

#SBATCH -t 8:00:00
#SBATCH --mem=30G
#SBATCH -J blastit
#SBATCH -p all
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o %x-%j.out

PDIR='/cluster/projects/mcgahalab/data/mcgahalab/wt_cutandrun/snakemake/results/alignment/blast'
#file='HPB67_Mono_Me3'
file=$1



cd ${PDIR}"/data"
mkdir map_unmap/${file}/

# Extract unmapped reads
module load samtools
samtools view -u  -f 4 -F 264 ${file}.bam  > ${file}.tmps1.bam
samtools view -u -f 8 -F 260 ${file}.bam  > ${file}.tmps2.bam
samtools view -u -f 12 -F 256 ${file}.bam > ${file}.tmps3.bam
samtools merge -u - ${file}.tmps[123].bam | samtools sort -n -o map_unmap/${file}/unmapped.bam -
rm ${file}.tmps1.bam ${file}.tmps2.bam ${file}.tmps3.bam

samtools view -b -F 4 ${file}.bam | samtools sort -n -o map_unmap/${file}/mapped.bam -





# Visualize the PHRED quality scores
cd ${PDIR}
mkdir -p qc/${file}/

module load picard/2.10.9
java -jar $picard_dir/picard.jar QualityScoreDistribution \
      I=data/map_unmap/${file}/unmapped.bam \
      O=qc/${file}/unmapped.qual_score_dist.txt \
      CHART=qc/${file}/unmapped.qual_score_dist.pdf

java -jar $picard_dir/picard.jar QualityScoreDistribution \
      I=data/map_unmap/${file}/mapped.bam \
      O=qc/${file}/mapped.qual_score_dist.txt \
      CHART=qc/${file}/mapped.qual_score_dist.pdf




cd ${PDIR}
mkdir -p blast/${file}/

module load samtools
module load blast+/2.2.30
# -s subsample 5% of the reads for blasting
samtools view -s 0.001 data/map_unmap/${file}/unmapped.bam | \
awk '{printf(">%s/%s\n%s\n",$1,(and(int($2),0x40)?1:2),$10);}' | \
blastn -db nt -out blast/${file}/unmapped_blast.xml -outfmt 5

# -s subsample 5% of the reads for blasting
samtools view -s 0.001 data/map_unmap/${file}/mapped.bam | \
awk '{printf(">%s/%s\n%s\n",$1,(and(int($2),0x40)?1:2),$10);}' | \
blastn -db nt -out blast/${file}/mapped_blast.xml -outfmt 5


# parse blast results
cd ${PDIR}/blast/${file}
python ../parseBlast.py
