cat mDC-367_L001/demux/mDC-367_GCGTAGTA+CGTCTAAT.L001.fq \
mDC-367_L002/demux/mDC-367_GCGTAGTA+CGTCTAAT.L002.fq > merge/unk1-367_merged_R1.fastq
gzip merge/unk1-367_merged_R1.fastq

cat mDC-367_L001/demux/mDC-367_GCGTAGTA+CTAAGCCT.L001.fq \
mDC-367_L002/demux/mDC-367_GCGTAGTA+CTAAGCCT.L002.fq > merge/unk2-367_merged_R1.fastq
gzip merge/unk2-367_merged_R1.fastq
