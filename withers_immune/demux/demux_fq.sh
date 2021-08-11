#file='mDC-367_S47_L001_R1_001.fastq.gz'
file=$1
bbmapdir='/cluster/projects/mcgahalab/bin/bbmap'
lane=${file/*_L/L}
lane=${lane/_*/}

${bbmapdir}/demuxbyname.sh \
-Xmx20g \
in=${file} \
out=demux/mDC-367_%.${lane}.fq \
outu=demux/unmatched.fq \
stats=demux/demux.stats \
barcode=t
