module load seqtk/git

id='Heni_28M_CTCF100'
target=$1
seed=1234
simpletarget=$(echo ${target} | sed "s/000000$/M/")

# Simplify the number to M or K
if [ "$simpletarget" = "$target" ]; then
  simpletarget=$(echo ${target} | sed "s/000$/K/")
fi

# Calculate the downsampling fraction
n=$(wc -l < ${id}.1.fastq.gz)
calc(){ awk "BEGIN { print "$*" }"; }
frac=$(calc ${target}/${n})
frac=$(printf "%.3g" "$frac")

# Downsample the fastqs
echo "Target: ${target}, SimpleTarget: ${simpletarget}, Seed=${seed}"
out_id=$(echo $id | sed "s/_[0-9]*M_/_${simpletarget}_/")
seqtk sample -s${seed} ${id}.1.fastq.gz ${frac} > ${out_id}.1.fastq
seqtk sample -s${seed} ${id}.2.fastq.gz ${frac} > ${out_id}.2.fastq
gzip ${out_id}.1.fastq
gzip ${out_id}.2.fastq
