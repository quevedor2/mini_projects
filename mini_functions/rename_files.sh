while read i; do
  id1=$(echo $i | cut -d$' ' -f1)
  id2=$(echo $i | cut -d$' ' -f2)
  cp -R ${id1} ../${id2}
done <mapping.tsv
