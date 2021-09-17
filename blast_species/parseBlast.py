from Bio.Blast import NCBIXML

result=open("unmapped_blast.xml","r")
records= NCBIXML.parse(result)
file = open("unmapped.blast_results.txt", "w")
for item in records:
  for alignment in item.alignments:
    file.write(">")
    for hsp in alignment.hsps:
      title=alignment.title
      title=title.replace(" ", "_")
      title=title.replace(",", "-")
      file.write(title + "," + str(hsp.score) + "," + str(hsp.expect) + "\n")
      break
    break

file.close()


result=open("mapped_blast.xml","r")
records= NCBIXML.parse(result)
file = open("mapped.blast_results.txt", "w")
for item in records:
  for alignment in item.alignments:
    file.write(">")
    for hsp in alignment.hsps:
      title=alignment.title
      title=title.replace(" ", "_")
      title=title.replace(",", "-")
      file.write(title + "," + str(hsp.score) + "," + str(hsp.expect) + "\n")
      break
    break

file.close()
