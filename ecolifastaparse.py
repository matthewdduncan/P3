from Bio import SeqIO
import gzip

output = open('parsedecolifasta.txt', 'w')

#parse the fasta file
for record in SeqIO.parse('ecoliprot.faa', 'fasta'):

	#take id and seq attributes, write them to a new file
	output.write(str(record.id) + '\t' + str(record.seq) + '\n')