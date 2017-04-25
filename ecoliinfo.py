from Bio import SeqIO
import gzip
import re

output = open('parsedecoli.txt', 'w')

taxid = ''

#genbank files canhandle multiple datasets, this iterates through them (though there is only one here, for e coli k12)
for record in SeqIO.parse('ecoli.gbff', 'gb'):
	
	#each dataset has multiple features
	for feature in record.features:
	
		#gets the taxid for e coli
		if feature.type == 'source':
			taxid = feature.qualifiers['db_xref'][0]
			
		#get information about each gene
		if feature.type == 'CDS':
			line= {}  #line is a dict where all information is stored.  this is necessary to beable to reorder everything.
			
			#can make a loop here instead of rewriting.  Extracts all relevant features from cds feature.qualifiers
			for x in ['protein_id', 'gene', 'locus_tag', 'gene_synonym', 'product', 'EC_number', 'db_xref']:
				first = True
				line[x] = ''
				
				#this loop adds multiple entries in cases where a feature has them.
				for ft in feature.qualifiers.get(x, '-'):
					if first:
						first = False
					else:
						line[x] += ','
					line[x] += ft
					
			#get features not in qualifiers
			line['coordinates'] = re.search('\[.*\]', str(feature.location)).group(0)
			line['strand'] = re.search('\(.*\)', str(feature.location)).group(0)
			line['taxid'] = taxid
			
			#print to output in the correct order, with tab spacing
			for x in ['protein_id', 'coordinates', 'strand', 'gene', 'locus_tag', 'gene_synonym', 'product', 'taxid', 'EC_number', 'db_xref']:
				output.write(line[x] + '\t')
				
			#start next line
			output.write('\n')