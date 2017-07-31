import sys
import os
from pybedtools import BedTool
import tabix
import re
import pysam
tb=tabix.open('/home/zjones/COSMIC/CosmicCoding_sort.bed.gz')
file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.strip().split()
	records=tb.query(m[0],int(m[1])-1 ,int(m[1]))
	index0=('@').join([m[0],m[1],m[3],m[4]])
	info=[]
	for record in records:
#		print(record)
		index=('@').join(record[-1].split('@')[0:4])
		#print(m[5],index)
		if index0 == index:
			found=1;
			all_ref=list(set(m[2]))
			all_alt=list(set(m[3]))
			if ( all_ref[0]=='O' )  & ( len(all_ref)==1  ):
				loc=m[0]+':'+str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)
				m[1]=str(int(m[1])-1)
				m[2]=base.strip()
				m[3]=base.strip()+m[3]
			elif ( all_alt[0]=='O' )  & ( len(all_alt)==1  ):
				loc=m[0]+':'+str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)
				m[1]=str(int(m[1])-1)
				m[3]=base.strip()
				m[2]=base.strip()+m[2]
			elif ('O' in all_ref ) & ('O' in all_alt )  :
				print('else')
			fields=record[-1].split('@')[-1].split(';')
			#print(fields)
			for field in fields:
				if 'CDS' in field:
					cds=field.replace('CDS=','')
				if 'AA' in field:
					aa=field.replace('AA=','')
				if 'STRAND' in field:
					strand=field.replace('STRAND=','')
				if 'COSM' in field:
					cosm=field
			info.append((',').join([cds,aa,strand,cosm]))
	entry=line.strip()+'\t'+';'+(';').join(info)
	print(entry)
