import re
import sys
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.strip().split('\t')
#	entry=('\t').join(m)
#	print(line.replace('\n','')+'\t'+'end')
	if 'O' in m[3]:
		m[3]=m[3].replace('O','')
		m[1]=str(int(m[1])-1)
		loc=m[0]+':'+str(m[1])+'-'+str(m[1])
		base=file.fetch(region=loc).upper()
		m[3]=base+m[3]
		m[4]=base+m[4]
		m.append('INS')
		entry=('\t').join(m)
		print(entry)
	elif 'O' in m[4]:
		m[4]=m[4].replace('O','')
		m[1]=str(int(m[1])-1)
		loc=m[0]+':'+str(m[1])+'-'+str(m[1])
		base=file.fetch(region=loc).upper()
		m[4]=base
		m[3]=base+m[3]
		m.append('DEL')
		entry=('\t').join(m)
		print(entry)
	elif len(m[3])>1:
		m.append('MNV')
		entry=('\t').join(m)
		print(entry)
	else:
		m.append('SNP')
		entry=('\t').join(m)
		print(entry)
