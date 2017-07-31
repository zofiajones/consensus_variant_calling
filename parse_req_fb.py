import re
import sys

for line in sys.stdin:
	m=line.split('\t')
	m[4]=m[4].replace('RKO@','')
	m[4]=m[4].replace('SW48@','')
	m[4]=m[4].replace('HCT116@','')
	index=m[0]+'@'+m[1]+'@'+m[2]+'@'+m[3]
	n=m[4].split(';')
	n[1]=n[0]
	n[0]=index
	m[4]=(';').join(n)
	col0=('\t').join(m)
	entry=index+'\t'+col0
	print(entry)
