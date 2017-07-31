import re
import sys
from pybedtools import BedTool

file=open('combined.bed','w')

for line in sys.stdin:
	m=line.split(',')
	print(m)
	n=m[1].split(';')
	p=n[1].split('@')
	if len(p)>3:
		entry=p[0]+'\t'+str(int(p[1])-1)+'\t'+p[1]+'\t' + p[2] + '@' + p[3]  + '\n'
		file.write(entry)
	else:
		print(m)
