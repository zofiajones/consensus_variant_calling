import re
import sys

for line in sys.stdin:
	#print(line)
	m=line.split('\t')
	info=m[1].split(';')
	#print(info,'info')
	n=info[1].split('@')
	#print(n,'n')
	if len(n)==4:
		chr=n[0]
		pos=n[1]
		ref=n[2]
		alt=n[3]
		id=m[1].split(';')[2]
		freq=m[9]
		error=m[5]
	else:
		continue
	print(('\t').join([chr,pos,id,ref,alt,'.',freq,error]))
