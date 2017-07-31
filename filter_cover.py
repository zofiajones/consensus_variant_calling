import re
import os
import sys

for line in sys.stdin:
	m=line.split('\t')
	if m[20]!='NS':
		if float(m[20]) < 0:
			m[20]='NS'
	elif m[25]!='NS':
		if float(m[25]) < 0:
			m[25]='NS'
	print(line.strip())
