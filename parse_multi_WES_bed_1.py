import re
import os
import sys
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.strip().split()
	if len(m)>5:
#		print('len',len(m))
#		print(line.strip())
		alt=m[4].split(',')
		n=m[7].split(';')
		for item in n:
			if 'DP=' in item:
				DP=item.replace('DP=','')
			if 'DP4=' in item:
				DP4=item.replace('DP4=','').split(',')
		if float(DP) == 0.0:
#			print('here')
			continue
#		print(DP,DP4)
		freq=[ ( float(DP4[2])+float(DP4[3]) ) /( float(DP4[0])+float(DP4[1])+float(DP4[2])+float(DP4[3])  )  ]
#		if int(DP)<100:
#			continue
#		if max(freq)<0.05:
#			continue
		j=freq.index(max(freq))
		info= str(freq[j]) + ',' + DP
		alt1=alt[j].replace('.','O')
		lref=len(m[3])
		lalt=len(alt1)
		trim=0
		mx=range(0)
		if lalt < lref:
			mx=range(lalt);
		elif lref < lalt:
			mx=range(lref);
		for i in mx:
			if m[3][lref-i-1] == alt1[lalt-i-1]:
				loc=m[0]+':'+ str(int(m[1])+lref-1)+'-'+str(int(m[1])+lref-1)
				base=file.fetch(region=loc).upper()
				#print( m[4],lref  )
				if base.strip() in m[3][lref-i-1].upper():
					trim=trim+1
				else:
					break
			else:
				break



		ref=m[3][0:lref-trim]
		alt=alt1[0:lalt-trim]
		while len(alt) < len(ref):
			alt=alt+'O'
		while len(ref) < len(alt):
			ref=ref+'O'

		lalt=len(alt)
		lref=len(ref)
		alt1=''
		ref1=''
		i=0
		while alt[i] == ref[i]:
			i=i+1
		alt1=alt[i:]
		ref1=ref[i:]
		m[1]=str(int(m[1])+i)


		ref=ref1;
		alt=alt1;
		lalt=len(alt)
		lref=len(ref)

		if alt=='O' and lalt ==1:
			loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
			base=file.fetch(region=loc)
			while base.strip() == ref:
				m[1]=str(int(m[1])-1)
				loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)

		if ref=='O' and lref == 1:
			loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
			base=file.fetch(region=loc)
			while base.strip() == alt:
				m[1]=str(int(m[1])-1)
				loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)


		ref0=''
		alt0=''
		for i in range(len(ref)):
			if alt[i] != ref[i]:
				alt0=alt0+alt[i]
				ref0=ref0+ref[i]
		info0=('@').join([m[0],str(int(m[1])),ref0.upper(),alt0.upper()])
		print(('\t').join([m[0],str(int(m[1])),ref0.upper(),alt0.upper(),info.strip(),info0]))
