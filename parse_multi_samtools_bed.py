import re
import os
import sys
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.split('\t')
	alt=m[4].split(',')
	tags=m[8].split(':')
	ind=tags.index('RD')
	RD=m[9].split(':')[ind]
	ind=tags.index('AD')
	AD=m[9].split(':')[ind].split(',')
	ind1=tags.index('DP')
	DP=m[9].split(':')[ind1]
	total_reads=0.0
	freq=[]
	for i in range(len(AD)):
		if int(DP) > 0:
			freq.append( float(AD[i])/float(DP)  )
#	if int(DP) < 100:
#		continue
	j=freq.index(max(freq))
#	if max(freq)<0.05:
#		continue
	info= str(freq[j]) + ',' + DP
#	print(info)
	alt1=alt[j].replace('.','O')
	lref=len(m[3])
	lalt=len(alt1)
#	print(m[4],alt1)
#	print(m)
	trim=0
	mx=range(0)
	if lalt < lref:
		mx=range(lalt);
	elif lref < lalt:
		mx=range(lref);
#                       print(range(lalt))
	for i in mx:
#                               print(i,m[4][lref-i-1],alt1[lalt-i-1])
		if m[3][lref-i-1] == alt1[lalt-i-1]:
			loc=m[0]+':'+ str(int(m[1])+lref-1)+'-'+str(int(m[1])+lref-1)
			base=file.fetch(region=loc)
#                                       print(loc,base)
                                        #print(type(base),type(m[4][lref-i-1]))
			if base.strip() in m[3][lref-i-1]:
				trim=trim+1
                                                #print(trim,str(lref-i-1))
			else:
				break
		else:
			break



	ref=m[3][0:lref-trim]
	alt=alt1[0:lalt-trim]
#	print('trim',trim)
	while len(alt) < len(ref):
		alt=alt+'O'
	while len(ref) < len(alt):
		ref=ref+'O'

	lalt=len(alt)
	lref=len(ref)
	alt1=''
	ref1=''
	i=0
#	print(alt,ref,'here1')
	while alt[i] == ref[i]:
		i=i+1
	alt1=alt[i:]
	ref1=ref[i:]
	m[1]=str(int(m[1])+i)


	ref=ref1;
	alt=alt1;
#	print(ref,alt,m[1],m[2])
	lalt=len(alt)
	lref=len(ref)

	if alt=='O' and lalt ==1:
		loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
		base=file.fetch(region=loc)
		while base.strip() == ref:
			m[1]=str(int(m[1])-1)
#			print(m[2],m[1])
			loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
			base=file.fetch(region=loc)

#	print(ref,alt,lref)
	if ref=='O' and lref == 1:
		loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
		base=file.fetch(region=loc)
#		print(loc,base,'hello')
		while base.strip() == alt:
			m[1]=str(int(m[1])-1)
#			print(loc,base)
			loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
			base=file.fetch(region=loc)


	for i in range(len(ref)):
		#print(i)
		altf=''
		reff=''
		if alt[i] != ref[i]:
			altf=altf+alt[i]
			reff=reff+ref[i]
	info0=('@').join([ m[0],str(int(m[1])),reff,altf  ])
	print(('\t').join([m[0],str(int(m[1])),reff,altf,info.strip(),info0]))

