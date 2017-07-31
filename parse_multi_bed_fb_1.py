import re
import os
import sys
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.split('\t')
#	print(m)
	alt=m[4].split(',')
	#print(alt)
	m2=re.search('AB=([\d\.,]+)',m[7]).group(1).split(',')
	m3=re.search('SAF=([,\d]+)',m[7]).group(1).split(',')
	m4=re.search('SAR=([,\d]+)',m[7]).group(1).split(',')
	m5=re.search('SRF=([,\d]+)',m[7]).group(1)
	m6=re.search('SRR=([,\d]+)',m[7]).group(1)
	tags=m[8].split(':')
	if len(m[9].split(':'))>1:
		ind=tags.index('AO')
#		print(tags,ind,m[8].split(':'),m[8])
		AO=m[9].split(':')[ind].split(',')
		ind1=tags.index('RO')
		RO=m[9].split(':')[ind1]
#		print(alt,m2,m3,m4,m5,m6)
#		print(line)
		freq0=[]
		total_reads0=float(m5)+float(m6)
		for i in range(len(alt)):
			total_reads0=float(m3[i])+float(m4[i])+total_reads0
		for i in range(len(alt)):
			if total_reads0 > 0:
				freq0.append(str( ( float(m3[i])+float(m4[i]) ) / total_reads0  ))
			else:
				freq0.append('0')
#		print(freq0,total_reads0)
		freq=[]
		total_reads=float(RO)
		for i in range(len(alt)):
			total_reads=total_reads+float(AO[i])
		for i in range(len(alt)):
			if total_reads > 0:
				freq.append( float(AO[i])/ total_reads  )
			else:
				freq.append(0.0)
		j=freq.index(max(freq))
		if total_reads < 100:
			continue
		if max(freq) < 0.05:
			continue
#		print(('\t').join([m[4],alt[j],m2[j],freq[j]]))
		info=m2[j] + ',' + str(freq[j]) + ',' + str(total_reads) + ',' + str(freq0[j]) + ',' + str(total_reads0)
	#print(info)
		alt1=alt[j].replace('.','O')
		lref=len(m[3])
		lalt=len(alt1)
#		print(m[4],alt1)
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
		lalt=len(alt)
		lref=len(ref)

		while len(alt) < len(ref):
			alt=alt+'O'
		while len(ref) < len(alt):
			ref=ref+'O'



		alt1=''
		ref1=''
		i=0
#		print(alt,ref)
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


		#print(trim,m[4],alt1,'o',ref,alt)
		while len(alt) < len(ref):
			alt=alt+'O'
		ref0=''
		alt0=''
		for i in range(len(ref)):
			#print(i)
			if alt[i] != ref[i]:
				alt0=alt0+alt[i]
				ref0=ref0+ref[i]
		info0=('@').join([m[0],str(int(m[1])),ref0,alt0])
		print(('\t').join([m[0],str(int(m[1])),ref0,alt0,info.strip(),info0  ]))

