import sys
import os
import re
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

for line in sys.stdin:
	m=line.split()
#	print(m)
	if not '#' in line :
#		print(line.strip())
		#print(('\t').join([m[0],m[1],m[2],m[3],m[4],m[5],m[6].strip()  ]))
		alt2=m[4].replace('.','O').split(',')
		m[3]=m[3].replace('.','O')
		lref=len(m[3])
#		info0=m[6].strip().split(';')
#		loc0=info0[0].split('@')
#		alt_base=loc0[3].split(',')
		alt_count=0
		for alt1 in alt2:
#			info=m[6] + ';' + ('@').join(loc0[0:3]) + '@' + alt_base[alt_count]
#			print('info',alt_base[alt_count],('@').join(loc0[0:3]),loc0[0:3])
			alt_count=alt_count+1
			lalt=len(alt1)
			lref=len(m[3])
			trim=0
			mx=range(0)
#			print(lalt,lref)
			if lalt < lref:
				mx=range(lalt);
				to_trim=m[3]
			elif lref < lalt:
				mx=range(lref);
				to_trim=alt1
#                       print(range(lalt))
#			print(lref,lalt,len(m[3]),len(alt1))
#			print(mx)
			for i in mx:
#                       	        print(i,m[4][lref-i-1],alt1[lalt-i-1])
#				print(i,lref,lalt,m[3],alt1)
				if m[3][lref-i-1] == alt1[lalt-i-1]:
					loc="chr"+m[0]+':'+ str(int(m[1])+lref-1)+'-'+str(int(m[1])+lref-1)
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
			while len(alt) < len(ref):
				alt=alt+'O'
			while len(ref) < len(alt):
				ref=ref+'O'

			lalt=len(alt)
			lref=len(ref)
			if alt=='O' and lalt == 1:
				loc="chr"+m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)
				while base.strip() == ref:
					m[1]=str(int(m[1])-1)
					loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
					base=file.fetch(region=loc)
			if ref=='O' and lref == 1:
				loc="chr"+m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
				base=file.fetch(region=loc)
				while base.strip() == alt:
					m[1]=str(int(m[1])-1)
					loc=m[0]+':'+ str(int(m[1])-1)+'-'+str(int(m[1])-1)
					base=file.fetch(region=loc)

			altf=''
			reff=''
			for i in range(len(ref)):
				if alt[i] != ref[i]:
					altf=altf+alt[i]
					reff=reff+ref[i]
			variant= m[0] + '@' + m[1] + '@' + reff + '@' + altf
			for i in range(len(reff)):
				print(('\t').join([m[0],str(int(m[1])+i),ref[i],alt[i],m[2].replace(',','@') ]))

