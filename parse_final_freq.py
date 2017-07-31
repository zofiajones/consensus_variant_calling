import re
import sys

for line in sys.stdin:
#	print(line)
	m=line.split('\t')
#	print(len(m))
	if len(m)!=29:
		#continue
		print('len',len(m))
		continue
	info=m[4]
	variant=info.split(';')[0]
	locs=variant.split('@')
	if len(locs)==4:
#		print(line)
		chr0=locs[0]
		pos0=locs[1]
		ref0=locs[2]
		alt0=locs[3]
	else:
		print('loc',len(locs))
		continue
	chr=m[0]
	pos=m[1]
	ref=m[3]
	alt=m[4]
	RKO_freq=float(m[5].replace('NA','0.0'))
	HCT116_freq=float(m[9].replace('NA','0.0'))
	SW48_freq=float(m[13].replace('NA','0.0'))
	HD701_freq=float(m[17].replace('NA','0.0'))
	#print(str(m[23]),'GTB')
	GTB_freq=float(m[23].replace('NA','0.0'))
	RKO_alt=m[7]
	HCT116_alt=m[11]
	SW48_alt=m[15]
	HD701_alt=m[19]
	GTB_alt=m[28]
	RKO_cov=float(m[6].replace('NA','0.0'))
	HCT116_cov=float(m[10].replace('NA','0.0'))
	SW48_cov=float(m[14].replace('NA','0.0'))
	HD701_cov=float(m[18].replace('NA','0.0'))
	GTB_cov=float(m[24].replace('NA','0.0'))
	alts=[RKO_alt,HCT116_alt,SW48_alt,HD701_alt]
	freqs_germ=[RKO_freq,HCT116_freq,SW48_freq]
	freqs_som=[HD701_freq,GTB_freq]
	freqmaxs=max(freqs_som)
	freqmaxg=max(freqs_germ)
	if freqmaxs < 0.01:
		print(freqmaxs)
		print('freqmaxs',freqmaxs,freqs_som)
		continue
#	if freqmaxg == 0.0:
		#print('freqmaxg',freqmaxg)
#		continue
	while 'NA' in alts:
		#print('alts',alts)
		alts.remove('NA')
	while '0' in alts:
		#print('alts',alts)
		alts.remove('0')
	if not alts:
		print('no alts',alts)
		continue
	alt1=list(set(alts))
#	print('alts',alts)
	if len(alt1)>1 or '0' in alt1:
		print('many alts',alts)
		continue
	else:
		print(line)
		z=1
