import pandas as pd

df=pd.read_csv('multi_freq.txt',sep='\t',header=None)
dfR=pd.read_csv('/home/zjones/exomes_1/RKO.vcf',sep='\t',header=None)
dfS=pd.read_csv('/home/zjones/exomes_1/SW48.vcf',sep='\t',header=None)
dfH=pd.read_csv('/home/zjones/exomes_1/HCT116.vcf',sep='\t',header=None)
dfGT=pd.read_csv('/home/zjones/exomes_1/GTB.vcf',sep='\t',header=None)
dfHD=pd.read_csv('/home/zjones/exomes_1/HD701.vcf',sep='\t',header=None)

variants00=list(set(df[0].tolist()))
#print(variants00,len(variants00))

variants0=[]
for var in variants00:
	#print('var',var,type(var))
	var0=var.split(';')[1]
	variants0.append(var0)

variants=list(set(variants0))
#print(variants)

file1=open('multi_freq_1.txt','w')

for vr in variants:
	#print(vr)
	R=list(set(dfR[dfR[5]==vr][4].tolist()))
	H=list(set(dfH[dfH[5]==vr][4].tolist()))
	S=list(set(dfS[dfS[5]==vr][4].tolist()))
	HD=list(set(dfHD[dfHD[5]==vr][4].tolist()))
	GB=list(set(dfGT[dfGT[5]==vr][4].tolist()))
	#print(R,H,S,HD,GB)
	if not R:
		R0='NS'
	else:
		R0=R[0].split(',')[0]

	if not H:
		H0='NS'
	else:
		H0=H[0].split(',')[0]
	if not S:
		S0='NS'
	else:
		S0=S[0].split(',')[0]
	if not HD:
		HD0='NS'
	else:
		HD0=HD[0].split(',')[0]
	if not GB:
		GB0='NS'
	else:
#		print(GB)
		GB0=GB[0].split(',')[1]

	items=[HD0,GB0,R0,H0,S0]
	#print(items)
	#item1=list(set(items))
	#while 'NS' in item1:
	#	item1.remove('NS')
	#if not item1:
	#	continue
	entry0=('\t').join(items)
	entry1=('\t').join(['NF','NF','NF','NF','NF','NF'])
	print(';' + vr + ';' + '\t' + entry0+'\t' +entry1)
	#file1.write(vr + '\t' + entry0+'\t' +entry1+'\t' + '\n')
