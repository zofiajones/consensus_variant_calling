import pandas as pd

df=pd.read_csv('multi.txt',sep='\t')
dfR=pd.read_csv('/home/zjones/exomes_1/RKO.vcf',sep='\t',header=None)
dfS=pd.read_csv('/home/zjones/exomes_1/SW48.vcf',sep='\t',header=None)
dfH=pd.read_csv('/home/zjones/exomes_1/HCT116.vcf',sep='\t',header=None)
dfGT=pd.read_csv('/home/zjones/exomes_1/GTB.vcf',sep='\t',header=None)
dfHD=pd.read_csv('/home/zjones/exomes_1/HD701.vcf',sep='\t',header=None)
dvar=pd.read_csv('req_var.txt',sep='\t')

variants=list(set(df['multi_var'].tolist()))
#print(dfS[5].tolist()[0])

print(variants)

file1=open('multi_freq_2.txt','w')

#multi_bases=[]
#multi_bases_found=[]
#multi_found=[]
for var in variants:
#	print(vr,'res0',df[df['name']==var]['res0'].tolist())
	multi_bases=str(df[df['multi_var']==var]['multi_bases'].tolist()[0])
	multi_bases_found=str(df[df['multi_var']==var]['multi_bases_found'].tolist()[0])
	multi_found=str(df[df['multi_var']==var]['multi_found'].tolist()[0])
	R=list(set(dfR[dfR[5]==var][4].tolist()))
	H=list(set(dfH[dfH[5]==var][4].tolist()))
	S=list(set(dfS[dfS[5]==var][4].tolist()))
	HD=list(set(dfHD[dfHD[5]==var][4].tolist()))
	GB=list(set(dfGT[dfGT[5]==var][4].tolist()))
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

	items=[H0,S0,R0,HD0,GB0]
	entry0=('\t').join(items)
	entry1=('\t').join(['NF','NF','NF','NF','NF','NF'])
	print(var + '\t' + entry0+'\t' +entry1)
#	print(multi_found,len(multi_found),type(multi_found))
	multi=multi_found.split(',')
	multi0=''
	for item in multi:
		multi0=multi0+','+item.replace('[','').replace(']','').replace('\'','').strip()
	print(multi0.replace('^,',''))
	file1.write(';' + var + ';' + '\t' + entry0+'\t' +entry1+'\t' + multi_bases + '\t' + multi_bases_found + '\t' + multi0 + '\n')
