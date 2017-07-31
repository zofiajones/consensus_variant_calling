import sys
import pandas as pd

df=pd.read_csv('final_in_uniq.txt',sep=',',header=None)
df1=pd.read_csv('got_genes.txt',sep='\t',header=None)

print(type(df),type(df1),type(len(df[0])))

print(len(df[0]))

loc=[]
for i in range(len(df[0].tolist())):
	info=df[1].tolist()[i]
	print(df[1].tolist()[i] , df[0].tolist()[i])
	info0=df[1].tolist()[i].split(';')
	#print(info0)
	loc=info.split(';')[1]
	print('loc',loc)
	print('gene',len( list(set(df1[df1[0]==loc][1].tolist())) ),df1[df1[0]==loc][1].tolist())
	gene=('+').join(list(set(df1[df1[0]==loc][1].tolist())))
	if len(info0)==3:
		info0[2]=gene
	else:
		info0.append(gene)
	new_info=(';').join(info0)
	df.set_value(i, 1, new_info)

df.to_csv('final_in_uniq_gene.txt',sep='\t',index=False,header=None)
