import re
import os
import sys
import pandas as pd
import pysam

file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')
df=pd.read_csv('variants_4.vcf',sep='\t',header=None)
df1=pd.read_csv('req_var.txt',sep='\t',header=None)

#df3=df[[3,4]]

#print(df.columns)
#print(df[6][:10])
#for line in sys.stdin:
#	m=line.split('\t')

found_store=[]
for i in range(len(df[0].tolist())):
	item=df[1][i]
	var=df[0][i]
	#print('var',var)
	bases=len(var.split('@')[2])
	m=item.split(',')[1:]
	find_lines=[]
	found=df1[df1[0]==item][0].tolist()
	if found:
		found_store.append(1)
	else:
		found_store.append(0)

df['found']=found_store



df2=df[df['found']==1]
print(df2.columns)
print(len(df2[0].tolist()),len(df2[1].tolist()),len(df2[2].tolist()),len(df2['found'].tolist()))
found_var=[]
found_var_name=[]
#print(df2[2].tolist()[4])
multi=dict()
multi['multi_var']=[]
multi['multi_found']=[]
multi['multi_bases']=[]
multi['multi_bases_found']=[]
for i in range(len(df2[0].tolist())):
	to_find=df2[2].tolist()[i]
	print(df2[0].tolist()[i])
	if to_find ==1:
		found_var.append(1)
		found_var_name.append(df2[0].tolist()[i])
		print('here',df2[1].tolist()[i],df2[0].tolist()[i])
	else:
		vars=df2[df2[0]==df2[0].tolist()[i]].index.tolist()
		#print('vars',vars,to_find,df2[0].tolist()[i])
		if len(vars)==to_find:
			if len(vars)==1:
				found_var.append(1)
				found_var_name.append(df2[0].tolist()[i])
				print('here',df2[0].tolist()[i],vars)
			elif len(vars)>1:
				found_var.append(2)
				found_var_name.append(df2[0].tolist()[i])
				print(len(vars),'vars')
		else:
			multi['multi_var'].append(df2[0].tolist()[i])
			varsf=df2[df2[0]==df2[0].tolist()[i]][1].tolist()
			multi['multi_found'].append(varsf)
			multi['multi_bases'].append(to_find)
			multi['multi_bases_found'].append(len(vars))
			found_var.append(0)
			found_var_name.append('')
			print(df2[0].tolist()[i])
		#	fv=[]
		#	for j in vars:
		#		fv.append(df2[1].tolist()[j])
		#	print(fv)
		#	if len(fv)==1:
		#		found_var.append(1)
		#		found_var_name.append(fv[0])
		#	else:
		#		found_var.append(0)
		#		print(fv)


dfmul=pd.DataFrame(multi)
dfmul.to_csv('multi.txt',sep='\t')

df2['confirmed']=found_var
df2['confirmed_var']=found_var_name

while 0 in found_var:
	found_var.remove(0)

print(len(found_var))

df3=df2[df2['confirmed']>0]

print(len(df3[0]),len(found_var),type(df3[0].tolist()))
indices=[]
selected=[]
variant=[]
for i in range(len(df1[0].tolist())):
	selected.append(0)
	variant.append('')
for i in range(len(df3[0].tolist())):
	index=df1[df3[1].tolist()[i]==df1[0]].index.tolist()
	print(index,'index')
	selected[index[0]]=df3['confirmed'].tolist()[i]
	variant[index[0]]=df3[0].tolist()[i]
	print(df3[0].tolist()[i],df1[df3[1].tolist()[i]==df1[0]][0].tolist(),'here1' ,index,df1.iloc[index[0],0],selected[index[0]],variant[index[0]])


df1['selected']=selected
df1['variant']=variant

info_list=[]
for i in range(len(df1[0].tolist())):
        info=df1[5].tolist()[i].split(';')
        info[1]=variant[i]
        info_list.append((';').join(info))
        print('variant',info_list[i])

df1[5]=info_list
df4=df1[df1['selected']==1]
#print(df1[df1[0]=='chr1@11317308@G@A'][0].index.tolist(),variant[13],'here2')
#print(df1[df1[0]=='chr1@11317308@G@A'][0].tolist(),df1[df1[0]=='chr1@11317308@G@A']['variant'].tolist(),'here3')
df4.to_csv('confirmed.csv',sep='\t')

df4=df1[df1['selected']==2]
df4.to_csv('confirmed_multi.csv',sep='\t')
