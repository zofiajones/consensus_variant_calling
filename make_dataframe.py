import re
import os
import sys
import pandas as pd

keys=['gtb','name','outlier0','ploidy','rel_error','res0','res1','test','test_alt','test_round', 'val', 'val0','washu','file','blend','RKOcov','HCT116cov','SW48cov','701cov','GTBcov']

df=pd.read_csv('final_in_uniq_gene.txt',sep='\t',names=keys,header=None)
#dfb=pd.read_csv('/home/zjones/scores/combined_scores_in.txt',sep=',',names=keys,header=None)
#df=df.append(dfb,ignore_index=True)
dfb=pd.read_csv('/home/zjones/scores_WES/final_in_uniq_gene.txt',sep='\t',names=keys,header=None)
df=df.append(dfb,ignore_index=True)
dfb=pd.read_csv('/home/zjones/scores_fil/final_in_uniq_gene.txt',sep='\t',names=keys,header=None)
df=df.append(dfb,ignore_index=True)
dfb=pd.read_csv('/home/zjones/scores_multiplex/final_in_uniq_gene.txt',sep='\t',header=None,names=keys)
df=df.append(dfb,ignore_index=True)
#dfb=pd.read_csv('/home/zjones/tru_0/combined_scores_in.txt',sep=',',header=None,names=keys)
#df=df.append(dfb,ignore_index=True)


df=df.drop_duplicates(keep='first')
df=df.reset_index(drop=True)

#df=pd.read_csv('uniq_all.txt',sep='\t',)

df=df[ df['blend'] < 1  ]
df=df.reset_index(drop=True)

df.to_csv('final_in_all.txt',sep='\t')

df_3=df[df['file']==0]
df_mul=df[df['file']==1]
print('mul',len(df_mul))
df_fil=df[df['file']==2]
df_WES=df[df['file']==3]
df_fb=df[df['file']==4]

var_list=[]
for var in df['name'].tolist():
	#print(var,'checkvar')
	var_list.append(var.split(';')[1])
variants=list(set(var_list))
df['variants']=var_list

count=[]
source=[]
for var in variants:
	print('var',var)
	count.append(len( df[ var == df['variants'] ]['variants'].tolist()))
	#print('count',len(df[ var == df[1] ][0].tolist()))
	group=''
	for item in df[ var == df['variants'] ]['file'].tolist():
		group=group+'+'+str(item)
	source.append( group  )

var_count=dict()
var_count['count']=count
var_count['var']=variants
var_count['source']=source

dfv=pd.DataFrame(var_count)
print(len(count),len(variants))
print( len(dfv[dfv['count']>1]['count'].tolist()) ,  dfv[dfv['count']>1]['source'].tolist()  )

print('lengths','0','fb','1','WES','2','mul','3','3')
print('lengths',len(df_fb['res0'].tolist()),len(df_WES['res0'].tolist()),len(df_mul['res0'].tolist()),len(df_3['res0'].tolist()))
df0=df[df['res0'] < 0.01]
df1=df[(df['res0'] >=  0.01) & (df['res0'] < 0.05)]
df2=df[(df['res0'] >=  0.05) & (df['res0'] < 0.1)]
df3=df[(df['res0'] >=  0.05) & (df['res0'] < 0.15)]
df4=df[(df['res0'] >= 0.15) & (df['res0'] < 1 ) ]

df_list=[df0,df1,df2,df3,df4]
print(len(df0.index),len(df1.index),len(df2.index),len(df3.index),len(df4.index))

for i in range(len(df_list)):
	df0=df_list[i]
	gene=[]
	print('var',len(df0.index))
	for j in range(len(df0.index)):
		gene.append(df0['name'].tolist()[j].split(';')[2])
	print('gene',i,len(list(set(gene))))

for i in range(len(df_list)):
	df0=df_list[i]
	for name,group in df0.groupby('res0'):
		print(i,name,len(group))

#print( list(dfv.groupby('source'))[0][0]  )
#print('groups')
#print(dfv.columns)
for name,group in list(dfv.groupby('source')):
	print( name,len(group) )

#print(dfv.columns)
#print(list(dfv.groupby('source')['var'])[3][0],list(list(dfv.groupby('source')['var'])[3][1]),'here0')
#print(len(list(dfv.groupby('source')['var'])[3]))
#multi_vars=list(list(dfv.groupby('source')['var'])[4][1])

d_uniq=dict()
for key in keys:
	d_uniq[key]=[]

keys.append('variants')
df_uniq=pd.DataFrame(d_uniq,columns=keys)
source=sorted(df['file'].unique().tolist())

print('source',type(source[0]),df['file'].unique() )
source=[0.0,1.0,2.0,3.0,4.0,5.0]
print('source',sorted(df['file'].unique().tolist()))

df.to_csv('test.csv','\t',index=False)

for var in variants:
#	print('var1',var)
	dfc=df[df['variants']==var]
	dfc=dfc.reset_index(drop=True)
	found=0
#	if 'chr17@41246245@C@A' in var:
#		print('dfc','chr17@41246245@C@A')
#		print(dfc,type(dfc),dfc['file'])
	for s in source:
		dfd=dfc[dfc['file']==s]
		if not dfd.empty and found ==0:
#			print('source',s,'var3',var)
			df_uniq=df_uniq.append(dfd,ignore_index=True)
			found=1


conf=['NA']*len(df_uniq['res0'])
df_uniq['conf']=conf
df_uniq.loc[df_uniq['rel_error'] <  0.1, 'conf' ]=0
df_uniq.loc[(df_uniq['rel_error'] >=  0.1) & (df_uniq['rel_error'] < 0.2),'conf']=1
df_uniq.loc[(df_uniq['rel_error'] >=  0.2) & (df_uniq['rel_error'] < 0.3),'conf']=2
df_uniq.loc[(df_uniq['rel_error'] >=  0.3) & (df_uniq['rel_error'] < 0.4),'conf']=3
df_uniq.loc[(df_uniq['rel_error'] >= 0.4) & (df_uniq['rel_error'] < 1 ) ,'conf']=4


df_uniq.to_csv('uniq_all.txt',sep='\t',index=False)
print(df_uniq.columns)
