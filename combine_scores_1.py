import re
import sys
import pandas as pd

keys=['gtb','name','outlier0','ploidy','rel_error','res0','res1','test','test_alt','test_round', 'val', 'val0','washu','file','blend']
df0=pd.read_csv('combined_scores.txt',sep='\t',names=keys)

df5=df0[df0['blend'] < 1 ]
df6=df0[df0['blend'] > 0 ]

print(df5['blend'].tolist()[:10],df6['blend'].tolist()[:10])

df7a=df5.append(df6,ignore_index=True)
df7=df7a.drop_duplicates(keep='first')
print(len(df5['gtb'].tolist()),len(df6['gtb'].tolist()),len(df7['gtb'].tolist()),len(df7a['gtb'].tolist()))

variants1=set(df5['name'].tolist())
variants2=set(df6['name'].tolist())
variants_unique_2=list(variants2.difference(variants1))
variants_unique_1=list(variants1.difference(variants2))
variants_int=list(variants1.intersection(variants2))

print(df7.columns)
df8=pd.DataFrame(columns=df7.columns,index=range(len(df7['gtb'].tolist())))
for i in range( len(df7['name'].tolist()) ):
	#print(df7[3].tolist()[i])
	ind=df7[df7['name'].tolist()[i]==df7['name']].index.tolist()
	blend=df7['blend'][ind].tolist()
	rel_error=df7['rel_error'][ind].tolist()
	#print(ind,df7[3][ind].tolist(),df7['blend'][ind].tolist(),df7[6][ind].tolist())
#	df2df2[3] == variants_unique_2[i]
#	print('i0',ind,i,len(ind),blend,rel_error)
	if len(ind)==1:
		#print(df7.loc[ind[0],:])
		df8.loc[i,:]=df7.loc[ind[0],:]
#		print('i',i,0)
		#print(df8.loc[i,:])
	else:
		if rel_error[blend.index(0)] > 1:
			if rel_error[blend.index(1)] < rel_error[blend.index(0)]:
				#df8.append(df7.loc[ind[blend.index(1)],:])
				df8.loc[i,:]=df7.loc[ind[blend.index(1)],:]
				#print(df8.loc[i,:])
			else:
				df8.loc[i,:]=df7.loc[ind[blend.index(0)],:]
				#print(df8.loc[i,:])
		else:
			df8.loc[i,:]=df7.loc[ind[blend.index(0)],:]

df8a=df8.drop_duplicates(keep='first')
df8=df8a.reset_index(drop=True)

for column in df8:
	for i in range(len(df8[column])):
		item=df8[column][i]
		if type(item) is str:
			if 'NS' not in item:
				if 'chr' not in item :
					if '@' not in item:
						df8.set_value(i,column,str(round(float(item),5)))
		else:
			df8.set_value(i,column,str(round(float(item),5)))
df8.to_csv('combined_scores_blend.txt',sep='\t',header=None,index=False)
print(df8.columns)

