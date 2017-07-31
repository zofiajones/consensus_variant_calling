import re
import pandas as pd

df=pd.read_csv('all_1.vcf',sep='\t',header=None)

i_store=[]
count=[]
for i in range(len(df[0].tolist())):
	chr=str(df[0].tolist()[i])
	pos=str(df[1].tolist()[i])
	ref=df[2].tolist()[i]
	alt=df[3].tolist()[i]
	index=('@').join([chr,pos,ref,alt])
	i_store.append(index)
	var=df[4].tolist()[i]
	count.append( len(df[df[4]==var][4].tolist()) )

df['index']=i_store
df['count']=count

df1=df[[ 4, 'index' , 'count'  ]]

print(len(df1[4].tolist()))
df1.to_csv('variants_4.vcf',sep='\t',header=None,index=False)
