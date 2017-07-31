import re
import sys
import pandas as pd

keys=['gtb','name','outlier0','ploidy','rel_error','res0','res1','test','test_alt','test_round', 'val', 'val0','washu','file','blend']
df0=pd.read_csv( 'combined_scores.txt' , '\t' , names=keys )
df0=df0[df0['blend'] < 1 ]
df1=pd.read_csv( 'confirmed_1.csv' , '\t' , header=None )
df1.rename(columns={6:'name'}, inplace=True)
df2 = pd.merge(df0, df1, how='left', on=['name'])

print(df2.columns)

df3=df2[[ 'name' ,  8 , 12 , 16 , 20 , 25  ]]

df4 = pd.merge(df0, df3, how='left', on=['name'])

print(df4.columns)

df4=df4.drop_duplicates()
df4=df4.reset_index(drop=True)

print(df4.columns)
print(df1.columns)
index = df4['name'][df4[8].isnull()].dropna().tolist()
print(index)
print(df1[32][0])
for item in index:
	print(item,'item')
	item0=item.replace(';','')
	print(df1[ item0 == df1[32] ])
	ind=df1[ item0 == df1[32] ][ [ 'name' , 8 , 12 , 16 , 20 , 25  ]].index.values[0]
	print('i',ind)
	for col in [8 , 12 , 16 , 20 , 25]:
		df4.set_value(  df4[df4['name']==item].index.values[0]  , col , df1.loc[ ind  , col ]  )
		print( item, 'z' , df4[df4['name']==item].index.values[0] ,  df1.loc[ ind  , col ] , df4.loc[ df4['name']==item  , col ].values ,'fill')

ind=df4.loc[ df4['name'].isnull() , 'name' ].index
df4=df4.drop( ind  )
df4=df4.reset_index(drop=True)


for column in df4:
        df4[column] = df4[column].astype(str)
        for i in range(len(df4[column])):
                item=df4[column][i]
                if type(item) is str:
                        if 'NS' not in item:
                                if 'chr' not in item :
                                        if '@' not in item:
                                                df4.set_value(i,column,str(round(float(item),5)))


df4.to_csv('combined_scores_DP.txt',sep=',', header=None , index=False )
