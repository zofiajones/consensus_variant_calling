import pandas as pd
import re
import sys
import os

dir=['/home/zjones/exomes_1']
dirs=os.listdir(os.getcwd())

d0=dict()
d1=dict()
d2=dict()

j=0;

for i in range(15):
	d0[i]=[]
	d1[i]=[]

for i in range(84):
	d2[i]=[]


df0=pd.DataFrame(d0)
df1=pd.DataFrame(d1)
df2=pd.DataFrame(d2)
for item in dir:
	df=pd.read_csv(item+'/scores.csv',sep='\t')
	if df0.empty:
		df0=df;
	else:
		df0=df0.append(df,ignore_index=True)
		#print('True')
	df=pd.read_csv(item+'/scores_2.csv',sep='\t')
	df0=df0.append(df,ignore_index=True)

	df=pd.read_csv(item+'/scores_1.csv',sep='\t')
	if df1.empty:
		df1=df;
	else:
		df1=df1.append(df,ignore_index=True)
	df=pd.read_csv(item+'/scores_4.csv',sep='\t')
	df1=df1.append(df,ignore_index=True)

	df=pd.read_csv(item+'/confirmed.csv',sep='\t')
	if df2.empty:
		df2=df;
	else:
		df2=df2.append(df,ignore_index=True)
	df=pd.read_csv(item+'/confirmed_multi.csv',sep='\t')
	df2=df2.append(df,ignore_index=True)
	print(item,'item')

print(df2.columns)

file=[str(j)]*len(df0['group'].tolist())
dfa=pd.DataFrame(file,columns=['file'])
df0=pd.concat([df0, dfa], axis=1)

file=[str(j)]*len(df1['group'].tolist())
dfa=pd.DataFrame(file,columns=['file'])
df1=pd.concat([df1, dfa], axis=1)

file=[str(j)]*len(df2['0'].tolist())
dfa=pd.DataFrame(file,columns=['file'])
df2=pd.concat([df2, dfa], axis=1)

#df0.to_csv('combined_scores_1.txt',sep='\t',header=None,index=False)
#df1.to_csv('combined_scores_3.txt',sep='\t',header=None,index=False)
df2.to_csv('confirmed_1.csv',sep='\t',header=None,index=False)

blend=['0']*len(df0['gtb'].tolist())
dfa=pd.DataFrame(blend,columns=['blend'])
df0=pd.concat([df0, dfa], axis=1)

blend=['1']*len(df1['gtb'].tolist())
dfa=pd.DataFrame(blend,columns=['blend'])
df1=pd.concat([df1, dfa], axis=1)
df0=df0.append(df1,ignore_index=True)

df0 = df0.drop(['Unnamed: 0','group'], 1)
df0=df0.drop_duplicates(keep='first')
df0=df0.reset_index(drop=True)


#df0.to_csv('combined_scores.txt',sep='\t',header=None,index=False)

for column in df0:
        df0[column] = df0[column].astype(str)
        for i in range(len(df0[column])):
                item=df0[column][i]
                if type(item) is str:
                        if 'NS' not in item:
                                if 'chr' not in item :
                                        if '@' not in item:
                                                df0.set_value(i,column,str(round(float(item),5)))

df0.to_csv('combined_scores.txt',sep='\t',header=None,index=False)
print(df0.columns)
