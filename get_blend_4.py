import re
import os
import sys
import numpy as np
import pandas as pd
import tabix
#3,4,11,12,13

url="/home/zjones/exomes/HD701_fil_washu_sam.vcf.gz"
tb = tabix.open(url)
gtb_info=[]

url1="/home/zjones/exomes/hd701_tru_fb_1.vcf.gz"
tb1 = tabix.open(url1)


test_ploidies=[2,3]
results=dict()
final_freq=pd.read_csv('request_variants_fb.txt',sep="\t", header=None)

file1=open('ploidies_4.txt','w')

keys=['val','res0','washu','gtb','test','test_round','test_alt','res1','group','name','val0','rel_error','ploidy']
print('keys',keys)

exceptions=dict()
exceptions['name']=[]
exceptions['reason']=[]

for item in keys:
	results[item]=[]

for line in sys.stdin:
	exp=0
	m=line.strip().split('\t')
	val=0.0
	info=m[0]
	print(m)
	if len(m) <  12:
		continue
	washu=m[1]
	gtb=m[2]
	dPCR_ploidy=list(set(m[6:]))
	if washu != "NS":
		washu=float(m[1])
	if gtb != "NS":
		gtb=float(m[2])
	if washu != "NS":
		val=washu
	elif gtb != "NS":
		val=gtb
	if washu == "NS" and gtb == "NS":
		print('exception','no somatic',washu,gtb)
		exceptions['name'].append(info)
		exceptions['reason'].append('no somatic')
		continue
	if m[3]=="NS" and m[4]=="NS" and m[5]=="NS":
		if len(dPCR_ploidy)==1:
			if dPCR_ploidy[0]=="NF":
				print('exception','no germline',m[3],m[4],m[5],dPCR_ploidy,info)
				exceptions['name'].append(info)
				exceptions['reason'].append('no germline')
				continue
	if float(val)==0.0:
		print('exception','freq 0',val)
		exceptions['name'].append(info)
		exceptions['reason'].append('0 val')
		continue
	if float(val) ==1:
		print('exception','freq 1',val)
		exceptions['name'].append(info)
		exceptions['reason'].append('1 val')
		continue
	freq_store=[m[3].replace('NS','0.0'),m[4].replace('NS','0.0'),m[5].replace('NS','0.0')]
	print(len(m))
	if len(m)> 9:
		if m[9]=="NF" and m[10]=="NF" and m[11]=="NF":
			freq_multiplex=[m[6],m[7],m[8]]
		else:
			freq_multiplex=[m[9],m[10],m[11]]
	else:
		freq_multiplex=[m[6],m[7],m[8]]
	print('freq_multiplex',freq_multiplex)
	test=float(0.15*float(freq_store[0])+0.2*float(freq_store[1])+0.65*float(freq_store[2]))
	k=-1
	found_freq0=[]
	found_ploidy=[]
	print('freq_store',freq_store)
	for freq in freq_store:
		k=k+1
		ploidies=[]
		freq_round0=[]
		residuals=[]
		print(freq,'freq')
		if float(freq) == 0.0 or float(freq) == 0:
			residuals.append(0.0)
			freq_round0.append(0.0)
			found_freq0.append(0.0)
			found_ploidy.append(str(0))
			print(k,'result0',2 , 'freq_round' , 0 , 'residuals' ,0 , 'freq', 0)
#			continue
		elif float(freq) > 0.9:
			residuals.append(abs(float(freq)-1.0))
			freq_round0.append(1.0)
			found_freq0.append(1.0)
			found_ploidy.append(str(0))
			print(k,'result0',2 , 'freq_round' , 1 , 'residuals' ,abs(float(freq)-1.0) , 'freq', 1)
#			continue
		else:
			print('test_ploidies',test_ploidies)
			for i in test_ploidies:
				residual_store=[]
#				print('ploidies0',i)
				for j in range(i):
					#print(freq)
					try:
#						print('residual',abs(float(freq)-float(j)/float(i)),i,j,float(j)/float(i))
						residual_store.append(abs( float(freq) - (float(j)/float(i)) ))
					except:
						#print('error',freq,'line',line)
						pass
                                	        #print(residual_store)
                                	        #print('residual',min(residual_store))
				imin=residual_store.index(min(residual_store))
                        	                #print('residual check',residual_store[imin],i,imin)
				residuals.append(residual_store[imin])
				freq_round0.append(float(imin)/float(i))
                        	        #print('ploidies',ploidies,'freq',freq,'res',ploidies[ploidies.index(min(ploidies))],'found ploidy',test_ploidies[ploidies.index(min(ploidies))])
			#gene_ploidy.append(test_ploidies[ploidies.index(min(ploidies))])
			ploidy_index=residuals.index(min(residuals))
			print('ploidy_index',ploidy_index)
	#		new_residuals=residuals.remove(residuals.index(min(residuals)))
	#		2_ploidy_index=new_residuals.index(min(new_residuals))
			print(k,'result',test_ploidies[ploidy_index] , 'freq_round' , freq_round0[ploidy_index] , 'residuals' ,residuals[ploidy_index] , 'freq', freq )
			found_freq0.append(freq_round0[ploidy_index])
			found_ploidy.append(str(test_ploidies[ploidy_index]))
			print('ploidy0',test_ploidies,ploidy_index)
	print('found0',found_freq0,freq_store)
	freq_default=["NA","NA","NA"]
	ploidy_default=["NA","NA","NA"]
	ploidy_test=["NA","NA","NA"]
	ploidy_store=[]
	test_i=[]
	for i in range(len(found_freq0)):
		if float(found_freq0[i]) == 0.0:
			freq_default[i]=float(found_freq0[i])
			ploidy_default[i]=0.0
		elif float(found_freq0[i]) == 1.0:
			freq_default[i]=float(found_freq0[i])
			ploidy_default[i]=0.0
		else:
			ploidy_test[i]=found_freq0[i]
			test_i.append(i)
	test_store={}
	test_store['0']=[]
	test_store['1']=[]
	test_store['2']=[]
	for i in range(3):
		if i not in test_i:
			test_store[str(i)].append(freq_default[i])
		else:
#			test_store[str(i)].append(0.5)
#			if ploidy_test[i] != 0.5:
#				test_store[str(i)].append(ploidy_test[i])
			test_store[str(i)].append(ploidy_test[i])

	test_default=[]
	test_default_store=[]
	for i in range(len(test_store['0'])):
		for j in range(len(test_store['1'])):
			for k in range(len(test_store['2'])):
				print( 'vals' , float(0.15*test_store['0'][i]+0.2*test_store['1'][j]+0.65*test_store['2'][k] )  )
				test_default.append( float(0.15*test_store['0'][i]+0.2*test_store['1'][j]+0.65*test_store['2'][k] )  )
				test_default_store.append([test_store['0'][i],test_store['1'][j],test_store['2'][k] ] )


#	index_min=test_default.index( min(test_default) )
#	freq_found=test_default_store[index_min]
#	test_round=test_default[index_min]
#	fit_ploidies['round_0'].append(freq_found[0])
#	fit_ploidies['round_1'].append(freq_found[1])
#	fit_ploidies['round_2'].append(freq_found[2])
	#print('index_min',index_min)
#	test_default=float(0.15*freq_default[0]+0.2*freq_default[1]+0.65*freq_default[2])
#	test_round=float(0.15*float(found_freq[0])+0.2*float(found_freq[1])+0.65*float(found_freq[2]))

	freq_alt=[]
	#print('found02',found_freq,freq_store,freq_alt,freq_multiplex,test_round)
	for i in range(3):
		print('i',i)
		if freq_multiplex[i]!= "NF":
			if round(float(found_freq0[i])*10) != round(float(freq_multiplex[i])*10):
				freq_alt.append(float(freq_multiplex[i]))
				print('found03',found_freq0,freq_store,freq_alt,freq_multiplex)
				exp=1
			else:
				freq_alt.append(float(found_freq0[i]))
				print('found03',found_freq0,freq_store,freq_alt,freq_multiplex)
		else:
			freq_alt.append(float(found_freq0[i]))
			print('found03',found_freq0,freq_store,freq_alt,freq_multiplex)
	print('found01',found_freq0,freq_store,freq_alt)
	test_alt=float(0.15*float(freq_alt[0])+0.2*float(freq_alt[1])+0.65*float(freq_alt[2]))
	if test_alt not in test_default:
		test_default.append(test_alt)
		test_default_store.append([freq_alt[0],freq_alt[1],freq_alt[2]])

	m=info.split(';')
	m[0]=''
	info=(';').replace(';','').join(m)
	m=info.split('@')
	loc=m[0]+':'+m[1]+'-'+m[1]
	records=tb.querys(loc)
	record_mark=0
	ind_samtools=-1
	test_samtools=''
	for record in records:
		if '0/1' in record[9]:
			record_mark=record_mark+1
			n=record[7].split(';')
			for i in range(len(n)):
				if 'DP4' in n[i]:
					DP4=n[i].replace('DP4=','')
			n=DP4.split(',')
			samtools=( float(n[2])+float(n[3]) )/( float(n[0])+float(n[1])+float(n[2])+float(n[3])  )
			samtools_default=[]
			for item in test_default:
				if item == 0:
					samtools_default.append(100.0)
				else:
					samtools_default.append( abs(samtools-item)  )
			ind_samtools=samtools_default.index(min(samtools_default))
			test_samtools=float(0.15*float(test_default_store[ind_samtools][0])+0.2*float(test_default_store[ind_samtools][1])+0.65*float(test_default_store[ind_samtools][2]))

	records=tb1.querys(loc)
	record_mark=0
	ind_fb=-1
	test_fb=''
	for record in records:
		if ':' in record[9]:
			record_mark=record_mark+1
			n=record[8].split(':')
			n0=record[9].split(':')
			for i in range(len(n)):
				if 'AO' in n[i]:
					AO=n0[i]
				elif 'RO' in n[i]:
					RO=n0[i]
			if ',' not in AO:
				fb=( float(AO) )/( float(AO) + float(RO)  )
				fb_default=[]
				for item in test_default:
					if item == 0:
						fb_default.append(100.0)
					else:
						fb_default.append( abs(fb-item)  )
				ind_fb=fb_default.index(min(fb_default))
				test_fb=float(0.15*float(test_default_store[ind_fb][0])+0.2*float(test_default_store[ind_fb][1])+0.65*float(test_default_store[ind_fb][2]))


	if washu != 'NS':
		washu_default=[]
		for item in test_default:
			if item == 0:
				washu_default.append(100.0)
			else:
				washu_default.append( abs(float(washu)-item)  )
			ind_washu=washu_default.index(min(washu_default))
			test_washu=float(0.15*float(test_default_store[ind_washu][0])+0.2*float(test_default_store[ind_washu][1])+0.65*float(test_default_store[ind_washu][2]))
	if gtb != 'NS':
		gtb_default=[]
		for item in test_default:
			if item == 0:
				gtb_default.append(100.0)
			else:
				gtb_default.append( abs(float(gtb)-item)  )
			ind_gtb=gtb_default.index(min(gtb_default))
			test_gtb=float(0.15*float(test_default_store[ind_gtb][0])+0.2*float(test_default_store[ind_gtb][1])+0.65*float(test_default_store[ind_gtb][2]))

	if washu!='NS':
		ind=ind_washu
	elif ind_fb > -1:
		ind=ind_fb
	elif ind_samtools > -1:
		ind=ind_samtools
	else:
		ind=ind_gtb
	test0=float(0.15*float(test_default_store[ind][0])+0.2*float(test_default_store[ind][1])+0.65*float(test_default_store[ind][2]))
	if test0==0:
		continue
	inds=[]
	tests=[]
	if washu!='NS':
		tests.append( abs( test0-float(washu)  ) )
		inds.append(washu)
	if ind_fb > -1:
		tests.append( abs( test0-float(fb)  )  )
		inds.append(fb)
	if ind_samtools > -1 :
		tests.append( abs( test0-float(samtools)  )  )
		inds.append(samtools)
	if gtb != 'NS':
		tests.append( abs( test0-float(gtb)  )  )
		inds.append(gtb)
	ind_min=tests.index(min(tests))
	val=inds[ind_min]
	print('tests',tests)
	res0=min(tests)
	if len(tests)>1:
		tests.remove(min(tests))
		ind_min0=tests.index(min(tests))
		val0=inds[ind_min0]
		res1=min(tests)
	else:
		val0=val
		res1=res0
	test_round=test0
	rel_error=res0/test0
	group=0
	ploidy=('@').join(found_ploidy)
	entry=[val,res0,washu,gtb,test,test_round,test_alt,res1,group,info,val0,rel_error,ploidy]
	for i in range(len(entry)):
		results[keys[i]].append(entry[i])
	entry=str(test_default_store[ind][0])+'\t'+str(test_default_store[ind][1])+'\t'+str(test_default_store[ind][2])+'\t'+info+'\n'
	file1.write(entry)


results['outlier0']=[0.0]*len(results['val'])

df=pd.DataFrame(results)
df.to_csv('scores_4.csv',sep="\t")

df1=pd.DataFrame(exceptions)
df1.to_csv('exceptions_4.csv',sep="\t")

