import subprocess
import os
import re
import MySQLdb

cmd_root="mysql --defaults-extra-file=~/config.cnf -D mysql"
cmd0='use mysql;'

db=MySQLdb.connect(db="mysql",read_default_file="~/config.cnf")
c=db.cursor()

#file=open('/home/zjones/parental_exomes/genes_test.txt','r')
#file=open('/home/zjones/parental_exomes/genes.txt','r')
file=open('/home/zjones/endogenous/genes.txt','r')
lines=file.readlines()
file1=open('output_ploidy.txt','w')

id=[]
par=['73','71','74']

par_dict=dict()
par_dict['73']='RKO'
par_dict['71']='HCT116'
par_dict['74']='SW48'
for line in lines:
#	print('line',line)
	m=line.strip().split('\t')
	gene=m[0]
#	print('gene',gene)
#	mutation=m[1]
#	mutation=m[8]
	freq=[]
	id=[]
	aa_change=[]
	aa=dict()
	#cmd0='select id,protein_change  from mutations where gene=\'' + gene + '\' and protein_change=\'' + mutation + '\';'
	cmd0='select id,protein_change from mutations where gene=\'' + gene + '\';'
#	print(cmd0)
	c.execute(cmd0)
	results=c.fetchall()
	if len(results)>0:
		for entry0 in results:
#			print(entry0)
			id.append(entry0[0])
			aa_change.append(entry0[1])
			aa[str(entry0[0])]=entry0[1]

	if len(id)>0:
		for id0 in id:
#			print('id',id0)
			freq=[]
			ploidy=[]
			for par_id in par:
				wt=''
				mut=''
				freq0=''
				ploidy0=0
				cmd0='select num_mutated_alleles,num_wildtype_alleles from mutation_cell_line where mutation_id=\'' + str(id0) + '\' and cell_line_id=\'' + par_id +'\';'
				c.execute(cmd0)
				results=c.fetchall()
				if len(results)>0:
					for entry0 in results:
						mut=entry0[0]
						wt=entry0[1]
#						print('mut',mut,wt,aa[str(id0)],par_dict[par_id])
						freq0=str(float(mut)/(float(mut)+float(wt)))
						ploidy0=str(int(mut)+int(wt))
				if not wt:
					freq0="NF"
					ploidy0="NF"
				freq.append(freq0)
				ploidy.append(ploidy0)

#			print('freq_vec',freq,aa[str(id0)])
			print(gene,aa[str(id0)],freq[0],freq[1],freq[2],str(id0))
			file1.write(("\t").join([gene,aa[str(id0)],freq[0],freq[1],freq[2],str(id0),ploidy[0],ploidy[1],ploidy[2],"\n" ]))
