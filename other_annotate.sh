sed  's/\tPASS//g' output_gene_hg19_PASS.vcf > output_gene_hg19_1.vcf

#awk 'BEGIN{OFS="@"}{print $1,$2,$3}' output_gene_hg19_1.vcf > loc.txt
#awk 'BEGIN{OFS="@"}{print $1,$2,$3}' confident_1.vcf > loc_1.txt
#paste -d '\t' missed_loc.txt confident_1.vcf > confident_2.vcf
#grep -vf loc.txt loc.txt > missed_loc.txt
#grep -f missed_loc.txt confident_2.vcf | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8}' > confident_3.vcf
#cat output_gene_hg19_1.vcf >> confident_3.vcf
#cat confident_3.vcf  > draft_1.vcf

cp output_gene_hg19_1.vcf draft_1.vcf

sort draft_1.vcf | uniq > draft.vcf
awk 'BEGIN{OFS="@"}{print $1,$3,$4,$5}' /hd1/HD701_set/request_end/all_parsed.vcf > cust_match.txt
rm draft_1.txt
rm rs.txt
rm cust_info.txt
rm rs_base.txt
while read p;do

array=($p)
echo -e ${array[0]}"@"${array[1]}"@"${array[3]}"@"${array[4]} | sed 's/\./\\./g' > match_base.txt
loc=$(cat match_base.txt | awk  -F '@' 'BEGIN{OFS=""}{print $1,":",$2-10,"-",$2+10}' | sed 's/chr//g')
#~/hb/tools/tabix ~/clinvar_20160831_1.vcf.gz $loc | awk 'BEGIN{OFS=""}{print "chr",$1,"@",$2,"@",$4,"@",$5,"@",$3}' > rs_base.txt
#rs_id=$(grep -f match_base.txt rs_base.txt | awk -F '@' '{print $5}')
echo $loc
~/hb/tools/tabix /home/zjones/All_20160601.vcf.gz $loc | python get_mnp_loc_1.py | awk 'BEGIN{OFS=""}{print "chr",$1,"@",$2,"@",$3,"@",$4,"@",$5}'  >> rs_base.txt
rs_id=$(grep -f match_base.txt rs_base.txt | sort | uniq | awk -F '@' '{print $5}' | awk '{printf "%s," ,$i}END{print ""}' | wc -l )
echo -e $rs_id"\trs_id"
cust_info=$(grep -f match_base.txt cust_match.txt | awk -F '@' '{printf "%s", $5}' | grep '[0-9]' | wc -l)
echo -e $cust_info"\tcust_info"
echo -e $p"\t"$rs_id"\t"$cust_info >> draft_1.txt

if [[ $rs_id == 0  ]];then
rs_id="NA"
else
rs_id=$(grep -f match_base.txt rs_base.txt | sort | uniq | awk -F '@' '{print $5}' | awk '{printf "%s," ,$i}END{print ""}')
fi

if [[  $cust_info == 0  ]];then
cust_info="NA"
else
cust_info=$(grep -f match_base.txt cust_match.txt | awk -F '@' '{printf "%s", $5}' )
fi

echo -e $rs_id >> rs.txt
echo -e $cust_info >> cust_info.txt

done < draft.vcf

#awk '{print $2}' final_gene_list.txt > genes.txt
#grep -of unique_gene.txt draft.vcf > genes_confident.txt
#cat draft_1.txt | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sed 's/\.;//g' | sed 's/;$//g' | sed 's/;$//g'  > confident.vcf
#paste -d '\t' draft.vcf rs.txt cust_info.txt genes_confident.txt > confident.vcf
#cat genes_confident.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}'  > genes_sorted.txt

#grep -f ~/endogenous/genes.txt genes.txt

awk '{print $3}' draft.vcf | sed 's/;//g' > genes_found.txt

rm where.txt
while read p;do

q=$(echo $p | grep -f ~/endogenous/genes.txt)

if [[ -z $q  ]];then

echo -e "0" >> where.txt

else

echo -e "1" >> where.txt

fi

done < genes_found.txt

paste -d '\t' draft.vcf rs.txt cust_info.txt where.txt > confident.txt
