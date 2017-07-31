#cat potential.vcf | python get_mnp_loc_3.py | sort | uniq > variants_4.vcf
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$6}' all.vcf | sort | uniq > all_1.vcf
python get_var_4.py

cat request_variants_fb.txt | python parse_final_freq.py | grep 'chr'  | python parse_req_fb.py | grep 'chr' | sort | uniq | sed 's/NA/NS/g'  > req_var.txt

awk -F '\t' '{print $2}' variants_4.vcf > locs.txt

rm found.txt
while read p;do

q=$(grep $p req_var.txt | wc -l)

if [[ $q == 1  ]];then

echo "1" >> found.txt

else

echo "0" >> found.txt

fi

done < locs.txt

paste -d '\t' variants_4.vcf found.txt > variants_5.vcf



