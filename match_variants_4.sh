cat all.vcf  | python get_mnp_loc_2.py  > potential.vcf

awk 'BEGIN{OFS=""}{print $1,"@",$2,"@",$3,"@",$4,"@",$6,";;"}' all.vcf | sort | uniq > potential.vcf
awk 'BEGIN{OFS=""}{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$6,";;"}' all.vcf | sort | uniq > potential_1.vcf

rm out_store.txt
cp potential.vcf out_variants_4_head.vcf
grep -v '#' out_variants_4_head.vcf | awk -F '@' 'BEGIN{OFS="@"}{print $1,$2,$3,$4}'  > loc.txt
sed 's/*/\*/g' loc.txt > loc_1.txt
grep -v '#' out_variants_4_head.vcf | awk -F '@' 'BEGIN{OFS="@"}{print $1,$2,$3,$4,$3}'  > alt_1.txt

#bash unique_loc_fb.sh hd701.vcf
#cp out.vcf hd7010.vcf

bash unique_loc.sh SW48.vcf
cp out.vcf SW480.vcf
bash unique_loc.sh RKO.vcf
cp out.vcf RKO0.vcf
bash unique_loc.sh HCT116.vcf
cp out.vcf HCT1160.vcf
bash unique_loc_fb.sh GTB.vcf
cp out.vcf GTB0.vcf
bash unique_loc_fb.sh HD701.vcf
cp out.vcf HD7010.vcf

cp SW480.vcf SW48.txt
cp RKO0.vcf RKO.txt
cp HCT1160.vcf HCT116.txt
cp GTB0.vcf GTB.txt
cp HD7010.vcf HD701.txt
#cp hd7010.vcf hd701.txt

rm tru.txt
rm out_store.txt
while read p;do

q=$(echo $p | awk -F '@' '{print $1,$2-1,$2,$3,$4}')
array=($q)
loc=$(echo -e ${array[0]}":"${array[1]}"-"${array[2]})


echo -e "loc\t"$loc

matchRKO=$(grep $p  RKO.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2}')
matchHCT116=$(grep $p  HCT116.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2}')
echo "HCT116 0"
echo $matchHCT116
matchSW48=$(grep $p  SW48.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2}')
matchGTB=$(grep $p   GTB.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5}')
#matchHD701=$(grep $p   HD701.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2}')
matchHD701=$(grep $p   hd701.txt | awk -F '@' '{print $5}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $2,$3}')


#if [[ -z $matchGTB  ]];then
#DP=$(~/hb/tools/samtools view /hd1/HD701_bams/GTB_348_IonXpress_028_.bam  $loc | wc -l)
#matchGTB="NA NA NA NA "$DP
#fi

matchRKO0=$(grep $p  RKO.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="\t"}{print $1,$2}')
matchHCT1160=$(grep $p  HCT116.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="\t"}{print $1,$2}')
echo "HCT116 1"
echo $matchHCT1160
matchSW480=$(grep $p  SW48.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="\t"}{print $1,$2}')
matchGTB0=$(grep $p   GTB.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="\t"}{print $1,$2}')
#matchHD7010=$(grep $p  HD701.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="\t"}{print $1,$2}')
matchHD7010=$(grep $p   hd701.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk 'BEGIN{OFS="\t"}{print $1,$2}')

r=$(echo $matchGTB0 | awk '{print NF}')
if [[ $r < 5  ]];then
matchGTB0="NA NA"
fi

r=$(echo $matchRKO0 | awk '{print NF}')
if [[ $r < 2  ]];then
matchRKO0="NA NA"
fi

r=$(echo $matchSW480 | awk '{print NF}')
if [[ $r < 2  ]];then
matchSW480="NA NA"
fi

r=$(echo $matchHCT1160 | awk '{print NF}')
if [[ $r < 2  ]];then
echo "HCT116 2"
echo $matchHCT1160
matchHCT1160="NA NA"
fi

r=$(echo $matchHD7010 | awk '{print NF}')
if [[ $r < 2  ]];then
matchHD7010="NA NA"
fi

matchRKO1=$(grep $p  RKO.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')
matchHCT1161=$(grep $p  HCT116.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')
matchSW481=$(grep $p  SW48.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')
matchGTB1=$(grep $p   GTB.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')
#matchHD7011=$(grep $p   HD701.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')
matchHD7011=$(grep $p   hd701.txt | awk -F '@' '{print $3,$4}'  | sort | uniq | awk  'BEGIN{OFS="@"}{print $1,$2}')

q=$(echo $p | awk -F '@' '{print $1,$2-10,$2+10,$3,$4}')
array=($q)
loc0=$(echo -e ${array[0]}":"${array[1]}"-"${array[2]})

~/hb/tools/tabix ~/clinvar_20160831.vcf.gz $loc0 | python get_mnp_loc_chr.py | awk 'BEGIN{OFS=""}{print "chr",$1,"@",$2,"@",$4,"@",$5,"@",$3}'  > rs_base.txt
rs_id=$(grep $p rs_base.txt | sort | uniq | awk -F '@' '{printf "%s,", $5}' )

echo -e $rs_id"\trs_id"

awk 'BEGIN{OFS="@"}{print $1,$3,$4,$5}' /hd1/HD701_set/request_end/all_parsed.vcf > cust_match.txt
cust_info=$(grep $p cust_match.txt | awk -F '@' '{printf "%s", $5}' )

echo -e $cust_info"\tcust_info"

r=$(echo $rs_id | awk '{print NF}')
if [[ $r < 1  ]];then
rs_id="NA"
fi

r=$(echo $cust_info | awk '{print NF}')
if [[ $r < 1  ]];then
cust_info="NA"
fi

r=$(echo $matchRKO | awk '{print NF}')
if [[ $r < 2  ]];then
DP=$(~/hb/tools/samtools view /hd1/HD701_set/RKO/RKO.alignment.bam  $loc | wc -l)
matchRKO="NA "$DP
fi

r=$(echo $matchHCT116 | awk '{print NF}')
if [[ $r < 2  ]];then
echo "HCT116 3"
echo $matchHCT116
DP=$(~/hb/tools/samtools view /hd1/HD701_set/HCT116/HCT116.alignment.bam  $loc | wc -l)
matchHCT116="NA "$DP
fi

r=$(echo $matchSW48 | awk '{print NF}')
if [[ $r < 2  ]];then
DP=$(~/hb/tools/samtools view /hd1/HD701_set/SW48/SW48.alignment.bam  $loc | wc -l)
matchSW48="NA "$DP
fi

r=$(echo $matchGTB | awk '{print NF}')
if [[ $r < 5  ]];then
DP=$(~/hb/tools/samtools view /hd1/HD701_bams/GTB_348_IonXpress_028_.bam  $loc | wc -l)
matchGTB="NA NA NA NA "$DP
fi

r=$(echo $matchHD701 | awk '{print NF}')
if [[ $r < 2  ]];then
DP=$(~/hb/tools/samtools view /home/zjones/WashU/CCATCCT/all_sorted_rmdup_sorted.bam  $loc | wc -l)
matchHD701="NA "$DP
fi


tru701=$matchHD701"\t"$matchHD7010
truRKO=$matchRKO"\t"$matchRKO0
truHCT116=$matchHCT116"\t"$matchHCT1160
truSW48=$matchSW48"\t"$matchSW480
truGTB=$matchGTB"\t"$matchGTB0

echo -e $truRKO"\t"$truHCT116"\t"$truSW48"\t"$tru701"\t"$truGTB

echo -e $p"\t"$truRKO"\t"$truHCT116"\t"$truSW48"\t"$tru701"\t"$truGTB >> tru.txt

echo -e "HCT116\t"$truHCT116 >> out_store.txt
echo -e "SW48\t"$truSW48 >> out_store.txt
echo -e "RKO\t"$truRKO >> out_store.txt
echo -e "HD701f\t"$tru701 >> out_store.txt
echo -e "GTB\t"$truGTB >> out_store.txt

echo -e "HCT116\t"$truHCT116"SW48\t"$truSW48"RKO\t"$truRKO"HD701f\t"$tru701"GTB\t"$truGTB"\t"$p

done < loc_1.txt
#done < loc_2.txt

awk 'BEGIN{OFS="\t"}{for(i = 2; i <= 24; i++) { printf "%s\t" , $i } printf "\n" }' tru.txt > tru1.txt
grep -v '#' potential_1.vcf > variants.vcf
paste -d '\t' variants.vcf tru1.txt > request_variants_fb.txt
#awk 'BEGIN{OFS="\t"}{for(i = 1; i <= NF; i++) { printf "%s\t" , $i } printf "\n" }' request_variants_fb_1.txt > request_variants_fb.txt

