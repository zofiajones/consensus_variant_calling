#python make_dataframe.py
#python make_dataframe_1.py

awk -F '\t' '{if ($22==0) print $0  }' uniq_all.txt > final_tier_00.txt
awk -F '\t' '{if ($22==1) print $0  }' uniq_all.txt > final_tier_10.txt
awk -F '\t' '{if ($22==2) print $0  }' uniq_all.txt > final_tier_20.txt
awk -F '\t' '{if ($22==3) print $0  }' uniq_all.txt > final_tier_30.txt
awk -F '\t' '{if ($22==4) print $0  }' uniq_all.txt > final_tier_40.txt

sort final_tier_00.txt | uniq > final_tier_0.txt
sort final_tier_10.txt | uniq > final_tier_1.txt
sort final_tier_20.txt | uniq > final_tier_2.txt
sort final_tier_30.txt | uniq > final_tier_3.txt
sort final_tier_40.txt | uniq > final_tier_4.txt

awk -F '\t' '{print $2}' final_tier_0.txt | awk -F ';' '{print $3}' | sort | uniq -c > tier_0_genes.txt
awk -F '\t' '{print $2}' final_tier_1.txt | awk -F ';' '{print $3}' | sort | uniq -c > tier_1_genes.txt
awk -F '\t' '{print $2}' final_tier_2.txt | awk -F ';' '{print $3}' | sort | uniq -c > tier_2_genes.txt
awk -F '\t' '{print $2}' final_tier_3.txt | awk -F ';' '{print $3}' | sort | uniq -c > tier_3_genes.txt
awk -F '\t' '{print $2}' final_tier_4.txt | awk -F ';' '{print $3}' | sort | uniq -c > tier_4_genes.txt

touch unique_genes.txt
awk '{print $2}' tier_0_genes.txt > unique_genes.txt
grep -vf unique_genes.txt tier_1_genes.txt | sort | uniq > 1_not_u.txt
awk '{print $2}' 1_not_u.txt >> unique_genes.txt
grep -vf tier_1_genes.txt tier_2_genes.txt | sort | uniq > 2_not_u.txt
awk '{print $2}' 2_not_u.txt >> unique_genes.txt
grep -vf unique_genes.txt tier_3_genes.txt | sort | uniq > 3_not_u.txt
awk '{print $2}' 3_not_u.txt >> unique_genes.txt
grep -vf unique_genes.txt tier_4_genes.txt | sort | uniq > 4_not_u.txt
awk '{print $2}' 4_not_u.txt >> unique_genes.txt

cp final_in_uniq_gene.txt unique_genes_3.txt
cp ~/scores_multiplex/final_in_uniq.txt unique_genes_fb.txt

awk -F '\t' '{print $2}' unique_genes_3.txt | awk -F ';' '{print $2}'  | sort | uniq | grep '[a-z]' > 3_var.txt
awk -F '\t' '{print $2}' unique_genes_fb.txt | awk -F ';' '{print $2}' | sort | uniq > var.txt

touch unique_var.txt
cat 3_var.txt > unique_var.txt
grep -vf 3_var.txt var.txt | sort | uniq > f_not_3.txt
cat f_not_3.txt >> unique_var.txt

cat ~/scores_multiplex/final_in_1.txt > final_in.txt
cat final_in_uniq_gene.txt >> final_in.txt

awk -F '\t' '{print $2}' final_in.txt | awk -F ';' '{print $2}' | sort | uniq -c | sort -rnk1 > vars.txt

bash annotate_COSMIC_new.sh final_tier_0.txt
#bash annotate_COSMIC_1_1.sh final_tier_0.txt
bash other_annotate.sh
cp confident.txt confident_0.txt

bash annotate_COSMIC_new.sh final_tier_1.txt
#bash annotate_COSMIC_1_1.sh final_tier_1.txt
bash other_annotate.sh
cp confident.txt confident_1.txt

bash annotate_COSMIC_new.sh final_tier_2.txt
#bash annotate_COSMIC_1_1.sh final_tier_2.txt
bash other_annotate.sh
cp confident.txt confident_2.txt

bash annotate_COSMIC_new.sh final_tier_3.txt
#bash annotate_COSMIC_1_1.sh final_tier_3.txt
bash other_annotate.sh
cp confident.txt confident_3.txt

bash annotate_COSMIC_new.sh final_tier_4.txt
#bash annotate_COSMIC_1_1.sh final_tier_4.txt
bash other_annotate.sh
cp confident.txt confident_4.txt

cat confident_0.txt | python reanchor.py > confident_0.vcf
cat confident_1.txt | python reanchor.py > confident_1.vcf
cat confident_2.txt | python reanchor.py > confident_2.vcf
cat confident_3.txt | python reanchor.py > confident_3.vcf
cat confident_4.txt | python reanchor.py > confident_4.vcf

awk -F '\t' '{print $3}' confident_0.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}' > tier0_genes.txt
awk -F '\t' '{print $3}' confident_1.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}' > tier1_genes.txt
awk -F '\t' '{print $3}' confident_2.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}' > tier2_genes.txt
awk -F '\t' '{print $3}' confident_3.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}' > tier3_genes.txt
awk -F '\t' '{print $3}' confident_4.vcf | sort | uniq -c | sort -rnk1 | awk 'BEGIN{OFS="\t"}{print $1,$2}' > tier4_genes.txt

awk -F '\t' '{print $3}' confident_0.vcf > genes_1.txt
awk -F '\t' '{print $3}' confident_1.vcf >> genes_1.txt
awk -F '\t' '{print $3}' confident_2.vcf >> genes_1.txt


sort genes_1.txt | uniq > uniq_genes.txt

sed -i 's/;//g' uniq_genes.txt
awk '{print $4}' gene.bed | grep -oP 'GENE=[A-Z 0-9]+' | grep -oP '[A-Z 0-9]+' | grep -v GENE > genes.txt
rm genes_in.txt
rm genes_missed.txt
while read p;do

awk -v x=$p '{if($1==x)print $0}' uniq_genes.txt >> genes_in.txt
q=$(awk -v x=$p '{if($1==x)print $0}' uniq_genes.txt | wc -l)

if [[ $q ==  0  ]];then
echo $p >> genes_missed.txt
fi

done < genes.txt

grep -f genes_in.txt uniq_genes.txt > genes_out.txt
