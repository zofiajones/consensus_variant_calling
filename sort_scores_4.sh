#DIR[0]="/home/zjones/BRCA_variants/"
#DIR[1]="/home/zjones/endogenous_4_WES/"
#DIR[2]="/home/zjones/endogenous_1_WES/"
#DIR[3]="/home/zjones/endogenous_no_genomic/"
#DIR[4]="/home/zjones/endogenous_2_WES/"
#DIR[5]="/home/zjones/endogenous/"
#DIR[6]="/home/zjones/endogenous_3_WES/"
#DIR[7]="/home/zjones/HD701_endogenous_validate/"

DIR[0]="/home/zjones/exomes_1/"
#DIR[1]="/home/zjones/endogenous_2/"
#DIR[2]="/home/zjones/endogenous_3/"
#DIR[3]="/home/zjones/endogenous_4/"


#cd ~/
#mkdir scores
#cd scores

head -n 1 ~/endogenous_1_WES/scores.csv > header.csv

sed -i 's/$/\tfile/g' header.csv
#cd ~/exomes_1/
python join_scores.py
python combine_scores_1.py

sort combined_scores.txt | uniq | awk '{print $2}'  | awk -F ';' '{print $3}' | grep '[A-Z]'  > genes.txt
grep -v outlier genes.txt | sort  | uniq -c | sort -nrk1 > gene_sort.txt
echo -e "total variants\t"$(grep -v group combined_scores.txt | awk '{print $2}'  | awk -F ';' '{print $2}' | sort | uniq | wc -l) >> summary.txt
echo -e "genes\t"$(cat gene_sort.txt | wc -l) >> summary.txt

python add_DP.py
awk -F ',' '{if ($5<1)print $0}' combined_scores_DP.txt  > combined_scores_in.txt
awk -F ',' '{if ($5>=1)print $0}' combined_scores_DP.txt > combined_scores_out.txt

#awk -F ',' '{if ( $61 < 70 || $66 < 70 || $71 < 70 || $43 < 200 ) print $0  }' combined_scores_out.txt  > combined_scores_out_coverage.txt
#awk -F ',' '{if ( $61 < 70 || $66 < 70 || $71 < 70 || $43 < 200 ) print $0  }' combined_scores_DP.txt   > combined_scores_coverage.txt

echo -e "total_variants\t"$(cat combined_scores_in.txt | awk -F ',' '{print $2}'  | awk -F ';' '{print $2}' | sort | uniq | wc -l) >> summary.txt
echo -e "genes\t"$(cat combined_scores_in.txt | awk -F ',' '{print $2}'  | awk -F ';' '{print $3}' | sort | uniq | wc -l) >> summary.txt

awk -F ',' '{print $2}' combined_scores_in.txt | awk -F ';' '{print $2}' | sort | uniq > loc.txt
cat combined_scores_in.txt | awk -F ',' '{print $2}'  | awk -F ';' '{print $3}' | sort | uniq -c | sort -rnk1 > final_genes.txt
#sed  's/,chr[@0-9a-zA-Z]\+;/,;/g'  combined_scores_in.txt | sort | uniq > final_in_uniq.txt

cat combined_scores_in.txt | sort | uniq > final_in_uniq.txt

python add_DP_1.py
awk -F ',' '{if ($5<1)print $0}' combined_scores_DP_blend.txt  > combined_scores_in_blend.txt
awk -F ',' '{if ($5>=1)print $0}' combined_scores_DP_blend.txt > combined_scores_out_blend.txt
sed  's/,chr[@0-9a-zA-Z]\+;/,;/g'  combined_scores_in_blend.txt | sort | uniq > final_in_uniq_blend.txt

