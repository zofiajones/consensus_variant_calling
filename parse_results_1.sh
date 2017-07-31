#grep -v ':' request_variants_sorted.txt | awk  '{for(i = 1; i <= NF; i++) { printf "%s\t", $i }printf "\n" }'  > tabbed.txt

#awk -F '\t' 'BEGIN{OFS="\t"}{print $3,$7,$42,$75,$58,$9,$12,$15}' tabbed.txt > results.txt


cat confirmed.csv | python filter_cover.py | awk -F '\t' 'BEGIN{OFS="\t"}{print $7,$20,$25,$8,$12,$16}' confirmed_multi.csv | grep -v variant | grep chr > results.txt
#awk -F '\t' 'BEGIN{OFS="\t"}{print $7,$20,$25,$8,$12,$16}'  confirmed_multi.csv | grep -v variant | grep chr > results.txt


cat results.txt | sort | uniq > final_freq.txt

sed -i 's/HD701,//g' final_freq.txt

sed -i 's/,\t/\t/g' final_freq.txt
