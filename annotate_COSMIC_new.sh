#bash annotate_genes.sh
cat $1 | python parse_final_in_1.py  | sort | uniq > confident_1.vcf
cat $1 | python parse_final_in_1.py | sort | uniq | awk -F '\t' 'BEGIN{OFS="@"}{print $1,$2}' | sort | uniq > loc.txt

paste -d '\t' loc.txt confident_1.vcf > confident_1a.vcf
rm genes_mul.txt
rm left.txt
rm right.txt
while read p;do

grep $p confident_1a.vcf | awk '{printf "%s;",$4}END{print ""}' >> genes_mul.txt
grep $p confident_1a.vcf | awk -F '\t' 'BEGIN{OFS="\t"}{print $2,$3}' | sort | uniq >> left.txt
grep $p confident_1a.vcf | awk -F '\t' 'BEGIN{OFS="\t"}{print $5,$6,$7,$8,$9}' | sort | uniq >> right.txt

done < loc.txt


paste -d '\t' left.txt genes_mul.txt right.txt > confident_2.vcf

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$7,$8}' confident_2.vcf > list.bed

cat list.bed | python get_COSMIC.py  > output_gene_hg19_PASS.vcf
