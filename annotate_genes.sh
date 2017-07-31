#grep 'GENE' gene.bed > gene_1.bed
#cp gene_1.bed gene.bed
#sed -i 's/^/chr/g' gene.bed
cat final_in_uniq.txt | grep chr > final_in_uniq_1.txt
cp final_in_uniq_1.txt final_in_uniq.txt
cat final_in_uniq.txt | grep chr | python find_genes.py
~/hb/tools/bedtools intersect -a combined.bed -b gene_1.bed -wb | awk '{print $1,$3,$4,$8}' > gene_find.txt
awk 'BEGIN{OFS=""}{print $1,"@",$2,"@",$3}' gene_find.txt > left.txt
grep -oP 'GENE=.*;' gene_find.txt | sed  's/;$//g' | awk -F '=' '{print $2}' > genes.txt
paste -d '\t' left.txt genes.txt > got_genes.txt
python add_genes.py
