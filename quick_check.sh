awk 'BEGIN{OFS=""}{print "chr",$1,":",$2+1,"-",$3+1}' gene.bed > gene_targets.txt
awk '{print $4}' gene.bed | grep -oP 'GENE=[A-Z 0-9]+' | grep -oP '[A-Z 0-9]+' | grep -v GENE > genes.txt

paste -d '\t' gene_targets.txt genes.txt > genes_targets.txt

#cat ~/Panel_beds/aml_comb.bed | awk 'BEGIN{OFS=""}{print $1,":",$2+1,"-",$3+1}' >   gene_targets.txt
#awk '{print $4}' ~/Panel_beds/aml_comb.bed > genes.txt

mkdir multiplex_endogenous
DIR="multiplex_endogenous/"

rm multiplex_endogenous/*


for i in $(ls *vcf.gz);do

#i=RKO.snp.vcf.gz

#~/hb/tools/tabix $i

k=$( basename $i)

j=$DIR${k%.vcf.gz}

#while read p;do
#echo "" > $j"_"${arr[1]}".txt"
#done < genes_targets.txt

while read p;do

arr=($p)
echo ${arr[0]}
echo ${arr[1]}

~/hb/tools/tabix $i ${arr[0]}  > found.txt

#while read q;do

#test0=$(echo $q | sed -n -e 's/.*SAF=\([0-9,]\+\).*SAR=\([0-9,]\+\);.*/\1 \2/p' | awk '{ if ($1>0 && $2>0){ print "1"  } else {print"0"}}' )
##test1=$(echo $q | sed -n -e 's/.*AB=\([0-9,\.]\+\).*/\1/p' | awk '{ if ($1<0.5){ print "1"  } else {print "0"}}' )
#test3=$(echo $q | sed -n -e 's/.*DP=\([0-9,\.]\+\).*/\1/p' | awk '{ if ($1>100){ print "1"  } else {print "0"}}' )
##test=$(( $test0 + $test1 + $test3))
#test=$(( $test0 + $test3 ))

#if [[ $test == 2 ]] ;then
cp found.txt $j"_"${arr[1]}".txt"
#fi

#done < found.txt

done < genes_targets.txt

done

rm GTB.vcf
rm RKO.vcf
rm SW48.vcf
rm HCT116.vcf
rm HD701.vcf
for i in $(ls mul*/GTB*);do

cat $i | python parse_multi_ion_bed.py >> GTB.vcf

done

for i in $(ls mul*/RKO*);do
echo $i
cat $i | python parse_multi_WES_bed.py >> RKO.vcf

done

for i in $(ls mul*/HCT116*);do
echo $i
cat $i | python parse_multi_WES_bed.py >> HCT116.vcf

done

for i in $(ls mul*/SW48*);do
#echo $i
cat $i | python parse_multi_WES_bed.py >> SW48.vcf

done

for i in $(ls mul*/hd701*varscan*);do

cat $i | python parse_multi_par_bed.py >> HD701.vcf

done

for i in $(ls mul*/hd701*fb*);do

cat $i | python parse_multi_bed_fb_1.py >> HD701_fb.vcf

done

for i in $(ls mul*/HD701*genes*);do

cat $i | python parse_multi_samtools_bed.py >> hd701_genes.vcf

done

for i in $(ls mul*/HD701*fil*);do

cat $i | python parse_multi_WES_bed_1.py >> hd701_fil.vcf

done


cat GTB.vcf > all.vcf
cat RKO.vcf >> all.vcf
cat HCT116.vcf >> all.vcf
cat SW48.vcf >> all.vcf
cat HD701.vcf >> all.vcf
