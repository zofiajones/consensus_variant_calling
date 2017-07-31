python get_ploidy.py | sed 's/p\.//g' | sed  's/ΔE746 - A750/ΔE746-A750/g' > variants_ploidy.txt

sed -i 's/ΔE746\t/ΔE746-A750\t/g' final_freq.txt

rm final_freq_ploidy_1.txt
while read p;do

array=($p)
q=$(echo ${array[0]} | awk -F ';' '{print $3,$4}' | sed 's/p\.//g' )

array1=($q)
echo -e ${array1[0]}"\t"${array1[1]}
r=$(grep -F ${array1[0]} variants_ploidy.txt | grep -F  ${array1[1]} | grep '[A-Z]' | awk '{print $3,$4,$5}')

if [[ -z $r   ]];then
r=$(echo -e "NF\tNF\tNF")
else

echo -e "found\t"$r

fi

echo -e $p"\t"$r >> final_freq_ploidy_1.txt

done < final_freq.txt

awk '{for(i = 1; i <= NF; i++) { printf "%s\t",$i } printf "\n"}' final_freq_ploidy_1.txt > final_freq_ploidy.txt
awk '{print $11/($10+$11),$7/($6+$7),$9/($8+$9)}' found_ploidy.txt > right.txt
awk '{print $1}' found_ploidy.txt | awk -F '_' '{print $1,$2}' > left.txt
paste -d '\t' left.txt right.txt > found_ploidy_1.txt

rm final_freq_ploidy_2.txt
while read p;do

r=$(echo -e "NF\tNF\tNF")
echo -e $p"\t"$r >> final_freq_ploidy_2.txt

done < final_freq_ploidy.txt

awk '{for(i = 1; i <= NF; i++) { printf "%s\t",$i } printf "\n"}' final_freq_ploidy_2.txt | sort | uniq > final_freq_ploidy.txt

