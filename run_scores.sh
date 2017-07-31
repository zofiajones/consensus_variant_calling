bash quick_check.sh
bash match_variants_4.sh
bash match_base.sh
python match_bases.py
bash parse_results_1.sh
bash match_ploidy.sh
cp final_freq_ploidy.txt multi_freq.txt
python multi_scores.py | sort | uniq > multi_freq_1.txt
python multi_scores_1.py
cat multi_freq_1.txt | python get_blend_3.py
cat multi_freq_1.txt | python get_blend_4.py
bash parse_results.sh
bash match_ploidy.sh
cat final_freq_ploidy.txt | python get_blend_1.py
cat final_freq_ploidy.txt | python get_blend_2.py
bash sort_scores_4.sh
bash annotate_genes.sh
python make_dataframe.py
bash run_annotate_tier.sh
