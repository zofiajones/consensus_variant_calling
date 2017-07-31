blends variant freq from parental cell lines and finds rel error to product data

README
 run_scores.sh – run script
quick_check.sh – collect raw data from vcf files
match_variants_4.sh – combine variant calls from parentals and products, one line per genomic location
match_base.sh – clean out useless lines, get backend file for checking indel partial matches
match_bases.py – get SNPS, and one base indels, and also match up longer indels
parse_results_1.sh – further tidying of longer indels
match_ploidy.sh – add info from in house data
multi_scores.py – get variant freq from output of quick_check.sh rather than combining multiple genomic locations
multi_scores_1.py – print out partial matches
get_blend_3.py – get blending scores for indels, allow some flexibility in ploidy
get_blend_4.py – get blending scores for indels, strictly take best fit to a ploidy of either 2 or 3
parse_results.sh – as parse_results_1.sh but for snps and one base indels
get_blend_1.py – get blending scores for SNPs, allow some flexibility in ploidy
get_blend_2.py – get blending scores for SNPs, strictly take best fit to a ploidy of either 2 or 3
sort_scores.sh – combine outputs of get_blend_1.py, get_blend_2.py, get_blend_3.py and get_blend_4.py
annotate_genes.sh – add gene info
make_dataframe.py – combine data from other variant callers
run_annotate_tier.sh – divide total dataset into confidence tiers.  Add annotations from COSMIC, dbSNP and info from customer data
