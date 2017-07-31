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
make_dataframe.py – combine data from other variant callers and finds final blending errors
run_annotate_tier.sh – divide total dataset into confidence tiers.  Add annotations from COSMIC, dbSNP and info from customer data

Example final output:
gtb     name    outlier0        ploidy  rel_error       res0    res1    test    test_alt        test_round      val     val0    washu   file    blend   RKOcov  HCT116cov       SW48cov 701cov       GTBcov  variants        conf
0.90625 ;chr19@3119184@T@G;GNA11        0.0     0@2@0   0.006940000000000001    0.00625 0.0157  0.90698 0.9     0.9     0.90625 0.8843  NS      0.0     0.0     221.0   257.0   136.0   749.0352.0   chr19@3119184@T@G       0
NS      ;chr16@89836507@T@C;FANCA       1.0     0@3@0   0.00541 0.00487 0.00487 0.93982 0.93333 0.9     0.89513 0.89513 0.89513 3.0     0.0     226.0   443.0   151.0   534.0   NS      chr16@89836507@T@C   0
0.11268 ;chr4@1805593@G@C;FGFR3 0.0     0@2@0   0.12675999999999998     0.01268 0.01707 0.1134  0.1     0.1     0.11268 0.11707 NS      0.0     0.0     77.0    97.0    54.0    643.0   213.0chr4@1805593@G@C        1
NS      ;chr5@176517292@A@G;FGFR4       0.0     0@2@0   0.007940000000000001    0.007140000000000001    0.007140000000000001    0.89333 0.9     0.9     0.89286 0.89286 0.89286 2.0     0.0 264.0    361.0   192.0   224.0   NS      chr5@176517292@A@G      0
NS      ;chr17@29667620@A@O;NF1;p.Q2340fs;@;NF1 0.0     2@0@0   0.09651 0.007240000000000002    0.007240000000000002    0.07037 0.075   0.075   0.06776 0.06776 0.06776 1.0     0.0     83.0100.0    69.0    974.0   NS      chr17@29667620@A@O      0
0.45205 ;chr17@58011498@G@T;RPS6KB1     0.0     0@0@2   0.048310000000000006    0.02295 0.02295 0.42576 0.475   0.475   0.45205 0.45205 NS      0.0     0.0     62.0    111.0   66.0    0.0 73.0     chr17@58011498@G@T      0
NS      ;chr7@128843169@C@T;SMO 0.0     3@2@0   0.02901 0.02466 0.02466 0.8593799999999999      0.85    0.85    0.87466 0.87466 0.87466 3.0     0.0     48.0    38.0    41.0    367.0   NS  chr7@128843169@C@T       0
0.29213 ;chr19@40743956@G@A;AKT2        0.0     0@0@2   0.05906 0.01919 0.02837 0.36507 0.325   0.325   0.30581 0.29663 NS      0.0     0.0     117.0   120.0   73.0    1046.0  89.0    chr19@40743956@G@A   0
0.20513 ;chr7@55228053@A@T;EGFR 0.0     3@0@0   0.17949 0.04487 0.06402000000000001     0.26111 0.25    0.25    0.20513 0.18598 NS      0.0     0.0     27.0    30.0    26.0    892.0   234.0chr7@55228053@A@T       1

