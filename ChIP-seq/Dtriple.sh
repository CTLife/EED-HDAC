
danpos.py     dtriple       \
1-BED/1_H3K27ac-EEDko-AdultCM_Rep1.bed,1-BED/1_H3K27ac-EEDko-AdultCM_Rep2.bed,1-BED/1_H3K27ac-WT-AdultCM_Rep1.bed,1-BED/1_H3K27ac-WT-AdultCM_Rep2.bed     \
--bg   1-BED/1_Input-H3K27ac-AdultCM_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 2-Dtriple   --save 1     --edge 1     \
--count 10000000    --span 20     --clonalcut 1e-10         >>2-Dtriple/2-Dtriple-1.runLog   2>&1




danpos.py     dtriple       \
1-BED/2_H3K27ac-AdultCM-P5aavEED_Rep1.bed,1-BED/2_H3K27ac-AdultCM-P5aavEED_Rep2.bed,1-BED/2_H3K27ac-AdultCM-P5aavGFP_Rep1.bed,1-BED/2_H3K27ac-AdultCM-P5aavGFP_Rep2.bed     \
--bg   1-BED/2_Input-AdultCM-P5aavBoth_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 2-Dtriple   --save 1     --edge 1     \
--count 10000000    --span 20     --clonalcut 1e-10         >>2-Dtriple/2-Dtriple-2.runLog   2>&1



danpos.py     dtriple       \
1-BED/3_H3K27ac-AdultCM-P14aavEED_Rep1.bed,1-BED/3_H3K27ac-AdultCM-P14aavEED_Rep2.bed      \
--bg   1-BED/3_Input-AdultCM-P14aavEED_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 2-Dtriple   --save 1     --edge 1     \
--count 10000000    --span 20     --clonalcut 1e-10         >>2-Dtriple/2-Dtriple-3.runLog   2>&1



danpos.py     dtriple       \
1-BED/4_H3K27ac-AdultCM-P25aavEED_Rep1.bed,1-BED/4_H3K27ac-AdultCM-P25aavEED_Rep2.bed,1-BED/4_H3K27ac-AdultCM-P25aavGFP_Rep1.bed,1-BED/4_H3K27ac-AdultCM-P25EEDheto_Rep1.bed     \
--bg   1-BED/4_Input-AdultCM-P25_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 2-Dtriple   --save 1     --edge 1     \
--count 10000000    --span 20     --clonalcut 1e-10         >>2-Dtriple/2-Dtriple-4.runLog   2>&1
