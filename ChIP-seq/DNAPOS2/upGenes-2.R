#############################################################################################################################
## Part  2:  Read in the input matrix data, and reduce the dimensions (rows or columns) of matrix, and save the information of input files.
##              We should get a n*m matrix that means it contains n rows (n DNA fragment) and m columns (m bin, 1 bin=20bp).
##              The element of matrix is reads density: ChIP-input (If ChIP-input <0, we let it = 0)
#############################################################################################################################





#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
A1_H2BGFP_WTCKO           <- read.table("3-heatmap/clusteredGenes.txt-F1_H2BGFP_CKO_out.txt",                      header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27me3_EEDrescue     <- read.table("3-heatmap/clusteredGenes.txt-F2_H3K27me3_rescue_out.txt",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_HDAC2_WTCKO            <- read.table("3-heatmap/clusteredGenes.txt-F5_HDAC_out.txt",                            header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27ac_WTCKO          <- read.table("3-heatmap/clusteredGenes.txt-F6_H3K27ac_out.txt",                         header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27me3_WTCKO         <- read.table("3-heatmap/clusteredGenes.txt-F7_H3K27me3_out.txt",                        header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27ac_HDACrescue     <- read.table("3-heatmap/clusteredGenes.txt-F8_H3K27ac_adult_rescue_out.txt",            header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27me3_HDACrescue    <- read.table("3-heatmap/clusteredGenes.txt-F9_H3K27me3_adult_rescue_out.txt",           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_H3K27ac_EEDrescue      <- read.table("3-heatmap/clusteredGenes.txt-F10_H3K27ac_rescue_out.txt",                 header=FALSE,   sep="\t",   quote = "",   comment.char = "")  
A1_MNase_WTCKO            <- read.table("3-heatmap/clusteredGenes.txt-G3_MNase_out.txt",                           header=FALSE,   sep="\t",   quote = "",   comment.char = "")  

dim(A1_H2BGFP_WTCKO )  
dim(A1_H3K27me3_EEDrescue)  
dim(A1_HDAC2_WTCKO)  
dim(A1_H3K27ac_WTCKO)  
dim(A1_H3K27me3_WTCKO)  
dim(A1_H3K27ac_HDACrescue)  
dim(A1_H3K27me3_HDACrescue)  
dim(A1_H3K27ac_EEDrescue)  
dim(A1_MNase_WTCKO)  





A1_H2BGFP_WTCKO        <-    A1_H2BGFP_WTCKO[-c(1,2), -c(1,2)]
A1_H3K27me3_EEDrescue  <-    A1_H3K27me3_EEDrescue[-c(1,2), -c(1,2)]
A1_HDAC2_WTCKO         <-    A1_HDAC2_WTCKO[-c(1,2), -c(1,2)] 
A1_H3K27ac_WTCKO       <-    A1_H3K27ac_WTCKO[-c(1,2), -c(1,2)]
A1_H3K27me3_WTCKO      <-    A1_H3K27me3_WTCKO[-c(1,2), -c(1,2)]
A1_H3K27ac_HDACrescue  <-    A1_H3K27ac_HDACrescue[-c(1,2), -c(1,2)]
A1_H3K27me3_HDACrescue <-    A1_H3K27me3_HDACrescue[-c(1,2), -c(1,2)]
A1_H3K27ac_EEDrescue   <-    A1_H3K27ac_EEDrescue[-c(1,2), -c(1,2)]
A1_MNase_WTCKO         <-    A1_MNase_WTCKO[-c(1,2), -c(1,2)]

dim(A1_H2BGFP_WTCKO )  
dim(A1_H3K27me3_EEDrescue)  
dim(A1_HDAC2_WTCKO)  
dim(A1_H3K27ac_WTCKO)  
dim(A1_H3K27me3_WTCKO)  
dim(A1_H3K27ac_HDACrescue)  
dim(A1_H3K27me3_HDACrescue)  
dim(A1_H3K27ac_EEDrescue)  
dim(A1_MNase_WTCKO)  





A2_H3K27me3_WT        <-    A1_H3K27me3_WTCKO[, c(1:500)]
A2_H3K27me3_CKO       <-    A1_H3K27me3_WTCKO[, c(511:1010)]
A2_H3K27ac_WT         <-    A1_H3K27ac_WTCKO[, c(1:500)]
A2_H3K27ac_CKO        <-    A1_H3K27ac_WTCKO[, c(511:1010)]
A2_HDAC2_WT           <-    A1_HDAC2_WTCKO[, c(1021:1520)] 
A2_HDAC2_CKO          <-    A1_HDAC2_WTCKO[, c(1531:2030)]  
A2_H2BGFP_WT_week0         <-    A1_H2BGFP_WTCKO[, c(1:500)]
A2_H2BGFP_CKO_week0        <-    A1_H2BGFP_WTCKO[, c(511:1010)]
A2_H2BGFP_WT_week4         <-    A1_H2BGFP_WTCKO[, c(1021:1520)]
A2_H2BGFP_CKO_week4        <-    A1_H2BGFP_WTCKO[, c(1531:2030)]
A2_MNase_WT          <-    A1_MNase_WTCKO[, c(511:1010)]
A2_MNase_CKO         <-    A1_MNase_WTCKO[, c(1531:2030)]
A2_H3K27me3_EEDrescue_P5EED    <-    A1_H3K27me3_EEDrescue[, c(1:500)]
A2_H3K27me3_EEDrescue_P5luci   <-    A1_H3K27me3_EEDrescue[, c(511:1010)]
A2_H3K27me3_EEDrescue_P14EED   <-    A1_H3K27me3_EEDrescue[, c(1021:1520)]
A2_H3K27me3_EEDrescue_P25EED   <-    A1_H3K27me3_EEDrescue[, c(1531:2030)]
A2_H3K27me3_EEDrescue_P25luci  <-    A1_H3K27me3_EEDrescue[, c(2041:2540)]
A2_H3K27me3_EEDrescue_P25hete  <-    A1_H3K27me3_EEDrescue[, c(2551:3050)]
A2_H3K27ac_EEDrescue_P5EED     <-    A1_H3K27ac_EEDrescue[, c(1:500)]
A2_H3K27ac_EEDrescue_P5luci    <-    A1_H3K27ac_EEDrescue[, c(511:1010)]
A2_H3K27ac_EEDrescue_P14EED    <-    A1_H3K27ac_EEDrescue[, c(1021:1520)]
A2_H3K27ac_EEDrescue_P25EED    <-    A1_H3K27ac_EEDrescue[, c(1531:2030)]
A2_H3K27ac_EEDrescue_P25luci   <-    A1_H3K27ac_EEDrescue[, c(2041:2540)]
A2_H3K27ac_EEDrescue_P25hete   <-    A1_H3K27ac_EEDrescue[, c(2551:3050)]
A2_H3K27me3_HDACrescue_EEDhete   <-    A1_H3K27me3_HDACrescue[, c(1:500)]
A2_H3K27me3_HDACrescue_AAVHDAC   <-    A1_H3K27me3_HDACrescue[, c(511:1010)]
A2_H3K27me3_HDACrescue_AAVluci   <-    A1_H3K27me3_HDACrescue[, c(1021:1520)]
A2_H3K27ac_HDACrescue_EEDhete    <-    A1_H3K27ac_HDACrescue[, c(1:500)]
A2_H3K27ac_HDACrescue_AAVHDAC    <-    A1_H3K27ac_HDACrescue[, c(511:1010)]
A2_H3K27ac_HDACrescue_AAVluci    <-    A1_H3K27ac_HDACrescue[, c(1021:1520)]




	
			
dim(A2_H3K27me3_WT)
dim(A2_H3K27me3_CKO)
dim(A2_H3K27ac_WT)
dim(A2_H3K27ac_CKO)
dim(A2_HDAC2_WT)
dim(A2_HDAC2_CKO)  
dim(A2_H2BGFP_WT_week0)
dim(A2_H2BGFP_CKO_week0)
dim(A2_H2BGFP_WT_week4)
dim(A2_H2BGFP_CKO_week4)
dim(A2_MNase_WT)
dim(A2_MNase_CKO)
dim(A2_H3K27me3_EEDrescue_P5EED)
dim(A2_H3K27me3_EEDrescue_P5luci)
dim(A2_H3K27me3_EEDrescue_P14EED)
dim(A2_H3K27me3_EEDrescue_P25EED)
dim(A2_H3K27me3_EEDrescue_P25luci)
dim(A2_H3K27me3_EEDrescue_P25hete)
dim(A2_H3K27ac_EEDrescue_P5EED)
dim(A2_H3K27ac_EEDrescue_P5luci)
dim(A2_H3K27ac_EEDrescue_P14EED)
dim(A2_H3K27ac_EEDrescue_P25EED)
dim(A2_H3K27ac_EEDrescue_P25luci)
dim(A2_H3K27ac_EEDrescue_P25hete)
dim(A2_H3K27me3_HDACrescue_EEDhete)
dim(A2_H3K27me3_HDACrescue_AAVHDAC)
dim(A2_H3K27me3_HDACrescue_AAVluci)
dim(A2_H3K27ac_HDACrescue_EEDhete)
dim(A2_H3K27ac_HDACrescue_AAVHDAC)
dim(A2_H3K27ac_HDACrescue_AAVluci)






cluster1_rows <- c(1:252)
cluster2_rows <- c(253:555)
cluster3_rows <- c(556:863)




C1_H3K27me3_WT        <-    A2_H3K27me3_WT[cluster1_rows, ]
C1_H3K27me3_CKO       <-    A2_H3K27me3_CKO[cluster1_rows, ]
C1_H3K27ac_WT         <-    A2_H3K27ac_WT[cluster1_rows, ]
C1_H3K27ac_CKO        <-    A2_H3K27ac_CKO[cluster1_rows, ]
C1_HDAC2_WT           <-    A2_HDAC2_WT[cluster1_rows, ]
C1_HDAC2_CKO          <-    A2_HDAC2_CKO[cluster1_rows, ]  
C1_H2BGFP_WT_week0         <-    A2_H2BGFP_WT_week0[cluster1_rows, ]
C1_H2BGFP_CKO_week0        <-    A2_H2BGFP_CKO_week0[cluster1_rows, ]
C1_H2BGFP_WT_week4         <-    A2_H2BGFP_WT_week4[cluster1_rows, ]
C1_H2BGFP_CKO_week4        <-    A2_H2BGFP_CKO_week4[cluster1_rows, ]
C1_MNase_WT          <-    A2_MNase_WT[cluster1_rows, ]
C1_MNase_CKO         <-    A2_MNase_CKO[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P5EED    <-    A2_H3K27me3_EEDrescue_P5EED[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P5luci   <-    A2_H3K27me3_EEDrescue_P5luci[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P14EED   <-    A2_H3K27me3_EEDrescue_P14EED[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P25EED   <-    A2_H3K27me3_EEDrescue_P25EED[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P25luci  <-    A2_H3K27me3_EEDrescue_P25luci[cluster1_rows, ]
C1_H3K27me3_EEDrescue_P25hete  <-    A2_H3K27me3_EEDrescue_P25hete[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P5EED     <-    A2_H3K27ac_EEDrescue_P5EED[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P5luci    <-    A2_H3K27ac_EEDrescue_P5luci[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P14EED    <-    A2_H3K27ac_EEDrescue_P14EED[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P25EED    <-    A2_H3K27ac_EEDrescue_P25EED[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P25luci   <-    A2_H3K27ac_EEDrescue_P25luci[cluster1_rows, ]
C1_H3K27ac_EEDrescue_P25hete   <-    A2_H3K27ac_EEDrescue_P25hete[cluster1_rows, ]
C1_H3K27me3_HDACrescue_EEDhete   <-    A2_H3K27me3_HDACrescue_EEDhete[cluster1_rows, ]
C1_H3K27me3_HDACrescue_AAVHDAC   <-    A2_H3K27me3_HDACrescue_AAVHDAC[cluster1_rows, ]
C1_H3K27me3_HDACrescue_AAVluci   <-    A2_H3K27me3_HDACrescue_AAVluci[cluster1_rows, ]
C1_H3K27ac_HDACrescue_EEDhete    <-    A2_H3K27ac_HDACrescue_EEDhete[cluster1_rows, ]
C1_H3K27ac_HDACrescue_AAVHDAC    <-    A2_H3K27ac_HDACrescue_AAVHDAC[cluster1_rows, ]
C1_H3K27ac_HDACrescue_AAVluci    <-    A2_H3K27ac_HDACrescue_AAVluci[cluster1_rows, ]

dim(C1_H3K27me3_WT)
dim(C1_H3K27me3_CKO)
dim(C1_H3K27ac_WT)
dim(C1_H3K27ac_CKO)
dim(C1_HDAC2_WT)
dim(C1_HDAC2_CKO)  
dim(C1_H2BGFP_WT_week0)
dim(C1_H2BGFP_CKO_week0)
dim(C1_H2BGFP_WT_week4)
dim(C1_H2BGFP_CKO_week4)
dim(C1_MNase_WT)
dim(C1_MNase_CKO)
dim(C1_H3K27me3_EEDrescue_P5EED)
dim(C1_H3K27me3_EEDrescue_P5luci)
dim(C1_H3K27me3_EEDrescue_P14EED)
dim(C1_H3K27me3_EEDrescue_P25EED)
dim(C1_H3K27me3_EEDrescue_P25luci)
dim(C1_H3K27me3_EEDrescue_P25hete)
dim(C1_H3K27ac_EEDrescue_P5EED)
dim(C1_H3K27ac_EEDrescue_P5luci)
dim(C1_H3K27ac_EEDrescue_P14EED)
dim(C1_H3K27ac_EEDrescue_P25EED)
dim(C1_H3K27ac_EEDrescue_P25luci)
dim(C1_H3K27ac_EEDrescue_P25hete)
dim(C1_H3K27me3_HDACrescue_EEDhete)
dim(C1_H3K27me3_HDACrescue_AAVHDAC)
dim(C1_H3K27me3_HDACrescue_AAVluci)
dim(C1_H3K27ac_HDACrescue_EEDhete)
dim(C1_H3K27ac_HDACrescue_AAVHDAC)
dim(C1_H3K27ac_HDACrescue_AAVluci)











C2_H3K27me3_WT        <-    A2_H3K27me3_WT[cluster2_rows, ]
C2_H3K27me3_CKO       <-    A2_H3K27me3_CKO[cluster2_rows, ]
C2_H3K27ac_WT         <-    A2_H3K27ac_WT[cluster2_rows, ]
C2_H3K27ac_CKO        <-    A2_H3K27ac_CKO[cluster2_rows, ]
C2_HDAC2_WT           <-    A2_HDAC2_WT[cluster2_rows, ]
C2_HDAC2_CKO          <-    A2_HDAC2_CKO[cluster2_rows, ]  
C2_H2BGFP_WT_week0         <-    A2_H2BGFP_WT_week0[cluster2_rows, ]
C2_H2BGFP_CKO_week0        <-    A2_H2BGFP_CKO_week0[cluster2_rows, ]
C2_H2BGFP_WT_week4         <-    A2_H2BGFP_WT_week4[cluster2_rows, ]
C2_H2BGFP_CKO_week4        <-    A2_H2BGFP_CKO_week4[cluster2_rows, ]
C2_MNase_WT          <-    A2_MNase_WT[cluster2_rows, ]
C2_MNase_CKO         <-    A2_MNase_CKO[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P5EED    <-    A2_H3K27me3_EEDrescue_P5EED[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P5luci   <-    A2_H3K27me3_EEDrescue_P5luci[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P14EED   <-    A2_H3K27me3_EEDrescue_P14EED[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P25EED   <-    A2_H3K27me3_EEDrescue_P25EED[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P25luci  <-    A2_H3K27me3_EEDrescue_P25luci[cluster2_rows, ]
C2_H3K27me3_EEDrescue_P25hete  <-    A2_H3K27me3_EEDrescue_P25hete[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P5EED     <-    A2_H3K27ac_EEDrescue_P5EED[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P5luci    <-    A2_H3K27ac_EEDrescue_P5luci[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P14EED    <-    A2_H3K27ac_EEDrescue_P14EED[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P25EED    <-    A2_H3K27ac_EEDrescue_P25EED[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P25luci   <-    A2_H3K27ac_EEDrescue_P25luci[cluster2_rows, ]
C2_H3K27ac_EEDrescue_P25hete   <-    A2_H3K27ac_EEDrescue_P25hete[cluster2_rows, ]
C2_H3K27me3_HDACrescue_EEDhete   <-    A2_H3K27me3_HDACrescue_EEDhete[cluster2_rows, ]
C2_H3K27me3_HDACrescue_AAVHDAC   <-    A2_H3K27me3_HDACrescue_AAVHDAC[cluster2_rows, ]
C2_H3K27me3_HDACrescue_AAVluci   <-    A2_H3K27me3_HDACrescue_AAVluci[cluster2_rows, ]
C2_H3K27ac_HDACrescue_EEDhete    <-    A2_H3K27ac_HDACrescue_EEDhete[cluster2_rows, ]
C2_H3K27ac_HDACrescue_AAVHDAC    <-    A2_H3K27ac_HDACrescue_AAVHDAC[cluster2_rows, ]
C2_H3K27ac_HDACrescue_AAVluci    <-    A2_H3K27ac_HDACrescue_AAVluci[cluster2_rows, ]

dim(C2_H3K27me3_WT)
dim(C2_H3K27me3_CKO)
dim(C2_H3K27ac_WT)
dim(C2_H3K27ac_CKO)
dim(C2_HDAC2_WT)
dim(C2_HDAC2_CKO)  
dim(C2_H2BGFP_WT_week0)
dim(C2_H2BGFP_CKO_week0)
dim(C2_H2BGFP_WT_week4)
dim(C2_H2BGFP_CKO_week4)
dim(C2_MNase_WT)
dim(C2_MNase_CKO)
dim(C2_H3K27me3_EEDrescue_P5EED)
dim(C2_H3K27me3_EEDrescue_P5luci)
dim(C2_H3K27me3_EEDrescue_P14EED)
dim(C2_H3K27me3_EEDrescue_P25EED)
dim(C2_H3K27me3_EEDrescue_P25luci)
dim(C2_H3K27me3_EEDrescue_P25hete)
dim(C2_H3K27ac_EEDrescue_P5EED)
dim(C2_H3K27ac_EEDrescue_P5luci)
dim(C2_H3K27ac_EEDrescue_P14EED)
dim(C2_H3K27ac_EEDrescue_P25EED)
dim(C2_H3K27ac_EEDrescue_P25luci)
dim(C2_H3K27ac_EEDrescue_P25hete)
dim(C2_H3K27me3_HDACrescue_EEDhete)
dim(C2_H3K27me3_HDACrescue_AAVHDAC)
dim(C2_H3K27me3_HDACrescue_AAVluci)
dim(C2_H3K27ac_HDACrescue_EEDhete)
dim(C2_H3K27ac_HDACrescue_AAVHDAC)
dim(C2_H3K27ac_HDACrescue_AAVluci)










C3_H3K27me3_WT        <-    A2_H3K27me3_WT[cluster3_rows, ]
C3_H3K27me3_CKO       <-    A2_H3K27me3_CKO[cluster3_rows, ]
C3_H3K27ac_WT         <-    A2_H3K27ac_WT[cluster3_rows, ]
C3_H3K27ac_CKO        <-    A2_H3K27ac_CKO[cluster3_rows, ]
C3_HDAC2_WT           <-    A2_HDAC2_WT[cluster3_rows, ]
C3_HDAC2_CKO          <-    A2_HDAC2_CKO[cluster3_rows, ]  
C3_H2BGFP_WT_week0         <-    A2_H2BGFP_WT_week0[cluster3_rows, ]
C3_H2BGFP_CKO_week0        <-    A2_H2BGFP_CKO_week0[cluster3_rows, ]
C3_H2BGFP_WT_week4         <-    A2_H2BGFP_WT_week4[cluster3_rows, ]
C3_H2BGFP_CKO_week4        <-    A2_H2BGFP_CKO_week4[cluster3_rows, ]
C3_MNase_WT          <-    A2_MNase_WT[cluster3_rows, ]
C3_MNase_CKO         <-    A2_MNase_CKO[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P5EED    <-    A2_H3K27me3_EEDrescue_P5EED[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P5luci   <-    A2_H3K27me3_EEDrescue_P5luci[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P14EED   <-    A2_H3K27me3_EEDrescue_P14EED[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P25EED   <-    A2_H3K27me3_EEDrescue_P25EED[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P25luci  <-    A2_H3K27me3_EEDrescue_P25luci[cluster3_rows, ]
C3_H3K27me3_EEDrescue_P25hete  <-    A2_H3K27me3_EEDrescue_P25hete[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P5EED     <-    A2_H3K27ac_EEDrescue_P5EED[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P5luci    <-    A2_H3K27ac_EEDrescue_P5luci[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P14EED    <-    A2_H3K27ac_EEDrescue_P14EED[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P25EED    <-    A2_H3K27ac_EEDrescue_P25EED[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P25luci   <-    A2_H3K27ac_EEDrescue_P25luci[cluster3_rows, ]
C3_H3K27ac_EEDrescue_P25hete   <-    A2_H3K27ac_EEDrescue_P25hete[cluster3_rows, ]
C3_H3K27me3_HDACrescue_EEDhete   <-    A2_H3K27me3_HDACrescue_EEDhete[cluster3_rows, ]
C3_H3K27me3_HDACrescue_AAVHDAC   <-    A2_H3K27me3_HDACrescue_AAVHDAC[cluster3_rows, ]
C3_H3K27me3_HDACrescue_AAVluci   <-    A2_H3K27me3_HDACrescue_AAVluci[cluster3_rows, ]
C3_H3K27ac_HDACrescue_EEDhete    <-    A2_H3K27ac_HDACrescue_EEDhete[cluster3_rows, ]
C3_H3K27ac_HDACrescue_AAVHDAC    <-    A2_H3K27ac_HDACrescue_AAVHDAC[cluster3_rows, ]
C3_H3K27ac_HDACrescue_AAVluci    <-    A2_H3K27ac_HDACrescue_AAVluci[cluster3_rows, ]

dim(C3_H3K27me3_WT)
dim(C3_H3K27me3_CKO)
dim(C3_H3K27ac_WT)
dim(C3_H3K27ac_CKO)
dim(C3_HDAC2_WT)
dim(C3_HDAC2_CKO)  
dim(C3_H2BGFP_WT_week0)
dim(C3_H2BGFP_CKO_week0)
dim(C3_H2BGFP_WT_week4)
dim(C3_H2BGFP_CKO_week4)
dim(C3_MNase_WT)
dim(C3_MNase_CKO)
dim(C3_H3K27me3_EEDrescue_P5EED)
dim(C3_H3K27me3_EEDrescue_P5luci)
dim(C3_H3K27me3_EEDrescue_P14EED)
dim(C3_H3K27me3_EEDrescue_P25EED)
dim(C3_H3K27me3_EEDrescue_P25luci)
dim(C3_H3K27me3_EEDrescue_P25hete)
dim(C3_H3K27ac_EEDrescue_P5EED)
dim(C3_H3K27ac_EEDrescue_P5luci)
dim(C3_H3K27ac_EEDrescue_P14EED)
dim(C3_H3K27ac_EEDrescue_P25EED)
dim(C3_H3K27ac_EEDrescue_P25luci)
dim(C3_H3K27ac_EEDrescue_P25hete)
dim(C3_H3K27me3_HDACrescue_EEDhete)
dim(C3_H3K27me3_HDACrescue_AAVHDAC)
dim(C3_H3K27me3_HDACrescue_AAVluci)
dim(C3_H3K27ac_HDACrescue_EEDhete)
dim(C3_H3K27ac_HDACrescue_AAVHDAC)
dim(C3_H3K27ac_HDACrescue_AAVluci)






####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################




