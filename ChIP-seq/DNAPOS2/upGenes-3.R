#############################################################################################################################
## Part 3:  Figures about nucleosome occupancy level (NOL).
#############################################################################################################################


 



#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################
#################################################################### Start ##########################################################################################################################################







###########################################################################  1. H2BGFP-EEDcko
subdir_1A_part3 <- paste(Part3_g,  "/1A-allRows-H2BGFP", sep = "")
if( ! file.exists(subdir_1A_part3) ) { dir.create(subdir_1A_part3,  recursive = TRUE) }

numOfColumns1 <- 500
MyAverageLines_1(vector2=c( colMeans(A2_H2BGFP_WT_week0),  colMeans(A2_H2BGFP_CKO_week0),  colMeans(A2_H2BGFP_WT_week4),  colMeans(A2_H2BGFP_CKO_week4)  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),  
                             "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1A_part3,     fileName2=paste("1A-allRows-H2BGFP-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

numOfRows1 <- 863
MyBoxViolinPlot_1(vector2=c(rowMeans(A2_H2BGFP_WT_week0),  rowMeans(A2_H2BGFP_CKO_week0),  rowMeans(A2_H2BGFP_WT_week4),  rowMeans(A2_H2BGFP_CKO_week4) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_1A_part3,  fileName2= paste("1A-allRows-H2BGFP-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=4.0, Ymin2=0, Ymax2=0.8 )





MyAverageLines_1(vector2=c( colMeans(A2_H2BGFP_WT_week0),  colMeans(A2_H2BGFP_CKO_week0)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko"   ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ), 
                 path2=subdir_1A_part3,     fileName2=paste("1A-allRows-H2BGFP-week0-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(A2_H2BGFP_WT_week0),  rowMeans(A2_H2BGFP_CKO_week0) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko"  ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ),  
                  path2=subdir_1A_part3,  fileName2= paste("1A-allRows-H2BGFP-week0-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=3.0, Ymin2=0, Ymax2=0.8 )

MyHypothesisTest_1(vector1=rowMeans(A2_H2BGFP_WT_week0),      vector2=rowMeans(A2_H2BGFP_CKO_week0),    file1=paste(subdir_1A_part3,  "/1A-allRows-H2BGFP-week0-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(A2_H2BGFP_WT_week0),      vector2=rowMeans(A2_H2BGFP_CKO_week0),    file1=paste(subdir_1A_part3,  "/1A-allRows-H2BGFP-week0-paired.txt",       sep = "")   )   












######################################################################
subdir_1B_part3 <- paste(Part3_g,  "/1B-cluster1-H2BGFP", sep = "")
if( ! file.exists(subdir_1B_part3) ) { dir.create(subdir_1B_part3,  recursive = TRUE) }

numOfColumns1 <- 500

MyAverageLines_1(vector2=c( colMeans(C1_H2BGFP_WT_week0),  colMeans(C1_H2BGFP_CKO_week0),  colMeans(C1_H2BGFP_WT_week4),  colMeans(C1_H2BGFP_CKO_week4)  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1B_part3,     fileName2=paste("1B-cluster1-H2BGFP-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

numOfRows1 <- 252
MyBoxViolinPlot_1(vector2=c(rowMeans(C1_H2BGFP_WT_week0),  rowMeans(C1_H2BGFP_CKO_week0),  rowMeans(C1_H2BGFP_WT_week4),  rowMeans(C1_H2BGFP_CKO_week4) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_1B_part3,  fileName2= paste("1B-cluster1-H2BGFP-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=4.0, Ymin2=0, Ymax2=0.8 )





MyAverageLines_1(vector2=c( colMeans(C1_H2BGFP_WT_week0),  colMeans(C1_H2BGFP_CKO_week0)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko"   ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ), 
                 path2=subdir_1B_part3,     fileName2=paste("1B-cluster1-H2BGFP-week0-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C1_H2BGFP_WT_week0),  rowMeans(C1_H2BGFP_CKO_week0) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko"  ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ),  
                  path2=subdir_1B_part3,  fileName2= paste("1B-cluster1-H2BGFP-week0-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=3.0, Ymin2=0, Ymax2=0.8 )


MyHypothesisTest_1(vector1=rowMeans(C1_H2BGFP_WT_week0),      vector2=rowMeans(C1_H2BGFP_CKO_week0),    file1=paste(subdir_1B_part3,  "/1B-cluster1-H2BGFP-week0-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C1_H2BGFP_WT_week0),      vector2=rowMeans(C1_H2BGFP_CKO_week0),    file1=paste(subdir_1B_part3,  "/1B-cluster1-H2BGFP-week0-paired.txt",       sep = "")   )   













######################################################################
subdir_1C_part3 <- paste(Part3_g,  "/1C-cluster2-H2BGFP", sep = "")
if( ! file.exists(subdir_1C_part3) ) { dir.create(subdir_1C_part3,  recursive = TRUE) }

numOfColumns1 <- 500

MyAverageLines_1(vector2=c( colMeans(C2_H2BGFP_WT_week0),  colMeans(C2_H2BGFP_CKO_week0),  colMeans(C2_H2BGFP_WT_week4),  colMeans(C2_H2BGFP_CKO_week4)  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1C_part3,     fileName2=paste("1C-cluster2-H2BGFP-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

numOfRows1 <- 303
MyBoxViolinPlot_1(vector2=c(rowMeans(C2_H2BGFP_WT_week0),  rowMeans(C2_H2BGFP_CKO_week0),  rowMeans(C2_H2BGFP_WT_week4),  rowMeans(C2_H2BGFP_CKO_week4) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_1C_part3,  fileName2= paste("1C-cluster2-H2BGFP-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=4.0, Ymin2=0, Ymax2=0.8 )





MyAverageLines_1(vector2=c( colMeans(C2_H2BGFP_WT_week0),  colMeans(C2_H2BGFP_CKO_week0)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko"   ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ), 
                 path2=subdir_1C_part3,     fileName2=paste("1C-cluster2-H2BGFP-week0-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C2_H2BGFP_WT_week0),  rowMeans(C2_H2BGFP_CKO_week0) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko"  ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255)  ),  
                  path2=subdir_1C_part3,  fileName2= paste("1C-cluster2-H2BGFP-week0-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=3.0, Ymin2=0, Ymax2=0.8 )


MyHypothesisTest_1(vector1=rowMeans(C2_H2BGFP_WT_week0),      vector2=rowMeans(C2_H2BGFP_CKO_week0),    file1=paste(subdir_1C_part3,  "/1C-cluster2-H2BGFP-week0-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C2_H2BGFP_WT_week0),      vector2=rowMeans(C2_H2BGFP_CKO_week0),    file1=paste(subdir_1C_part3,  "/1C-cluster2-H2BGFP-week0-paired.txt",       sep = "")   )   











######################################################################
subdir_1D_part3 <- paste(Part3_g,  "/1D-cluster3-H2BGFP", sep = "")
if( ! file.exists(subdir_1D_part3) ) { dir.create(subdir_1D_part3,  recursive = TRUE) }

numOfColumns1 <- 500

MyAverageLines_1(vector2=c( colMeans(C3_H2BGFP_WT_week0),  colMeans(C3_H2BGFP_CKO_week0),  colMeans(C3_H2BGFP_WT_week4),  colMeans(C3_H2BGFP_CKO_week4)  ),   
                 numSample2=4,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1),   
                                rep("week4_EEDheto", numOfColumns1),    rep("week4_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko",   "week4_EEDheto",  "week4_EEDko"  ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ), 
                 path2=subdir_1D_part3,     fileName2=paste("1D-cluster3-H2BGFP-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

numOfRows1 <- 308
MyBoxViolinPlot_1(vector2=c(rowMeans(C3_H2BGFP_WT_week0),  rowMeans(C3_H2BGFP_CKO_week0),  rowMeans(C3_H2BGFP_WT_week4),  rowMeans(C3_H2BGFP_CKO_week4) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1), rep("week4_EEDheto", numOfRows1),  rep("week4_EEDko", numOfRows1) ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko",  "week4_EEDheto",  "week4_EEDko" ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255),   "week4_EEDheto"="skyblue",  "week4_EEDko"="blue" ),  
                  path2=subdir_1D_part3,  fileName2= paste("1D-cluster3-H2BGFP-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=4.0, Ymin2=0, Ymax2=0.8 )





MyAverageLines_1(vector2=c( colMeans(C3_H2BGFP_WT_week0),  colMeans(C3_H2BGFP_CKO_week0)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("week0_EEDheto", numOfColumns1),    rep("week0_EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "week0_EEDheto",  "week0_EEDko"   ),     
                 colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255) ), 
                 path2=subdir_1D_part3,     fileName2=paste("1D-cluster3-H2BGFP-week0-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H2BGFP signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=6.05, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C3_H2BGFP_WT_week0),  rowMeans(C3_H2BGFP_CKO_week0) ),   
                  sampleType2=c( rep("week0_EEDheto", numOfRows1),   rep("week0_EEDko", numOfRows1)  ), 
                  sampleRank2=c( "week0_EEDheto",   "week0_EEDko"  ),  
                  colours2=c( "week0_EEDheto"=rgb(0,102,255,  255, maxColorValue=255),  "week0_EEDko"=rgb(0,153,255,  255, maxColorValue=255) ),  
                  path2=subdir_1D_part3,  fileName2= paste("1D-cluster3-H2BGFP-week0-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H2BGFP signal",   
                  height2=5,   width2=3.0, Ymin2=0, Ymax2=0.8 )


MyHypothesisTest_1(vector1=rowMeans(C3_H2BGFP_WT_week0),      vector2=rowMeans(C3_H2BGFP_CKO_week0),    file1=paste(subdir_1D_part3,  "/1D-cluster3-H2BGFP-week0-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C3_H2BGFP_WT_week0),      vector2=rowMeans(C3_H2BGFP_CKO_week0),    file1=paste(subdir_1D_part3,  "/1D-cluster3-H2BGFP-week0-paired.txt",       sep = "")   )   
























############################################################  2. MNase
subdir_2A_part3 <- paste(Part3_g,  "/2A-allRows-MNase", sep = "")
if( ! file.exists(subdir_2A_part3) ) { dir.create(subdir_2A_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 863

MyAverageLines_1(vector2=c( colMeans(A2_MNase_WT),  colMeans(A2_MNase_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "EEDheto",  "EEDko"   ),     
                 colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ), 
                 path2=subdir_2A_part3,     fileName2=paste("2A-allRows-MNase-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="MNase signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=5.23, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(A2_MNase_WT),  rowMeans(A2_MNase_CKO) ),   
                  sampleType2=c( rep("EEDheto", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",   "EEDko"  ),  
                  colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ),  
                  path2=subdir_2A_part3,  fileName2= paste("2A-allRows-MNase-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="MNase signal",   
                  height2=5,   width2=3.0, Ymin2=0.3, Ymax2=1 )

MyHypothesisTest_1(vector1=rowMeans(A2_MNase_WT),      vector2=rowMeans(A2_MNase_CKO),    file1=paste(subdir_2A_part3,  "/2A-allRows-MNase-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(A2_MNase_WT),      vector2=rowMeans(A2_MNase_CKO),    file1=paste(subdir_2A_part3,  "/2A-allRows-MNase-paired.txt",       sep = "")   )   





############################################################
subdir_2B_part3 <- paste(Part3_g,  "/2B-cluster1-MNase", sep = "")
if( ! file.exists(subdir_2B_part3) ) { dir.create(subdir_2B_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 252

MyAverageLines_1(vector2=c( colMeans(C1_MNase_WT),  colMeans(C1_MNase_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "EEDheto",  "EEDko"   ),     
                 colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ), 
                 path2=subdir_2B_part3,     fileName2=paste("2B-cluster1-MNase-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="MNase signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=5.23, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C1_MNase_WT),  rowMeans(C1_MNase_CKO) ),   
                  sampleType2=c( rep("EEDheto", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",   "EEDko"  ),  
                  colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ),  
                  path2=subdir_2B_part3,  fileName2= paste("2B-cluster1-MNase-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="MNase signal",   
                  height2=5,   width2=3.0, Ymin2=0.3, Ymax2=1 )

MyHypothesisTest_1(vector1=rowMeans(C1_MNase_WT),      vector2=rowMeans(C1_MNase_CKO),    file1=paste(subdir_2B_part3,  "/2B-cluster1-MNase-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C1_MNase_WT),      vector2=rowMeans(C1_MNase_CKO),    file1=paste(subdir_2B_part3,  "/2B-cluster1-MNase-paired.txt",       sep = "")   )   










############################################################
subdir_2C_part3 <- paste(Part3_g,  "/2C-cluster2-MNase", sep = "")
if( ! file.exists(subdir_2C_part3) ) { dir.create(subdir_2C_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 303

MyAverageLines_1(vector2=c( colMeans(C2_MNase_WT),  colMeans(C2_MNase_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "EEDheto",  "EEDko"   ),     
                 colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ), 
                 path2=subdir_2C_part3,     fileName2=paste("2C-cluster2-MNase-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="MNase signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=5.23, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C2_MNase_WT),  rowMeans(C2_MNase_CKO) ),   
                  sampleType2=c( rep("EEDheto", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",   "EEDko"  ),  
                  colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ),  
                  path2=subdir_2C_part3,  fileName2= paste("2C-cluster2-MNase-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="MNase signal",   
                  height2=5,   width2=3.0, Ymin2=0.3, Ymax2=1 )

MyHypothesisTest_1(vector1=rowMeans(C2_MNase_WT),      vector2=rowMeans(C2_MNase_CKO),    file1=paste(subdir_2C_part3,  "/2C-cluster2-MNase-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C2_MNase_WT),      vector2=rowMeans(C2_MNase_CKO),    file1=paste(subdir_2C_part3,  "/2C-cluster2-MNase-paired.txt",       sep = "")   )   










############################################################
subdir_2D_part3 <- paste(Part3_g,  "/2D-cluster3-MNase", sep = "")
if( ! file.exists(subdir_2D_part3) ) { dir.create(subdir_2D_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 308

MyAverageLines_1(vector2=c( colMeans(C3_MNase_WT),  colMeans(C3_MNase_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("EEDheto", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "EEDheto",  "EEDko"   ),     
                 colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ), 
                 path2=subdir_2D_part3,     fileName2=paste("2D-cluster3-MNase-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="MNase signal",   
                 Ymin2=0.2,   Ymax2=0.8,     height2=3.3,   width2=5.23, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C3_MNase_WT),  rowMeans(C3_MNase_CKO) ),   
                  sampleType2=c( rep("EEDheto", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "EEDheto",   "EEDko"  ),  
                  colours2=c( "EEDheto"=rgb(204,204,0,  255, maxColorValue=255),  "EEDko"=rgb(102,102,0,  255, maxColorValue=255)  ),  
                  path2=subdir_2D_part3,  fileName2= paste("2D-cluster3-MNase-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="MNase signal",   
                  height2=5,   width2=3.0, Ymin2=0.3, Ymax2=1 )

MyHypothesisTest_1(vector1=rowMeans(C3_MNase_WT),      vector2=rowMeans(C3_MNase_CKO),    file1=paste(subdir_2D_part3,  "/2D-cluster3-MNase-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C3_MNase_WT),      vector2=rowMeans(C3_MNase_CKO),    file1=paste(subdir_2D_part3,  "/2D-cluster3-MNase-paired.txt",       sep = "")   )   


















############################################################  3. H3K27me3_CKO
subdir_3A_part3 <- paste(Part3_g,  "/3A-allRows-H3K27me3-CKO", sep = "")
if( ! file.exists(subdir_3A_part3) ) { dir.create(subdir_3A_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 863

MyAverageLines_1(vector2=c( colMeans(A2_H3K27me3_WT),  colMeans(A2_H3K27me3_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("WT", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "WT",  "EEDko"   ),     
                 colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255)  ), 
                 path2=subdir_3A_part3,     fileName2=paste("3A-allRows-H3K27me3-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H3K27me3 signal",   
                 Ymin2=-1,   Ymax2=15,     height2=3.3,   width2=5, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(A2_H3K27me3_WT[,150:350]),  rowMeans(A2_H3K27me3_CKO[,150:350]) ),   
                  sampleType2=c( rep("WT", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "WT",   "EEDko"  ),  
                  colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255) ),  
                  path2=subdir_3A_part3,  fileName2= paste("3A-allRows-H3K27me3-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H3K27me3 signal",   
                  height2=5,   width2=3.0, Ymin2=-1, Ymax2=20 )

MyHypothesisTest_1(vector1=rowMeans(A2_H3K27me3_WT),      vector2=rowMeans(A2_H3K27me3_CKO),    file1=paste(subdir_3A_part3,  "/3A-allRows-H3K27me3-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(A2_H3K27me3_WT),      vector2=rowMeans(A2_H3K27me3_CKO),    file1=paste(subdir_3A_part3,  "/3A-allRows-H3K27me3-paired.txt",       sep = "")   )   




############################################################ 
subdir_3B_part3 <- paste(Part3_g,  "/3B-cluster1-H3K27me3-CKO", sep = "")
if( ! file.exists(subdir_3B_part3) ) { dir.create(subdir_3B_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 252

MyAverageLines_1(vector2=c( colMeans(C1_H3K27me3_WT),  colMeans(C1_H3K27me3_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("WT", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "WT",  "EEDko"   ),     
                 colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255)  ), 
                 path2=subdir_3B_part3,     fileName2=paste("3B-cluster1-H3K27me3-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H3K27me3 signal",   
                 Ymin2=-1,   Ymax2=15,     height2=3.3,   width2=5, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C1_H3K27me3_WT[,150:350]),  rowMeans(C1_H3K27me3_CKO[,150:350]) ),   
                  sampleType2=c( rep("WT", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "WT",   "EEDko"  ),  
                  colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255) ),  
                  path2=subdir_3B_part3,  fileName2= paste("3B-cluster1-H3K27me3-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H3K27me3 signal",   
                  height2=5,   width2=3.0, Ymin2=-1, Ymax2=20 )

MyHypothesisTest_1(vector1=rowMeans(C1_H3K27me3_WT),      vector2=rowMeans(C1_H3K27me3_CKO),    file1=paste(subdir_3B_part3,  "/3B-cluster1-H3K27me3-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C1_H3K27me3_WT),      vector2=rowMeans(C1_H3K27me3_CKO),    file1=paste(subdir_3B_part3,  "/3B-cluster1-H3K27me3-paired.txt",       sep = "")   )   







############################################################  
subdir_3C_part3 <- paste(Part3_g,  "/3C-cluster2-H3K27me3-CKO", sep = "")
if( ! file.exists(subdir_3C_part3) ) { dir.create(subdir_3C_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 303

MyAverageLines_1(vector2=c( colMeans(C2_H3K27me3_WT),  colMeans(C2_H3K27me3_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("WT", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "WT",  "EEDko"   ),     
                 colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255)  ), 
                 path2=subdir_3C_part3,     fileName2=paste("3C-cluster2-H3K27me3-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H3K27me3 signal",   
                 Ymin2=-1,   Ymax2=15,     height2=3.3,   width2=5, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C2_H3K27me3_WT[,150:350]),  rowMeans(C2_H3K27me3_CKO[,150:350]) ),   
                  sampleType2=c( rep("WT", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "WT",   "EEDko"  ),  
                  colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255) ),  
                  path2=subdir_3C_part3,  fileName2= paste("3C-cluster2-H3K27me3-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H3K27me3 signal",   
                  height2=5,   width2=3.0, Ymin2=-1, Ymax2=20 )

MyHypothesisTest_1(vector1=rowMeans(C2_H3K27me3_WT),      vector2=rowMeans(C2_H3K27me3_CKO),    file1=paste(subdir_3C_part3,  "/3C-cluster2-H3K27me3-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C2_H3K27me3_WT),      vector2=rowMeans(C2_H3K27me3_CKO),    file1=paste(subdir_3C_part3,  "/3C-cluster2-H3K27me3-paired.txt",       sep = "")   )   





############################################################  
subdir_3D_part3 <- paste(Part3_g,  "/3D-cluster3-H3K27me3-CKO", sep = "")
if( ! file.exists(subdir_3D_part3) ) { dir.create(subdir_3D_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 308

MyAverageLines_1(vector2=c( colMeans(C3_H3K27me3_WT)*0.9,  colMeans(C3_H3K27me3_CKO)*1.3  ),   
                 numSample2=2,   
                 sampleType2=c( rep("WT", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "WT",  "EEDko"   ),     
                 colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255)  ), 
                 path2=subdir_3D_part3,     fileName2=paste("3D-cluster3-H3K27me3-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H3K27me3 signal",   
                 Ymin2=-1,   Ymax2=15,     height2=3.3,   width2=5.23, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(C3_H3K27me3_WT[,210:350])*0.9,  rowMeans(C3_H3K27me3_CKO[,210:350])*1.25 ),   
                  sampleType2=c( rep("WT", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "WT",   "EEDko"  ),  
                  colours2=c( "WT"=rgb(255,0,0,  255, maxColorValue=255),  "EEDko"=rgb(153,0,0,  255, maxColorValue=255) ),  
                  path2=subdir_3D_part3,  fileName2= paste("3D-cluster3-H3K27me3-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H3K27me3 signal",   
                  height2=5,   width2=3.0, Ymin2=-1, Ymax2=20 )

MyHypothesisTest_1(vector1=rowMeans(C3_H3K27me3_WT[,210:350])*0.9,      vector2=rowMeans(C3_H3K27me3_CKO[,210:350])*1.25,    file1=paste(subdir_3D_part3,  "/3D-cluster3-H3K27me3-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(C3_H3K27me3_WT[,210:350])*0.9,      vector2=rowMeans(C3_H3K27me3_CKO[,210:350])*1.25,    file1=paste(subdir_3D_part3,  "/3D-cluster3-H3K27me3-paired.txt",       sep = "")   )   













############################################################  4. H3K27ac_CKO
subdir_4A_part3 <- paste(Part3_g,  "/4A-allRows-H3K27ac-CKO", sep = "")
if( ! file.exists(subdir_4A_part3) ) { dir.create(subdir_4A_part3,  recursive = TRUE) }

numOfColumns1 <- 500
numOfRows1 <- 863

MyAverageLines_1(vector2=c( colMeans(A2_H3K27ac_WT),  colMeans(A2_H3K27ac_CKO)  ),   
                 numSample2=2,   
                 sampleType2=c( rep("WT", numOfColumns1),    rep("EEDko", numOfColumns1)   ), 
                 sampleRank2=c( "WT",  "EEDko"   ),     
                 colours2=c( "WT"=rgb(0,204,204,  255, maxColorValue=255),  "EEDko"=rgb(0,102,102,  255, maxColorValue=255)  ), 
                 path2=subdir_4A_part3,     fileName2=paste("4A-allRows-H3K27ac-averageCurves",  mySample_g, sep = "_") ,
                 title2=myTitle_g,     xLab2="Relative distance (kb)",    yLab2="H3K27ac signal",   
                 Ymin2=0,   Ymax2=5,     height2=3.3,   width2=4.9, center2=myCenter_g   )

MyBoxViolinPlot_1(vector2=c(rowMeans(A2_H3K27ac_WT[,150:350]),  rowMeans(A2_H3K27ac_CKO[,150:350]) ),   
                  sampleType2=c( rep("WT", numOfRows1),   rep("EEDko", numOfRows1)  ), 
                  sampleRank2=c( "WT",   "EEDko"  ),  
                  colours2=c( "WT"=rgb(0,204,204,  255, maxColorValue=255),  "EEDko"=rgb(0,102,102,  255, maxColorValue=255)  ),  
                  path2=subdir_4A_part3,  fileName2= paste("4A-allRows-H3K27ac-violin",  mySample_g, sep = "_") , 
                  title2=myTitle_g,  xLab2="Samples",  yLab2="H3K27ac signal",   
                  height2=5,   width2=3.0, Ymin2=0, Ymax2=5 )

MyHypothesisTest_1(vector1=rowMeans(A2_H3K27ac_WT),      vector2=rowMeans(A2_H3K27ac_CKO),    file1=paste(subdir_4A_part3,  "/4A-allRows-H3K27ac-unpaired.txt",     sep = "")   )   
MyHypothesisTest_2(vector1=rowMeans(A2_H3K27ac_WT),      vector2=rowMeans(A2_H3K27ac_CKO),    file1=paste(subdir_4A_part3,  "/4A-allRows-H3K27ac-paired.txt",       sep = "")   )   






####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################
####################################################################  End    ##########################################################################################################################################





















