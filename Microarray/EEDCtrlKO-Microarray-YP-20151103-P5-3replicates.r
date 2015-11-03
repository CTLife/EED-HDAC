

######################################################################################################################### （一）Normalize and get gene expression level
library(oligo)

celFiles <- list.celfiles('0-CEL-Files', full.names=TRUE)
celFiles
rawData <- read.celfiles(celFiles)
rawData
class(rawData)
dim(rawData)

Part1_g <- paste( "1-ExpressinoLevel-noFilter",     sep = "")
if( ! file.exists(Part1_g ) ) { dir.create(Part1_g ) }

sink( paste(Part1_g,  "/1-rawData-info.txt",     sep = "") )
rawData
sink()

pdf(file= paste(Part1_g,  "/2-images.pdf",     sep = "") ) 
image( rawData[,1], col = gray((64:0)/64) ) 
image( rawData[,2], col = gray((64:0)/64) ) 
image( rawData[,3], col = gray((64:0)/64) ) 
image( rawData[,4], col = gray((64:0)/64) ) 
image( rawData[,5], col = gray((64:0)/64) ) 
image( rawData[,6], col = gray((64:0)/64) ) 
dev.off()

svg(filename = paste(Part1_g,  "/3-hist.svg", sep = ""),  width = 6, height = 5 )
hist(rawData, col = rainbow(6), lty = 1, xlim = c(3, 15))
dev.off()

svg(filename = paste(Part1_g,  "/4-boxplot.svg", sep = ""),  width = 6, height = 5 )
boxplot(rawData, col = rainbow(6), names = 1:6)
dev.off()

pdf(file = paste(Part1_g,  "/5-MAplot.pdf", sep = "")  )
MAplot( rawData, arrays = 1, lowessPlot = TRUE, ylim = c(-1, 1) )
dev.off()



pmSeq <- pmSequence(rawData)
pmSeq[1:5]
write.table(pmSeq, file = paste(Part1_g,  "/6-pmSequence.txt", sep = "") ,  append = FALSE,  quote = TRUE,  sep = " ",
            eol = "\n",  na = "NA",  dec = ".",  row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


pmsLog2 <- log2(pm(rawData))
dim(pmsLog2)
counts <- Biostrings::alphabetFrequency(pmSeq, baseOnly = TRUE)
dim(counts)
counts[1:10, ]
GCcontent <- ordered(counts[, "G"] + counts[, "C"])
length(GCcontent)
xL <- "GC Frequency in 25-mers"
yL <- expression(log[2] ~ intensity)

pdf(file = paste(Part1_g,  "/7-GC-boxplot-1.pdf", sep = "")  )
boxplot(pmsLog2[,1] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()
pdf(file = paste(Part1_g,  "/7-GC-boxplot-2.pdf", sep = "")  )
boxplot(pmsLog2[,2] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()
pdf(file = paste(Part1_g,  "/7-GC-boxplot-3.pdf", sep = "")  )
boxplot(pmsLog2[,3] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()
pdf(file = paste(Part1_g,  "/7-GC-boxplot-4.pdf", sep = "")  )
boxplot(pmsLog2[,4] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()
pdf(file = paste(Part1_g,  "/7-GC-boxplot-5.pdf", sep = "")  )
boxplot(pmsLog2[,5] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()
pdf(file = paste(Part1_g,  "/7-GC-boxplot-6.pdf", sep = "")  )
boxplot(pmsLog2[,6] ~ GCcontent, xlab = xL, ylab = yL, range = 0, col =rainbow(10))
dev.off()




library(mogene10sttranscriptcluster.db)
library(RColorBrewer)

eset <- rma(rawData )
eset
write.exprs(eset, file=paste(Part1_g,  "/7-RMA.xls", sep = "")) 
annoPackage <- paste(annotation(eset) )
annoPackage   #"pd.mogene.1.0.st.v1"
affy_IDs <- featureNames(eset)
an <- "mogene10sttranscriptcluster"

symbols <- as.character( unlist( mget(affy_IDs, get(paste(an, "SYMBOL", sep="") ), ifnotfound=NA  ) ) )
gene_names <- as.character( unlist( mget(affy_IDs, get(paste(an, "GENENAME", sep="") ), ifnotfound=NA  ) ) )
entrez <- unlist( mget(affy_IDs, get(paste(an, "ENTREZID", sep="") ), ifnotfound=NA  ) ) 
unigenes <- as.character( unlist( lapply(mget(affy_IDs, get(paste(an, "UNIGENE", sep="") ),   ifnotfound=NA), paste, collapse="//"  ) ) )
refseqs <- as.character( unlist( lapply(mget(affy_IDs, get(paste(an, "REFSEQ", sep="") ),   ifnotfound=NA), paste, collapse="//"  ) ) )

expression_matrix <- exprs(eset)
dim(expression_matrix)
colnames(expression_matrix)
sink(paste(Part1_g,  "/8-summary-npFiltered-Expression-Level.txt", sep = ""))
summary(expression_matrix)
sink()


ctrlMeans <- rowMeans(2^expression_matrix[,1:3])
KOMeans <- rowMeans(2^expression_matrix[,4:6])
FC <- KOMeans/ctrlMeans
length(ctrlMeans)
length(KOMeans)
length(FC)

out <- data.frame(ProbeID=affy_IDs, Symbol=symbols, Name=gene_names, EntrezGene=entrez,  UniGene=unigenes, RefSeq=refseqs, expression_matrix, CtrlMean=ctrlMeans,  KOMean=KOMeans,  FoldChange=FC, stringsAsFactors = FALSE)
row.names(out) <- 1:length(affy_IDs)
write.table(out, paste(Part1_g,  "/9-ExpressionLevel_P30_RMA_exprs.xls", sep = "") , sep="\t", row.names=FALSE)

out1 <- out[!(is.na(out$Symbol)), ]
dim(out1)
write.table(out1, paste(Part1_g,  "/10-ExpressionLevel_P30_noNA.xls", sep = ""),   sep="\t",   row.names=FALSE)



svg(file=paste(Part1_g,  "/11-normalize-and-unnormalize.svg", sep = "") , width=10, height=5) 
op <- par(mfrow=c(1,2))     #1 row, 2 columns
cols <- brewer.pal(6, "Set3")    #表示使用Set3的10种颜色
boxplot(rawData, col=cols, names=1:6, main="Unnormalized Data", xlab="Samples", ylab="Expression Level")
boxplot( data.frame(expression_matrix) ,  names=1:6, main="Normalized Data", xlab="Samples",  ylab="Expression Level", col="blue", border="brown" )
par(op)
dev.off()



distance1 <- dist(t(expression_matrix), method="euclidean")
clusters <- hclust(distance1)
svg(file=paste(Part1_g,  "/12-cluster1-euclidean.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance2 <- dist(t(expression_matrix), method="maximum")
clusters <- hclust(distance2)
svg(file=paste(Part1_g,  "/12-cluster2-max.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance3 <- dist(t(expression_matrix), method="manhattan")
clusters <- hclust(distance3)
svg(file=paste(Part1_g,  "/12-cluster3-manhattan.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance4 <- dist(t(expression_matrix), method="canberra")
clusters <- hclust(distance4)
svg(file=paste(Part1_g,  "/12-cluster4-canberra.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance5 <- dist(t(expression_matrix), method="binary")
clusters <- hclust(distance5)
svg(file=paste(Part1_g,  "/12-cluster5-binary.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance6 <- dist(t(expression_matrix), method="minkowski")
clusters <- hclust(distance6)
svg(file=paste(Part1_g,  "/12-cluster5-minkowski.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()








################################################################################################# (二) 数据过滤 Quality control checks and Filtering data
# 原始数据读入，经AffyBatch目标转成ExpressionSet目标后，为提高后续分析（如差异表达基因的检测）的统计功效，
# 往往需要进一步经过Detection Call Filter和IQR filter等过滤（“基因芯片数据的特异性过滤与非特异性过滤”）。

Part2_g <- paste( "2-ExpressinoLevel-Filtered",     sep = "")
if( ! file.exists(Part2_g ) ) { dir.create(Part2_g ) }

library(genefilter)
dim(eset)
eset2 <-  varFilter(eset, var.func=IQR, var.cutoff=0.2, filterByQuantile=TRUE)
dim(eset2)
write.exprs(eset2, file=paste(Part2_g,  "/1-filtered-RMA-varFilter.xls", sep = "")) 


affy_IDs2 <- featureNames(eset2)
an2 <- "mogene10sttranscriptcluster"
symbols2 <- as.character( unlist( mget(affy_IDs2, get(paste(an2, "SYMBOL", sep="") ), ifnotfound=NA  ) ) )
gene_names2 <- as.character( unlist( mget(affy_IDs2, get(paste(an2, "GENENAME", sep="") ), ifnotfound=NA  ) ) )
entrez2 <- unlist( mget(affy_IDs2, get(paste(an2, "ENTREZID", sep="") ), ifnotfound=NA  ) ) 
unigenes2 <- as.character( unlist( lapply(mget(affy_IDs2, get(paste(an2, "UNIGENE", sep="") ),   ifnotfound=NA), paste, collapse="//"  ) ) )
refseqs2 <- as.character( unlist( lapply(mget(affy_IDs2, get(paste(an2, "REFSEQ", sep="") ),   ifnotfound=NA), paste, collapse="//"  ) ) )

expMatrix3 <- exprs(eset2)
expMatrix3
dim(expMatrix3)

ctrlMeans3 <- rowMeans(2^expMatrix3[,1:3])
KOMeans3 <- rowMeans(2^expMatrix3[,4:6])
FC3 <- KOMeans3/ctrlMeans3
length(ctrlMeans3)
length(KOMeans3)
length(FC3)

sink(paste(Part2_g,  "/2-summary-filtered-Expression-Level.txt", sep = ""))
summary(expMatrix3)
sink()

out3 <- data.frame(ProbeID=affy_IDs2, Symbol=symbols2, Name=gene_names2, EntrezGene=entrez2,  UniGene=unigenes2, RefSeq=refseqs2, expMatrix3, CtrlMean=ctrlMeans3,  KOMean=KOMeans3,  FoldChange=FC3, stringsAsFactors = FALSE)
row.names(out3) <- 1:length(affy_IDs2)
write.table(out3, paste(Part2_g,  "/3-filtered-all-Expression-Level.xls", sep = ""),  sep="\t",  row.names=FALSE)

out4 <- out3[!(is.na(out3$Symbol)), ]
dim(out4)
write.table(out4, paste(Part2_g,  "/3-filtered-noNA-Expression-Level.xls", sep = ""),   sep="\t",   row.names=FALSE)




distance1 <- dist(t(expMatrix3), method="euclidean")
clusters <- hclust(distance1)
svg(file=paste(Part2_g,  "/4-cluster1-euclidean.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance2 <- dist(t(expMatrix3), method="maximum")
clusters <- hclust(distance2)
svg(file=paste(Part2_g,  "/4-cluster2-max.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance3 <- dist(t(expMatrix3), method="manhattan")
clusters <- hclust(distance3)
svg(file=paste(Part2_g,  "/4-cluster3-manhattan.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance4 <- dist(t(expMatrix3), method="canberra")
clusters <- hclust(distance4)
svg(file=paste(Part2_g,  "/4-cluster4-canberra.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance5 <- dist(t(expMatrix3), method="binary")
clusters <- hclust(distance5)
svg(file=paste(Part2_g,  "/4-cluster5-binary.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()


distance6 <- dist(t(expMatrix3), method="minkowski")
clusters <- hclust(distance6)
svg(file=paste(Part2_g,  "/4-cluster5-minkowski.svg", sep = ""), width=5, height=5) 
plot(clusters)
dev.off()








######################################################################################################## （三）聚类分析 (找到差异表达基因以前，找到差异表达基因以后还会做一次。)
# 常规做法是先筛选出差异表达基因，然后只用差异表达基因进行聚类分析

Part3_g <- paste( "3-cluster",     sep = "")
if( ! file.exists(Part3_g ) ) { dir.create(Part3_g ) }

library(latticeExtra)
library(gplots) 


#（1）样本聚类
dd <- dist2(log2(exprs(eset)))
diag(dd) <- 0
dd.row <- as.dendrogram(hclust(as.dist(dd)))
row.ord <- order.dendrogram(dd.row)
legend <- list(top = list(fun = dendrogramGrob, args = list(x = dd.row, side = "top")))
lp <- levelplot(dd[row.ord, row.ord], scales = list(x = list(rot = 90)), xlab = "", ylab = "", legend = legend)
svg(file=paste(Part3_g,  "/1-cluster-noFilter.svg", sep = ""), width=12, height=6) 
plot(lp)
dev.off()

dd2 <- dist2(log2(exprs(eset2)))
diag(dd2) <- 0
dd2.row <- as.dendrogram(hclust(as.dist(dd2)))
row.ord <- order.dendrogram(dd2.row)
legend <- list(top = list(fun = dendrogramGrob, args = list(x = dd2.row, side = "top")))
lp <- levelplot(dd2[row.ord, row.ord], scales = list(x = list(rot = 90)), xlab = "", ylab = "", legend = legend)
svg(file=paste(Part3_g,  "/1-cluster-Filtered.svg", sep = ""), width=12, height=6) 
plot(lp)
dev.off()




# （2）二维聚类
my.colorFct <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
  if (n < 2) stop("n must be greater than 2")
  n1 <- n%/%2
  n2 <- n - n1
  c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}
mydata<-exprs(eset)
mydatascale <- t(scale(t(mydata)))
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete")  
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete") 
svg(file=paste(Part3_g,  "/2-cluster-noFilter.svg", sep = ""), width=12, height=6) 
heatmap.2(mydata, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=redgreen(75), scale="row", trace="none", key=T)
dev.off()


#  上述聚类图一般和论文里的聚类图有点不同，聚类的模式不太直观，也可以用下面的语句进行更直观的作图：
mycl <- cutree(hr, h=max(hr$height)/1.5);
mycolhc <- sample(rainbow(256)); mycolhc <- mycolhc[as.vector(mycl)]
myc2 <- cutree(hc, h=max(hc$height)/1.5); mycolhr <- sample(rainbow(256)); mycolhr <- mycolhr[as.vector(myc2)]
svg(file=paste(Part3_g,  "/3-heatmap-cluster-noFilter.svg", sep = ""), width=12, height=6)
heatmap(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", ColSideColors=mycolhr, RowSideColors=mycolhc, margins = c(20,0))
dev.off()






mydata2<-exprs(eset2)
mydatascale2 <- t(scale(t(mydata2)))
hr <- hclust(as.dist(1-cor(t(mydatascale2), method="pearson")), method="complete")  
hc <- hclust(as.dist(1-cor(mydatascale2, method="spearman")), method="complete") 
svg(file=paste(Part3_g,  "/2-cluster-Filtered.svg", sep = ""), width=12, height=6) 
heatmap.2(mydata, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=redgreen(75), scale="row", trace="none", key=T)
dev.off()

#  上述聚类图一般和论文里的聚类图有点不同，聚类的模式不太直观，也可以用下面的语句进行更直观的作图：
mycl <- cutree(hr, h=max(hr$height)/1.5);
mycolhc <- sample(rainbow(256)); mycolhc <- mycolhc[as.vector(mycl)]
myc2 <- cutree(hc, h=max(hc$height)/1.5); mycolhr <- sample(rainbow(256)); mycolhr <- mycolhr[as.vector(myc2)]
svg(file=paste(Part3_g,  "/3-heatmap-cluster-Filtered.svg", sep = ""), width=12, height=6)
heatmap(mydatascale2, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", ColSideColors=mycolhr, RowSideColors=mycolhc, margins = c(20,0))
dev.off()







library(affycoretools)
library(rgl)
library(plot3D)

svg(file=paste(Part3_g,  "/4-2D-PCA1-noFilter.svg", sep = ""), width=4, height=4.5)
affycoretools::plotPCA(eset, groups = NULL,   groupnames = NULL,   
                       addtext = c(  "Ctrl3",    "Ctrl4",   "Ctrl5",     "EEDKO3",  "EEDKO4",   "EEDKO5"),
                       x.coord = "right",   y.coord = "center",    screeplot = FALSE,   squarepca = FALSE,
                       pch = c(19, 19,  19, 19, 19, 19 ), 
                       col = c("purple",  "red", "red4",  "yellowgreen",  "green",   "skyblue4"), 
                       pcs = c(1, 2), legend = FALSE,
                       main = "Principal Components Plot", plot3d = FALSE, outside = FALSE )
dev.off()




pFigure <- affycoretools::plotPCA(eset, groups = NULL,   groupnames = NULL,   
                       addtext = c(  "Ctrl3",    "Ctrl4",   "Ctrl5",      "EEDKO3",  "EEDKO4",   "EEDKO5"),
                       x.coord = "right",   y.coord = "center",    screeplot = FALSE,   squarepca = FALSE,
                       pch = c(19, 19,  19, 19, 19, 19 ), 
                       col = c("red", "red", "red",   "green",    "green", "green"), 
                       pcs = c(1, 2, 3), legend = TRUE,
                       main = "Principal Components Plot", plot3d = TRUE, outside = FALSE )
print(pFigure)
rgl.postscript(paste(Part3_g,  "/4-3D-PCA2-noFilter.ps", sep = ""),   "ps") 










svg(file=paste(Part3_g,  "/5-2D-PCA1-Filtered.svg", sep = ""), width=4, height=4.5)
affycoretools::plotPCA(eset2, groups = NULL,   groupnames = NULL,   
                       addtext = c(  "Ctrl3",    "Ctrl4",   "Ctrl5",      "EEDKO3",  "EEDKO4",   "EEDKO5"),
                       x.coord = "right",   y.coord = "center",    screeplot = FALSE,   squarepca = FALSE,
                       pch = c(19, 19,  19, 19,    19, 19), 
                       col = c("purple",  "red", "red4",  "yellowgreen",  "green",   "skyblue4"), 
                       pcs = c(1, 2), legend = FALSE,
                       main = "Principal Components Plot", plot3d = FALSE, outside = FALSE )
dev.off()




pFigure <- affycoretools::plotPCA(eset2, groups = NULL,   groupnames = NULL,   
                                  addtext = c(  "Ctrl3",    "Ctrl4",   "Ctrl5",    "EEDKO3",  "EEDKO4",   "EEDKO5"),
                                  x.coord = "right",   y.coord = "center",    screeplot = FALSE,   squarepca = FALSE,
                                  pch = c(19, 19,  19, 19,  19, 19), 
                                  col = c("red", "red",  "red",  "green",    "green", "green"), 
                                  pcs = c(1, 2, 3), legend = TRUE,
                                  main = "Principal Components Plot", plot3d = TRUE, outside = FALSE )
print(pFigure)
rgl.postscript(paste(Part3_g,  "/5-3D-PCA2-Filtered.ps", sep = ""),   "ps") 










############################################################################################### (四)差异表达基因的分析  Finding differentially expressed probesets
Part4_g <- paste( "4-DiffGenes",     sep = "")
if( ! file.exists(Part4_g ) ) { dir.create(Part4_g ) }

library(limma)  # two-color pre-processing; differential expression

samples <- as.factor( c("Ctrl",     "Ctrl", "Ctrl",       "EEDKO",  "EEDKO",   "EEDKO") )       
design <- model.matrix(~ -1 + samples) 
colnames(design) <- c("Ctrl", "EEDKO")  
design
contrast.matrix <- makeContrasts(contrasts = "EEDKO - Ctrl",  levels=design)
fit1 <- lmFit(exprs(eset), design)
head(fit1$coefficients)

fit2  <- contrasts.fit(fit1, contrast.matrix)
fit3 <- eBayes(fit2)
probeset.list <- topTable(fit3, number=1000000)
probes=row.names(probeset.list) 
Symbols = unlist(mget(probes, mogene10sttranscriptclusterSYMBOL, ifnotfound=NA))
length(Symbols)
results1 = cbind(probes,Symbols, probeset.list)

dim(results1)  
write.table(results1, file = paste(Part4_g,  "/1-nofilter-all-diffGeneExpression_RMA_eBayes.xls", sep = "") ,   quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = TRUE)

results2 <- results1[!(is.na(results1[,2])), ]   
dim(results2)  
write.table(results2, file = paste(Part4_g,  "/1-nofilter-noNA-diffGeneExpression_RMA_eBayes.xls", sep = "") ,   quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = TRUE)




samples <- as.factor( c("Ctrl",     "Ctrl",       "Ctrl",      "EEDKO",  "EEDKO",   "EEDKO") )       
design <- model.matrix(~ -1 + samples) 
colnames(design) <- c("Ctrl", "EEDKO")  
design
contrast.matrix <- makeContrasts(contrasts = "EEDKO - Ctrl",  levels=design)
fit1 <- lmFit(exprs(eset2), design)
head(fit1$coefficients)

fit2  <- contrasts.fit(fit1, contrast.matrix)
fit3 <- eBayes(fit2)
probeset.list <- topTable(fit3, number=1000000)
probes=row.names(probeset.list) 
Symbols = unlist(mget(probes, mogene10sttranscriptclusterSYMBOL, ifnotfound=NA))
length(Symbols)
results1 = cbind(probes,Symbols, probeset.list)

dim(results1)  
write.table(results1, file = paste(Part4_g,  "/2-filtered-all-diffGeneExpression_RMA_eBayes.xls", sep = ""),   quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = TRUE)

results2 <- results1[!(is.na(results1[,2])), ]   
dim(results2)  
write.table(results2, file = paste(Part4_g,  "/2-filtered-noNA-diffGeneExpression_RMA_eBayes.xls", sep = ""),   quote = FALSE,   sep = "\t",   row.names = FALSE,   col.names = TRUE)



sink(paste(Part4_g,  "/3-version-number.txt", sep = ""))
print(.libPaths())
cat("################################\n\n\n")
print(sessionInfo())
cat("################################\n\n\n")
print(version)
cat("################################\n\n\n")
packinfo <- installed.packages (fields = c ("Package", "Version"))
packinfo[,c("Package", "Version")]
cat("################################\n\n\n")
package_version(R.version)
cat("################################\n\n\n")
getRversion()
sink() 













