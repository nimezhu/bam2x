# 
# contact: Wei Gufieng,<guifengwei@gmail.com>
#
print("########################################")
print("                                        ")
print("      Dependency: gplots                ")
print("                                        ")
print("########################################")

library(gplots)

# read the data into dataframe
BamCorr <- read.table("cor_matrix.txt",header=T,sep="\t")

row.names(BamCorr) <- BamCorr$title
# remove the first column
BamCorr <- BamCorr[,-1]
matrix<-as.matrix(BamCorr)

# draw the heatmap
cr <- colorRampPalette(colors = c("#00FF00","#FFFFFF","#FF0000"), bias=1)
heatmap.2(matrix,col=cr,trace="none",margin=c(12,12),keysize=1.2,density.info="none")

