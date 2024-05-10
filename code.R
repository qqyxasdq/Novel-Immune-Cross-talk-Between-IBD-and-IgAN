####Data acquisition####
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
###
chooseBioCmirror()
gset = getGEO('GSE93798', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###
gset[[1]]

#grouping
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#Set reference level
group_list <- ifelse(str_detect(pdata$description, "Control"), "Control",
                     "IgAN")
#Factor type
group_list = factor(group_list,
                    levels = c("Control","IgAN"))
## The expression matrix was obtained and corrected
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###data correction
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
write.csv(exp ,file = "Expression matrix of GSE93798 after normalized.csv")




####GSVA####
library(readxl)
library(dplyr)
dat <- read.csv("Expression matrix of GSE66407  after normalized.csv") 

row.names(dat) <- dat$gene
dat <- dat[,-1]
head(dat)

library(GSVA)
geneSets <- getGmt('h.all.v7.5.1.symbols.gmt')    
GSVA_hall <- gsva(expr=as.matrix(dat), 
                  gset.idx.list=geneSets, 
                  mx.diff=T, 
                  kcdf="Gaussian", 
                  parallel.sz=4) #
head(GSVA_hall)
group <- factor(c(rep("IgAN", 20), rep("Control", 22)), levels = c('IgAN', 'Control'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design

compare <- makeContrasts(IgAN - Control, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)
## barplot
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)



####Merge datasets####

setwd("batch")

library(sva)
library(tidyverse)

load("Expression matrix of GSE93798  after normalized.Rda")
load("Expression matrix of GSE66409  after normalized.Rda")

merge_eset=inner_join(exprSet_GSE93798,exprSet_GSE66407,by="symbol")
rownames(merge_eset) <- merge_eset$symbol
merge_eset <- merge_eset[,-1]
dim(merge_eset)
exp <- as.matrix(merge_eset)
dimnames <- list(rownames(exp),colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
dim(data)
class(data)
batchType <- c(rep(1,42),rep(2,363))
modType <- c(rep("Control",22),rep("IgAN",20),rep("IBD",264),rep("Control",99))
mod  <-  model.matrix(~as.factor(modType))
outTab <- data.frame(ComBat(data, batchType,mod, par.prior=TRUE))

write.table(outTab,file="Expression files after removing batch effects for both diseases.csv",sep="\t",quote=F,col.names=F)

####Cibersort####
setwd("cibersort")   
install.packages('e1071')
install.packages('parallel')
#install.packages("BiocManager")

library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'Expression files after removing batch effects for both diseases.txt'   
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
 
write.csv(res_cibersort,file = "IBD IgAN datasets Cibersort algorithm.csv")


res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   
#
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) 
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, 
        names.arg = rep("",nrow(ciber.res)), 
        yaxt = "n", 
        ylab = "Relative percentage", 
        col = mycol) 
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   


####X cell#####

system("g++ -v")
system("where make")

chooseBioCmirror()

library(xCell)
library(ggpubr)
library(tidyverse)
setwd("Xcell")
exp <- read.table("Expression files after removing batch effects for both diseases.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# cell types
celltypeuse<-xCell.data$spill$K
rs<-xCellAnalysis(exp,parallel.sz=10) #
#

library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=xCell,value = Expression,-c(group,sample))


write.csv(b ,file = "IBD IgAN datasets  X cell  algorithm.csv")
ggboxplot(b, x = "xCell", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####ssGSEA####
## ssGSEA
setwd("ssGSEA")
BiocManager::install('GSVA')
library(tidyverse)
library(data.table)
library(GSVA)
#1.2 
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"
#
type <- split(cellMarker,cellMarker$celltype)

cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})

##
expr <- data.table::fread("Expression files after removing batch effects for both diseases.txt",data.table = F)   #读取表达文件
rownames(expr) <- expr[,1]   
expr <- expr[,-1]   
expr <- as.matrix(expr)   

#2. 
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")

a <- gsva_data %>% t() %>% as.data.frame()
identical(rownames(a),rownames(group))
a$group <- group$group
a <- a %>% rownames_to_column("sample")
write.table(a,"IBD IgAN datasets ssGSEA algorithm.csv",sep = "\t",row.names = T,col.names = NA,quote = F)
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))

ggboxplot(b, x = "ssGSEA", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()


####WGCNA####
library("tidyverse")
library("WGCNA")

library(DESeq2)
powerEstimates <- sft$powerEstimate # 
##
exp_mt = as.data.frame(t(expGSE193677))
exp_mt[1:4,1:4]
##(2)
gsg = goodSamplesGenes(exp_mt)
gsg$allOK


##(3)
sampleTree = hclust(dist(exp_mt), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 15000, col = "red") #

##
clust = cutreeStatic(sampleTree, cutHeight = 15000, minSize = 10)
table(clust)

keepSamples = (clust==1)
exp_mt_f = exp_mt[keepSamples, ]

##Sample clinical data were sorted out
trait_dat = read.csv("IBD IgAN datasets Cibersort algorithm.csv",row.names = 2) %>% 
  .[,setdiff(11:37,c(15,30))]
trait_dat_f = trait_dat[rownames(exp_mt_f),]
     6.96

identical(rownames(exp_mt_f), rownames(trait_dat_f))


##Summarizing the final data
exp_dat = exp_mt_f
dim(exp_dat)

exp_dat[1:4,1:4]

trait_dat = trait_dat_f
dim(trait_dat)

trait_dat[1:4,1:4]

##
#
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 
sft = pickSoftThreshold(exp_dat, powerVector = powers, verbose = 5)
##
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
###
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
### 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

par(mfrow = c(1,1))
#Establish the network and identify the module
##(1) According to the selected soft threshold, the adjacency matrix is obtained
softPower=sft$powerEstimate
adjacency = adjacency(exp_dat, power = 7)
adjacency[1:4,1:4]



##(2)
TOM = TOMsimilarity(adjacency)
TOM[1:4,1:4]

##(3)
dissTOM = 1-TOM
dissTOM[1:4,1:4]

##(4) hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")
##(5) Dynamic cutting tree, identification module
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)
# dynamicMods

## unassigned genes

##(6)Map module names to color names and visualize them
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors) #module0 
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and primary module colors")
#
MEList = moduleEigengenes(exp_dat, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.4, col = "red")  #
##(8)Based on the chosen threshold, the modules are merged
merge = mergeCloseModules(exp_dat, dynamicColors, cutHeight = 0.4)
#Color mapping of the new partition module
mergedColors = merge$colors
#Visualize the feature values of the newly divided modules before and after merging
mergedMEs = merge$newMEs
#Visualize the modules before and after merging
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#
net = blockwiseModules(exp_dat, power = 7,
                       corType = "pearson",
                       networkType="unsigned",
                       TOMType = "unsigned", 
                       minModuleSize = 30,
                       mergeCutHeight = 0.4,
                       verbose = 3)

names(net)

#(1) The resulting network module
unique(net$colors)
table(net$colors)
# 

##Merge the previous network modules
unique(net$unmergedColors)
table(net$unmergedColors)

#(2) 
dim(net$MEs)
net$MEs[1:4,1:4]
net$MEsOK

#(3)
plotDendroAndColors(dendro = net$dendrograms[[1]],
                    colors = net$colors)

# (1)Calculate module eigenvalues

# 
MEs = net$MEs
MEs = orderMEs(MEs)

# (2)
moduleTraitCor = cor(MEs, trait_dat, use = "p")
moduleTraitCor[1:2,1:2]


moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(MEs))
moduleTraitPvalue[1:2,1:2]


# (3)Visual correlation and P value

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(trait_dat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




##Correlation analysis of a specified phenotype calculated by visual blue module gene signature
weight = trait_dat[,"IgAN",drop=F]
colnames(weight) = "weight"
# (1) Gene significance，GS：
GS_weight = as.data.frame(cor(exp_dat, weight, use = "p"))
colnames(GS_weight) = "GS_weight"
head(GS_weight)


GS.p_weight = as.data.frame(corPvalueStudent(as.matrix(GS_weight), nrow(exp_dat)))



# 
modNames = substring(names(MEs), 3)
# 
MM = as.data.frame(cor(exp_dat, MEs, use = "p"))
colnames(MM) = paste("MM", modNames, sep="")
MM[1:4,1:4]

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nrow(exp_dat)))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nrow(exp_dat)));
colnames(MMPvalue) = paste("p.MM", modNames, sep="");
MMPvalue[1:4,1:4]

# (3) Visualized brown module gene characteristics
# identical(rownames(MM), names(net$colors))
# TRUE
module = "brown"
moduleGenes = names(net$colors)[net$colors=="brown"]
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(MM[moduleGenes, "MMbrown"]),
                   abs(GS_weight[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

for(mod in 1:nrow(table(net$colors)))
{
  modules=names(table(net$colors))[mod]
  probes=colnames(exp_dat)
  inModule=(net$colors==modules)
  eigengenes=probes[inModule]
  write.table(eigengenes,file=paste0(modules,".txt"),sep="\t",row.names = F,col.names = F,quote = F)
}  


##
datExpr <- exp_dat
datTraits <- trait_dat

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 
dissTOM <- 1-TOMsimilarityFromExpr(datExpr, power = 6);

# 
nSelect = 400

select <- sample(nGenes, size = nSelect);
selectTOM <- dissTOM[select, select];
# Reunion class
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- net$colors[select]
# visualization
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")











####Radar map####
library(fmsb)

exp <- read.table("The 10 crosstalk genes' expression in GSE93798.txt")
exam_scores <- exp[,comgene$gene]
#
exam_scores <- exp


max_min <- data.frame(
  MAPK3 = c(7.8, 7.0 ), RIPK2 = c(6.6, 5.7),NFKB1=c(8.4,7.5), HSDL2 = c(8.4, 7.4),
  EPHX2 = c(8.0, 6.7), RIDA = c(9.0, 7.2), FDX1 = c(8.1, 7.0),
  SYNPO = c(9.8, 9.3), KDF1 = c(4.6, 3.7), METTL7A= c(7.9, 6.6)
)


rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, exam_scores)
df
library(fmsb)
MAPK3_data <- df[c("Max", "Min", "MAPK3"), ]
radarchart(MAPK3_data)

radarchart(
  MAPK3_data, axistype = 1,
  # Customize the polygon
  pcol = "#00AFBB", pfcol = scales::alpha("#00AFBB", 0.5), plwd = 2, plty = 1,
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "grey", 
  # Variable labels
  vlcex = 0.7, vlabels = colnames(MAPK3_data),
  caxislabels = c(
    4.91e-40
    , 1.73e-23
    , 1.58e-11
    , 5.22e-05
    , 3.52e-14
  ))


radarchart(
  df, axistype = 1,
  # Customize the polygon
  pcol = c("#00AFBB", "#E7B800", "#FC4E07"), pfcol = scales::alpha(c("#00AFBB", "#E7B800", "#FC4E07"),0.5), plwd = 2, plty = 1,
  # Customize the grid
  cglcol = "grey", cglty = 1, cglwd = 0.8,
  # Customize the axis
  axislabcol = "grey", 
  # Variable labels
  vlcex = 0.7, vlabels = colnames(GSM1621630_data),
  caxislabels = c(0, 5, 10, 15, 20))
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)


library("ggradar")
library(tidyverse)
# Change the row names to group columns
df <- exam_scores %>% rownames_to_column("group")
df
ggradar(
  df[1, ], 
  values.radar = c("0", "10", "20"),# 
  grid.min = 0, # 
  grid.mid = 10, # 
  grid.max = 20 # 
)

ggradar(
  df[1, ], 
  values.radar = c("0", "10", "20"),
  grid.min = 0, grid.mid = 10, grid.max = 20,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = "#00AFBB",
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey"
)


ggradar(
  df[1:20,], 
  values.radar = c("0", "15", "20"),
  grid.min = 0, grid.mid = 10, grid.max = 20,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  
)







