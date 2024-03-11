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



####

# 
MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes

# Add module characteristics
MET <- orderMEs(cbind(MEs, weight))
# visualization
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "",
                      marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = 0.8, xLabelsAngle
                      = 90)

####Heat map####

library(pheatmap)
group <- read.csv("The 40 common genes' expression in GSE66407.csv")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(exp_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,
         fontsize_col=3)
pheatmap(c,annotation = annotation,
         annotation_colors = list(group = c("1" ="#01468b","2"= "#ee0000")),
         cluster_cols = F,fontsize=5,fontsize_row=5,
         scale="row",show_colnames=F,cluster_row = F,
         fontsize_col=3)

dev.off()



####LASSO####

library(glmnet)
library(readxl)
library(plyr)
library(caret)
library(corrplot)
library(ggplot2)
library(Hmisc)
library(openxlsx)
library(broom)

data <- read.csv("The 40 common genes' expression in GSE93798.csv")
x <- as.matrix(data[,-41])
y <- as.matrix(data$group)
fit <- glmnet(x,y,family="binomial",nlambda = 1000,alpha=1)
print(fit)
plot(fit,xvar="lambda")
lasso_fit <- cv.glmnet(x,y,family="binomial",alpha=1,type.measure="auc",nlambda=1000)
plot(lasso_fit)
print(lasso_fit)

lasso_best <- glmnet(x=x,y=y,alpha=1,lambda = lasso_fit$lambda.min)
coef(lasso_best)


cvfit=cv.glmnet(x,y,type.measure = "mse",nfolds = 5,alpha=1)
plot(cvfit)
cvfit$lambda.min

lasso <- glmnet(x,y,family="binomial",alpha=1)
print(lasso)
plot(lasso,label=TRUE)
plot(lasso,xvar="lambda",label=TRUE)
lasso.coef <- predict(lasso,s=0.01864089,type="coefficients")
lasso.coef
plot(lasso,xvar="dev",label=TRUE)
lasso.y <- predict(lasso,newx=x,type = "response",s=0.005)
plot(lasso.y,y,xlab="Predicted",ylab="Actural",main="Lasso REgression")



####RFE####
library(caret)
subsets <- generateTestVariableSet(ncol(train_data))
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, repeats=5)

train_data <- read.csv("The 40 common genes' expression in GSE93798.csv")
rfe <- rfe(x=train_data[,(1:40)], y=train_data$group, size=subsets, rfeControl=control)
print(rfe, top=16)
plot(rfe, type=c("g", "o"))
#
caretRfe_variables <- data.frame(Item=rfe$optVariables, Type="Caret_RFE")
vairables <- rbind(boruta.finalVars, boruta.finalVarsWithTentative, caretRfe_variables)
library(VennDiagram)
library(ImageGP)
sp_vennDiagram2(vairables, item_variable = "Item", set_variable = "Type", manual_color_vector ="Set1")




####Model selection####
library(breakDown)
wine <- read.csv("The 40 common genes' expression in GSE93798.csv")

wine$group <- factor(wine$group)
trainIndex <- createDataPartition(wine$group,p=0.6,list=FALSE,times=1)
wineTrain <- wine[trainIndex,]
wineTest <- wine[-trainIndex,]
classif_rf <- train(group~,data=wineTrain,method="rf",ntree=100,tuneLength=1)
classif_glm <- train(group~,data=wineTrain,method="glm",family="binomial")
classif_svm <- train(group~,data=wineTrain,method="svmRadial",prob.model=TRUE,tuneLength=1)

explainer_classif_rf <- DALEX::explain(classif_rf,data = wineTest, y = yTest,predict_function = p_fun)

explainer_classif_glm <- DALEX::explain(classif_glm, label = "glm", 
                                        data = wineTest, y = yTest,
                                        predict_function = p_fun)
explainer_classif_svm <- DALEX::explain(classif_svm,  label = "svm", 
                                        data = wineTest, y = yTest,
                                        predict_function = p_fun)

mp_classif_rf <- model_performance(explainer_classif_rf)
mp_classif_glm<- model_performance(explainer_classif_glm)
mp_classif_svm <- model_performance(explainer_classif_svm)
plot(mp_classif_rf, mp_classif_glm, mp_classif_svm)
plot(mp_classif_rf, mp_classif_glm, mp_classif_svm, geom = "boxplot")
vi_classif_rf <- variable_importance(explainer_classif_rf, loss_function = loss_root_mean_square)
vi_classif_glm <- variable_importance(explainer_classif_glm, loss_function = loss_root_mean_square)
vi_classif_svm <- variable_importance(explainer_classif_svm, loss_function = loss_root_mean_square)
plot(vi_classif_rf, vi_classif_glm, vi_classif_svm)
pdp_classif_rf  <- variable_response(explainer_classif_rf, variable = "pH", type = "pdp")     



##########Boruta########
# install.packages("Boruta")
library(Boruta)
exp1 <- read.csv("The 40 common genes' expression in GSE93798.csv")
boruta <- Boruta(x=exp1[,(1:40)], y=exp1$group, pValue=0.01, mcAdj=T, 
                 maxRuns=300)
boruta



table(boruta$finalDecision)

boruta$finalDecision[boruta$finalDecision=="Confirmed"]

Boruta::plotImpHistory(boruta)


library(dplyr)

boruta.variable.imp <- boruta.imp(boruta)

head(boruta.variable.imp)



library(YSX)

sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
           legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
           xtics_angle = 90)

boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")
caret::featurePlot(exp1[,boruta.finalVarsWithTentative$Item], exp1$group, plot="box")


#
boruta_train_data <- exp1[, boruta.finalVarsWithTentative$Item]
boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
#
library(caret)
# Create model with default parameters
trControl <- trainControl(method="repeatedcv", number=10, repeats=5)


# 
tuneGrid <- expand.grid(mtry=boruta_mtry)

borutaConfirmed_rf_default <- train(x=boruta_train_data, y=exp1$group, method="rf", 
                                    tuneGrid = tuneGrid, # 
                                    metric="Accuracy", #metric='Kappa'
                                     trControl=trControl)
borutaConfirmed_rf_default

plot(borutaConfirmed_rf_default)
dotPlot(varImp(borutaConfirmed_rf_default))




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




####RF####
library(foreach)
library(randomForest)
wineTest <- read.csv("The 40 common genes' expression in GSE93798.csv")

trains <- createDataPartition(
  y=wineTest$group,
  p=0.75,
  list=F
)
traindata <- wineTest[trains,]
testdata <- wineTest[-trains,]
colnames(wineTest)
form_cls <- as.formula(
  paste0(
    "group~",
    paste(colnames(traindata)[1:10],collapse = "+")
  )
)
form_cls
traindata$group <- as.factor(traindata$group)
testdata$group <- as.factor(testdata$group)
##model building###

fit_rf_cls <- randomForest(
  form_cls,
  data=traindata,
  htree=500,
  mtry=6,
  importance=T
)


##
plot(fit_rf_cls,main = "ERROR&TREES")
legend("top",
       legend = colnames(fit_rf_cls$err.rate),      
       lty=1:3,
       col = 1:3,
       horiz = T)

###
importance(fit_rf_cls)

varImpPlot(fit_rf_cls,main="VarImpPlot")
varImpPlot(fit_rf_cls,main="VarImpPlot",type=1)


###
partialPlot(x=fit_rf_cls,
            pred.data = traindata,
            x.var=ChestPain,
            which.class = "Yes",
            ylab="Yes")



trainpredprob <- predict(fit_rf_cls,newdata = traindata,type="prob")
##
trainroc <- roc(response=traindata$group,
                predictor=trainpredprob[,2])
##
plot(trainroc,
     print.auc=TRUE,
     auc.polygon=TRUE,
     grid=T,
     max.auc.polygon.col="skyblue",
     print.thres=T,
     legacy.axes=T,
     bty="l")


##
bestp <- trainroc$thresholds[
  which.max(trainroc$sensitivities+trainroc$specificities-1)
]
###
trainpredlab <- as.factor(
  ifelse(trainpredprob[,2]>bestp,"1","0"))

confusionMatrix(data=trainpredlab,
                reference = traindata$group,
                mode="everything")



##
testpredprob <- predict(fit_rf_cls,newdata = testdata,type="prob")
testpredlab <- as.factor(
  ifelse(testpredprob[,2]>bestp,"1","0")
)

confusionMatrix(data=testpredlab,
                reference = testdata$group,
                positive = "1",
                mode="everything")

testroc <- roc(response=testdata$group,
               predictor=testpredprob[,2])

plot(testsetroc,
     print.auc=TRUE,
     grid=c(0.1,0.2),
     auc.polygon=F,
     max.auc.polygon=T,
     main="RF--ROC",
     grif.col=c("green","red"))
plot(testroc,
     print.auc=TRUE,
     print.auc.y=0.4,
     add=T,
     col="red")
legend("bottomright",
       legend=c("traindata","testdata"),
       col=c(par("fg"),"red"),
       lwd=2,
       cex=0.9)


####Lollipop Chart####

dat = read.csv("The 10 corsstalk genes and cibersort algorithm correlation.csv")
head(dat)
# 
dat$cor <- cut(abs(dat$cor),# 
                breaks = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                labels = c("< 0.3","0.3 - 0.5","0.5 - 0.7","0.7 - 0.9","> 0.9"),
                right=FALSE) 
dat$pvalue <- cut(dat$pvalue,
                   breaks = c(0, 0.001, 0.01, 0.05, 1),
                   labels = c("< 0.001","< 0.01","< 0.05","> 0.05"),
                   right=FALSE) 
# 
dat = dat[order(dat$cor),]
dat$Cell = factor(dat$Cell, levels = dat$Cell)
p = ggplot(dat, aes(x = cor, y = Cell, color = pvalue1)) +
  scale_color_manual(name="pvalue",
                     values = c("#E69F00", "#56B4E9", "#009E73", "gray"))+
  geom_segment(aes(x = 0, y = Cell, xend = cor, yend = Cell),size = 1) +
  geom_point(aes(size = cor1))+
  theme_bw()+
  labs(size = "Cor")
p

## 
dat$pvalue <- cut(dat$pvalue,
                   breaks = c(0, 0.05,1),
                   labels = c("< 0.05","> 0.05"),
                   right=FALSE) 
p1 = ggplot()+
  geom_text(dat,mapping = aes(x = 0, y = Cell, color = pvalue2, 
                              label = round(pvalue,3)))+
  scale_color_manual(name="",
                     values = c("red", "black"))+
  theme_void()+
  guides(color=F)

p1
library(patchwork)
p|p1

ggsave("NFKB1 and 22 immune cell correlations.pdf",width = 8,height = 5)





####Correlation heat map####
setwd("cor")
install.packages("corrplot")
library(corrplot)
library(tidyverse)
expr <- read.table("Expression files after removing batch effects for both diseases.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gene <- c("MAPK3","RIPK2","NFKB1","HSDL2",'EPHX2',
          'RIDA','FDX1','SYNPO','KDF1',"METTL7A")
exp <- expr[gene,]
exp <- exp %>% t() %>% as.data.frame()
ciber <- read.table("IBD IgAN datasets Cibersort algorithm.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber[,1:22]
identical(rownames(ciber),rownames(exp))
class(exp$MAPK3)
class(ciber$B.cells.naive)
cor<-sapply(ciber,function(x,y) cor(x,y,method="spearman"),exp)
rownames(cor)<-colnames(exp)
cor_res <- cor.mtest(cor,#
                     conf.level = 0.95)#
corrplot(cor,
         method = "color",#
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         addCoef.col = "black",#
         tl.col="black",#
         number.cex = 0.5,
         tl.cex = 0.7,
         cl.align = "l")
dev.off()



####GO enrichment####
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

##
rt=read.table("Expression files after removing batch effects for both diseases.txt",sep="\t",check.names=F,header=T)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="Enrichment analysis of 40 common genes.csv",sep="\t",quote=F,row.names=F)

ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db, 
                pvalueCutoff =0.05, 
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)

##
geneList <- rt$avg_log2FC
names(geneList) <- rt$entrezID
cnetplot(ego, 
         foldChange = geneList, 
         #foldChange = NULL, #
         circular = TRUE,
         #node_label = FALSE, #
         showCategory = 4, #
         colorEdge = TRUE)
ggsave("clusterProfiler_circle.pdf", width = 13, height = 11)


####mantel test####

library(dplyr)
library(linkET)
library(ggplot2)

env <- read.csv("The 10 crosstalk genes.csv")
speciese <- read.csv("IBD IgAN datasets Cibersort algorithm.csv")



mantel01 <- mantel_test(env, speciese,
                        spec_select = list(ITPR1=1,FAM13A=2,DSG2=3,ANP32B=4,XPO1=5,HNRNPD=6,PROS1=7,ITGA6=8,TLN2=9,MTRNR2L8=10,NCL=11)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


qcorrplot(correlate(speciese), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = mantel01, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  labs(title = "Mantel test")



plantcharacter <- ssGSEA
qcorrplot(correlate(cibersort., ssGSEA)) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), data = mantel01, curvature = 0.1) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  expand_limits(x = 20) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))



####single cell ssGSEA#####

surv <- read.csv("Expression matrix of GSE93798 after normalized.csv")
group_list <- ifelse(str_detect(surv$group, "Low"), "Low",
                     "High")
#
group_list = factor(group_list,
                    levels = c("Low","High"))
group_list

#
library(limma)
design=model.matrix(~group_list)
exp <- surv[,-17207]
exp <- as.data.frame(t(exp))
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)


####GSEA####
setwd("GSEA")

library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
deg <- read.csv("Expression matrix of GSE66407 after normalized.csv")
DEG <- as.data.frame(deg)%>% 
  arrange(P.Value) %>% 
  dplyr::filter(abs(logFC) > 0.5, P.Value < 0.05)

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))


msigdb_GMTs <- "msigdb_v7.0_GMTs"
msigdb <- "c5.all.v7.5.1.entrez.gmt"  

kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))

geneList = DEG[,2]#
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)


KEGG<-GSEA(geneList,TERM2GENE = kegmt) 
#
KEGG_result_df <- as.data.frame(KEGG)


#
library(enrichplot)
gseaplot2(KEGG,1,color="red")
gseaplot2(KEGG,3,color="red",pvalue_table = T)


gseaplot2(KEGG, geneSetID = c(8:11,51), subplots = 1:3)
gseaplot2(KEGG, geneSetID = c(3:7), subplots = 1:3)
gseaplot2(KEGG, geneSetID = 1:4, subplots = 1)
gseaplot2(KEGG, geneSetID = 1:10, subplots = 1:3)

dev.off()



####Box plot####

library(tidyverse)
library(ggsci)
library(tidyr)
library(ggpubr)
a <- read.table("The 10 crosstalk genes' expression in GSE66407.txt", sep = "\t",row.names = 1,check.names = F,header = T)
a <- a %>% rownames_to_column("genes") 
group=a[,19]  #
group= as.matrix(group)
b <- gather(a,key=m6A, value = Value,-c(group,m6A))
pdf("boxplot.pdf",height=8,width=15)   
ggboxplot(b, x = "m6A", y = "Value",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "t.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values=c("#56B4E9", "#E69F00"))+
  theme_bw()
dev.off()

####Mountain map####

library(tidyverse)
data <- wineTest
data2 <- data %>% 
  
  add_column(class = rep(c("Normal","Cancer"), c(2029,461)), .before = 1) %>% 
  
  gather("cell_types","percentages", -class) %>% 
  
  mutate(cell_types = factor(cell_types, levels = colnames(data))) %>% 
  mutate(class = factor(class, levels = c("Normal","Cancer")))


library(ggridges)
library(ggplot2)
ggplot(data2, aes(x=percentages, y= cell_types, fill=..x..))+
  geom_density_ridges_gradient(scale=1.5, rel_min_height=0.01, gradient_lwd = 1)+
  scale_fill_gradient2(low = "blue", mid = "orange", high = "red", midpoint = 0.3)+
  scale_x_continuous(expand = c(0.01, 0))+
  scale_y_discrete(expand = c(0.01,0))+
  labs(title="")+ 
  theme_ridges(font_size = 12, grid = FALSE)+
  theme(axis.title.y = element_blank())+
  facet_grid(~class )

