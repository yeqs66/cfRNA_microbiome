library(WGCNA)
library(tidyverse)
library(readxl)
library(DESeq2)
library(sjmisc)
library(sva)

args <- commandArgs(trailingOnly = TRUE)
tomtype <- as.character(args[1])
networktype <- as.character(args[2])
cor_method <- as.character(args[3])

tomtype=gsub("_"," ",tomtype)
networktype<-gsub("_"," ",networktype)

print(tomtype)
print(networktype)

#load("__PATH_001__")
#load("__PATH_002__")
#load("__PATH_003__")
#load("__PATH_004__")
load("__PATH_005__")


coldata_wgcna<-coldata_wgcna_ko
coldata_wgcna<-as.data.frame(lapply(coldata_wgcna, as.numeric))


sft = as.numeric(4)
cor <- WGCNA::cor

net = blockwiseModules(vst_cfRNA, power = sft,
                       TOMType = tomtype, minModuleSize = 30,
                       maxBlockSize = 13000,
                       networkType = networktype,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste(networktype,tomtype,sep="-"),
                       verbose = 3)
cor<-stats::cor

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;

tomtype=gsub(" ","_",tomtype)
networktype<-gsub(" ","_",networktype)
#聚类看模块间区分
for (i in 1:length(net$blockGenes)) {
pdf(file=paste("Gene_dendrogram_and_module_colors-",tomtype,"-",networktype,"-",cor_method,"-",i,".pdf",sep=""))

plotDendroAndColors(net$dendrograms[[i]], moduleColors[net$blockGenes[[i]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = paste("Gene dendrogram and module colors",i,sep="-"))

dev.off()
}

##绘制模块之间相关性图
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
pdf(file=paste("modules_cor-",tomtype,"-",networktype,"-",cor_method,".pdf",sep=""),height=40,width=25)
if (ncol(MEs_col)>3) {
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
} else {
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = F,
                      xLabelsAngle = 90)
}
dev.off()

##模块和性状的关系
nGenes = ncol(vst_cfRNA)
nSamples = nrow(vst_cfRNA)
# 用颜色标签重新计算MEs
MEs0 = moduleEigengenes(vst_cfRNA, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, coldata_wgcna, use = "p", method=cor_method)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 查看每个模块的颜色及包含的基因数目
table(moduleColors)

#绘制基因模块和性状的相关性热图
pdf(file=paste("ori-Module-trait_associations-",tomtype,"-",networktype,"-",cor_method,".pdf",sep=""),width = 100, height=7.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(20, 8.5, 3, 3))
# 通过热图显示相关性
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(coldata_wgcna), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
# 想改变热图配色可以更改 colors = greenWhiteRed(50)
dev.off()

##挑选显著性状并画图
tmp<-as.data.frame(moduleTraitCor) %>% mutate_all(~tidyr::replace_na(.,0))
wgcna_sign<-data.frame(module=rownames(tmp))
for (col in names(tmp)[1:ncol(tmp)]) {
        if (any(tmp[[col]] > 0.2 | tmp[[col]] < -0.2)) {
        wgcna_sign[[col]] <- tmp[[col]]
        }
}

#coldata_wgcna_select<-coldata_wgcna[,which(colnames(coldata_wgcna) %in% colnames(wgcna_sign))]
if (ncol(wgcna_sign)>1) {
coldata_wgcna_select<-as.data.frame(coldata_wgcna[,which(colnames(coldata_wgcna) %in% colnames(wgcna_sign))])
colnames(coldata_wgcna_select)<-colnames(wgcna_sign)[-1]
rownames(coldata_wgcna_select)<-rownames(coldata_wgcna)

moduleTraitCor_flt = cor(MEs, coldata_wgcna_select, use = "p", method=cor_method)
moduleTraitPvalue_flt = corPvalueStudent(moduleTraitCor_flt, nSamples)

pdf(file=paste("flt-Module-trait_associations-",tomtype,"-",networktype,"-",cor_method,".pdf",sep=""),width = 8, height=7.5)
textMatrix_flt = paste(signif(moduleTraitCor_flt, 2), "\n(", signif(moduleTraitPvalue_flt, 1), ")", sep = "")
dim(textMatrix_flt) = dim(moduleTraitCor_flt)
par(mar = c(20, 8.5, 3, 3))
# 通过热图显示相关性
labeledHeatmap(Matrix = moduleTraitCor_flt, xLabels = names(coldata_wgcna_select), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix_flt, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
# 想改变热图配色可以更改 colors = greenWhiteRed(50)
dev.off()

write.csv(t(wgcna_sign),file = paste("wgcna_sign_ori-",tomtype,"-",networktype,"-",cor_method,".csv",sep=""))
}
#write.table(colnames(wgcna_sign)[-1],file = "select_col.txt")
write.csv(t(moduleTraitCor),file=paste("moduleTraitCor-",tomtype,"-",networktype,"-",cor_method,".csv",sep=""))
write.csv(t(moduleTraitPvalue),file=paste("moduleTraitCorPvalue-",tomtype,"-",networktype,"-",cor_method,".csv",sep=""))

##输出module内gene
module_gene <- as.data.frame(moduleLabels)
module_gene <- cbind(module_gene,moduleColors)
write.csv(module_gene,file = paste("module_gene-",tomtype,"-",networktype,"-",cor_method,".csv",sep=""))

###输出到cytoscape
TOM=TOMsimilarityFromExpr(vst_cfRNA,power=sft,TOMType = tomtype,networkType=networktype)
for (i in unique(moduleColors)){
module=i
probes = colnames(vst_cfRNA)
inModule = (moduleColors==module)
modProbes=probes[inModule]
head(modProbes)
modTOM = TOM[inModule,inModule]
dimnames(modTOM)=list(modProbes,modProbes)
### 这里只选了top100的基因
nTop=100
IMConn = softConnectivity(vst_cfRNA[,modProbes])
top=(rank(-IMConn)<=nTop)
filterTOM=modTOM[top,top]
# for visANT
vis = exportNetworkToVisANT(filterTOM,file = paste("visANTinput-",module,".txt",sep = ""),
weighted = T,threshold = 0)
# for cytoscape
cyt = exportNetworkToCytoscape(filterTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes[top],
nodeAttr = moduleColors[inModule][top])
}

## 构建完整的保存路径
current_path <- getwd()
save_file_name <- paste("Result-",tomtype,"-",networktype,"-",cor_method,".RData",sep="")
save_path <- file.path(current_path, save_file_name)
save.image(file = save_path)

