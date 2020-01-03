#############################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
#############################
BAT_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_M_Y/")
dense.size <- object.size(x = as.matrix(x = BAT_M_Y.data))
sparse.size <- object.size(x = BAT_M_Y.data)
BAT_M_Y <- CreateSeuratObject(raw.data = BAT_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_BAT_M_Y")
BAT_M_Y@meta.data$gender<- "M"
BAT_M_Y@meta.data$age<- "Y"
BAT_M_Y@meta.data$sample<- "M.Y"

mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_M_Y@raw.data[mito.gens, ])/Matrix::colSums(BAT_M_Y@raw.data)
BAT_M_Y <- AddMetaData(object = BAT_M_Y, metadata = percent.mito, col.name = "percent.mito")
BAT_M_Y <- FilterCells(object = BAT_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_M_Y <- NormalizeData(object = BAT_M_Y, normalization.method = "LogNormalize")
BAT_M_Y <- ScaleData(object = BAT_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_Y <- FindVariableGenes(object = BAT_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

BAT_M_Y <- RunPCA(BAT_M_Y, verbose = FALSE)
BAT_M_Y <- FindClusters(BAT_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_M_Y <- RunTSNE(BAT_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BAT_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_M_Y <- doubletFinder(BAT_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_M_Y <- doubletFinder(BAT_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)

BAT_M_Y@meta.data$Doublet <- BAT_M_Y@meta.data$DF.classifications_0.25_0.04_95

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_M_Y.pdf',width=6,height=6)
DimPlot(BAT_M_Y,reduction.use='tsne',do.label=T)
DimPlot(BAT_M_Y,reduction.use='tsne',group.by='Doublet')
dev.off()

########## Male,O-AL

BAT_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_M_O/")
dense.size <- object.size(x = as.matrix(x = BAT_M_O.data))
sparse.size <- object.size(x = BAT_M_O.data)
BAT_M_O <- CreateSeuratObject(raw.data = BAT_M_O.data, min.cells = 5, min.genes = 200, project = "10X_BAT_M_O")
BAT_M_O@meta.data$gender<- "M"
BAT_M_O@meta.data$age<- "O"
BAT_M_O@meta.data$sample<- "M.O"

mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_M_O@raw.data[mito.gens, ])/Matrix::colSums(BAT_M_O@raw.data)
BAT_M_O <- AddMetaData(object = BAT_M_O, metadata = percent.mito, col.name = "percent.mito")
BAT_M_O <- FilterCells(object = BAT_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_M_O <- NormalizeData(object = BAT_M_O, normalization.method = "LogNormalize")
BAT_M_O <- ScaleData(object = BAT_M_O, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_O <- FindVariableGenes(object = BAT_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
length(x = BAT_M_O@var.genes)

BAT_M_O <- RunPCA(BAT_M_O, verbose = FALSE)
BAT_M_O <- FindClusters(BAT_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_M_O <- RunTSNE(BAT_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BAT_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_M_O <- doubletFinder(BAT_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_M_O <- doubletFinder(BAT_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)

table(BAT_M_O@meta.data$DF.classifications_0.25_0.02_115)

BAT_M_O@meta.data$Doublet <- BAT_M_O@meta.data$DF.classifications_0.25_0.02_115

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_M_O.pdf',width=6,height=6)
DimPlot(BAT_M_O,reduction.use='tsne',do.label=T)
DimPlot(BAT_M_O,reduction.use='tsne',group.by='Doublet')
dev.off()

########## Male,O-CR

BAT_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_M_CR/")
dense.size <- object.size(x = as.matrix(x = BAT_M_CR.data))
sparse.size <- object.size(x = BAT_M_CR.data)
BAT_M_CR <- CreateSeuratObject(raw.data = BAT_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_BAT_M_CR")
BAT_M_CR@meta.data$gender<- "M"
BAT_M_CR@meta.data$age<- "CR"
BAT_M_CR@meta.data$sample<- "M.CR"

mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_M_CR@raw.data[mito.gens, ])/Matrix::colSums(BAT_M_CR@raw.data)
BAT_M_CR <- AddMetaData(object = BAT_M_CR, metadata = percent.mito, col.name = "percent.mito")
BAT_M_CR <- FilterCells(object = BAT_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_M_CR <- NormalizeData(object = BAT_M_CR, normalization.method = "LogNormalize")
BAT_M_CR <- ScaleData(object = BAT_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_CR <- FindVariableGenes(object = BAT_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
length(x = BAT_M_CR@var.genes)

BAT_M_CR <- RunPCA(BAT_M_CR, verbose = FALSE)
BAT_M_CR <- FindClusters(BAT_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_M_CR <- RunTSNE(BAT_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(BAT_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_M_CR <- doubletFinder(BAT_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_M_CR <- doubletFinder(BAT_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)


BAT_M_CR@meta.data$Doublet <- BAT_M_CR@meta.data$DF.classifications_0.25_0.03_160

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_M_CR.pdf',width=6,height=6)
DimPlot(BAT_M_CR,reduction.use='tsne',do.label=T)
DimPlot(BAT_M_CR,reduction.use='tsne',group.by='Doublet')
dev.off()

########## Female,Y-AL

BAT_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_F_Y/")
dense.size <- object.size(x = as.matrix(x = BAT_F_Y.data))
sparse.size <- object.size(x = BAT_F_Y.data)
BAT_F_Y <- CreateSeuratObject(raw.data = BAT_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_BAT_F_Y")
BAT_F_Y@meta.data$gender<- "F"
BAT_F_Y@meta.data$age<- "Y"
BAT_F_Y@meta.data$sample<- "F.Y"

mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_F_Y@raw.data[mito.gens, ])/Matrix::colSums(BAT_F_Y@raw.data)
BAT_F_Y <- AddMetaData(object = BAT_F_Y, metadata = percent.mito, col.name = "percent.mito")
BAT_F_Y <- FilterCells(object = BAT_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_F_Y <- NormalizeData(object = BAT_F_Y, normalization.method = "LogNormalize")
BAT_F_Y <- ScaleData(object = BAT_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_Y <- FindVariableGenes(object = BAT_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
length(x = BAT_F_Y@var.genes)

BAT_F_Y <- RunPCA(BAT_F_Y, verbose = FALSE)
BAT_F_Y <- FindClusters(BAT_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_F_Y <- RunTSNE(BAT_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BAT_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_F_Y <- doubletFinder(BAT_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_F_Y <- doubletFinder(BAT_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
	 
BAT_F_Y@meta.data$Doublet <- BAT_F_Y@meta.data$DF.classifications_0.25_0.005_90

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_F_Y.pdf',width=6,height=6)
DimPlot(BAT_F_Y,reduction.use='tsne',do.label=T)
DimPlot(BAT_F_Y,reduction.use='tsne',group.by='Doublet')
dev.off()


########## Female,O-AL

BAT_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_F_O/")
dense.size <- object.size(x = as.matrix(x = BAT_F_O.data))
sparse.size <- object.size(x = BAT_F_O.data)
BAT_F_O <- CreateSeuratObject(raw.data = BAT_F_O.data, min.cells = 5, min.genes = 200, project = "10X_BAT_F_O")
BAT_F_O@meta.data$gender<- "F"
BAT_F_O@meta.data$age<- "O"
BAT_F_O@meta.data$sample<- "F.O"

mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_F_O@raw.data[mito.gens, ])/Matrix::colSums(BAT_F_O@raw.data)
BAT_F_O <- AddMetaData(object = BAT_F_O, metadata = percent.mito, col.name = "percent.mito")
BAT_F_O <- FilterCells(object = BAT_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_F_O <- NormalizeData(object = BAT_F_O, normalization.method = "LogNormalize")
BAT_F_O <- ScaleData(object = BAT_F_O, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_O <- FindVariableGenes(object = BAT_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
length(x = BAT_F_O@var.genes)

BAT_F_O <- RunPCA(BAT_F_O, verbose = FALSE)
BAT_F_O <- FindClusters(BAT_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_F_O <- RunTSNE(BAT_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.061*length(BAT_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_F_O <- doubletFinder(BAT_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_F_O <- doubletFinder(BAT_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)

table(BAT_F_O@meta.data$DF.classifications_0.25_0.005_400)
	 
BAT_F_O@meta.data$Doublet <- BAT_F_O@meta.data$DF.classifications_0.25_0.005_400

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_F_O.pdf',width=6,height=6)
DimPlot(BAT_F_O,reduction.use='tsne',do.label=T)
DimPlot(BAT_F_O,reduction.use='tsne',group.by='Doublet')
dev.off()

########## Female,O-CR

BAT_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/zongzhifang-7/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = BAT_F_CR.data))
sparse.size <- object.size(x = BAT_F_CR.data)
BAT_F_CR <- CreateSeuratObject(raw.data = BAT_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_BAT_F_CR")
BAT_F_CR@meta.data$gender<- "F"
BAT_F_CR@meta.data$age<- "CR"
BAT_F_CR@meta.data$sample<- "F.CR"


mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BAT_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(BAT_F_CR@raw.data[mito.gens, ])/Matrix::colSums(BAT_F_CR@raw.data)
BAT_F_CR <- AddMetaData(object = BAT_F_CR, metadata = percent.mito, col.name = "percent.mito")
BAT_F_CR <- FilterCells(object = BAT_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))

BAT_F_CR <- NormalizeData(object = BAT_F_CR, normalization.method = "LogNormalize")
BAT_F_CR <- ScaleData(object = BAT_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_CR <- FindVariableGenes(object = BAT_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
length(x = BAT_F_CR@var.genes)


BAT_F_CR <- RunPCA(BAT_F_CR, verbose = FALSE)
BAT_F_CR <- FindClusters(BAT_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BAT_F_CR <- RunTSNE(BAT_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 

sweep.res.list_SM <- paramSweep(BAT_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BAT_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.046*length(BAT_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BAT_F_CR <- doubletFinder(BAT_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BAT_F_CR <- doubletFinder(BAT_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)

table(BAT_F_CR@meta.data$DF.classifications_0.25_0.24_224)
	 
BAT_F_CR@meta.data$Doublet <- BAT_F_CR@meta.data$DF.classifications_0.25_0.24_224

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_F_CR.pdf',width=6,height=6)
DimPlot(BAT_F_CR,reduction.use='tsne',do.label=T)
DimPlot(BAT_F_CR,reduction.use='tsne',group.by='Doublet')
dev.off()

#####################################
BAT_M_Y <- SubsetData(BAT_M_Y,subset.name='Doublet',accept.value='Singlet')
BAT_M_Y <- NormalizeData(object = BAT_M_Y, normalization.method = "LogNormalize")
BAT_M_Y <- ScaleData(object = BAT_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_Y <- FindVariableGenes(object = BAT_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BAT_M_O <- SubsetData(BAT_M_O,subset.name='Doublet',accept.value='Singlet')
BAT_M_O <- NormalizeData(object = BAT_M_O, normalization.method = "LogNormalize")
BAT_M_O <- ScaleData(object = BAT_M_O, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_O <- FindVariableGenes(object = BAT_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BAT_M_CR <- SubsetData(BAT_M_CR,subset.name='Doublet',accept.value='Singlet')
BAT_M_CR <- NormalizeData(object = BAT_M_CR, normalization.method = "LogNormalize")
BAT_M_CR <- ScaleData(object = BAT_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
BAT_M_CR <- FindVariableGenes(object = BAT_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BAT_F_Y <- SubsetData(BAT_F_Y,subset.name='Doublet',accept.value='Singlet')
BAT_F_Y <- NormalizeData(object = BAT_F_Y, normalization.method = "LogNormalize")
BAT_F_Y <- ScaleData(object = BAT_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_Y <- FindVariableGenes(object = BAT_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BAT_F_O <- SubsetData(BAT_F_O,subset.name='Doublet',accept.value='Singlet')
BAT_F_O <- NormalizeData(object = BAT_F_O, normalization.method = "LogNormalize")
BAT_F_O <- ScaleData(object = BAT_F_O, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_O <- FindVariableGenes(object = BAT_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BAT_F_CR <- SubsetData(BAT_F_CR,subset.name='Doublet',accept.value='Singlet')
BAT_F_CR <- NormalizeData(object = BAT_F_CR, normalization.method = "LogNormalize")
BAT_F_CR <- ScaleData(object = BAT_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
BAT_F_CR <- FindVariableGenes(object = BAT_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)


g.1 <- head(rownames(BAT_M_O@hvg.info), 1200) 
g.2 <- head(rownames(BAT_M_Y@hvg.info), 1200) 
g.3 <- head(rownames(BAT_M_CR@hvg.info), 1200) 
g.4 <- head(rownames(BAT_F_O@hvg.info), 1200) 
g.5 <- head(rownames(BAT_F_Y@hvg.info), 1200) 
g.6 <- head(rownames(BAT_F_CR@hvg.info), 1200) 

genes.use <- unique(c(g.1, g.2,g.3,g.4,g.5,g.6)) 
genes.use <- intersect(genes.use, rownames(BAT_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(BAT_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(BAT_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(BAT_F_CR@scale.data))

###########CCA

BAT.list <- list(BAT_M_Y, BAT_M_O, BAT_M_CR, BAT_F_Y ,BAT_F_O, BAT_F_CR)
BAT_cca <- RunMultiCCA(object.list = BAT.list, genes.use = genes.use, num.ccs = 30)
BAT_cca <- AlignSubspace(BAT_cca, reduction.type = "cca", grouping.var = "sample", dims.align = 1:20)

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_DimHeatmap.pdf',width=6,height=6)
DimHeatmap(object = BAT_cca, reduction.type = "cca", cells.use = 500, dim.use = 1:12, do.balanced = TRUE) 
DimHeatmap(object = BAT_cca, reduction.type = "cca", cells.use = 500, dim.use = 13:24, do.balanced = TRUE) 
dev.off()

BAT_cca <- FindClusters(BAT_cca, reduction.type =  "cca.aligned", resolution = 1.2, dims.use = 1:7)
BAT_cca <- RunTSNE(BAT_cca, reduction.use =  "cca.aligned", dims.use = 1:7, do.fast = T) 

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_tSNE_sample.pdf',width=7,height=6)
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("seagreen1","orange1", "dodgerblue1", "seagreen4", "orange4", "dodgerblue4")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("seagreen1","NA", "NA", "NA", "NA", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","orange1", "NA", "NA", "NA", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "dodgerblue1", "NA", "NA", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "seagreen4", "NA", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "NA", "orange4", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "NA", "NA", "dodgerblue4")) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_tSNE_age.pdf',width=7,height=6)
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "age" ,colors.use= c("seagreen2","orange2", "dodgerblue2")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("seagreen2","NA", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("NA","orange2", "NA")) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("NA","NA", "dodgerblue2")) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_tSNE_gender.pdf',width=7,height=6)
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('plum2','steelblue1')) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('plum2','NA')) 
TSNEPlot(BAT_cca, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('NA','steelblue1')) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_tSNE_cluster.pdf',width=7,height=6)
TSNEPlot(BAT_cca, do.label = T, do.return = T, pt.size = 0.4,label.size=5) 
dev.off()

png('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_FeaturePlot.png',width=1000,height=1000)
FeaturePlot(object = BAT_cca, features.plot = c("Lum","Dcn","Pdgfra","Anpep","Pecam1","Kdr", "Acta2","Ptprc"), cols.use = c("grey95","red"), reduction.use = "tsne", pt.size=0.5,nCol=3, no.axes =T)
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_VlnPlot.pdf',width=24,height=20)
VlnPlot(object = BAT_cca, features.plot = c("Lum","Dcn","Pdgfra","Anpep","Pecam1","Kdr", "Acta2","Ptprc"), point.size.use = 0,nCol=1)
dev.off()

markers <- FindAllMarkers(object = BAT_cca, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_markers.csv')	
saveRDS(BAT_cca, file = "/data1/mashuai/data/rat_10X/RDS/New/BAT_Y_O_CR.rds")

BAT_cca_CT <- BAT_cca

new.ident <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24")

new.ident <- c("IC","IC","IC","IC","IC","IC","ASC","IC","Fib","IC","ASC","Fib","Fib","EC","IC","IC","IC","EC","Fib","ASC","IC","IC","EC","IC","IC")

for (i in 0:24) {
    BAT_cca_CT <- RenameIdent(object = BAT_cca_CT, old.ident.name = i, new.ident.name = new.ident[i + 1])
}


pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_tSNE_celltype.pdf',width=7,height=6)
TSNEPlot(BAT_cca_CT, do.label = T, do.return = T, pt.size = 0.4,label.size=5,colors.use=c("#f46d43","#fee08b","#addd8e","#3288bd"))
TSNEPlot(BAT_cca_CT, do.label = F, do.return = T, pt.size = 0.4,colors.use=c("#f46d43","#fee08b","#addd8e","#3288bd"))  
dev.off()

########### DEGs

BAT_cca_new <- BAT_cca_CT
BAT_cca_new@meta.data$celltype.age <- paste0(BAT_cca_new@ident, "_", BAT_cca_new@meta.data$age)
BAT_cca_new <- StashIdent(BAT_cca_new, save.name = "celltype")
BAT_cca_new <- SetAllIdent(BAT_cca_new, id = "celltype.age")
cell.type <- BAT_cca_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]

identity <- paste(type,"_O",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident){
OvsY <- FindMarkers(BAT_cca_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_cca_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}


markers <- FindAllMarkers(object = BAT_cca_CT, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_markers_CT.csv')

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_CellTree.pdf") 
BuildClusterTree(object= BAT_cca_CT)
dev.off()

BAT_cca_CT_new <- BAT_cca_CT
BAT_cca_CT_new@meta.data$celltype.age <- paste0(BAT_cca_CT_new@ident, "_", BAT_cca_CT_new@meta.data$age)
cell_age <- table(BAT_cca_CT_new@meta.data$celltype.age)
write.csv(cell_age,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_CT_num.csv')

############ Immune cell

BAT_IC <- SubsetData(object= BAT_cca_CT, ident.use ="IC" )
BAT_IC <- RunPCA(object = BAT_IC, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
BAT_IC <- JackStraw(object = BAT_IC, num.replicate = 1000)

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_IC_PCHeatmap_JackStrawPlot.pdf',width=6,height=12)
PCHeatmap(object = BAT_IC, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
JackStrawPlot(object = BAT_IC, PCs = 1:20)
dev.off()

BAT_IC <- FindClusters(object = BAT_IC, reduction.type = "pca", dims.use = c(1,2,5:11),     resolution = 0.8, save.SNN = TRUE)
BAT_IC <- RunTSNE(BAT_IC, reduction.use =  "pca", dims.use = c(1,2,5:11), do.fast = T) 

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_tSNE_sample.pdf',width=7,height=6)
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("seagreen1","orange1", "dodgerblue1", "seagreen4", "orange4", "dodgerblue4")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("seagreen1","NA", "NA", "NA", "NA", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","orange1", "NA", "NA", "NA", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "dodgerblue1", "NA", "NA", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "seagreen4", "NA", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "NA", "orange4", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "sample", colors.use= c("NA","NA", "NA", "NA", "NA", "dodgerblue4")) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_tSNE_age.pdf',width=7,height=6)
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "age" ,colors.use= c("seagreen2","orange2", "dodgerblue2")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("seagreen2","NA", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("NA","orange2", "NA")) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "age", colors.use= c("NA","NA", "dodgerblue2")) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_tSNE_gender.pdf',width=7,height=6)
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('plum2','steelblue1')) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('plum2','NA')) 
TSNEPlot(BAT_IC, do.return = T, pt.size = 0.4, group.by = "gender", colors.use= c('NA','steelblue1')) 
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_tSNE_cluster.pdf',width=7,height=6)
TSNEPlot(BAT_IC, do.label = T, do.return = T, pt.size = 0.4,label.size=5) 
dev.off()

png('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_FeaturePlot.png',width=1500,height=1500)
FeaturePlot(object = BAT_IC, features.plot = c("Ptprc","Cd3e","Cd4","Cd8a","Cd19","Cd79b","Klrb1a","Cd69","Cd68","Cd14","Cd163","S100a8","S100a9","Cd74","Irf5","Irf8","Mzb1","Igh-6","Nkg7","Mrc1","Fcgr3a"), cols.use = c("grey95","red"), reduction.use = "tsne", pt.size=0.5,nCol=5, no.axes =T)
dev.off()

pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_VlnPlot.pdf',width=21,height=52.5)
VlnPlot(object = BAT_IC, features.plot = c("Ptprc","Cd3e","Cd4","Cd8a","Cd19","Cd79b","Klrb1a","Cd69","Cd68","Cd14","Cd163","S100a8","S100a9","Cd74","Irf5","Irf8","Mzb1","Igh-6","Nkg7","Mrc1","Fcgr3a"), point.size.use = 0,nCol=1)
dev.off()

saveRDS(BAT_IC, file = "/data1/mashuai/data/rat_10X/RDS/New/BAT_IC_Y_O_CR.rds")

BAT_IC_CT <- BAT_IC

new.ident <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")

new.ident <- c("T2","Neu","NK","M1","Neu","DC1","M1","M1","M2","Neu","Neu","M2","M2","T2","B","NKT","DC2","DC1","Neu","NKT","Pla")

for (i in 0:20) {
    BAT_IC_CT <- RenameIdent(object = BAT_IC_CT, old.ident.name = i, new.ident.name = new.ident[i + 1])
}


pdf('/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_tSNE_celltype.pdf',width=7,height=6)
TSNEPlot(BAT_IC_CT, do.label = T, do.return = T, pt.size = 0.4,label.size=5,colors.use=c("#9e0142","#f46d43","#fee08b","#ffffbf","#e6f598","#abdda4","#addd8e","#66c2a5","#6baed6","#5e4fa2","#dadaeb"))
TSNEPlot(BAT_IC_CT, do.label = F, do.return = T, pt.size = 0.4,colors.use=c("#9e0142","#f46d43","#fee08b","#ffffbf","#e6f598","#abdda4","#addd8e","#66c2a5","#6baed6","#5e4fa2","#dadaeb"))  
dev.off()

########### IC DEGs

BAT_IC_new <- BAT_IC_CT
BAT_IC_new@meta.data$celltype.age <- paste0(BAT_IC_new@ident, "_", BAT_IC_new@meta.data$age)
BAT_IC_new <- StashIdent(BAT_IC_new, save.name = "celltype")
BAT_IC_new <- SetAllIdent(BAT_IC_new, id = "celltype.age")
cell.type <- BAT_IC_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident)
{
if ((length(which(BAT_IC_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_IC_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_IC_new@ident %in% identity_CR)) >= 3 ))
{
OvsY <- FindMarkers(BAT_IC_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_IC_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}
}

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_IC_CellTree.pdf") 
BuildClusterTree(object= BAT_IC_CT)
dev.off()

BAT_IC_CT_new <- BAT_IC_CT
BAT_IC_CT_new@meta.data$celltype.age <- paste0(BAT_IC_CT_new@ident, "_", BAT_IC_CT_new@meta.data$age)
cell_age <- table(BAT_IC_CT_new@meta.data$celltype.age)
write.csv(cell_age,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_IC_CT_num.csv')

markers <- FindAllMarkers(object = BAT_IC_CT, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/BAT_IC_markers_CT.csv')
