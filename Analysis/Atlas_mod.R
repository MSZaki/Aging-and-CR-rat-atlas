##########################################################

############################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)

###################################################  Tissue: BAT
############################BAT_M_Y
BAT_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/BAT/BAT_M_Y/")
dense.size <- object.size(x = as.matrix(x = BAT_M_Y.data))
sparse.size <- object.size(x = BAT_M_Y.data)
BAT_M_Y <- CreateSeuratObject(raw.data = BAT_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_BAT_M_Y")
BAT_M_Y@meta.data$gender<- "M"
BAT_M_Y@meta.data$age<- "Y"
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

############################BAT_M_O
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
BAT_M_O@meta.data$Doublet <- BAT_M_O@meta.data$DF.classifications_0.25_0.02_115

############################BAT_M_CR
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

############################BAT_F_Y
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

############################BAT_F_O
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
BAT_F_O@meta.data$Doublet <- BAT_F_O@meta.data$DF.classifications_0.25_0.005_400

############################BAT_F_CR
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
BAT_F_CR@meta.data$Doublet <- BAT_F_CR@meta.data$DF.classifications_0.25_0.24_224

###################################################  Tissue: WAT
############################WAT_M_Y
WAT_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/WAT/WAT_M_Y/")
dense.size <- object.size(x = as.matrix(x = WAT_M_Y.data))
sparse.size <- object.size(x = WAT_M_Y.data)
WAT_M_Y <- CreateSeuratObject(raw.data = WAT_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_WAT_M_Y")
WAT_M_Y@meta.data$gender<- "M"
WAT_M_Y@meta.data$age<- "Y"
WAT_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_M_Y@raw.data[mito.gens, ])/Matrix::colSums(WAT_M_Y@raw.data)
WAT_M_Y <- AddMetaData(object = WAT_M_Y, metadata = percent.mito, col.name = "percent.mito")
WAT_M_Y <- FilterCells(object = WAT_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_M_Y <- NormalizeData(object = WAT_M_Y, normalization.method = "LogNormalize")
WAT_M_Y <- ScaleData(object = WAT_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_Y <- FindVariableGenes(object = WAT_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_M_Y <- RunPCA(WAT_M_Y, verbose = FALSE)
WAT_M_Y <- FindClusters(WAT_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_M_Y <- RunTSNE(WAT_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.046*length(WAT_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_M_Y <- doubletFinder(WAT_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_M_Y <- doubletFinder(WAT_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
WAT_M_Y@meta.data$Doublet <- WAT_M_Y@meta.data$DF.classifications_0.25_0.27_240

############################WAT_M_O
WAT_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/WAT/WAT_M_O/")
dense.size <- object.size(x = as.matrix(x = WAT_M_O.data))
sparse.size <- object.size(x = WAT_M_O.data)
WAT_M_O <- CreateSeuratObject(raw.data = WAT_M_O.data, min.cells = 5, min.genes = 200, project = "10X_WAT_M_O")
WAT_M_O@meta.data$gender<- "M"
WAT_M_O@meta.data$age<- "O"
WAT_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_M_O@raw.data[mito.gens, ])/Matrix::colSums(WAT_M_O@raw.data)
WAT_M_O <- AddMetaData(object = WAT_M_O, metadata = percent.mito, col.name = "percent.mito")
WAT_M_O <- FilterCells(object = WAT_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_M_O <- NormalizeData(object = WAT_M_O, normalization.method = "LogNormalize")
WAT_M_O <- ScaleData(object = WAT_M_O, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_O <- FindVariableGenes(object = WAT_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_M_O <- RunPCA(WAT_M_O, verbose = FALSE)
WAT_M_O <- FindClusters(WAT_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_M_O <- RunTSNE(WAT_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(WAT_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_M_O <- doubletFinder(WAT_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_M_O <- doubletFinder(WAT_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
WAT_M_O@meta.data$Doublet <- WAT_M_O@meta.data$DF.classifications_0.25_0.005_157

############################WAT_M_CR
WAT_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/WAT/WAT_M_CR/")
dense.size <- object.size(x = as.matrix(x = WAT_M_CR.data))
sparse.size <- object.size(x = WAT_M_CR.data)
WAT_M_CR <- CreateSeuratObject(raw.data = WAT_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_WAT_M_CR")
WAT_M_CR@meta.data$gender<- "M"
WAT_M_CR@meta.data$age<- "CR"
WAT_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_M_CR@raw.data[mito.gens, ])/Matrix::colSums(WAT_M_CR@raw.data)
WAT_M_CR <- AddMetaData(object = WAT_M_CR, metadata = percent.mito, col.name = "percent.mito")
WAT_M_CR <- FilterCells(object = WAT_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_M_CR <- NormalizeData(object = WAT_M_CR, normalization.method = "LogNormalize")
WAT_M_CR <- ScaleData(object = WAT_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_CR <- FindVariableGenes(object = WAT_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_M_CR <- RunPCA(WAT_M_CR, verbose = FALSE)
WAT_M_CR <- FindClusters(WAT_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_M_CR <- RunTSNE(WAT_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(WAT_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_M_CR <- doubletFinder(WAT_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_M_CR <- doubletFinder(WAT_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
WAT_M_CR@meta.data$Doublet <- WAT_M_CR@meta.data$DF.classifications_0.25_0.07_113

############################WAT_F_Y
WAT_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/WAT/WAT_F_Y/")
dense.size <- object.size(x = as.matrix(x = WAT_F_Y.data))
sparse.size <- object.size(x = WAT_F_Y.data)
WAT_F_Y <- CreateSeuratObject(raw.data = WAT_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_WAT_F_Y")
WAT_F_Y@meta.data$gender<- "F"
WAT_F_Y@meta.data$age<- "Y"
WAT_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_F_Y@raw.data[mito.gens, ])/Matrix::colSums(WAT_F_Y@raw.data)
WAT_F_Y <- AddMetaData(object = WAT_F_Y, metadata = percent.mito, col.name = "percent.mito")
WAT_F_Y <- FilterCells(object = WAT_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_F_Y <- NormalizeData(object = WAT_F_Y, normalization.method = "LogNormalize")
WAT_F_Y <- ScaleData(object = WAT_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_Y <- FindVariableGenes(object = WAT_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_Y <- RunPCA(WAT_F_Y, verbose = FALSE)
WAT_F_Y <- FindClusters(WAT_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_F_Y <- RunTSNE(WAT_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.046*length(WAT_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_F_Y <- doubletFinder(WAT_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_F_Y <- doubletFinder(WAT_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
WAT_F_Y@meta.data$Doublet <- WAT_F_Y@meta.data$DF.classifications_0.25_0.24_241

############################WAT_F_O
WAT_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/baizhifang-6/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = WAT_F_O.data))
sparse.size <- object.size(x = WAT_F_O.data)
WAT_F_O <- CreateSeuratObject(raw.data = WAT_F_O.data, min.cells = 5, min.genes = 200, project = "10X_WAT_F_O")
WAT_F_O@meta.data$gender<- "F"
WAT_F_O@meta.data$age<- "O"
WAT_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_F_O@raw.data[mito.gens, ])/Matrix::colSums(WAT_F_O@raw.data)
WAT_F_O <- AddMetaData(object = WAT_F_O, metadata = percent.mito, col.name = "percent.mito")
WAT_F_O <- FilterCells(object = WAT_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_F_O <- NormalizeData(object = WAT_F_O, normalization.method = "LogNormalize")
WAT_F_O <- ScaleData(object = WAT_F_O, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_O <- FindVariableGenes(object = WAT_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_O <- RunPCA(WAT_F_O, verbose = FALSE)
WAT_F_O <- FindClusters(WAT_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_F_O <- RunTSNE(WAT_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(WAT_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_F_O <- doubletFinder(WAT_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_F_O <- doubletFinder(WAT_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
WAT_F_O@meta.data$Doublet <- WAT_F_O@meta.data$DF.classifications_0.25_0.25_173

############################WAT_F_CR
WAT_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/baizhifang-7/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = WAT_F_CR.data))
sparse.size <- object.size(x = WAT_F_CR.data)
WAT_F_CR <- CreateSeuratObject(raw.data = WAT_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_WAT_F_CR")
WAT_F_CR@meta.data$gender<- "F"
WAT_F_CR@meta.data$age<- "CR"
WAT_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = WAT_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(WAT_F_CR@raw.data[mito.gens, ])/Matrix::colSums(WAT_F_CR@raw.data)
WAT_F_CR <- AddMetaData(object = WAT_F_CR, metadata = percent.mito, col.name = "percent.mito")
WAT_F_CR <- FilterCells(object = WAT_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
WAT_F_CR <- NormalizeData(object = WAT_F_CR, normalization.method = "LogNormalize")
WAT_F_CR <- ScaleData(object = WAT_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_CR <- FindVariableGenes(object = WAT_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_CR <- RunPCA(WAT_F_CR, verbose = FALSE)
WAT_F_CR <- FindClusters(WAT_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
WAT_F_CR <- RunTSNE(WAT_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(WAT_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- WAT_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.046*length(WAT_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
WAT_F_CR <- doubletFinder(WAT_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
WAT_F_CR <- doubletFinder(WAT_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
WAT_F_CR@meta.data$Doublet <- WAT_F_CR@meta.data$DF.classifications_0.25_0.005_222

###################################################  Tissue: Liver
############################Liver_M_Y
Liver_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_M_Y/")
dense.size <- object.size(x = as.matrix(x = Liver_M_Y.data))
sparse.size <- object.size(x = Liver_M_Y.data)
Liver_M_Y <- CreateSeuratObject(raw.data = Liver_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_Liver_M_Y")
Liver_M_Y@meta.data$gender<- "M"
Liver_M_Y@meta.data$age<- "Y"
Liver_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_M_Y@raw.data[mito.gens, ])/Matrix::colSums(Liver_M_Y@raw.data)
Liver_M_Y <- AddMetaData(object = Liver_M_Y, metadata = percent.mito, col.name = "percent.mito")
Liver_M_Y <- FilterCells(object = Liver_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_M_Y <- NormalizeData(object = Liver_M_Y, normalization.method = "LogNormalize")
Liver_M_Y <- ScaleData(object = Liver_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_Y <- FindVariableGenes(object = Liver_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_M_Y <- RunPCA(Liver_M_Y, verbose = FALSE)
Liver_M_Y <- FindClusters(Liver_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_M_Y <- RunTSNE(Liver_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Liver_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_M_Y <- doubletFinder(Liver_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_M_Y <- doubletFinder(Liver_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Liver_M_Y@meta.data$Doublet <- Liver_M_Y@meta.data$DF.classifications_0.25_0.005_107

############################Liver_M_O
Liver_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_M_O/")
dense.size <- object.size(x = as.matrix(x = Liver_M_O.data))
sparse.size <- object.size(x = Liver_M_O.data)
Liver_M_O <- CreateSeuratObject(raw.data = Liver_M_O.data, min.cells = 5, min.genes = 200, project = "10X_Liver_M_O")
Liver_M_O@meta.data$gender<- "M"
Liver_M_O@meta.data$age<- "O"
Liver_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_M_O@raw.data[mito.gens, ])/Matrix::colSums(Liver_M_O@raw.data)
Liver_M_O <- AddMetaData(object = Liver_M_O, metadata = percent.mito, col.name = "percent.mito")
Liver_M_O <- FilterCells(object = Liver_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_M_O <- NormalizeData(object = Liver_M_O, normalization.method = "LogNormalize")
Liver_M_O <- ScaleData(object = Liver_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_O <- FindVariableGenes(object = Liver_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_M_O <- RunPCA(Liver_M_O, verbose = FALSE)
Liver_M_O <- FindClusters(Liver_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_M_O <- RunTSNE(Liver_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.016*length(Liver_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_M_O <- doubletFinder(Liver_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_M_O <- doubletFinder(Liver_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Liver_M_O@meta.data$Doublet <- Liver_M_O@meta.data$DF.classifications_0.25_0.01_34

############################Liver_M_CR
Liver_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_M_CR/")
dense.size <- object.size(x = as.matrix(x = Liver_M_CR.data))
sparse.size <- object.size(x = Liver_M_CR.data)
Liver_M_CR <- CreateSeuratObject(raw.data = Liver_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_Liver_M_CR")
Liver_M_CR@meta.data$gender<- "M"
Liver_M_CR@meta.data$age<- "CR"
Liver_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_M_CR@raw.data[mito.gens, ])/Matrix::colSums(Liver_M_CR@raw.data)
Liver_M_CR <- AddMetaData(object = Liver_M_CR, metadata = percent.mito, col.name = "percent.mito")
Liver_M_CR <- FilterCells(object = Liver_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_M_CR <- NormalizeData(object = Liver_M_CR, normalization.method = "LogNormalize")
Liver_M_CR <- ScaleData(object = Liver_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_CR <- FindVariableGenes(object = Liver_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_M_CR <- RunPCA(Liver_M_CR, verbose = FALSE)
Liver_M_CR <- FindClusters(Liver_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_M_CR <- RunTSNE(Liver_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Liver_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_M_CR <- doubletFinder(Liver_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_M_CR <- doubletFinder(Liver_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Liver_M_CR@meta.data$Doublet <- Liver_M_CR@meta.data$DF.classifications_0.25_0.005_107

############################Liver_F_Y
Liver_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_F_Y/")
dense.size <- object.size(x = as.matrix(x = Liver_F_Y.data))
sparse.size <- object.size(x = Liver_F_Y.data)
Liver_F_Y <- CreateSeuratObject(raw.data = Liver_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_Liver_F_Y")
Liver_F_Y@meta.data$gender<- "F"
Liver_F_Y@meta.data$age<- "Y"
Liver_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_F_Y@raw.data[mito.gens, ])/Matrix::colSums(Liver_F_Y@raw.data)
Liver_F_Y <- AddMetaData(object = Liver_F_Y, metadata = percent.mito, col.name = "percent.mito")
Liver_F_Y <- FilterCells(object = Liver_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_F_Y <- NormalizeData(object = Liver_F_Y, normalization.method = "LogNormalize")
Liver_F_Y <- ScaleData(object = Liver_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_Y <- FindVariableGenes(object = Liver_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_Y <- RunPCA(Liver_F_Y, verbose = FALSE)
Liver_F_Y <- FindClusters(Liver_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_F_Y <- RunTSNE(Liver_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Liver_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_F_Y <- doubletFinder(Liver_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_F_Y <- doubletFinder(Liver_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Liver_F_Y@meta.data$Doublet <- Liver_F_Y@meta.data$DF.classifications_0.25_0.29_51

############################Liver_F_O
Liver_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_F_O/")
dense.size <- object.size(x = as.matrix(x = Liver_F_O.data))
sparse.size <- object.size(x = Liver_F_O.data)
Liver_F_O <- CreateSeuratObject(raw.data = Liver_F_O.data, min.cells = 5, min.genes = 200, project = "10X_Liver_F_O")
Liver_F_O@meta.data$gender<- "F"
Liver_F_O@meta.data$age<- "O"
Liver_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_F_O@raw.data[mito.gens, ])/Matrix::colSums(Liver_F_O@raw.data)
Liver_F_O <- AddMetaData(object = Liver_F_O, metadata = percent.mito, col.name = "percent.mito")
Liver_F_O <- FilterCells(object = Liver_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_F_O <- NormalizeData(object = Liver_F_O, normalization.method = "LogNormalize")
Liver_F_O <- ScaleData(object = Liver_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_O <- FindVariableGenes(object = Liver_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_O <- RunPCA(Liver_F_O, verbose = FALSE)
Liver_F_O <- FindClusters(Liver_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_F_O <- RunTSNE(Liver_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Liver_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_F_O <- doubletFinder(Liver_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_F_O <- doubletFinder(Liver_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
Liver_F_O@meta.data$Doublet <- Liver_F_O@meta.data$DF.classifications_0.25_0.29_66

############################Liver_F_cR
Liver_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Liver/Liver_F_CR/")
dense.size <- object.size(x = as.matrix(x = Liver_F_CR.data))
sparse.size <- object.size(x = Liver_F_CR.data)
Liver_F_CR <- CreateSeuratObject(raw.data = Liver_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_Liver_F_CR")
Liver_F_CR@meta.data$gender<- "F"
Liver_F_CR@meta.data$age<- "CR"
Liver_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Liver_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Liver_F_CR@raw.data[mito.gens, ])/Matrix::colSums(Liver_F_CR@raw.data)
Liver_F_CR <- AddMetaData(object = Liver_F_CR, metadata = percent.mito, col.name = "percent.mito")
Liver_F_CR <- FilterCells(object = Liver_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Liver_F_CR <- NormalizeData(object = Liver_F_CR, normalization.method = "LogNormalize")
Liver_F_CR <- ScaleData(object = Liver_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_CR <- FindVariableGenes(object = Liver_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_CR <- RunPCA(Liver_F_CR, verbose = FALSE)
Liver_F_CR <- FindClusters(Liver_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Liver_F_CR <- RunTSNE(Liver_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Liver_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Liver_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Liver_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Liver_F_CR <- doubletFinder(Liver_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Liver_F_CR <- doubletFinder(Liver_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Liver_F_CR@meta.data$Doublet <- Liver_F_CR@meta.data$DF.classifications_0.25_0.2_95

###################################################  Tissue: Kidney
############################Kidney_M_Y
Kidney_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_M_Y/")
dense.size <- object.size(x = as.matrix(x = Kidney_M_Y.data))
sparse.size <- object.size(x = Kidney_M_Y.data)
Kidney_M_Y <- CreateSeuratObject(raw.data = Kidney_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_M_Y")
Kidney_M_Y@meta.data$gender<- "M"
Kidney_M_Y@meta.data$age<- "Y"
Kidney_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_M_Y@raw.data[mito.gens, ])/Matrix::colSums(Kidney_M_Y@raw.data)
Kidney_M_Y <- AddMetaData(object = Kidney_M_Y, metadata = percent.mito, col.name = "percent.mito")
Kidney_M_Y <- FilterCells(object = Kidney_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_M_Y <- NormalizeData(object = Kidney_M_Y, normalization.method = "LogNormalize")
Kidney_M_Y <- ScaleData(object = Kidney_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_Y <- FindVariableGenes(object = Kidney_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_M_Y <- RunPCA(Kidney_M_Y, verbose = FALSE)
Kidney_M_Y <- FindClusters(Kidney_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_M_Y <- RunTSNE(Kidney_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Kidney_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_M_Y <- doubletFinder(Kidney_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_M_Y <- doubletFinder(Kidney_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Kidney_M_Y@meta.data$Doublet <- Kidney_M_Y@meta.data$DF.classifications_0.25_0.005_77

############################Kidney_M_O
Kidney_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_M_O/")
dense.size <- object.size(x = as.matrix(x = Kidney_M_O.data))
sparse.size <- object.size(x = Kidney_M_O.data)
Kidney_M_O <- CreateSeuratObject(raw.data = Kidney_M_O.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_M_O")
Kidney_M_O@meta.data$gender<- "M"
Kidney_M_O@meta.data$age<- "O"
Kidney_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_M_O@raw.data[mito.gens, ])/Matrix::colSums(Kidney_M_O@raw.data)
Kidney_M_O <- AddMetaData(object = Kidney_M_O, metadata = percent.mito, col.name = "percent.mito")
Kidney_M_O <- FilterCells(object = Kidney_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_M_O <- NormalizeData(object = Kidney_M_O, normalization.method = "LogNormalize")
Kidney_M_O <- ScaleData(object = Kidney_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_O <- FindVariableGenes(object = Kidney_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_M_O <- RunPCA(Kidney_M_O, verbose = FALSE)
Kidney_M_O <- FindClusters(Kidney_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_M_O <- RunTSNE(Kidney_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Kidney_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_M_O <- doubletFinder(Kidney_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_M_O <- doubletFinder(Kidney_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Kidney_M_O@meta.data$Doublet <- Kidney_M_O@meta.data$DF.classifications_0.25_0.01_90

############################Kidney_M_CR
Kidney_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_M_CR/")
dense.size <- object.size(x = as.matrix(x = Kidney_M_CR.data))
sparse.size <- object.size(x = Kidney_M_CR.data)
Kidney_M_CR <- CreateSeuratObject(raw.data = Kidney_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_M_CR")
Kidney_M_CR@meta.data$gender<- "M"
Kidney_M_CR@meta.data$age<- "CR"
Kidney_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_M_CR@raw.data[mito.gens, ])/Matrix::colSums(Kidney_M_CR@raw.data)
Kidney_M_CR <- AddMetaData(object = Kidney_M_CR, metadata = percent.mito, col.name = "percent.mito")
Kidney_M_CR <- FilterCells(object = Kidney_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_M_CR <- NormalizeData(object = Kidney_M_CR, normalization.method = "LogNormalize")
Kidney_M_CR <- ScaleData(object = Kidney_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_CR <- FindVariableGenes(object = Kidney_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_M_CR <- RunPCA(Kidney_M_CR, verbose = FALSE)
Kidney_M_CR <- FindClusters(Kidney_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_M_CR <- RunTSNE(Kidney_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(Kidney_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_M_CR <- doubletFinder(Kidney_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_M_CR <- doubletFinder(Kidney_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Kidney_M_CR@meta.data$Doublet <- Kidney_M_CR@meta.data$DF.classifications_0.25_0.005_140

############################Kidney_F_Y
Kidney_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_F_Y/")
dense.size <- object.size(x = as.matrix(x = Kidney_F_Y.data))
sparse.size <- object.size(x = Kidney_F_Y.data)
Kidney_F_Y <- CreateSeuratObject(raw.data = Kidney_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_F_Y")
Kidney_F_Y@meta.data$gender<- "F"
Kidney_F_Y@meta.data$age<- "Y"
Kidney_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_F_Y@raw.data[mito.gens, ])/Matrix::colSums(Kidney_F_Y@raw.data)
Kidney_F_Y <- AddMetaData(object = Kidney_F_Y, metadata = percent.mito, col.name = "percent.mito")
Kidney_F_Y <- FilterCells(object = Kidney_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_F_Y <- NormalizeData(object = Kidney_F_Y, normalization.method = "LogNormalize")
Kidney_F_Y <- ScaleData(object = Kidney_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_Y <- FindVariableGenes(object = Kidney_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_Y <- RunPCA(Kidney_F_Y, verbose = FALSE)
Kidney_F_Y <- FindClusters(Kidney_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_F_Y <- RunTSNE(Kidney_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Kidney_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_F_Y <- doubletFinder(Kidney_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_F_Y <- doubletFinder(Kidney_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
Kidney_F_Y@meta.data$Doublet <- Kidney_F_Y@meta.data$DF.classifications_0.25_0.29_49

############################Kidney_F_O
Kidney_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_F_O/")
dense.size <- object.size(x = as.matrix(x = Kidney_F_O.data))
sparse.size <- object.size(x = Kidney_F_O.data)
Kidney_F_O <- CreateSeuratObject(raw.data = Kidney_F_O.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_F_O")
Kidney_F_O@meta.data$gender<- "F"
Kidney_F_O@meta.data$age<- "O"
Kidney_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_F_O@raw.data[mito.gens, ])/Matrix::colSums(Kidney_F_O@raw.data)
Kidney_F_O <- AddMetaData(object = Kidney_F_O, metadata = percent.mito, col.name = "percent.mito")
Kidney_F_O <- FilterCells(object = Kidney_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_F_O <- NormalizeData(object = Kidney_F_O, normalization.method = "LogNormalize")
Kidney_F_O <- ScaleData(object = Kidney_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_O <- FindVariableGenes(object = Kidney_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_O <- RunPCA(Kidney_F_O, verbose = FALSE)
Kidney_F_O <- FindClusters(Kidney_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_F_O <- RunTSNE(Kidney_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.016*length(Kidney_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_F_O <- doubletFinder(Kidney_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_F_O <- doubletFinder(Kidney_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Kidney_F_O@meta.data$Doublet <- Kidney_F_O@meta.data$DF.classifications_0.25_0.29_13

############################Kidney_F_CR
Kidney_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Kidney/Kidney_F_CR/")
dense.size <- object.size(x = as.matrix(x = Kidney_F_CR.data))
sparse.size <- object.size(x = Kidney_F_CR.data)
Kidney_F_CR <- CreateSeuratObject(raw.data = Kidney_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_Kidney_F_CR")
Kidney_F_CR@meta.data$gender<- "F"
Kidney_F_CR@meta.data$age<- "CR"
Kidney_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Kidney_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Kidney_F_CR@raw.data[mito.gens, ])/Matrix::colSums(Kidney_F_CR@raw.data)
Kidney_F_CR <- AddMetaData(object = Kidney_F_CR, metadata = percent.mito, col.name = "percent.mito")
Kidney_F_CR <- FilterCells(object = Kidney_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.5))
Kidney_F_CR <- NormalizeData(object = Kidney_F_CR, normalization.method = "LogNormalize")
Kidney_F_CR <- ScaleData(object = Kidney_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_CR <- FindVariableGenes(object = Kidney_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_CR <- RunPCA(Kidney_F_CR, verbose = FALSE)
Kidney_F_CR <- FindClusters(Kidney_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Kidney_F_CR <- RunTSNE(Kidney_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Kidney_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Kidney_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Kidney_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Kidney_F_CR <- doubletFinder(Kidney_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Kidney_F_CR <- doubletFinder(Kidney_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Kidney_F_CR@meta.data$Doublet <- Kidney_F_CR@meta.data$DF.classifications_0.25_0.005_83

###################################################  Tissue: Aorta
############################Aorta_M_Y
Aorta_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_M_Y/")
dense.size <- object.size(x = as.matrix(x = Aorta_M_Y.data))
sparse.size <- object.size(x = Aorta_M_Y.data)
Aorta_M_Y <- CreateSeuratObject(raw.data = Aorta_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_M_Y")
Aorta_M_Y@meta.data$gender<- "M"
Aorta_M_Y@meta.data$age<- "Y"
Aorta_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_M_Y@raw.data[mito.gens, ])/Matrix::colSums(Aorta_M_Y@raw.data)
Aorta_M_Y <- AddMetaData(object = Aorta_M_Y, metadata = percent.mito, col.name = "percent.mito")
Aorta_M_Y <- FilterCells(object = Aorta_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_M_Y <- NormalizeData(object = Aorta_M_Y, normalization.method = "LogNormalize")
Aorta_M_Y <- ScaleData(object = Aorta_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_Y <- FindVariableGenes(object = Aorta_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_M_Y <- RunPCA(Aorta_M_Y, verbose = FALSE)
Aorta_M_Y <- FindClusters(Aorta_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_M_Y <- RunTSNE(Aorta_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Aorta_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_M_Y <- doubletFinder(Aorta_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_M_Y <- doubletFinder(Aorta_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Aorta_M_Y@meta.data$Doublet <- Aorta_M_Y@meta.data$DF.classifications_0.25_0.3_55

############################Aorta_M_O
Aorta_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_M_O/")
dense.size <- object.size(x = as.matrix(x = Aorta_M_O.data))
sparse.size <- object.size(x = Aorta_M_O.data)
Aorta_M_O <- CreateSeuratObject(raw.data = Aorta_M_O.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_M_O")
Aorta_M_O@meta.data$gender<- "M"
Aorta_M_O@meta.data$age<- "O"
Aorta_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_M_O@raw.data[mito.gens, ])/Matrix::colSums(Aorta_M_O@raw.data)
Aorta_M_O <- AddMetaData(object = Aorta_M_O, metadata = percent.mito, col.name = "percent.mito")
Aorta_M_O <- FilterCells(object = Aorta_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_M_O <- NormalizeData(object = Aorta_M_O, normalization.method = "LogNormalize")
Aorta_M_O <- ScaleData(object = Aorta_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_O <- FindVariableGenes(object = Aorta_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_M_O <- RunPCA(Aorta_M_O, verbose = FALSE)
Aorta_M_O <- FindClusters(Aorta_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_M_O <- RunTSNE(Aorta_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Aorta_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_M_O <- doubletFinder(Aorta_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_M_O <- doubletFinder(Aorta_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Aorta_M_O@meta.data$Doublet <- Aorta_M_O@meta.data$DF.classifications_0.25_0.12_98

############################Aorta_M_CR
Aorta_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_M_CR/")
dense.size <- object.size(x = as.matrix(x = Aorta_M_CR.data))
sparse.size <- object.size(x = Aorta_M_CR.data)
Aorta_M_CR <- CreateSeuratObject(raw.data = Aorta_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_M_CR")
Aorta_M_CR@meta.data$gender<- "M"
Aorta_M_CR@meta.data$age<- "CR"
Aorta_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_M_CR@raw.data[mito.gens, ])/Matrix::colSums(Aorta_M_CR@raw.data)
Aorta_M_CR <- AddMetaData(object = Aorta_M_CR, metadata = percent.mito, col.name = "percent.mito")
Aorta_M_CR <- FilterCells(object = Aorta_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_M_CR <- NormalizeData(object = Aorta_M_CR, normalization.method = "LogNormalize")
Aorta_M_CR <- ScaleData(object = Aorta_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_CR <- FindVariableGenes(object = Aorta_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_M_CR <- RunPCA(Aorta_M_CR, verbose = FALSE)
Aorta_M_CR <- FindClusters(Aorta_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_M_CR <- RunTSNE(Aorta_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Aorta_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_M_CR <- doubletFinder(Aorta_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_M_CR <- doubletFinder(Aorta_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
Aorta_M_CR@meta.data$Doublet <- Aorta_M_CR@meta.data$DF.classifications_0.25_0.1_110

############################Aorta_F_Y
Aorta_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_F_Y/")
dense.size <- object.size(x = as.matrix(x = Aorta_F_Y.data))
sparse.size <- object.size(x = Aorta_F_Y.data)
Aorta_F_Y <- CreateSeuratObject(raw.data = Aorta_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_F_Y")
Aorta_F_Y@meta.data$gender<- "F"
Aorta_F_Y@meta.data$age<- "Y"
Aorta_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_F_Y@raw.data[mito.gens, ])/Matrix::colSums(Aorta_F_Y@raw.data)
Aorta_F_Y <- AddMetaData(object = Aorta_F_Y, metadata = percent.mito, col.name = "percent.mito")
Aorta_F_Y <- FilterCells(object = Aorta_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_F_Y <- NormalizeData(object = Aorta_F_Y, normalization.method = "LogNormalize")
Aorta_F_Y <- ScaleData(object = Aorta_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_Y <- FindVariableGenes(object = Aorta_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_Y <- RunPCA(Aorta_F_Y, verbose = FALSE)
Aorta_F_Y <- FindClusters(Aorta_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_F_Y <- RunTSNE(Aorta_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.054*length(Aorta_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_F_Y <- doubletFinder(Aorta_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_F_Y <- doubletFinder(Aorta_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Aorta_F_Y@meta.data$Doublet <- Aorta_F_Y@meta.data$DF.classifications_0.25_0.09_325

############################Aorta_F_O
Aorta_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_F_O/")
dense.size <- object.size(x = as.matrix(x = Aorta_F_O.data))
sparse.size <- object.size(x = Aorta_F_O.data)
Aorta_F_O <- CreateSeuratObject(raw.data = Aorta_F_O.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_F_O")
Aorta_F_O@meta.data$gender<- "F"
Aorta_F_O@meta.data$age<- "O"
Aorta_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_F_O@raw.data[mito.gens, ])/Matrix::colSums(Aorta_F_O@raw.data)
Aorta_F_O <- AddMetaData(object = Aorta_F_O, metadata = percent.mito, col.name = "percent.mito")
Aorta_F_O <- FilterCells(object = Aorta_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_F_O <- NormalizeData(object = Aorta_F_O, normalization.method = "LogNormalize")
Aorta_F_O <- ScaleData(object = Aorta_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_O <- FindVariableGenes(object = Aorta_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_O <- RunPCA(Aorta_F_O, verbose = FALSE)
Aorta_F_O <- FindClusters(Aorta_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_F_O <- RunTSNE(Aorta_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Aorta_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_F_O <- doubletFinder(Aorta_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_F_O <- doubletFinder(Aorta_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Aorta_F_O@meta.data$Doublet <- Aorta_F_O@meta.data$DF.classifications_0.25_0.06_103

############################Aorta_F_CR
Aorta_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Vessel/Vessel_F_CR/")
dense.size <- object.size(x = as.matrix(x = Aorta_F_CR.data))
sparse.size <- object.size(x = Aorta_F_CR.data)
Aorta_F_CR <- CreateSeuratObject(raw.data = Aorta_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_Aorta_F_CR")
Aorta_F_CR@meta.data$gender<- "F"
Aorta_F_CR@meta.data$age<- "CR"
Aorta_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Aorta_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Aorta_F_CR@raw.data[mito.gens, ])/Matrix::colSums(Aorta_F_CR@raw.data)
Aorta_F_CR <- AddMetaData(object = Aorta_F_CR, metadata = percent.mito, col.name = "percent.mito")
Aorta_F_CR <- FilterCells(object = Aorta_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Aorta_F_CR <- NormalizeData(object = Aorta_F_CR, normalization.method = "LogNormalize")
Aorta_F_CR <- ScaleData(object = Aorta_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_CR <- FindVariableGenes(object = Aorta_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_CR <- RunPCA(Aorta_F_CR, verbose = FALSE)
Aorta_F_CR <- FindClusters(Aorta_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Aorta_F_CR <- RunTSNE(Aorta_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Aorta_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Aorta_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Aorta_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Aorta_F_CR <- doubletFinder(Aorta_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Aorta_F_CR <- doubletFinder(Aorta_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Aorta_F_CR@meta.data$Doublet <- Aorta_F_CR@meta.data$DF.classifications_0.25_0.07_52

###################################################  Tissue: Skin
############################Skin_M_Y
Skin_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Skin/Skin_M_Y/")
dense.size <- object.size(x = as.matrix(x = Skin_M_Y.data))
sparse.size <- object.size(x = Skin_M_Y.data)
Skin_M_Y <- CreateSeuratObject(raw.data = Skin_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_Skin_M_Y")
Skin_M_Y@meta.data$gender<- "M"
Skin_M_Y@meta.data$age<- "Y"
Skin_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_M_Y@raw.data[mito.gens, ])/Matrix::colSums(Skin_M_Y@raw.data)
Skin_M_Y <- AddMetaData(object = Skin_M_Y, metadata = percent.mito, col.name = "percent.mito")
Skin_M_Y <- FilterCells(object = Skin_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_M_Y <- NormalizeData(object = Skin_M_Y, normalization.method = "LogNormalize")
Skin_M_Y <- ScaleData(object = Skin_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_Y <- FindVariableGenes(object = Skin_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_M_Y <- RunPCA(Skin_M_Y, verbose = FALSE)
Skin_M_Y <- FindClusters(Skin_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_M_Y <- RunTSNE(Skin_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Skin_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_M_Y <- doubletFinder(Skin_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_M_Y <- doubletFinder(Skin_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Skin_M_Y@meta.data$Doublet <- Skin_M_Y@meta.data$DF.classifications_0.25_0.3_64

############################Skin_M_O
Skin_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/pifu_L3_X41/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = Skin_M_O.data))
sparse.size <- object.size(x = Skin_M_O.data)
Skin_M_O <- CreateSeuratObject(raw.data = Skin_M_O.data, min.cells = 5, min.genes = 200, project = "10X_Skin_M_O")
Skin_M_O@meta.data$gender<- "M"
Skin_M_O@meta.data$age<- "O"
Skin_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_M_O@raw.data[mito.gens, ])/Matrix::colSums(Skin_M_O@raw.data)
Skin_M_O <- AddMetaData(object = Skin_M_O, metadata = percent.mito, col.name = "percent.mito")
Skin_M_O <- FilterCells(object = Skin_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_M_O <- NormalizeData(object = Skin_M_O, normalization.method = "LogNormalize")
Skin_M_O <- ScaleData(object = Skin_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_O <- FindVariableGenes(object = Skin_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_M_O <- RunPCA(Skin_M_O, verbose = FALSE)
Skin_M_O <- FindClusters(Skin_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_M_O <- RunTSNE(Skin_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(Skin_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_M_O <- doubletFinder(Skin_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_M_O <- doubletFinder(Skin_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
Skin_M_O@meta.data$Doublet <- Skin_M_O@meta.data$DF.classifications_0.25_0.005_60

############################Skin_M_CR
Skin_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Skin/Skin_M_CR/")
dense.size <- object.size(x = as.matrix(x = Skin_M_CR.data))
sparse.size <- object.size(x = Skin_M_CR.data)
Skin_M_CR <- CreateSeuratObject(raw.data = Skin_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_Skin_M_CR")
Skin_M_CR@meta.data$gender<- "M"
Skin_M_CR@meta.data$age<- "CR"
Skin_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_M_CR@raw.data[mito.gens, ])/Matrix::colSums(Skin_M_CR@raw.data)
Skin_M_CR <- AddMetaData(object = Skin_M_CR, metadata = percent.mito, col.name = "percent.mito")
Skin_M_CR <- FilterCells(object = Skin_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_M_CR <- NormalizeData(object = Skin_M_CR, normalization.method = "LogNormalize")
Skin_M_CR <- ScaleData(object = Skin_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_CR <- FindVariableGenes(object = Skin_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_M_CR <- RunPCA(Skin_M_CR, verbose = FALSE)
Skin_M_CR <- FindClusters(Skin_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_M_CR <- RunTSNE(Skin_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.054*length(Skin_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_M_CR <- doubletFinder(Skin_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_M_CR <- doubletFinder(Skin_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Skin_M_CR@meta.data$Doublet <- Skin_M_CR@meta.data$DF.classifications_0.25_0.005_319

############################Skin_F_Y
Skin_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/pifu-5/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = Skin_F_Y.data))
sparse.size <- object.size(x = Skin_F_Y.data)
Skin_F_Y <- CreateSeuratObject(raw.data = Skin_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_Skin_F_Y")
Skin_F_Y@meta.data$gender<- "F"
Skin_F_Y@meta.data$age<- "Y"
Skin_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_F_Y@raw.data[mito.gens, ])/Matrix::colSums(Skin_F_Y@raw.data)
Skin_F_Y <- AddMetaData(object = Skin_F_Y, metadata = percent.mito, col.name = "percent.mito")
Skin_F_Y <- FilterCells(object = Skin_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_F_Y <- NormalizeData(object = Skin_F_Y, normalization.method = "LogNormalize")
Skin_F_Y <- ScaleData(object = Skin_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_Y <- FindVariableGenes(object = Skin_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_Y <- RunPCA(Skin_F_Y, verbose = FALSE)
Skin_F_Y <- FindClusters(Skin_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_F_Y <- RunTSNE(Skin_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(Skin_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_F_Y <- doubletFinder(Skin_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_F_Y <- doubletFinder(Skin_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Skin_F_Y@meta.data$Doublet <- Skin_F_Y@meta.data$DF.classifications_0.25_0.005_108

############################Skin_F_O
Skin_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Skin/Skin_F_O/")
dense.size <- object.size(x = as.matrix(x = Skin_F_O.data))
sparse.size <- object.size(x = Skin_F_O.data)
Skin_F_O <- CreateSeuratObject(raw.data = Skin_F_O.data, min.cells = 5, min.genes = 200, project = "10X_Skin_F_O")
Skin_F_O@meta.data$gender<- "F"
Skin_F_O@meta.data$age<- "O"
Skin_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_F_O@raw.data[mito.gens, ])/Matrix::colSums(Skin_F_O@raw.data)
Skin_F_O <- AddMetaData(object = Skin_F_O, metadata = percent.mito, col.name = "percent.mito")
Skin_F_O <- FilterCells(object = Skin_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_F_O <- NormalizeData(object = Skin_F_O, normalization.method = "LogNormalize")
Skin_F_O <- ScaleData(object = Skin_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_O <- FindVariableGenes(object = Skin_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_O <- RunPCA(Skin_F_O, verbose = FALSE)
Skin_F_O <- FindClusters(Skin_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_F_O <- RunTSNE(Skin_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(Skin_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_F_O <- doubletFinder(Skin_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_F_O <- doubletFinder(Skin_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
Skin_F_O@meta.data$Doublet <- Skin_F_O@meta.data$DF.classifications_0.25_0.01_181

############################Skin_F_CR
Skin_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Skin/Skin_F_CR/")
dense.size <- object.size(x = as.matrix(x = Skin_F_CR.data))
sparse.size <- object.size(x = Skin_F_CR.data)
Skin_F_CR <- CreateSeuratObject(raw.data = Skin_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_Skin_F_CR")
Skin_F_CR@meta.data$gender<- "F"
Skin_F_CR@meta.data$age<- "CR"
Skin_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = Skin_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(Skin_F_CR@raw.data[mito.gens, ])/Matrix::colSums(Skin_F_CR@raw.data)
Skin_F_CR <- AddMetaData(object = Skin_F_CR, metadata = percent.mito, col.name = "percent.mito")
Skin_F_CR <- FilterCells(object = Skin_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
Skin_F_CR <- NormalizeData(object = Skin_F_CR, normalization.method = "LogNormalize")
Skin_F_CR <- ScaleData(object = Skin_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_CR <- FindVariableGenes(object = Skin_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_CR <- RunPCA(Skin_F_CR, verbose = FALSE)
Skin_F_CR <- FindClusters(Skin_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
Skin_F_CR <- RunTSNE(Skin_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(Skin_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- Skin_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(Skin_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
Skin_F_CR <- doubletFinder(Skin_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
Skin_F_CR <- doubletFinder(Skin_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
Skin_F_CR@meta.data$Doublet <- Skin_F_CR@meta.data$DF.classifications_0.25_0.07_158

###################################################  Tissue: BM
############################BM_M_Y
BM_M_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Marrow/Marrow_M_Y/")
dense.size <- object.size(x = as.matrix(x = BM_M_Y.data))
sparse.size <- object.size(x = BM_M_Y.data)
BM_M_Y <- CreateSeuratObject(raw.data = BM_M_Y.data, min.cells = 5, min.genes = 200, project = "10X_BM_M_Y")
BM_M_Y@meta.data$gender<- "M"
BM_M_Y@meta.data$age<- "Y"
BM_M_Y@meta.data$sample<- "M.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_M_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_M_Y@raw.data[mito.gens, ])/Matrix::colSums(BM_M_Y@raw.data)
BM_M_Y <- AddMetaData(object = BM_M_Y, metadata = percent.mito, col.name = "percent.mito")
BM_M_Y <- FilterCells(object = BM_M_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_M_Y <- NormalizeData(object = BM_M_Y, normalization.method = "LogNormalize")
BM_M_Y <- ScaleData(object = BM_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_Y <- FindVariableGenes(object = BM_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_M_Y <- RunPCA(BM_M_Y, verbose = FALSE)
#DimHeatmap(object = BM_M_Y, reduction.type = "pca", cells.use = 500, dim.use = 1:12, do.balanced = TRUE) 
BM_M_Y <- FindClusters(BM_M_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_M_Y <- RunTSNE(BM_M_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_M_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_M_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.023*length(BM_M_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_M_Y <- doubletFinder(BM_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_M_Y <- doubletFinder(BM_M_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
BM_M_Y@meta.data$Doublet <- BM_M_Y@meta.data$DF.classifications_0.25_0.05_68

############################BM_M_O
BM_M_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/buce/result/gusui_L3_X53/outs/filtered_gene_bc_matrices/refdata-cellranger-rat6/")
dense.size <- object.size(x = as.matrix(x = BM_M_O.data))
sparse.size <- object.size(x = BM_M_O.data)
BM_M_O <- CreateSeuratObject(raw.data = BM_M_O.data, min.cells = 5, min.genes = 200, project = "10X_BM_M_O")
BM_M_O@meta.data$gender<- "M"
BM_M_O@meta.data$age<- "O"
BM_M_O@meta.data$sample<- "M.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_M_O@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_M_O@raw.data[mito.gens, ])/Matrix::colSums(BM_M_O@raw.data)
BM_M_O <- AddMetaData(object = BM_M_O, metadata = percent.mito, col.name = "percent.mito")
BM_M_O <- FilterCells(object = BM_M_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_M_O <- NormalizeData(object = BM_M_O, normalization.method = "LogNormalize")
BM_M_O <- ScaleData(object = BM_M_O, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_O <- FindVariableGenes(object = BM_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_M_O <- RunPCA(BM_M_O, verbose = FALSE)
BM_M_O <- FindClusters(BM_M_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_M_O <- RunTSNE(BM_M_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_M_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_M_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.039*length(BM_M_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_M_O <- doubletFinder(BM_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_M_O <- doubletFinder(BM_M_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
BM_M_O@meta.data$Doublet <- BM_M_O@meta.data$DF.classifications_0.25_0.03_156

############################BM_M_CR
BM_M_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Marrow/Marrow_M_CR/")
dense.size <- object.size(x = as.matrix(x = BM_M_CR.data))
sparse.size <- object.size(x = BM_M_CR.data)
BM_M_CR <- CreateSeuratObject(raw.data = BM_M_CR.data, min.cells = 5, min.genes = 200, project = "10X_BM_M_CR")
BM_M_CR@meta.data$gender<- "M"
BM_M_CR@meta.data$age<- "CR"
BM_M_CR@meta.data$sample<- "M.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_M_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_M_CR@raw.data[mito.gens, ])/Matrix::colSums(BM_M_CR@raw.data)
BM_M_CR <- AddMetaData(object = BM_M_CR, metadata = percent.mito, col.name = "percent.mito")
BM_M_CR <- FilterCells(object = BM_M_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_M_CR <- NormalizeData(object = BM_M_CR, normalization.method = "LogNormalize")
BM_M_CR <- ScaleData(object = BM_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_CR <- FindVariableGenes(object = BM_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_M_CR <- RunPCA(BM_M_CR, verbose = FALSE)
BM_M_CR <- FindClusters(BM_M_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_M_CR <- RunTSNE(BM_M_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_M_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_M_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BM_M_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_M_CR <- doubletFinder(BM_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_M_CR <- doubletFinder(BM_M_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
BM_M_CR@meta.data$Doublet <- BM_M_CR@meta.data$DF.classifications_0.25_0.005_108

############################BM_F_Y
BM_F_Y.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Marrow/Marrow_F_Y/")
dense.size <- object.size(x = as.matrix(x = BM_F_Y.data))
sparse.size <- object.size(x = BM_F_Y.data)
BM_F_Y <- CreateSeuratObject(raw.data = BM_F_Y.data, min.cells = 5, min.genes = 200, project = "10X_BM_F_Y")
BM_F_Y@meta.data$gender<- "F"
BM_F_Y@meta.data$age<- "Y"
BM_F_Y@meta.data$sample<- "F.Y"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_F_Y@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_F_Y@raw.data[mito.gens, ])/Matrix::colSums(BM_F_Y@raw.data)
BM_F_Y <- AddMetaData(object = BM_F_Y, metadata = percent.mito, col.name = "percent.mito")
BM_F_Y <- FilterCells(object = BM_F_Y, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_F_Y <- NormalizeData(object = BM_F_Y, normalization.method = "LogNormalize")
BM_F_Y <- ScaleData(object = BM_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_Y <- FindVariableGenes(object = BM_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_Y <- RunPCA(BM_F_Y, verbose = FALSE)
BM_F_Y <- FindClusters(BM_F_Y, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_F_Y <- RunTSNE(BM_F_Y, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_F_Y, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_F_Y@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BM_F_Y@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_F_Y <- doubletFinder(BM_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_F_Y <- doubletFinder(BM_F_Y, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)
BM_F_Y@meta.data$Doublet <- BM_F_Y@meta.data$DF.classifications_0.25_0.24_241

############################BM_F_O
BM_F_O.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Marrow/Marrow_F_O/")
dense.size <- object.size(x = as.matrix(x = BM_F_O.data))
sparse.size <- object.size(x = BM_F_O.data)
BM_F_O <- CreateSeuratObject(raw.data = BM_F_O.data, min.cells = 5, min.genes = 200, project = "10X_BM_F_O")
BM_F_O@meta.data$gender<- "F"
BM_F_O@meta.data$age<- "O"
BM_F_O@meta.data$sample<- "F.O"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_F_O@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_F_O@raw.data[mito.gens, ])/Matrix::colSums(BM_F_O@raw.data)
BM_F_O <- AddMetaData(object = BM_F_O, metadata = percent.mito, col.name = "percent.mito")
BM_F_O <- FilterCells(object = BM_F_O, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_F_O <- NormalizeData(object = BM_F_O, normalization.method = "LogNormalize")
BM_F_O <- ScaleData(object = BM_F_O, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_O <- FindVariableGenes(object = BM_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_O <- RunPCA(BM_F_O, verbose = FALSE)
BM_F_O <- FindClusters(BM_F_O, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_F_O <- RunTSNE(BM_F_O, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_F_O, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_F_O@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BM_F_O@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_F_O <- doubletFinder(BM_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_F_O <- doubletFinder(BM_F_O, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
BM_F_O@meta.data$Doublet <- BM_F_O@meta.data$DF.classifications_0.25_0.3_116

############################BM_F_CR
BM_F_CR.data <- Read10X(data.dir = "/data1/mashuai/data/rat_10X/matrix/Marrow/Marrow_F_CR/")
dense.size <- object.size(x = as.matrix(x = BM_F_CR.data))
sparse.size <- object.size(x = BM_F_CR.data)
BM_F_CR <- CreateSeuratObject(raw.data = BM_F_CR.data, min.cells = 5, min.genes = 200, project = "10X_BM_F_CR")
BM_F_CR@meta.data$gender<- "F"
BM_F_CR@meta.data$age<- "CR"
BM_F_CR@meta.data$sample<- "F.CR"
mito.gens <- grep(pattern = "^Mt-", x = rownames(x = BM_F_CR@data), value = TRUE)
percent.mito <- Matrix::colSums(BM_F_CR@raw.data[mito.gens, ])/Matrix::colSums(BM_F_CR@raw.data)
BM_F_CR <- AddMetaData(object = BM_F_CR, metadata = percent.mito, col.name = "percent.mito")
BM_F_CR <- FilterCells(object = BM_F_CR, subset.names = c("nGene", "percent.mito"),low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.1))
BM_F_CR <- NormalizeData(object = BM_F_CR, normalization.method = "LogNormalize")
BM_F_CR <- ScaleData(object = BM_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_CR <- FindVariableGenes(object = BM_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_CR <- RunPCA(BM_F_CR, verbose = FALSE)
BM_F_CR <- FindClusters(BM_F_CR, reduction.type =  "pca", resolution = 0.6, dims.use = 1:12)
BM_F_CR <- RunTSNE(BM_F_CR, reduction.use =  "pca", dims.use = 1:12, do.fast = T) 
sweep.res.list_SM <- paramSweep(BM_F_CR, PCs = 1:12)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- BM_F_CR@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(BM_F_CR@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
BM_F_CR <- doubletFinder(BM_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
BM_F_CR <- doubletFinder(BM_F_CR, PCs = 1:12, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value) 
BM_F_CR@meta.data$Doublet <- BM_F_CR@meta.data$DF.classifications_0.25_0.01_115

#############################  SaveRDS
saveRDS(BAT_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_Y.rds")
saveRDS(BAT_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_O.rds")
saveRDS(BAT_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_CR.rds")
saveRDS(WAT_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_Y.rds")
saveRDS(WAT_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_O.rds")
saveRDS(WAT_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_CR.rds")
saveRDS(Liver_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_Y.rds")
saveRDS(Liver_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_O.rds")
saveRDS(Liver_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_CR.rds")
saveRDS(Kidney_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_Y.rds")
saveRDS(Kidney_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_O.rds")
saveRDS(Kidney_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_CR.rds")
saveRDS(Aorta_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_Y.rds")
saveRDS(Aorta_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_O.rds")
saveRDS(Aorta_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_CR.rds")
saveRDS(Skin_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_Y.rds")
saveRDS(Skin_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_O.rds")
saveRDS(Skin_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_CR.rds")
saveRDS(BM_M_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_Y.rds")
saveRDS(BM_M_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_O.rds")
saveRDS(BM_M_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_CR.rds")

saveRDS(BAT_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_Y.rds")
saveRDS(BAT_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_O.rds")
saveRDS(BAT_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_CR.rds")
saveRDS(WAT_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_Y.rds")
saveRDS(WAT_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_O.rds")
saveRDS(WAT_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_CR.rds")
saveRDS(Liver_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_Y.rds")
saveRDS(Liver_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_O.rds")
saveRDS(Liver_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_CR.rds")
saveRDS(Kidney_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_Y.rds")
saveRDS(Kidney_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_O.rds")
saveRDS(Kidney_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_CR.rds")
saveRDS(Aorta_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_Y.rds")
saveRDS(Aorta_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_O.rds")
saveRDS(Aorta_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_CR.rds")
saveRDS(Skin_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_Y.rds")
saveRDS(Skin_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_O.rds")
saveRDS(Skin_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_CR.rds")
saveRDS(BM_F_Y, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_Y.rds")
saveRDS(BM_F_O, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_O.rds")
saveRDS(BM_F_CR, file = "/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_CR.rds")

BAT_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_Y.rds")
BAT_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_O.rds")
BAT_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_M_CR.rds")
WAT_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_Y.rds")
WAT_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_O.rds")
WAT_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_M_CR.rds")
Liver_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_Y.rds")
Liver_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_O.rds")
Liver_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_M_CR.rds")
Kidney_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_Y.rds")
Kidney_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_O.rds")
Kidney_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_M_CR.rds")
Aorta_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_Y.rds")
Aorta_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_O.rds")
Aorta_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_M_CR.rds")
Skin_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_Y.rds")
Skin_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_O.rds")
Skin_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_M_CR.rds")
BM_M_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_Y.rds")
BM_M_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_O.rds")
BM_M_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_M_CR.rds")

BAT_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_Y.rds")
BAT_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_O.rds")
BAT_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BAT_F_CR.rds")
WAT_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_Y.rds")
WAT_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_O.rds")
WAT_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/WAT_F_CR.rds")
Liver_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_Y.rds")
Liver_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_O.rds")
Liver_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Liver_F_CR.rds")
Kidney_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_Y.rds")
Kidney_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_O.rds")
Kidney_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Kidney_F_CR.rds")
Aorta_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_Y.rds")
Aorta_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_O.rds")
Aorta_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Aorta_F_CR.rds")
Skin_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_Y.rds")
Skin_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_O.rds")
Skin_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/Skin_F_CR.rds")
BM_F_Y <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_Y.rds")
BM_F_O <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_O.rds")
BM_F_CR <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/samples/BM_F_CR.rds")

######################### Subset  
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

WAT_M_Y <- SubsetData(WAT_M_Y,subset.name='Doublet',accept.value='Singlet')
WAT_M_Y <- NormalizeData(object = WAT_M_Y, normalization.method = "LogNormalize")
WAT_M_Y <- ScaleData(object = WAT_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_Y <- FindVariableGenes(object = WAT_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_M_O <- SubsetData(WAT_M_O,subset.name='Doublet',accept.value='Singlet')
WAT_M_O <- NormalizeData(object = WAT_M_O, normalization.method = "LogNormalize")
WAT_M_O <- ScaleData(object = WAT_M_O, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_O <- FindVariableGenes(object = WAT_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_M_CR <- SubsetData(WAT_M_CR,subset.name='Doublet',accept.value='Singlet')
WAT_M_CR <- NormalizeData(object = WAT_M_CR, normalization.method = "LogNormalize")
WAT_M_CR <- ScaleData(object = WAT_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
WAT_M_CR <- FindVariableGenes(object = WAT_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_Y <- SubsetData(WAT_F_Y,subset.name='Doublet',accept.value='Singlet')
WAT_F_Y <- NormalizeData(object = WAT_F_Y, normalization.method = "LogNormalize")
WAT_F_Y <- ScaleData(object = WAT_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_Y <- FindVariableGenes(object = WAT_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_O <- SubsetData(WAT_F_O,subset.name='Doublet',accept.value='Singlet')
WAT_F_O <- NormalizeData(object = WAT_F_O, normalization.method = "LogNormalize")
WAT_F_O <- ScaleData(object = WAT_F_O, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_O <- FindVariableGenes(object = WAT_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
WAT_F_CR <- SubsetData(WAT_F_CR,subset.name='Doublet',accept.value='Singlet')
WAT_F_CR <- NormalizeData(object = WAT_F_CR, normalization.method = "LogNormalize")
WAT_F_CR <- ScaleData(object = WAT_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
WAT_F_CR <- FindVariableGenes(object = WAT_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

Liver_M_Y <- SubsetData(Liver_M_Y,subset.name='Doublet',accept.value='Singlet')
Liver_M_Y <- NormalizeData(object = Liver_M_Y, normalization.method = "LogNormalize")
Liver_M_Y <- ScaleData(object = Liver_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_Y <- FindVariableGenes(object = Liver_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_M_O <- SubsetData(Liver_M_O,subset.name='Doublet',accept.value='Singlet')
Liver_M_O <- NormalizeData(object = Liver_M_O, normalization.method = "LogNormalize")
Liver_M_O <- ScaleData(object = Liver_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_O <- FindVariableGenes(object = Liver_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_M_CR <- SubsetData(Liver_M_CR,subset.name='Doublet',accept.value='Singlet')
Liver_M_CR <- NormalizeData(object = Liver_M_CR, normalization.method = "LogNormalize")
Liver_M_CR <- ScaleData(object = Liver_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Liver_M_CR <- FindVariableGenes(object = Liver_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_Y <- SubsetData(Liver_F_Y,subset.name='Doublet',accept.value='Singlet')
Liver_F_Y <- NormalizeData(object = Liver_F_Y, normalization.method = "LogNormalize")
Liver_F_Y <- ScaleData(object = Liver_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_Y <- FindVariableGenes(object = Liver_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_O <- SubsetData(Liver_F_O,subset.name='Doublet',accept.value='Singlet')
Liver_F_O <- NormalizeData(object = Liver_F_O, normalization.method = "LogNormalize")
Liver_F_O <- ScaleData(object = Liver_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_O <- FindVariableGenes(object = Liver_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Liver_F_CR <- SubsetData(Liver_F_CR,subset.name='Doublet',accept.value='Singlet')
Liver_F_CR <- NormalizeData(object = Liver_F_CR, normalization.method = "LogNormalize")
Liver_F_CR <- ScaleData(object = Liver_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Liver_F_CR <- FindVariableGenes(object = Liver_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

Kidney_M_Y <- SubsetData(Kidney_M_Y,subset.name='Doublet',accept.value='Singlet')
Kidney_M_Y <- NormalizeData(object = Kidney_M_Y, normalization.method = "LogNormalize")
Kidney_M_Y <- ScaleData(object = Kidney_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_Y <- FindVariableGenes(object = Kidney_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_M_O <- SubsetData(Kidney_M_O,subset.name='Doublet',accept.value='Singlet')
Kidney_M_O <- NormalizeData(object = Kidney_M_O, normalization.method = "LogNormalize")
Kidney_M_O <- ScaleData(object = Kidney_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_O <- FindVariableGenes(object = Kidney_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_M_CR <- SubsetData(Kidney_M_CR,subset.name='Doublet',accept.value='Singlet')
Kidney_M_CR <- NormalizeData(object = Kidney_M_CR, normalization.method = "LogNormalize")
Kidney_M_CR <- ScaleData(object = Kidney_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_M_CR <- FindVariableGenes(object = Kidney_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_Y <- SubsetData(Kidney_F_Y,subset.name='Doublet',accept.value='Singlet')
Kidney_F_Y <- NormalizeData(object = Kidney_F_Y, normalization.method = "LogNormalize")
Kidney_F_Y <- ScaleData(object = Kidney_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_Y <- FindVariableGenes(object = Kidney_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_O <- SubsetData(Kidney_F_O,subset.name='Doublet',accept.value='Singlet')
Kidney_F_O <- NormalizeData(object = Kidney_F_O, normalization.method = "LogNormalize")
Kidney_F_O <- ScaleData(object = Kidney_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_O <- FindVariableGenes(object = Kidney_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Kidney_F_CR <- SubsetData(Kidney_F_CR,subset.name='Doublet',accept.value='Singlet')
Kidney_F_CR <- NormalizeData(object = Kidney_F_CR, normalization.method = "LogNormalize")
Kidney_F_CR <- ScaleData(object = Kidney_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Kidney_F_CR <- FindVariableGenes(object = Kidney_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

Aorta_M_Y <- SubsetData(Aorta_M_Y,subset.name='Doublet',accept.value='Singlet')
Aorta_M_Y <- NormalizeData(object = Aorta_M_Y, normalization.method = "LogNormalize")
Aorta_M_Y <- ScaleData(object = Aorta_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_Y <- FindVariableGenes(object = Aorta_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_M_O <- SubsetData(Aorta_M_O,subset.name='Doublet',accept.value='Singlet')
Aorta_M_O <- NormalizeData(object = Aorta_M_O, normalization.method = "LogNormalize")
Aorta_M_O <- ScaleData(object = Aorta_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_O <- FindVariableGenes(object = Aorta_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_M_CR <- SubsetData(Aorta_M_CR,subset.name='Doublet',accept.value='Singlet')
Aorta_M_CR <- NormalizeData(object = Aorta_M_CR, normalization.method = "LogNormalize")
Aorta_M_CR <- ScaleData(object = Aorta_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_M_CR <- FindVariableGenes(object = Aorta_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_Y <- SubsetData(Aorta_F_Y,subset.name='Doublet',accept.value='Singlet')
Aorta_F_Y <- NormalizeData(object = Aorta_F_Y, normalization.method = "LogNormalize")
Aorta_F_Y <- ScaleData(object = Aorta_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_Y <- FindVariableGenes(object = Aorta_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_O <- SubsetData(Aorta_F_O,subset.name='Doublet',accept.value='Singlet')
Aorta_F_O <- NormalizeData(object = Aorta_F_O, normalization.method = "LogNormalize")
Aorta_F_O <- ScaleData(object = Aorta_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_O <- FindVariableGenes(object = Aorta_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Aorta_F_CR <- SubsetData(Aorta_F_CR,subset.name='Doublet',accept.value='Singlet')
Aorta_F_CR <- NormalizeData(object = Aorta_F_CR, normalization.method = "LogNormalize")
Aorta_F_CR <- ScaleData(object = Aorta_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Aorta_F_CR <- FindVariableGenes(object = Aorta_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

Skin_M_Y <- SubsetData(Skin_M_Y,subset.name='Doublet',accept.value='Singlet')
Skin_M_Y <- NormalizeData(object = Skin_M_Y, normalization.method = "LogNormalize")
Skin_M_Y <- ScaleData(object = Skin_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_Y <- FindVariableGenes(object = Skin_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_M_O <- SubsetData(Skin_M_O,subset.name='Doublet',accept.value='Singlet')
Skin_M_O <- NormalizeData(object = Skin_M_O, normalization.method = "LogNormalize")
Skin_M_O <- ScaleData(object = Skin_M_O, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_O <- FindVariableGenes(object = Skin_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_M_CR <- SubsetData(Skin_M_CR,subset.name='Doublet',accept.value='Singlet')
Skin_M_CR <- NormalizeData(object = Skin_M_CR, normalization.method = "LogNormalize")
Skin_M_CR <- ScaleData(object = Skin_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
Skin_M_CR <- FindVariableGenes(object = Skin_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_Y <- SubsetData(Skin_F_Y,subset.name='Doublet',accept.value='Singlet')
Skin_F_Y <- NormalizeData(object = Skin_F_Y, normalization.method = "LogNormalize")
Skin_F_Y <- ScaleData(object = Skin_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_Y <- FindVariableGenes(object = Skin_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_O <- SubsetData(Skin_F_O,subset.name='Doublet',accept.value='Singlet')
Skin_F_O <- NormalizeData(object = Skin_F_O, normalization.method = "LogNormalize")
Skin_F_O <- ScaleData(object = Skin_F_O, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_O <- FindVariableGenes(object = Skin_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
Skin_F_CR <- SubsetData(Skin_F_CR,subset.name='Doublet',accept.value='Singlet')
Skin_F_CR <- NormalizeData(object = Skin_F_CR, normalization.method = "LogNormalize")
Skin_F_CR <- ScaleData(object = Skin_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
Skin_F_CR <- FindVariableGenes(object = Skin_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

BM_M_Y <- SubsetData(BM_M_Y,subset.name='Doublet',accept.value='Singlet')
BM_M_Y <- NormalizeData(object = BM_M_Y, normalization.method = "LogNormalize")
BM_M_Y <- ScaleData(object = BM_M_Y, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_Y <- FindVariableGenes(object = BM_M_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_M_O <- SubsetData(BM_M_O,subset.name='Doublet',accept.value='Singlet')
BM_M_O <- NormalizeData(object = BM_M_O, normalization.method = "LogNormalize")
BM_M_O <- ScaleData(object = BM_M_O, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_O <- FindVariableGenes(object = BM_M_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_M_CR <- SubsetData(BM_M_CR,subset.name='Doublet',accept.value='Singlet')
BM_M_CR <- NormalizeData(object = BM_M_CR, normalization.method = "LogNormalize")
BM_M_CR <- ScaleData(object = BM_M_CR, vars.to.regress = c("nUMI", "percent.mito"))
BM_M_CR <- FindVariableGenes(object = BM_M_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_Y <- SubsetData(BM_F_Y,subset.name='Doublet',accept.value='Singlet')
BM_F_Y <- NormalizeData(object = BM_F_Y, normalization.method = "LogNormalize")
BM_F_Y <- ScaleData(object = BM_F_Y, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_Y <- FindVariableGenes(object = BM_F_Y, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_O <- SubsetData(BM_F_O,subset.name='Doublet',accept.value='Singlet')
BM_F_O <- NormalizeData(object = BM_F_O, normalization.method = "LogNormalize")
BM_F_O <- ScaleData(object = BM_F_O, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_O <- FindVariableGenes(object = BM_F_O, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)
BM_F_CR <- SubsetData(BM_F_CR,subset.name='Doublet',accept.value='Singlet')
BM_F_CR <- NormalizeData(object = BM_F_CR, normalization.method = "LogNormalize")
BM_F_CR <- ScaleData(object = BM_F_CR, vars.to.regress = c("nUMI", "percent.mito"))
BM_F_CR <- FindVariableGenes(object = BM_F_CR, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125 , x.high.cutoff = 3, y.cutoff = 0.5)

##########################################

g.1 <- head(rownames(BAT_M_Y@hvg.info), 1000) 
g.2 <- head(rownames(BAT_M_O@hvg.info), 1000) 
g.3 <- head(rownames(BAT_M_CR@hvg.info), 1000) 
g.4 <- head(rownames(WAT_M_Y@hvg.info), 1000) 
g.5 <- head(rownames(WAT_M_O@hvg.info), 1000) 
g.6 <- head(rownames(WAT_M_CR@hvg.info), 1000)
g.7 <- head(rownames(Skin_M_Y@hvg.info), 1000) 
g.8 <- head(rownames(Skin_M_O@hvg.info), 1000) 
g.9 <- head(rownames(Skin_M_CR@hvg.info), 1000)
g.10 <- head(rownames(Kidney_M_Y@hvg.info), 1000) 
g.11 <- head(rownames(Kidney_M_O@hvg.info), 1000) 
g.12 <- head(rownames(Kidney_M_CR@hvg.info), 1000)
g.13 <- head(rownames(Aorta_M_Y@hvg.info), 1000) 
g.14 <- head(rownames(Aorta_M_O@hvg.info), 1000) 
g.15 <- head(rownames(Aorta_M_CR@hvg.info), 1000)
g.16 <- head(rownames(Liver_M_Y@hvg.info), 1000) 
g.17 <- head(rownames(Liver_M_O@hvg.info), 1000) 
g.18 <- head(rownames(Liver_M_CR@hvg.info), 1000)
g.19 <- head(rownames(BM_M_Y@hvg.info), 1000) 
g.20 <- head(rownames(BM_M_O@hvg.info), 1000) 
g.21 <- head(rownames(BM_M_CR@hvg.info), 1000)
g.22 <- head(rownames(BAT_F_Y@hvg.info), 1000) 
g.23 <- head(rownames(BAT_F_O@hvg.info), 1000) 
g.24 <- head(rownames(BAT_F_CR@hvg.info), 1000) 
g.25 <- head(rownames(WAT_F_Y@hvg.info), 1000) 
g.26 <- head(rownames(WAT_F_O@hvg.info), 1000) 
g.27 <- head(rownames(WAT_F_CR@hvg.info), 1000)
g.28 <- head(rownames(Skin_F_Y@hvg.info), 1000) 
g.29 <- head(rownames(Skin_F_O@hvg.info), 1000) 
g.30 <- head(rownames(Skin_F_CR@hvg.info), 1000)
g.31 <- head(rownames(Kidney_F_Y@hvg.info), 1000)
g.32 <- head(rownames(Kidney_F_O@hvg.info), 1000)
g.33 <- head(rownames(Kidney_F_CR@hvg.info), 1000)
g.34 <- head(rownames(Aorta_F_Y@hvg.info), 1000) 
g.35 <- head(rownames(Aorta_F_O@hvg.info), 1000) 
g.36 <- head(rownames(Aorta_F_CR@hvg.info), 1000)
g.37 <- head(rownames(Liver_F_Y@hvg.info), 1000) 
g.38 <- head(rownames(Liver_F_O@hvg.info), 1000) 
g.39 <- head(rownames(Liver_F_CR@hvg.info), 1000)
g.40<- head(rownames(BM_F_Y@hvg.info), 1000) 
g.41 <- head(rownames(BM_F_O@hvg.info), 1000) 
g.42 <- head(rownames(BM_F_CR@hvg.info), 1000)

genes.use <- unique(c(g.1,g.2,g.3,g.4,g.5,g.6,g.7,g.8,g.9,g.10,g.11,g.12,g.13,g.14,g.15,g.16,g.17,g.18,g.19,g.20,g.21,g.22,g.23,g.24,g.25,g.26,g.27,g.28,g.29,g.30,g.31,g.32,g.33,g.34,g.35,g.36,g.37,g.38,g.39,g.40,g.41,g.42)) 

####################################

genes.use <- intersect(genes.use, rownames(BAT_M_Y@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(WAT_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(WAT_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(WAT_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Skin_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Skin_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Skin_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Kidney_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Kidney_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Kidney_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Aorta_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Aorta_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Aorta_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Liver_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Liver_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Liver_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(BM_M_Y@scale.data))
genes.use <- intersect(genes.use, rownames(BM_M_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BM_M_CR@scale.data))
genes.use <- intersect(genes.use, rownames(BAT_F_Y@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BAT_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(WAT_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(WAT_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(WAT_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Skin_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Skin_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Skin_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Kidney_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Kidney_F_O@scale.data))
genes.use <- intersect(genes.use, rownames(Kidney_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Aorta_F_Y@scale.data)) 
genes.use <- intersect(genes.use, rownames(Aorta_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Aorta_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(Liver_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(Liver_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(Liver_F_CR@scale.data))
genes.use <- intersect(genes.use, rownames(BM_F_Y@scale.data))
genes.use <- intersect(genes.use, rownames(BM_F_O@scale.data)) 
genes.use <- intersect(genes.use, rownames(BM_F_CR@scale.data))

################################################### Atlas

Atlas.list <- list(BAT_M_Y,BAT_M_O, BAT_M_CR, WAT_M_Y, WAT_M_O, WAT_M_CR, Skin_M_Y, Skin_M_O, Skin_M_CR, Kidney_M_Y, Kidney_M_O, Kidney_M_CR, Aorta_M_Y, Aorta_M_O, Aorta_M_CR, Liver_M_Y, Liver_M_O, Liver_M_CR, BM_M_O, BM_M_CR, BM_M_Y, BAT_F_Y,BAT_F_O, BAT_F_CR, WAT_F_Y, WAT_F_O, WAT_F_CR, Skin_F_Y, Skin_F_O, Skin_F_CR, Kidney_F_Y,Kidney_F_O, Kidney_F_CR, Aorta_F_Y, Aorta_F_O, Aorta_F_CR, Liver_F_Y, Liver_F_O, Liver_F_CR, BM_F_O, BM_F_CR, BM_F_Y)
Atlas_cca <- RunMultiCCA(object.list = Atlas.list, genes.use = genes.use, num.ccs = 30)
Atlas_cca <- AlignSubspace(Atlas_cca, reduction.type = "cca", grouping.var = "sample", dims.align = 1:30)
Atlas_cca <- RunPCA(object = Atlas_cca, pc.genes = genes.use, pcs.compute = 100, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

saveRDS(Atlas_cca, file = "/data1/mashuai/data/rat_10X/RDS/New/RatCellAtlas_cca.rds")

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/PCHeatmap.pdf")
PCElbowPlot(object = Atlas_cca, num.pc = 100)
PCHeatmap(Atlas_cca, pc.use = c(1:3, 61:63,75:80,95:100), cells.use = 500, do.balanced = TRUE)
dev.off()

Atlas_cca <- FindClusters(object = Atlas_cca, reduction.type = "pca", dims.use = 1:100, resolution = 3.5, save.SNN = TRUE, n.start = 100, nn.eps = 0)
Atlas_cca <- RunTSNE(object = Atlas_cca, reduction.use = "pca", dims.use = 1:100, min_dist = 0.75)

Atlas_cca@meta.data$tissue <- factor(Atlas_cca@meta.data$tissue,levels=c('BAT','WAT','Liver','Kidney','Aorta','Skin','BM'))

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-tissue.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne",do.return = TRUE, pt.size = 1,group.by="tissue", cols.use = c("#AC8226","#9A979B","#78B976","#309ACD","#D82B24","#E9B42C","#906AA3")) 
dev.off()
 
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-age.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne", do.return = TRUE, pt.size = 1,group.by="age", cols.use = c("seagreen2","orange2", "dodgerblue2")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-Y.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne", do.return = TRUE, pt.size = 1,group.by="age", cols.use = c("NA","NA", "dodgerblue2")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-O.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne", do.return = TRUE, pt.size = 1,group.by="age", cols.use = c("NA","orange2", "NA")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-CR.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne", do.return = TRUE, pt.size = 1,group.by="age", cols.use = c("seagreen2","NA", "NA")) 
dev.off()

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-gender.png",height=1500,width=1550)
DimPlot(object = Atlas_cca, reduction.use = "tsne", do.return = TRUE, pt.size = 1, group.by = "gender", cols.use = c('plum2','steelblue1')) 
dev.off()

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-cluster.png",height=1500,width=1700)
DimPlot(object = Atlas_cca, reduction.use = "tsne",do.return = TRUE, pt.size = 1, do.label =T, label.size=8) 
dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Genes_Vln_Plot-1.pdf",height=8,width=36)
VlnPlot(object=Atlas_cca,features.plot=c("nGene","nUMI","percent.mito"),point.size.use=0,nCol=1,group.by="SAMPLE")
dev.off()

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/FeaturePlot-Ptprc.png",height=1500,width=1550)
FeaturePlot(object = Atlas_cca, features.plot=c("Ptprc"),cols.use = c("grey90","red")) 
dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Marker_VlnPlot.pdf",height=88,width=35)
VlnPlot(object = Atlas_cca, features.plot = c("Pecam1","Kdr","Tagln","Lum","Dcn","Krt15","Krt5","Krt14","Krt1","Krt10","Krt79","Sox9","Gjb2","Mlana","Stmn1","Mki67","Cyp2e1","Krt18","Onecut1","Sox9","Hnf1b","Slc27a2","Lrp2","Slc5a2","Slc5a12","Slc27a2","Lrp2","Rida","Slc27a2","Lrp2","Aadat","Kap","Slc12a1","Umod","Slc12a3","Aqp1","Aqp2","Hsd11b2","Atp6v1g3","Atp6v0d2","Slc4a1","Atp6v1g3","Atp6v0d2","Slc26a4","Hmx2","Upk3bl1","Krt18","Krt8","Dcn","Gsn","Dcn","Pdgfra","Anpep","Hbb","Ptprc"), size.title.use=14,point.size.use = 0 ,nCol=1)
dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Marker_VlnPlot2.pdf",height=66,width=35)
VlnPlot(object = Atlas_cca, features.plot = c("Pf4","Meis1","Cd3d","Cd4","Cd8a","Klrb1a","Irf5","Irf8","Cd14","Fcgr3a","Cd68","Cd207","Cd163","Hbb","Tfrc","Cd74","Prtn3","Ctsg","Jchain","Mzb1","Ms4a2","Cd63","Cd19","Cd79b","Dntt","Rag1","Pglyrp1","Lcn2","S100a8","Mki67","Ezh2"), point.size.use = 0,nCol=1)
dev.off()

Atlas_cca@meta.data$age <- factor(Atlas_cca@meta.data$age,levels = c('Y','O','CR'))
pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/RidgePlot.pdf",height=5,width=10)
RidgePlot(Atlas_cca,features.plot = 'nGene',group.by="age",cols.use = c("dodgerblue2","orange2","seagreen2"))
RidgePlot(Atlas_cca,features.plot = 'nUMI',group.by="age",cols.use = c("dodgerblue2","orange2","seagreen2"))
dev.off()

######################################### Atlas_cca new

new.ident <- c("BC1","IC","IC","BC2","HFC","IC","PT","IC","Fib","IC","IC","ASC","IC","Fib","IC","IC","IC","LOH","ASC","IC","IC","Fib","SMC","IC","IC","IC","IC","IC","IC","EC","IC","EC","BC1","IC","ASC","IC","EC","IC","EC","IC","IC","IC","IC","ASC","IC","ASC","IC","Fib","IC","Fib","IC","ASC","Fib","Mit2","IC","IC","PT","Fib","SMC","IC","Fib","IC","Fib","IC","Fib","DLOH","LOH","IC","IC","CD2","Spi","ASC","Fib","IC","CD3","IC","Mit1","IC","SMC","EC","IC","IC","Fib","Fib","IC","IC","CC","IC","EC","CD1","IC","Hep","EC","ProE","EC","Ery","IC","IC","Fib","IC","IC","Epi","IC","IC","IC","IC","IC","Cho","IC")

new.ident <- c("BC1","Neu","T1","BC2","HFC","M1","PT","Neu","Fib","Neu(PT)","DC1","ASC","M2","Fib","M2","ProN","NK","LOH","ASC","Pla","M2","Fib","SMC","T2","Neu(PT)","Neu","T2","PB","M1","EC","Mon","EC","BC1","M2","ASC","DC1","EC","M2","EC","DC1","Neu","NK","KC","ASC","M1","ASC","KC","Fib","M1","Fib","DC1","ASC","Fib","Mit2","KC","M1","PT","Fib","SMC","LC","Fib","DC2","Fib","M1","Fib","DLOH","LOH","DC2","M1","CD2","Spi","ASC","Fib","B","CD3","M2","Mit1","Neu(PT)","SMC","EC","KC","NKT","Fib","Fib","Bas","NKT","CC","M1","EC","CD1","IB","Hep","EC","ProE","EC","Ery","M2","LPB","Fib","DC1","T1","Epi","B","DC1","NKT","M2","Fib","Cho","Epi")

Atlas_cca_new <- Atlas_cca

for (i in 0:108) {
    Atlas_cca_new <- RenameIdent(object = Atlas_cca_new, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

Atlas_cca_new@ident <- factor(Atlas_cca_new@ident,levels=c("Fib","M2","ASC","M1","Neu","EC","DC1","BC1","Neu(PT)","KC","T2","T1","PT","SMC","NK","BC2","LOH","HFC","ProN","Pla","DC2","PB","Mon","NKT","Mit2","LC","B","DLOH","CD2","Spi","CD3","Mit1","Bas","CC","CD1","IB","Hep","ProE","Ery","Epi","LPB","Cho"))

Atlas_cca_new <- StashIdent(Atlas_cca_new, save.name = "celltype")
Atlas_cca_new@meta.data$celltype <- factor(Atlas_cca_new@meta.data$celltype,levels=c("Fib","M2","ASC","M1","Neu","EC","DC1","BC1","Neu(PT)","KC","T2","T1","PT","SMC","NK","BC2","LOH","HFC","ProN","Pla","DC2","PB","Mon","NKT","Mit2","LC","B","DLOH","CD2","Spi","CD3","Mit1","Bas","CC","CD1","IB","Hep","ProE","Ery","Epi","LPB","Cho"))

Atlas_cca_new@meta.data$celltype.tissue <- paste0(Atlas_cca_new@meta.data$celltype, "_", Atlas_cca_new@meta.data$tissue)
Atlas_celltype_num <- table(Atlas_cca_new@meta.data$celltype)
Atlas_celltype_tissue_num <- table(Atlas_cca_new@meta.data$celltype.tissue)
write.csv(Atlas_celltype_num,'/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Atlas_celltype_num.csv')
write.csv(Atlas_celltype_tissue_num,'/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Atlas_celltype_tissue_num.csv')

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-celltype1.png",height=1500,width=1600)
DimPlot(object = Atlas_cca_new, reduction.use = "tsne",do.return = TRUE, pt.size = 1,cols.use =c("paleturquoise","turquoise","turquoise3","skyblue2","deepskyblue2","grey30","dodgerblue2","plum2","plum4","orchid2","orchid4","mediumpurple2","mediumslateblue","purple","purple2","purple4","brown","brown2","brown4","tomato2","coral2","coral4","chocolate3","darkorange4","darkorange2","darkorange","darkgoldenrod2","goldenrod1","goldenrod4","gold3","khaki2","khaki4","yellow","yellow3","yellowgreen","olivedrab2","olivedrab4","palegreen","palegreen4","seagreen2","seagreen4","grey70")) 
dev.off()

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-celltype2.png",height=1500,width=1600)
DimPlot(object = Atlas_cca_new, reduction.use = "tsne",do.return = TRUE, pt.size = 1,cols.use =c("paleturquoise","turquoise","turquoise3","skyblue2","deepskyblue2","grey30","dodgerblue2","plum2","plum4","orchid2","orchid4","mediumpurple2","mediumslateblue","purple","purple2","purple4","brown","brown2","brown4","tomato2","coral2","coral4","chocolate3","darkorange4","darkorange2","darkorange","darkgoldenrod2","goldenrod1","goldenrod4","gold3","khaki2","khaki4","yellow","yellow3","yellowgreen","olivedrab2","olivedrab4","palegreen","palegreen4","seagreen2","seagreen4","grey70"),do.label = T,label.size=10) 
dev.off()

Atlas_IC <- SubsetData(object=Atlas_cca_new,ident.use=c('IC'))
Atlas_IC <- SubsetData(Atlas_IC,subset.name='tissue',accept.value=c('WAT','BAT','Liver','Kidney','Aorta','Skin'))


png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_age.png",height=1500,width=1550)
DimPlot(object = Atlas_IC, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("seagreen2","orange2", "dodgerblue2")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_Y.png",height=1500,width=1550)
DimPlot(object = Atlas_IC, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("NA","NA", "dodgerblue2")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_O.png",height=1500,width=1550)
DimPlot(object = Atlas_IC, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("NA","orange2", "NA")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_CR.png",height=1500,width=1550)
DimPlot(object = Atlas_IC, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("seagreen2","NA", "NA")) 
dev.off()

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/FeaturePlot-Ptprc.png",height=1500,width=1550)
FeaturePlot(object = Atlas_cca, features.plot=c("Ptprc"),cols.use = c("grey90","red"),pt.size = 2) 
dev.off()

######################################### Atlas gender

Atlas_IC_F <- SubsetData(Atlas_IC,subset.name='gender',accept.value='F')
Atlas_IC_M <- SubsetData(Atlas_IC,subset.name='gender',accept.value='M')

png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_F_age.png",height=1500,width=1550)
DimPlot(object = Atlas_IC_F, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("seagreen2","orange2", "dodgerblue2")) 
dev.off()
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/tSNE_Plot-IC_M_age.png",height=1500,width=1550)
DimPlot(object = Atlas_IC_M, reduction.use = "tsne", do.return = TRUE, pt.size = 2,group.by="age", cols.use = c("seagreen2","orange2", "dodgerblue2")) 
dev.off()

Atlas_cca_new.cell.80000.list <-  sample(length(Atlas_cca_new@cell.names),size=80000)
Atlas_cca_new.cell.80000 <- SubsetData(object= Atlas_cca_new,cells.use= Atlas_cca_new@cell.names[Atlas_cca_new.cell.80000.list])
table(Atlas_cca_new.cell.80000@meta.data$age)
Atlas.markers <- FindAllMarkers(object = Atlas_cca_new.cell.80000,logfc.threshold=0.5,only.pos=T)
write.csv(Atlas.markers,'/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Atlas_markers_CT.csv')
top15 <- Atlas.markers %>% group_by(cluster) %>% top_n(15, avg_logFC)
png("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Atlas_heatmap.png",height=600,width=2400)
DoHeatmap(object= Atlas_cca_new.cell.80000,genes.use = top15$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Dotplot_markers.pdf",height=12,width=40)
sdp <- DotPlot(Atlas_cca_new, genes.plot = rev(c("Dcn","Lum","Cd68","Cd14","Cd163","Anpep","Pdgfra","Pglyrp1","Lcn2","Slc6a13","Pecam1","Kdr","Cd74","Irf5","Krt5","Krt15","S100a8","S100a9","Clec4f","Cd3d","Cd8a","Slc27a2","Lrp2","Tagln","Myh11","Klrb1a","Krt14","Slc12a1","Umod","Krt79","Sox9","Mki67","Mzb1","Jchain","Irf8","Cd19","Cd79b","Ezh2","Prtn3","Ctsg","Stmn1","Cd207","Aqp1","Atp6v1g3","Slc4a1","Krt1","Krt10","Slc26a4","Ms4a2","Gjb2","Aqp2","Hsd11b2","Cyp2e1","Hbb","Tfrc","Krt8","Krt18","Dntt","Rag1","Onecut1")), cols.use = c("royalblue4","orangered2"), x.lab.rot = T, plot.legend = T, dot.scale = 4, do.return = T)
dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Dotplot_age_markers.pdf",height=36,width=40)
sdp <- SplitDotPlotGG(Atlas_cca_new, genes.plot = rev(c("Dcn","Lum","Cd68","Cd14","Cd163","Anpep","Pdgfra","Pglyrp1","Lcn2","Slc6a13","Pecam1","Kdr","Cd74","Irf5","Krt5","Krt15","S100a8","S100a9","Clec4f","Cd3d","Cd8a","Slc27a2","Lrp2","Tagln","Myh11","Klrb1a","Krt14","Slc12a1","Umod","Krt79","Sox9","Mki67","Mzb1","Jchain","Irf8","Cd19","Cd79b","Ezh2","Prtn3","Ctsg","Stmn1","Cd207","Aqp1","Atp6v1g3","Slc4a1","Krt1","Krt10","Slc26a4","Ms4a2","Gjb2","Aqp2","Hsd11b2","Cyp2e1","Hbb","Tfrc","Krt8","Krt18","Dntt","Rag1","Onecut1")), x.lab.rot = T, plot.legend = T, dot.scale = 4, do.return = T,cols.use = c("seagreen2","orange2", "dodgerblue2"),grouping.var = "age")
dev.off()


