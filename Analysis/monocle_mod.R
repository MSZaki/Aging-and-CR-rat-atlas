##############################################  monocle
#################
#################

###################
library(monocle)
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)

Neu_merge <- readRDS('/data1/mashuai/data/rat_10X/RDS/New/RatCellAtlas_Neu.rds')
Neu_merge <- StashIdent(Neu_merge, save.name = "celltype")
Neu_merge

Neu_merge<- SubsetData(object=Neu_merge, subset.name='tissue',accept.value=c('BAT','WAT','Liver','Kidney','BM'))
Neu_merge
table(Neu_merge@meta.data$celltype)
table(Neu_merge@meta.data$tissue)

Neu_merge <- StashIdent(Neu_merge, save.name = "celltype")

Neu_m <- importCDS(Neu_merge,import_all=T)
Neu_m <- estimateSizeFactors(Neu_m)
Neu_m <- estimateDispersions(Neu_m)
pData(Neu_m)$Total_mRNAs <- Matrix::colSums(exprs(Neu_m))
Neu_m <- Neu_m[,pData(Neu_m)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(Neu_m)$Total_mRNAs)) +
            2*sd(log10(pData(Neu_m)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(Neu_m)$Total_mRNAs)) -
            2*sd(log10(pData(Neu_m)$Total_mRNAs)))
Neu_m<- Neu_m[,pData(Neu_m)$Total_mRNAs > lower_bound &
      pData(Neu_m)$Total_mRNAs < upper_bound]
Neu_m<- detectGenes(Neu_m, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(Neu_m), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(Neu_m[expressed_genes,],fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
Neu_m <- setOrderingFilter(Neu_m, ordering_genes)
Neu_m <- reduceDimension(Neu_m, max_components = 2,
method = 'DDRTree')

GM_state <- function(Neu_m){
  if (length(unique(pData(Neu_m)$State)) > 1){
    T0_counts <- table(pData(Neu_m)$State, pData(Neu_m)$celltype)[,'ProN']
    return(as.numeric(names(T0_counts)[which
          (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

Neu_m <- orderCells(Neu_m)

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Neu_trajectory.pdf",width=10,height=6)
plot_cell_trajectory(Neu_m, color_by = "Pseudotime") 
plot_cell_trajectory(Neu_m, color_by = "celltype")+scale_color_manual(values= c('pink2','slateblue1','seagreen2' ))
plot_cell_trajectory(Neu_m, color_by = "tissue")+scale_color_manual(values= c('darkgoldenrod4','grey80','seagreen3','skyblue2','plum3'))
plot_cell_trajectory(Neu_m, color_by = "age")+scale_color_manual(values= c('seagreen2','orange2','dodgerblue2'))

dev.off()

pdf("/data1/mashuai/data/rat_10X/DoubletFinder/Atlas/Neu_pseudotime_gene.pdf",width=12,height=12)
Neu_expressed_genes <-  row.names(subset(fData(Neu_m), num_cells_expressed >= 10))
Neu_filtered <- Neu_m[Neu_expressed_genes,]
my_genes <- row.names(subset(fData(Neu_filtered), gene_short_name %in% c("Mki67","Lcn2","Cxcl2", "Ccl2")))
cds_subset <- Neu_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "celltype")+scale_color_manual(values= c('pink2','slateblue1','seagreen2'))

dev.off()

