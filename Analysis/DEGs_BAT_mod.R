#####################################################
######### Function

###########################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)

###########################

BAT_cca <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/BAT_Y_O_CR.rds")
BAT_IC <- readRDS("/data1/mashuai/data/rat_10X/RDS/New/BAT_IC_Y_O_CR.rds")

BAT_cca_CT <- BAT_cca
new.ident <- c("IC","IC","IC","IC","IC","IC","ASC","IC","Fib","IC","ASC","Fib","Fib","EC","IC","IC","IC","EC","Fib","ASC","IC","IC","EC","IC","IC")
for (i in 0:24) {
    BAT_cca_CT <- RenameIdent(object = BAT_cca_CT, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

BAT_IC_CT <- BAT_IC
new.ident <- c("T2","Neu","NK","M1","Neu","DC1","M1","M1","M2","Neu","Neu","M2","M2","T2","B","NKT","DC2","DC1","Neu","NKT","Pla")
for (i in 0:20) {
    BAT_IC_CT <- RenameIdent(object = BAT_IC_CT, old.ident.name = i, new.ident.name = new.ident[i + 1])
}

BAT_cca_CT.F <- SubsetData(BAT_cca_CT,subset.name='gender',accept.value='F')
BAT_cca_CT.M <- SubsetData(BAT_cca_CT,subset.name='gender',accept.value='M')
BAT_IC_CT.F <- SubsetData(BAT_IC_CT,subset.name='gender',accept.value='F')
BAT_IC_CT.M <- SubsetData(BAT_IC_CT,subset.name='gender',accept.value='M')

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
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_cca_new@ident)$BAT_cca_new.ident)
{
if ((length(which(BAT_cca_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_cca_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_cca_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_cca_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsY <- FindMarkers(BAT_cca_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsY.up <- CRvsY[which(CRvsY$avg_logFC>=0.5&CRvsY$p_val_adj<=0.05),]
CRvsY.down <- CRvsY[which(CRvsY$avg_logFC<=-0.5&CRvsY$p_val_adj<=0.05),]
nonres.up <- CRvsY.up[which(rownames(CRvsY.up) %in% rownames(OvsY.up)),]
nonres.down <- CRvsY.down[which(rownames(CRvsY.down) %in% rownames(OvsY.down)),]
write.table(CRvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_CRvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_CRvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(nonres.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_nonres.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(nonres.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/DEGs/",type,"_nonres.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")

}
}
}
############################# BAT_IC_new
BAT_IC_new <- BAT_IC_CT
BAT_IC_new@meta.data$celltype.age <- paste0(BAT_IC_new@ident, "_", BAT_IC_new@meta.data$age)
BAT_IC_new <- StashIdent(BAT_IC_new, save.name = "celltype")
BAT_IC_new <- SetAllIdent(BAT_IC_new, id = "celltype.age")
cell.type <- BAT_IC_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_IC_new@ident)$BAT_IC_new.ident)
{
if ((length(which(BAT_IC_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_IC_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_IC_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_IC_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsY <- FindMarkers(BAT_IC_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsY.up <- CRvsY[which(CRvsY$avg_logFC>=0.5&CRvsY$p_val_adj<=0.05),]
CRvsY.down <- CRvsY[which(CRvsY$avg_logFC<=-0.5&CRvsY$p_val_adj<=0.05),]
nonres.up <- CRvsY.up[which(rownames(CRvsY.up) %in% rownames(OvsY.up)),]
nonres.down <- CRvsY.down[which(rownames(CRvsY.down) %in% rownames(OvsY.down)),]
write.table(CRvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_CRvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_CRvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(nonres.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_nonres.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(nonres.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/BAT_IC/DEGs/",type,"_nonres.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")

}
}
}

########### Gender DEGs

BAT_cca_F_new <- BAT_cca_CT.F
BAT_cca_F_new@meta.data$celltype.age <- paste0(BAT_cca_F_new@ident, "_", BAT_cca_F_new@meta.data$age)
BAT_cca_F_new <- StashIdent(BAT_cca_F_new, save.name = "celltype")
BAT_cca_F_new <- SetAllIdent(BAT_cca_F_new, id = "celltype.age")
cell.num <- table(BAT_cca_F_new@meta.data$celltype.age)
write.csv(cell.num,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/BAT_F_CT_num.csv')
cell.type <- BAT_cca_F_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_cca_F_new@ident)$BAT_cca_F_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_cca_F_new@ident)$BAT_cca_F_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_cca_F_new@ident)$BAT_cca_F_new.ident)
{
if ((length(which(BAT_cca_F_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_cca_F_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_cca_F_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_cca_F_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_cca_F_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}
}

######################### BAT_IC_F_new
BAT_IC_F_new <- BAT_IC_CT.F
BAT_IC_F_new@meta.data$celltype.age <- paste0(BAT_IC_F_new@ident, "_", BAT_IC_F_new@meta.data$age)
BAT_IC_F_new <- StashIdent(BAT_IC_F_new, save.name = "celltype")
BAT_IC_F_new <- SetAllIdent(BAT_IC_F_new, id = "celltype.age")
cell.num <- table(BAT_IC_F_new@meta.data$celltype.age)
write.csv(cell.num,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/BAT_IC_F_CT_num.csv')
cell.type <- BAT_IC_F_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_IC_F_new@ident)$BAT_IC_F_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_IC_F_new@ident)$BAT_IC_F_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_IC_F_new@ident)$BAT_IC_F_new.ident)
{
if ((length(which(BAT_IC_F_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_IC_F_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_IC_F_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_IC_F_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_IC_F_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Female/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}
}

###################### BAT_cca_M_new

BAT_cca_M_new <- BAT_cca_CT.M
BAT_cca_M_new@meta.data$celltype.age <- paste0(BAT_cca_M_new@ident, "_", BAT_cca_M_new@meta.data$age)
BAT_cca_M_new <- StashIdent(BAT_cca_M_new, save.name = "celltype")
BAT_cca_M_new <- SetAllIdent(BAT_cca_M_new, id = "celltype.age")
cell.num <- table(BAT_cca_M_new@meta.data$celltype.age)
write.csv(cell.num,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/BAT_M_CT_num.csv')
cell.type <- BAT_cca_M_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_cca_M_new@ident)$BAT_cca_M_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_cca_M_new@ident)$BAT_cca_M_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_cca_M_new@ident)$BAT_cca_M_new.ident)
{
if ((length(which(BAT_cca_M_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_cca_M_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_cca_M_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_cca_M_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_cca_M_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}
}
######################## BAT_IC_M_new
BAT_IC_M_new <- BAT_IC_CT.M
BAT_IC_M_new@meta.data$celltype.age <- paste0(BAT_IC_M_new@ident, "_", BAT_IC_M_new@meta.data$age)
BAT_IC_M_new <- StashIdent(BAT_IC_M_new, save.name = "celltype")
BAT_IC_M_new <- SetAllIdent(BAT_IC_M_new, id = "celltype.age")
cell.num <- table(BAT_IC_M_new@meta.data$celltype.age)
write.csv(cell.num,'/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/BAT_IC_M_CT_num.csv')
cell.type <- BAT_IC_M_new@meta.data$celltype
cell.type <- cell.type[!duplicated(cell.type)]
n<-length(cell.type)
for (i in 1:n){
type <- cell.type[i]
#cell number
identity_Y <- paste(type,"_Y",sep="")
identity_O <- paste(type,"_O",sep="")
identity_CR <- paste(type,"_CR",sep="")
if (paste(type,"_O",sep="")%in%data.frame(BAT_IC_M_new@ident)$BAT_IC_M_new.ident & paste(type,"_Y",sep="")%in%data.frame(BAT_IC_M_new@ident)$BAT_IC_M_new.ident & paste(type,"_CR",sep="")%in%data.frame(BAT_IC_M_new@ident)$BAT_IC_M_new.ident)
{
if ((length(which(BAT_IC_M_new@ident %in% identity_Y)) >= 3 )&(length(which(BAT_IC_M_new@ident %in% identity_O)) >= 3 )&(length(which(BAT_IC_M_new@ident %in% identity_CR)) >= 3 ))
{
#dif gene
OvsY <- FindMarkers(BAT_IC_M_new,ident.1=paste(type,"_O",sep=""),ident.2=paste(type,"_Y",sep=""))
CRvsO <- FindMarkers(BAT_IC_M_new,ident.1=paste(type,"_CR",sep=""),ident.2=paste(type,"_O",sep=""))
CRvsO <- CRvsO[which(abs(CRvsO$avg_logFC)>=0.5&CRvsO$p_val_adj<=0.05),]
OvsY.up <- OvsY[which(OvsY$avg_logFC>=0.5&OvsY$p_val_adj<=0.05),]
OvsY.down <- OvsY[which(OvsY$avg_logFC<=-0.5&OvsY$p_val_adj<=0.05),]
CRvsO.up <- CRvsO[which(CRvsO$avg_logFC>=0.5&CRvsO$p_val_adj<=0.05),]
CRvsO.down <- CRvsO[which(CRvsO$avg_logFC<=-0.5&CRvsO$p_val_adj<=0.05),]
res.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.down)),]
res.down <- CRvsO.down[which(rownames(CRvsO.down)%in%rownames(OvsY.up)),]
ng.up <- CRvsO.up[which(rownames(CRvsO.up) %in% rownames(OvsY.up)),]
ng.down <- CRvsO.down[which(rownames(CRvsO.down) %in% rownames(OvsY.down)),]
write.table(OvsY.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_OvsY.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(OvsY.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_OvsY.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_CRvsO.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(CRvsO.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_CRvsO.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_res.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(res.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_res.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.up,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_ng.up.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
write.table(ng.down,file=paste("/data1/mashuai/data/rat_10X/DoubletFinder/BAT/Male/DEGs/",type,"_ng.down.list",sep=""),quote = FALSE,row.names = T, col.names = T,sep="\t")
}
}
}
