################################### ageR O_vs_Y
data.path="G:\\result\\result\\09_sc_new\\pathway\\"
gene.df=read.csv(paste(data.path,"list\\ageR_pathway.csv",sep = ""))
gene.list=as.character(gene.df$V1) 
all.stat=data.frame(Gene_name=gene.list)
gene.df=all.stat

###############
library(data.table)
library(pheatmap)
library(xlsx)

tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
name=tissue.list[1]
for (name in tissue.list) {

all.DEG=read.csv(paste("G:\\result\\result\\09_sc_new\\pathway\\",name,".DEG.csv",sep=""))
all.DEG=all.DEG[grep("OvsY",all.DEG$class),]
tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
tmp.DEG=tmp.DEG[,-1]  

type.list=unique(tmp.DEG$Cell_type)
type=type.list[2]
for (type in type.list) {
 
  type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
  stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
  colnames(stat.df)[2]=type
  if(type=="BAT:ASC"){
    all.stat=stat.df
  }else{
    all.stat=merge(all.stat,stat.df,by="Gene_name",all.x=T,sort=F)
  }
  
}

}

plot.df=all.stat

plot.df=plot.df[,-1]
plot.df[is.na(plot.df)]=0

plot.df=plot.df[rowSums(abs(plot.df) )>0,]
plot.df=plot.df[,colSums(abs(plot.df) )>0]

plot.df$pos=rowSums(plot.df>0)
plot.df$neg=rowSums(plot.df<0)
plot.df$diff=plot.df$pos - plot.df$neg
plot.df=plot.df[order(plot.df$diff,decreasing = T),]
plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]
plot.df2=t(plot.df)
write.csv(plot.df2,"G:\\result\\result\\09_sc_new\\pathway\\ageR_OvsY.csv",row.names = T)

type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 

pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(table(type.anno) )[-7],
         #annotation_row =row_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T,
         filename = "G:\\result\\result\\09_sc_new\\pathway\\ageR_OvsY.pdf"
)

########################################## other pathway  O_vs_Y
way.list=c("AMPK","ILSR","mTOR","PPAR")

for (way in way.list) {
  
gene.df=data.frame(fread(paste("G:\\result\\result\\09_sc_new\\pathway\\list\\",way,"_pathway.list",sep=""),header=F)) 
colnames(gene.df)="Gene_name"
gene.list=gene.df$Gene_name
all.stat=gene.df
tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
name=tissue.list[1]

for (name in tissue.list) {
  
  all.DEG=read.csv(paste("G:\\result\\result\\09_sc_new\\pathway\\",name,".DEG.csv",sep=""))
  all.DEG=all.DEG[grep("OvsY",all.DEG$class),]
  tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
  tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
  tmp.DEG=tmp.DEG[,-1]  
  
  type.list=paste(name,unique(all.DEG$Cell_type),sep = ":") 
  type=type.list[9]
  
  for (type in type.list) {
    
    type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
    stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
    colnames(stat.df)[2]=type
    if(type=="BAT:ASC"){
      all.stat=stat.df
    }else{
      all.stat=merge(all.stat,stat.df,by="Gene_name",sort=F)
    }
    
  }
  
}

all.stat$Gene_name=paste(way,all.stat$Gene_name,sep = ':')
rownames(all.stat)=all.stat$Gene_name

plot.df=all.stat

plot.df=plot.df[,-1]
plot.df=plot.df[,-ncol(plot.df)]
plot.df[is.na(plot.df)]=0

plot.df$pos=rowSums(plot.df>0)
plot.df$neg=rowSums(plot.df<0)
plot.df$diff=plot.df$pos - plot.df$neg
plot.df=plot.df[order(plot.df$diff,decreasing = T),]
plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]

if(way==way.list[1]){
  plot.all=plot.df
}else{
  plot.all=rbind(plot.all,plot.df)
}

}

plot.all=plot.all[rowSums(abs(plot.all) )>0,]
plot.all=plot.all[,colSums(abs(plot.all) )>0]

type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 
plot.df2=t(plot.all)
dim(plot.df2)
row.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 
row_anno=data.frame(row.names = rownames(plot.df2),anno=row.anno)
col.anno=factor(unlist(strsplit(colnames(plot.df2),split = ':') )[seq(1,ncol(plot.df2)*2,2)],levels = unique(unlist(strsplit(colnames(plot.df2),split = ':') )[seq(1,ncol(plot.df2)*2,2)])) 
col_anno=data.frame(row.names = colnames(plot.df2),anno=col.anno)
row_gap=data.frame(table(row.anno) ) 
col_gap=data.frame(table(col.anno) ) 
pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(row_gap$Freq )[-7],
         gaps_col = cumsum(col_gap$Freq )[-4],
         annotation_col =col_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T
         #filename = "G:\\result\\result\\09_sc_new\\pathway\\ageR_OvsY.pdf"
)

dev.off()

########################################## rescue AgeR
###################################

data.path="G:\\result\\result\\09_sc_new\\pathway\\"
gene.df=read.csv(paste(data.path,"list\\ageR_pathway.csv",sep = ""))
gene.list=as.character(gene.df$V1) 
head(gene.list)
all.stat=data.frame(Gene_name=gene.list)

gene.df=all.stat

###############
library(data.table)
library(pheatmap)
library(xlsx)

tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
name=tissue.list[1]
for (name in tissue.list) {
  
  all.DEG=read.xlsx(paste("G:\\result\\result\\09_sc_new\\pathway\\all.DEG2.xls",sep=""),sheetName = name)
  all.DEG=all.DEG[grep("res",all.DEG$class),]
  tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
  tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
  tmp.DEG=tmp.DEG[,-1]  
  
  type.list=unique(tmp.DEG$Cell_type)
  type=type.list[2]
  for (type in type.list) {
    
    type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
    stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
    colnames(stat.df)[2]=type
    if(type=="BAT:ASC"){
      all.stat=stat.df
    }else{
      all.stat=merge(all.stat,stat.df,by="Gene_name",all.x=T,sort=F)
    }
    
  }
  
}

rownames(all.stat)=all.stat$Gene_name
plot.df=all.stat

plot.df=plot.df[,-1]
plot.df[is.na(plot.df)]=0

plot.df=plot.df[rowSums(abs(plot.df) )>0,]
plot.df=plot.df[,colSums(abs(plot.df) )>0]

plot.df$pos=rowSums(plot.df>0)
plot.df$neg=rowSums(plot.df<0)
plot.df$diff=plot.df$pos - plot.df$neg
plot.df=plot.df[order(plot.df$diff,decreasing = T),]
plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]
plot.df2=t(plot.df)
write.csv(plot.df2,"G:\\result\\result\\09_sc_new\\pathway\\ageR_RSE.csv",row.names = T)
type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 

pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(table(type.anno) )[-7],
         #annotation_row =row_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T,
         filename = "G:\\result\\result\\09_sc_new\\pathway\\ageR_CRvsO.pdf"
)

########################################## other pathway Rescue

way.list=c("AMPK","ILSR","mTOR","PPAR")
way=way.list[1]


for (way in way.list) {
  
  print(way)
  gene.df=data.frame(fread(paste("G:\\result\\result\\09_sc_new\\pathway\\list\\",way,"_pathway.list",sep=""),header=F)) 
  colnames(gene.df)="Gene_name"
  gene.list=gene.df$Gene_name
  all.stat=gene.df
  tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
  name=tissue.list[1]
  
  for (name in tissue.list) {
    
    all.DEG=read.xlsx(paste("G:\\result\\result\\09_sc_new\\pathway\\all.DEG2.xls",sep=""),sheetName = name)
    all.DEG=all.DEG[grep("res",all.DEG$class),]
    tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
    tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
    tmp.DEG=tmp.DEG[,-1]  
    
    type.list=paste(name,unique(all.DEG$Cell_type),sep = ":") 
    type=type.list[9]
    
    for (type in type.list) {
      
      type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
      stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
      colnames(stat.df)[2]=type
      if(type=="BAT:ASC"){
        all.stat=stat.df
      }else{
        all.stat=merge(all.stat,stat.df,by="Gene_name",sort=F)
      }
      
    }
    
  }
  
  all.stat$Gene_name=paste(way,all.stat$Gene_name,sep = ':')
  #all.stat$Gene_name=paste(all.stat$anno,all.stat$Gene_name,sep = ':')
  rownames(all.stat)=all.stat$Gene_name
  
  plot.df=all.stat
  
  plot.df=plot.df[,-1]
  plot.df=plot.df[,-ncol(plot.df)]
  plot.df[is.na(plot.df)]=0
  
  #plot.df=plot.df[rowSums(abs(plot.df) )>0,]
  #plot.df=plot.df[,colSums(abs(plot.df) )>0]
  
  
  plot.df$pos=rowSums(plot.df>0)
  plot.df$neg=rowSums(plot.df<0)
  plot.df$diff=plot.df$pos - plot.df$neg
  plot.df=plot.df[order(plot.df$diff,decreasing = T),]
  plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]
  
  if(way==way.list[1]){
    plot.all=plot.df
  }else{
    plot.all=rbind(plot.all,plot.df)
  }
  
}


plot.all=plot.all[rowSums(abs(plot.all) )>0,]
plot.all=plot.all[,colSums(abs(plot.all) )>0]

type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 
plot.df2=t(plot.all)
row.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 
row_anno=data.frame(row.names = rownames(plot.df2),anno=row.anno)
col.anno=factor(unlist(strsplit(colnames(plot.df2),split = ':') )[seq(1,ncol(plot.df2)*2,2)],levels = unique(unlist(strsplit(colnames(plot.df2),split = ':') )[seq(1,ncol(plot.df2)*2,2)])) 
col_anno=data.frame(row.names = colnames(plot.df2),anno=col.anno)
row_gap=data.frame(table(row.anno) ) 
col_gap=data.frame(table(col.anno) ) 
pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(row_gap$Freq )[-7],
         gaps_col = cumsum(col_gap$Freq )[-4],
         annotation_col =col_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T
         #filename = "G:\\result\\result\\09_sc_new\\pathway\\ageR_OvsY.pdf"
)

dev.off()

######################################### SASP related gene list  O_vs_Y
###################################
data.path="G:\\result\\result\\09_sc_new\\pathway\\"
gene.df=read.csv(paste(data.path,"list\\SASP_pathway.csv",sep = ""),header = F)
gene.list=as.character(gene.df$V1) 
head(gene.list)
all.stat=data.frame(Gene_name=gene.list)

gene.df=all.stat

###############
library(data.table)
library(pheatmap)
library(xlsx)

tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
name=tissue.list[1]
for (name in tissue.list) {
  
  all.DEG=read.csv(paste("G:\\result\\result\\09_sc_new\\pathway\\",name,".DEG.csv",sep=""))
  all.DEG=all.DEG[grep("OvsY",all.DEG$class),]
  tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
  tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
  tmp.DEG=tmp.DEG[,-1]  
  
  type.list=unique(tmp.DEG$Cell_type)
  type=type.list[2]
  for (type in type.list) {
    
    type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
    stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
    colnames(stat.df)[2]=type
    if(type=="BAT:ASC"){
      all.stat=stat.df
    }else{
      all.stat=merge(all.stat,stat.df,by="Gene_name",all.x=T,sort=F)
    }
    
  }
  
}

rownames(all.stat)=all.stat$Gene_name
plot.df=all.stat

plot.df=plot.df[,-1]
plot.df[is.na(plot.df)]=0

plot.df=plot.df[rowSums(abs(plot.df) )>0,]
plot.df=plot.df[,colSums(abs(plot.df) )>0]

plot.df
plot.df$pos=rowSums(plot.df>0)
plot.df$neg=rowSums(plot.df<0)
plot.df$diff=plot.df$pos - plot.df$neg
plot.df=plot.df[order(plot.df$diff,decreasing = T),]
plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]
plot.df2=t(plot.df)
rownames(plot.df2)
write.csv(plot.df2,"G:\\result\\result\\09_sc_new\\pathway\\SASP_OvsY.csv",row.names = T)
type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 

pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(table(type.anno) )[-7],
         #annotation_row =row_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T,
         filename = "G:\\result\\result\\09_sc_new\\pathway\\SASP_OvsY.pdf"
)

######################################### SASP related gene list  rescue
###################################
data.path="G:\\result\\result\\09_sc_new\\pathway\\"
gene.df=read.csv(paste(data.path,"list\\SASP_pathway.csv",sep = ""),header = F)
gene.list=as.character(gene.df$V1) 
head(gene.list)
all.stat=data.frame(Gene_name=gene.list)

gene.df=all.stat

###############
library(data.table)
library(pheatmap)
library(xlsx)

tissue.list=c("BAT",'WAT',"Aorta","Kidney","Liver","Skin","BM")
name=tissue.list[1]
for (name in tissue.list) {
  
  all.DEG=read.xlsx(paste("G:\\result\\result\\09_sc_new\\pathway\\all.DEG2.xls",sep=""),sheetName = name)
  all.DEG=all.DEG[grep("res",all.DEG$class),]
  tmp.DEG=all.DEG[all.DEG$Gene_name %in% gene.list,c(1,2,4,5)]
  tmp.DEG$Cell_type=paste(tmp.DEG$Tissue,tmp.DEG$Cell_type,sep = ":")
  tmp.DEG=tmp.DEG[,-1]  
  
  type.list=unique(tmp.DEG$Cell_type)
  type=type.list[2]
  for (type in type.list) {
    
    type.df=tmp.DEG[tmp.DEG$Cell_type==type,2:3]
    stat.df=merge(gene.df,type.df,by="Gene_name",all.x=T,sort=F)
    colnames(stat.df)[2]=type
    if(type=="BAT:ASC"){
      all.stat=stat.df
    }else{
      all.stat=merge(all.stat,stat.df,by="Gene_name",all.x=T,sort=F)
    }
    
  }
  
}

rownames(all.stat)=all.stat$Gene_name
plot.df=all.stat

plot.df=plot.df[,-1]
plot.df[is.na(plot.df)]=0

plot.df=plot.df[rowSums(abs(plot.df) )>0,]
plot.df=plot.df[,colSums(abs(plot.df) )>0]

plot.df
plot.df$pos=rowSums(plot.df>0)
plot.df$neg=rowSums(plot.df<0)
plot.df$diff=plot.df$pos - plot.df$neg
plot.df=plot.df[order(plot.df$diff,decreasing = T),]
plot.df=plot.df[,-c((ncol(plot.df)-2):ncol(plot.df))]
plot.df2=t(plot.df)
rownames(plot.df2)
write.csv(plot.df2,"G:\\result\\result\\09_sc_new\\pathway\\SASP_RES.csv",row.names = T)
type.anno=factor(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)],levels = unique(unlist(strsplit(rownames(plot.df2),split = ':') )[seq(1,nrow(plot.df2)*2,2)])) 

pheatmap(plot.df2  ,
         main =  '',
         cluster_rows = F,
         cluster_cols = F,
         scale = "none",
         color = c(colorRampPalette(c('#4575b4', "grey95"))(500),colorRampPalette(c("grey95",'#d73027'))(500)) ,
         breaks=seq(-1,1,length.out = 1000),
         border_color = NA, #"grey60",
         fontsize = 10,
         fontsize_row = 5, #fontsize,
         fontsize_col = 5, #fontsize,
         fontsize_number = 0.8* fontsize,
         kmeans_k = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         show_rownames = T,
         show_colnames = T,
         legend = TRUE,
         gaps_row = cumsum(table(type.anno) )[-7],
         #annotation_row =row_anno,
         #annotation_colors =anno_color  ,
         annotation_legend = T,
         # angle_col = "45",
         drop_levels = TRUE,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = 5,
         cellheight = 5,
         revC=T,
         filename = "G:\\result\\result\\09_sc_new\\pathway\\SASP_RES.pdf"
)
