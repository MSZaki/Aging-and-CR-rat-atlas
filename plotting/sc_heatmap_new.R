############################ 

################# generate all_DEG.csv
library(xlsx)
data.path='G:\\result\\result\\09_sc_new\\heatmap\\'
name.list=list.files(data.path)

for (sample.name in name.list) {
  
  tmp.path=paste(data.path,sample.name,"\\DEGs_all\\",sep = '')
 # tmp.path=paste(data.path,sample.name,"\\",sep = '')
  type.list=list.files(tmp.path)[grep('CRvsO.down.list',list.files(tmp.path))] 
  type.list=substr(type.list,1,nchar(type.list)-16)
  all.deg=data.frame(Tissue=character(),Cell_type=character(),Gene_name=character(),avg_logFC=numeric(),p_val=numeric(),p_val_adj=numeric())
  for (type in type.list) {
    
    CO.down=data.frame(fread(paste(tmp.path,'\\',type,'_CRvsO.down.list',sep=''))) 
    CO.up=data.frame(fread(paste(tmp.path,'\\',type,'_CRvsO.up.list',sep=''))) 
    
    OY.down=data.frame(fread(paste(tmp.path,'\\',type,'_OvsY.down.list',sep=''))) 
    OY.up=data.frame(fread(paste(tmp.path,'\\',type,'_OvsY.up.list',sep=''))) 
    
    if(nrow(CO.down)+nrow(CO.up)+nrow(OY.down)+nrow(OY.up)==0){
      type.df=data.frame(Tissue=character(),Cell_type=character(),Gene_name=character(),avg_logFC=numeric(),p_val=numeric(),p_val_adj=numeric())
    }else{
      type.df=data.frame(Tissue=sample.name,Cell_type=type,
                         class=rep(c('CRvsO.down','CRvsO.up','OvsY.down','OvsY.up'),c(nrow(CO.down),nrow(CO.up),nrow(OY.down),nrow(OY.up))),
                         Gene_name=c(as.character(CO.down[,1]),as.character(CO.up[,1]),as.character(OY.down[,1]),as.character(OY.up[,1])),
                         avg_logFC=c(as.numeric(CO.down[,'avg_logFC']),as.numeric(CO.up[,'avg_logFC']),as.numeric(OY.down[,'avg_logFC']),as.numeric(OY.up[,'avg_logFC'])),
                         p_val=c(as.numeric(CO.down[,'p_val']),as.numeric(CO.up[,'p_val']),as.numeric(OY.down[,'p_val']),as.numeric(OY.up[,'p_val'])),
                         p_val_adj=c(as.numeric(CO.down[,'p_val_adj']),as.numeric(CO.up[,'p_val_adj']),as.numeric(OY.down[,'p_val_adj']),as.numeric(OY.up[,'p_val_adj'])))
    }
    
    all.deg=rbind(all.deg,type.df)
    
  }
  
  write.xlsx(all.deg,paste(data.path,sample.name,"\\",sample.name,'.all.DEG.xlsx',sep = ''),sheetName = sample.name,row.names = F,append = T)
  
}

################# generate heatmap.matrix
data.path="G:\\result\\result\\09_sc_new\\heatmap\\"
name.list=list.files(data.path)
name.list=c("Brain","Muscle")
sample.name=name.list[1]

for (sample.name in name.list) {
  
all.deg=read.xlsx(paste(data.path,sample.name,"\\",sample.name,'.all.DEG.xlsx',sep = ''),sheetIndex=1)
all.deg$class2= unlist(strsplit(as.character(all.deg$class),split = "[.]"))[seq(1,nrow(all.deg)*2,by=2)]
all.deg$class2=paste(all.deg$Cell_type,all.deg$class2,sep = '_')
all.gene=as.character(unique(all.deg$Gene_name)) 
all.cell=unique(as.character(all.deg$Cell_type) ) 

heat.df=matrix(0,ncol = length(all.cell)*2,nrow = length(all.gene))
colnames(heat.df)=c(paste(all.cell,"_OvsY",sep = ""),paste(all.cell,"_CRvsO",sep = ""))

for (i in 1:length(all.gene) ) {
gene=all.gene[i]
all.deg[all.deg$Gene_name==gene,]
tmp.class=all.deg[all.deg$Gene_name==gene,]$class2
heat.df[i,tmp.class]=all.deg[all.deg$Gene_name==gene,]$avg_logFC
}
rownames(heat.df)=all.gene
write.csv(heat.df,paste(data.path,sample.name,"\\",sample.name,'.heatmap.csv',sep = ''),row.names=T)
}

############################################# heatmap plotting

data.path="G:\\result\\result\\09_sc_new\\heatmap\\"
plot.path='G:\\result\\result\\09_sc_new\\heatmap\\'

data.list=list.files(data.path)
library(pheatmap)

check.pos= function(df,n){
  tmp.df=df
  or=c()
  
  for (i in 1:n) {
    tmp.or= isTRUE(tmp.df[i]>0 & tmp.df[i+n]<0) 
    or=c(or,tmp.or)
  }
  or.final= any(or)
  return(or.final)
}

check.neg= function(df,n){
  tmp.df=df
  or=c()
  
  for (i in 1:n) {
    tmp.or= isTRUE(tmp.df[i]<0 & tmp.df[i+n]>0) 
    or=c(or,tmp.or)
  }
  or.final= any(or)
  return(or.final)
}


data.list=name.list
 
for (l in 1:2 ) {
  
  #plot.data= plot.data[ which(rowSums(plot.data)==0) ,]
  plot.data=read.csv(paste(data.path,data.list[l],"\\",data.list[l],'.heatmap.csv',sep = '') )
  rownames(plot.data)=plot.data$X
  
  plot.data=plot.data[,-1]
  
  perc0.df=data.frame(perc=apply(plot.data,2,perc0)) 
  perc0.df$type=unlist(strsplit(rownames(perc0.df),split = '_'))[seq(1,nrow(perc0.df)*2,2)] 
  perc0.df=data.frame(perc=tapply(perc0.df$perc,perc0.df$type,mean) ) 
  perc0.df$name=rownames(perc0.df)
  
  colnames(plot.data) = gsub('[+]',replacement = '', colnames(plot.data))
  
  plot.data=data.frame(plot.data)
  
  all.type=colnames(plot.data)
  
  n.type=length(all.type)/2
  
  #name.list=substr(all.type[1:n.type],1,nchar(all.type)-5)
  
  name.list=rev(perc0.df[order(perc0.df$perc),]$name) 
  
  plot.data=data.frame(plot.data)
  
  plot.data=plot.data[,paste(rep(name.list,2),rep(c('CRvsO','OvsY'),each=n.type),sep = '_' )] 
  
  plot.data$rank=seq(1,nrow(plot.data),1)
  
  com.df=data.frame(combn(name.list,1)) 

  
  for (o in 1:ncol(com.df)) {
    
    sample.list=as.character(com.df[,o]) 
    
    #other.list=setdiff(name.all,paste(sample.list,c('OvsY','CRvsO'),sep = '_'))
    
    sample.ab=data.frame(plot.data[,c(paste(sample.list,c('OvsY','CRvsO'),sep = '_'),'rank') ]) 
    
    #other.ab=data.frame(plot.data[,other.list]) 
    
    #tar.df= sample.ab[sample.ab[,1] * sample.ab[,2] <= 0 & (sample.ab[,1]!=0 | sample.ab[,2]!=0) ,]
    tar.df= sample.ab[sample.ab[,1] * sample.ab[,2] < 0 ,]
    assign(paste(sample.list,'rank',sep = '.'),tar.df$rank)
    
    tmp.rank1=tar.df$rank
    
    if(o==1){
      all.rank1=tmp.rank1
    }else{
      all.rank1=c(all.rank1,tmp.rank1)
    }
    
  }
  
  
  get_merge=function(sample){
    
    merge.rank = intersect(get(paste(sample[1],'rank',sep = '.')),get(paste(sample[2],'rank',sep = '.')))
    
    if(length(sample)>2){
      
      for (n in 3:length(sample)) {
        
        tmp.rank= get(paste(sample[n],'rank',sep = '.'))
        merge.rank= intersect(merge.rank,tmp.rank) 
        
      } 
      
    }
    
    return(merge.rank)
    
  }
  
  other.rank=c()
  
  for (i in 2:n.type) {
    
    com.df=data.frame(combn(name.list,i)) 
    for (o in 1:ncol(com.df)) {
      
      sample.list=as.character(com.df[,o]) 
      tmp.other.rank=get_merge(sample=sample.list)
      
      other.rank=c(other.rank,tmp.other.rank)
      
    }
    
  }
  
  all.rank=c(all.rank1,other.rank)
  #all.rank.rev=rev(all.rank)
  all.rank.rev=data.frame(rank=unique(rev(all.rank)))
  
  plot.df= merge(all.rank.rev,plot.data,by='rank',sort=F)
  
  rownames(plot.df)=plot.df$rank
  plot.df=plot.df[,-1]
  
  plot.df1=plot.df[apply(plot.df[,1:n.type],1, function(x) any(x>0)) &apply(plot.df[,(n.type+1):(n.type*2)],1, function(x) any(x<0)) ,] 
  plot.df1=plot.df1[apply(plot.df1, 1, check.pos,n=n.type),]
  plot.df2=plot.df[apply(plot.df[,1:n.type],1, function(x) any(x<0)) &apply(plot.df[,(n.type+1):(n.type*2)],1, function(x) any(x>0)) ,]
  plot.df2=plot.df2[apply(plot.df2, 1, check.neg,n=n.type),]
  
  
  no.r.data= plot.data[!(plot.data$rank %in% all.rank.rev$rank),]
  
  no.r.data$pos=rowSums(as.matrix(no.r.data[,(n.type+1):(n.type*2)]) > 0)
  no.r.pos=no.r.data[order(no.r.data$pos,decreasing = T),]
  no.r.pos=no.r.pos[no.r.pos$pos>0,]  
  
  no.r.data$neg=rowSums(as.matrix(no.r.data[,(n.type+1):(n.type*2)]) < 0)
  no.r.neg=no.r.data[order(no.r.data$neg,decreasing = T),]
  no.r.neg=no.r.neg[no.r.neg$neg>0,]  
  
  #############################################
  
  no.r.data[,'CR.pos']=rowSums(as.matrix(no.r.data[,1:n.type]) > 0)
  no.r.data[,'CR.neg']=rowSums(as.matrix(no.r.data[,1:n.type]) < 0)
  
  
  no.r.cr.pos=no.r.data[(no.r.data$pos+no.r.data$neg)==0 & (no.r.data$CR.pos > 0) ,]
  no.r.cr.pos=no.r.cr.pos[order(no.r.cr.pos$CR.pos,decreasing = T),]
  no.r.cr.neg=no.r.data[(no.r.data$pos+no.r.data$neg)==0 & (no.r.data$CR.neg > 0) ,]
  no.r.cr.neg=no.r.cr.neg[order(no.r.cr.neg$CR.neg,decreasing = T),]
  
  ##############################################
  
  if(length(intersect(no.r.pos$rank,no.r.neg$rank))==0){
    
    plot.df3=no.r.pos[,1:(n.type*2)]
    plot.df4=no.r.neg[,1:(n.type*2)]
    plot.df5=plot.df3[F,]
  }else{
    
    plot.df3=no.r.pos[- which(no.r.pos$rank %in% intersect(no.r.pos$rank,no.r.neg$rank)),1:(n.type*2)]
    plot.df4=no.r.neg[- which(no.r.neg$rank %in% intersect(no.r.pos$rank,no.r.neg$rank)),1:(n.type*2)]
    
    plot.df5= no.r.pos[(no.r.pos$rank %in% intersect(no.r.pos$rank,no.r.neg$rank)),1:(n.type*2)]
  }
  
  ##########################################
  
  if(length(intersect(no.r.cr.pos$rank,no.r.cr.neg$rank))==0){
    
    plot.df6=no.r.cr.pos[,1:(n.type*2)]
    plot.df7=no.r.cr.neg[,1:(n.type*2)]
    plot.df8=plot.df6[F,]
  }else{
    
    plot.df6=no.r.cr.pos[- which(no.r.cr.pos$rank %in% intersect(no.r.cr.pos$rank,no.r.cr.neg$rank) ),1:(n.type*2)]
    plot.df7=no.r.cr.neg[- which(no.r.cr.neg$rank %in% intersect(no.r.cr.pos$rank,no.r.cr.neg$rank)),1:(n.type*2)]
    plot.df8= no.r.cr.pos[(no.r.cr.pos$rank %in% intersect(no.r.cr.pos$rank,no.r.cr.neg$rank)),1:(n.type*2)]
  }
  
  ##########################################
  rescue.df=rbind(plot.df1,plot.df2)
  no.res.df=rbind(plot.df3,plot.df5)
  no.res.df=rbind(no.res.df,plot.df4)
  
  no.res.cr.df=rbind(plot.df6,plot.df8)
  no.res.cr.df=rbind(no.res.cr.df,plot.df7)
  
  stat.df=data.frame(row.names = c('rescue_up','rescue_down','no_rescue_up','no_rescue_down','no_rescue_overlap','no_rescue_CR_up','no_rescue_CR_down','no_rescue_CR_overlap'),
                     Counts=c(nrow(plot.df1),nrow(plot.df2),nrow(plot.df3),nrow(plot.df4),nrow(plot.df5),nrow(plot.df6),nrow(plot.df7),nrow(plot.df8)  ) )
  
  p.df=rbind(rescue.df,no.res.df)
  p.df=rbind(p.df,no.res.cr.df)
  
  col_anno=data.frame(row.names = all.type, type= rep(c('OvsY','CRvsO'),each=n.type))
  anno_color=list(type=c('OvsY'='goldenrod2','CRvsO'='seagreen3'))
  
  pheatmap(p.df,
           main = data.list[l],
           cluster_rows = F,
           cluster_cols = FALSE,
           scale = "none",
           color =  c(colorRampPalette(c('#4575b4', "grey90"))(50),colorRampPalette(c("grey90",'#d73027'))(50)),
           #breaks=seq(-2,2,length.out = 1000),
           breaks=seq(-0.01,0.01,length.out = 100),
           border_color = NA, #"grey60",
           fontsize = 10,
           fontsize_row = 12, #fontsize,
           fontsize_col = 12, #fontsize,
           fontsize_number = 0.8* fontsize,
           show_rownames = FALSE,
           show_colnames = T,
           legend = TRUE,
           gaps_row = cumsum(c(nrow(rescue.df),nrow(no.res.df)) ),
           gaps_col = c(n.type),
           annotation_col=col_anno,
           annotation_colors =anno_color  ,
           annotation_legend = T,
           revC=T,
           height = 18,
           width = 4.5,
           filename = paste(plot.path,data.list[l],'.pdf',sep = '')
           
  )

  
  pheatmap(p.df,
           main = data.list[l],
           cluster_rows = F,
           cluster_cols = FALSE,
           scale = "none",
           color =   c(colorRampPalette(c('#4575b4', "grey90"))(750),colorRampPalette("grey90")(500),colorRampPalette(c("grey90",'#d73027'))(750)),
           breaks=seq(-2,2,length.out = 2000),
           legend_breaks=seq(-2,2,by=0.5),
           legend_labels = seq(-2,2,by=0.5),
           border_color = NA, #"grey60",
           fontsize = 10,
           fontsize_row = 12, #fontsize,
           fontsize_col = 12, #fontsize,
           fontsize_number = 0.8* fontsize,
           show_rownames = FALSE,
           show_colnames = T,
           legend = TRUE,
           gaps_row = cumsum(c(nrow(rescue.df),nrow(no.res.df)) ),
           gaps_col = c(n.type),
           annotation_col=col_anno,
           annotation_colors =anno_color  ,
           annotation_legend = T,
           revC=T,
           height = 18,
           width = 4.5,
           filename = paste(plot.path,data.list[l],'2.pdf',sep = '')
           
  )
  
  write.csv(stat.df,paste(plot.path,data.list[l],'.stat.csv',sep = ''),row.names = T)
  write.csv(p.df,paste(plot.path,data.list[l],'.plot.csv',sep = ''),row.names = T)
}
