###############################
library(xlsx)
data.path='G:\\result\\result\\09_sc_new\\heatmap\\'
name.list=list.files(data.path)

name.list=name.list[-1]

sample.name=name.list[1]

for (sample.name in name.list) {

tmp.path=paste(data.path,sample.name,"\\DEGs_all\\",sep = '')
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

write.xlsx(all.deg,paste(data.path,'all.DEG.xlsx',sep = ''),sheetName = sample.name,row.names = F,append = T)

}

####################################################################################
library(xlsxjars)
options(java.parameters="-Xmx8g")
sample.name=name.list[1]
data.path='G:\\result\\result\\09_sc_new\\heatmap\\'

class="OvsY.down"

sample.name=name.list[1]
all.list=c()
for (sample.name in name.list) {

tmp.df=read.xlsx(paste(data.path,'all.DEG.xlsx',sep = ''),sheetName = sample.name)
tmp.list=as.character(tmp.df[grep("CRvsO.down",tmp.df$class),]$Gene_name) 
all.list=c(all.list,tmp.list)
    
}

stat.df=matrix(ncol = length(name.list),nrow = length(unique(all.list) ))
stat.df[is.na(stat.df)]=0
rownames(stat.df)=unique(all.list)
colnames(stat.df)=name.list
stat.df=data.frame(stat.df)

sample.name=name.list[1]
for (sample.name in name.list) {
  
  tmp.df=read.xlsx(paste(data.path,'all.DEG.xlsx',sep = ''),sheetName = sample.name)
  tmp.list=as.character(tmp.df[grep("CRvsO.down",tmp.df$class),]$Gene_name) 
  stat.df[rownames(stat.df) %in% unique(tmp.list),sample.name]=1
  
}
write.csv(stat.df,paste(data.path,'DEG.',class,".csv",sep = ''),row.names = T)


#######  rescue

library(xlsx)
library(data.table)
data.path='G:\\result\\result\\09_sc_new\\heatmap\\'
name.list=list.files(data.path)

name.list=name.list[-1]

sample.name=name.list[1]

for (sample.name in name.list) {
  
  tmp.path=paste(data.path,sample.name,"\\DEGs_all\\",sep = '')
  type.list=list.files(tmp.path)[grep('CRvsO.down.list',list.files(tmp.path))] 
  type.list=substr(type.list,1,nchar(type.list)-16)
  all.deg=data.frame(Tissue=character(),Cell_type=character(),Gene_name=character(),avg_logFC=numeric(),p_val=numeric(),p_val_adj=numeric())
  for (type in type.list) {
    
    CO.down=data.frame(fread(paste(tmp.path,'\\',type,'_CRvsO.down.list',sep=''))) 
    CO.up=data.frame(fread(paste(tmp.path,'\\',type,'_CRvsO.up.list',sep=''))) 
    
    OY.down=data.frame(fread(paste(tmp.path,'\\',type,'_OvsY.down.list',sep=''))) 
    OY.up=data.frame(fread(paste(tmp.path,'\\',type,'_OvsY.up.list',sep=''))) 
    
    res.down=data.frame(fread(paste(tmp.path,'\\',type,'_res.down.list',sep=''))) 
    res.up=data.frame(fread(paste(tmp.path,'\\',type,'_res.up.list',sep=''))) 
    
    ng.down=data.frame(fread(paste(tmp.path,'\\',type,'_ng.down.list',sep='')))
    ng.up=data.frame(fread(paste(tmp.path,'\\',type,'_ng.up.list',sep=''))) 
    
    if(nrow(CO.down)+nrow(CO.up)+nrow(OY.down)+nrow(OY.up)==0){
      type.df=data.frame(Tissue=character(),Cell_type=character(),Gene_name=character(),avg_logFC=numeric(),p_val=numeric(),p_val_adj=numeric())
    }else{
      type.df=data.frame(Tissue=sample.name,Cell_type=type,
                         class=rep(c('CRvsO.down','CRvsO.up','OvsY.down','OvsY.up',"res.down","res.up","ng.down","ng.up"),c(nrow(CO.down),nrow(CO.up),nrow(OY.down),nrow(OY.up),nrow(res.down),nrow(res.up),nrow(ng.down),nrow(ng.up)   ) ),
                         Gene_name=c(as.character(CO.down[,1]),as.character(CO.up[,1]),as.character(OY.down[,1]),as.character(OY.up[,1]),as.character(res.down[,1]),as.character(res.up[,1]),as.character(ng.down[,1]),as.character(ng.up[,1]) ),
                         avg_logFC=c(as.numeric(CO.down[,'avg_logFC']),as.numeric(CO.up[,'avg_logFC']),as.numeric(OY.down[,'avg_logFC']),as.numeric(OY.up[,'avg_logFC']),as.numeric(res.down[,'avg_logFC']),as.numeric(res.up[,'avg_logFC']),as.numeric(ng.down[,'avg_logFC']),as.numeric(ng.up[,'avg_logFC'])   ),
                         p_val=c(as.numeric(CO.down[,'p_val']),as.numeric(CO.up[,'p_val']),as.numeric(OY.down[,'p_val']),as.numeric(OY.up[,'p_val']),as.numeric(res.down[,'p_val']),as.numeric(res.up[,'p_val']),as.numeric(ng.down[,'p_val']),as.numeric(ng.up[,'p_val'])),
                         p_val_adj=c(as.numeric(CO.down[,'p_val_adj']),as.numeric(CO.up[,'p_val_adj']),as.numeric(OY.down[,'p_val_adj']),as.numeric(OY.up[,'p_val_adj']),as.numeric(res.down[,'p_val_adj']),as.numeric(res.up[,'p_val_adj']),as.numeric(ng.down[,'p_val_adj']),as.numeric(ng.up[,'p_val_adj']))
                         )
    }
    
    all.deg=rbind(all.deg,type.df)
    
  }
  
  write.xlsx(all.deg,paste(data.path,'all.DEG3.xlsx',sep = ''),sheetName = sample.name,row.names = F,append = T)
  
}

################################
library(xlsxjars)
options(java.parameters="-Xmx8g")
sample.name=name.list[1]
data.path='G:\\result\\result\\09_sc_new\\heatmap\\'

class.list=c("res.up","res.down")

class="res.down"

for (class in class.list) {
  all.list=c()
for (sample.name in name.list) {
  
  tmp.df=read.xlsx(paste(data.path,'all.DEG2.xlsx',sep = ''),sheetName = sample.name)
  tmp.list=as.character(tmp.df[grep(class,tmp.df$class),]$Gene_name) 
  all.list=c(all.list,tmp.list)
  
}

stat.df=matrix(ncol = length(name.list),nrow = length(unique(all.list) ))
stat.df[is.na(stat.df)]=0
rownames(stat.df)=unique(all.list)
colnames(stat.df)=name.list
stat.df=data.frame(stat.df)


for (sample.name in name.list) {
  
  tmp.df=read.xlsx(paste(data.path,'all.DEG2.xlsx',sep = ''),sheetName = sample.name)
  tmp.list=as.character(tmp.df[grep(class,tmp.df$class),]$Gene_name) 
  stat.df[rownames(stat.df) %in% unique(tmp.list),sample.name]=1
  
}
write.csv(stat.df,paste(data.path,'DEG.',class,".csv",sep = ''),row.names = T)
}

############################################
class.list=c("res.up","res.down","CRvsO.down","CRvsO.up","OvsY.down","OvsY.up")
class=class.list[1]
class
for (class in class.list) {
tmp.df=read.csv(paste(data.path,"DEG.",class,".csv",sep = ''))
tmp.df$freq=rowSums(tmp.df[,2:8])
write.csv(tmp.df,paste(data.path,'DEG.freq.',class,".csv",sep = ''),row.names = F)
}


##############################################
class.list=c("CRvsO.down","CRvsO.up","OvsY.down","OvsY.up","res.up","res.down","ng.up","ng.down")
class=class.list[1]
sample.name=name.list[1]
all.freq=data.frame(Var1=c("CRvsO.down","CRvsO.up","OvsY.down","OvsY.up","res.up","res.down","ng.up","ng.down"))
sample.name=name.list[1]
for (sample.name in name.list) {
  
  sample.df=read.xlsx(paste(data.path,'all.DEG3.xlsx',sep = ''),sheetName = sample.name)
  all.count=c()
  for (class in class.list) {
  tmp.df=sample.df[sample.df$class==class,]
  n.count=length(unique(tmp.df$Gene_name))
  all.count=c(all.count,n.count)
  }
  
  tmp.freq=data.frame(Var1=class.list,count=all.count)
  colnames(tmp.freq)[2]=sample.name
  
  if(sample.name==name.list[1]){
  all.freq=tmp.freq
  }else{
  all.freq=merge(all.freq,tmp.freq,by="Var1")
  }
  
}


all.freq[,2:8][is.na(all.freq[,2:8])]=0
write.csv(all.freq,paste(data.path,'DEG.stat.csv',sep = ''),row.names = F)
