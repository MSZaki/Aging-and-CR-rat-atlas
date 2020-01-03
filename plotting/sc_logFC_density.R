##############################
library(xlsx)

class.list=c("Skin","BM","WAT","BAT","BAT","Liver","Kidney","Aorta")

for (class in class.list) {
 sc.df=read.xlsx("G:\\result\\result\\09_sc_new\\heatmap\\all.DEG2.xlsx",sheetName = class) 
 rescue.up=sc.df[grep("res.up",sc.df$class),]
 for (i in 1:nrow(rescue.up)) {
  
name=as.character(rescue.up$Gene_name[i]) 
tissue=as.character(rescue.up$Tissue[i]) 
type=as.character( rescue.up$Cell_type[i])
fc=as.numeric( rescue.up$avg_logFC[i])

co.up=sc.df[sc.df$Gene_name==name & sc.df$Tissue==tissue&sc.df$Cell_type==type&sc.df$class=="CRvsO.up",]$avg_logFC
oy.down=sc.df[sc.df$Gene_name==name & sc.df$Tissue==tissue&sc.df$Cell_type==type&sc.df$class=="OvsY.down",]$avg_logFC
stat.df=data.frame(Tissue=tissue,Cell_type=type,Gene_name=name,CRvsO_up=co.up,OvsY_down=oy.down,rescue_up=fc)

if(i==1 & class==class.list[1]){
  all.stat=stat.df
}else{
  all.stat=rbind(all.stat,stat.df)
}
}
}

write.csv(all.stat,"G:\\result\\result\\09_sc_new\\heatmap\\all.rescue.up.stat.csv",row.names = F)


for (class in class.list) {
  sc.df=read.xlsx("G:\\result\\result\\09_sc_new\\heatmap\\all.DEG2.xlsx",sheetName = class) 
  rescue.up=sc.df[grep("res.down",sc.df$class),]
  for (i in 1:nrow(rescue.up)) {
    
    name=as.character(rescue.up$Gene_name[i]) 
    tissue=as.character(rescue.up$Tissue[i]) 
    type=as.character( rescue.up$Cell_type[i])
    fc=as.numeric( rescue.up$avg_logFC[i])
    
    co.up=sc.df[sc.df$Gene_name==name & sc.df$Tissue==tissue&sc.df$Cell_type==type&sc.df$class=="CRvsO.down",]$avg_logFC
    oy.down=sc.df[sc.df$Gene_name==name & sc.df$Tissue==tissue&sc.df$Cell_type==type&sc.df$class=="OvsY.up",]$avg_logFC
    stat.df=data.frame(Tissue=tissue,Cell_type=type,Gene_name=name,CRvsO_down=co.up,OvsY_up=oy.down,rescue_down=fc)
    
    if(i==1 & class==class.list[1]){
      all.stat=stat.df
    }else{
      all.stat=rbind(all.stat,stat.df)
    }
  }
  
}

write.csv(all.stat,"G:\\result\\result\\09_sc_new\\heatmap\\all.rescue.down.stat.csv",row.names = F)

all.stat=read.csv("G:\\result\\result\\09_sc_new\\heatmap\\all.rescue.down.stat.csv")

all.stat$ratio=abs(all.stat$CRvsO_down/all.stat$OvsY_up) 

#nrow(all.stat[all.stat$ratio>0.5,])/nrow(all.stat)
#nrow(all.stat[all.stat$ratio>1,])/nrow(all.stat)

pdf(paste("G:\\result\\result\\09_sc_new\\heatmap\\rescue_down.pdf",sep = ""),width=6,height = 4)
ggplot(all.stat, aes(x = all.stat$ratio))+geom_density(color = "black", fill = "gray")+
  theme_few()+expand_limits(x=0)+
  geom_vline(xintercept =as.numeric(quantile(all.stat$ratio,0.05)))+
  geom_vline(xintercept =as.numeric(quantile(all.stat$ratio,0.95)))+coord_cartesian(ylim = c(0,2))
dev.off()

all.stat=read.csv("G:\\result\\result\\09_sc_new\\heatmap\\all.rescue.up.stat.csv")
all.stat$ratio=abs(all.stat$CRvsO_up/all.stat$OvsY_down) 
nrow(all.stat[all.stat$ratio>1,])/nrow(all.stat)

pdf(paste("G:\\result\\result\\09_sc_new\\heatmap\\rescue_up.pdf",sep = ""),width=6,height = 4)
ggplot(all.stat, aes(x = all.stat$ratio))+geom_density(color = "black", fill = "gray")+
  theme_few()+expand_limits(x=0)+
  geom_vline(xintercept =as.numeric(quantile(all.stat$ratio,0.05)))+
  geom_vline(xintercept =as.numeric(quantile(all.stat$ratio,0.95)))+coord_cartesian(ylim = c(0,2))
dev.off()
