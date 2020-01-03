#####################################
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
data.path="G:\\result\\result\\09_sc_new\\Goplot\\GO-terms\\"

#tissue.list=c("BAT","WAT","Liver", "Kidney", "Aorta", "Skin", "BM")
class.list[1]="Fig7"
go.df=read.csv(paste(data.path,class.list[1],"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
heat.df=read.csv(paste(data.path,class.list[1],"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0050673","R-RNO-1474244","R-RNO-381426","GO:0031099","GO:0097435",
             "GO:0032963","GO:0071363","GO:0001944","GO:0048545","GO:0061061")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1,10)]
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))
i=1
for (i in 1:nrow(heat.df)) {
  
  if(heat.df[i,]$value != 0){
    
    
    des=as.character(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}
library(SpectralTAD)

heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
#order.list[9]="cellular response to growth factor stimulusX"
#heat.df[grep("cellular response to growth factor stimulus",heat.df$Description),]$Description="cellular response to growth factor stimulusX"

heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=factor(heat.df$name,levels = tissue.list)

plot.df=heat.df
pdf(paste(data.path,class.list[1],'.pdf',sep = ''),width = 10,height = 7)
ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(10,40,70),range = c(-1,9),name='Size')+
  scale_color_gradientn(colors = brewer.pal(7,'YlGnBu') ,breaks=c(0,8,16,24),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black'),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
dev.off()

#####################################
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
data.path="G:\\result\\result\\09_sc_new\\Goplot\\GO-terms\\"
class.list=list.files(data.path)
class.list=class.list[-2]

go.df=read.csv(paste(data.path,class.list[2],"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
heat.df=read.csv(paste(data.path,class.list[2],"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0010942","GO:0009611","GO:0070555","GO:0009991","GO:0006954",
          "GO:0002237","R-RNO-168249","GO:0019882","GO:0008285","GO:0007568")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1,10)]
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))
i=1
for (i in 1:nrow(heat.df)) {
  
  if(heat.df[i,]$value != 0){
    
    
    des=as.character(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=factor(heat.df$name,levels = tissue.list)

plot.df=heat.df
range(plot.df$logP)
range(plot.df$count)
pdf(paste(data.path,class.list[2],'.pdf',sep = ''),width = 10,height = 7)
ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(10,40,70),range = c(-1,9),name='Size')+
  scale_color_gradientn(colors = brewer.pal(7,'YlOrRd') ,breaks=c(0,10,20,30),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black'),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
dev.off()

#####################################3
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
data.path="G:\\result\\result\\09_sc_new\\Goplot\\GO-terms\\"
#class.list=list.files(data.path)
#class.list=class.list[-3]

go.df=read.csv(paste(data.path,class.list[3],"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
heat.df=read.csv(paste(data.path,class.list[3],"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=c("GO:0010942","GO:0006954","GO:0032496","R-RNO-168249","GO:0051384",
          "rno04657","GO:0009636","GO:0070555","GO:0072593","GO:0009408")
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1,10)]
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))
i=1
for (i in 1:nrow(heat.df)) {
  
  if(heat.df[i,]$value != 0){
    
    
    des=as.character(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=factor(heat.df$name,levels = tissue.list)

plot.df=heat.df
pdf(paste(data.path,class.list[3],'.pdf',sep = ''),width = 10,height = 7)
ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(10,20,40),range = c(-1,9),name='Size')+
  scale_color_gradientn(colors = brewer.pal(7,'YlGnBu') ,breaks=c(0,5,10,20),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black'),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
dev.off()

#####################################
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
data.path="G:\\result\\result\\09_sc_new\\Goplot\\GO-terms\\"
#class.list=list.files(data.path)
#class.list=class.list[-4]

go.df=read.csv(paste(data.path,class.list[4],"\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
heat.df=read.csv(paste(data.path,class.list[4],"\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=as.character(heat.df$GO[c(1,3,7,9,10,16,20,21,27)])
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1,10)]
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))

for (i in 1:nrow(heat.df)) {
  
  if(heat.df[i,]$value != 0){
    
    
    des=as.character(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=factor(heat.df$name,levels = tissue.list)

plot.df=heat.df
range(plot.df$logP)
range(plot.df$count)
pdf(paste(data.path,class.list[4],'.pdf',sep = ''),width = 10,height = 7)
ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(5,10,20),range = c(-1,9),name='Size')+
  scale_color_gradientn(colors = brewer.pal(7,'YlOrRd') ,breaks=c(0,3,6,9),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black'),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
dev.off()

####################################

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
data.path="G:\\result\\result\\09_sc_new\\Goplot\\GO-terms\\"
#class.list=list.files(data.path)
#class.list=class.list[-4]

go.df=read.csv(paste(data.path,"Fig7\\Enrichment_GO\\GO_AllLists.csv",sep = "") )
heat.df=read.csv(paste(data.path,"Fig7\\Enrichment_heatmap\\HeatmapSelectedGOTop100.csv",sep = "") )
go.list=as.character(heat.df$GO[c(1,2,3,4,6,13,16,17,20,29)])
heat.df=heat.df[heat.df$GO %in% go.list,]
heat.df=heat.df[,-c(1,10)]
heat.df$count=rowSums(as.matrix(heat.df[,2:8]) != 0)
heat.df$sumP=abs(rowSums(as.matrix(heat.df[,2:8]))) 
order.list=as.character(heat.df[order(heat.df$count,heat.df$sumP),'Description']) 
heat.df=heat.df[,-c(9,10)]
heat.df= melt(heat.df,id.vars = c('Description'))

for (i in 1:nrow(heat.df)) {
  
  if(heat.df[i,]$value != 0){
    
    
    des=as.character(heat.df[i,]$Description) 
    gro=as.character(heat.df[i,]$variable)
    #gro=substr(gro,8,nchar(gro)) 
    
    count=go.df[go.df$Description==des & go.df$GeneList== gro, ]$X.GeneInGOAndHitList
    
  }else{
    count=0
  }
  
  if(i==1){
    all.count=count
  }else{
    all.count=c(all.count,count)
  }
  
}

heat.df$count=all.count
heat.df$name=heat.df$variable
heat.df$Description=as.character(heat.df$Description)
#order.list[5]="IGF"
#heat.df[grep("IGF",heat.df$Description),]$Description="IGF"
heat.df$logP=abs(heat.df$value) 
heat.df$Description=factor(heat.df$Description,levels = order.list)
heat.df$name=factor(heat.df$name,levels = tissue.list)

plot.df=heat.df
range(plot.df$logP)
range(plot.df$count)
pdf(paste(data.path,'F7.pdf',sep = ''),width = 10,height = 7)
ggplot(plot.df, aes(x=name,y=Description,size=count,colour=logP))+
  geom_point()+
  scale_size_continuous(breaks = c(5,10,20),range = c(-1,9),name='Size')+
  scale_color_gradientn(colors = brewer.pal(7,'YlGnBu') ,breaks=c(0,5,10,20),name='LogP')+  theme_classic()+ylab('Go terms')+xlab('')+
  theme(legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.y  = element_text(size = 15,colour = 'black'),
        axis.text.x = element_text(size = 15,colour = 'black'),
        legend.title = element_text(size = 15),
        #axis.ticks = element_blank(),
        legend.position ="right",legend.direction = "vertical")
dev.off()
