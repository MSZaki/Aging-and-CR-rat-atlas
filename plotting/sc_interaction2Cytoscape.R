ibrary(optparse)

option_list <- list(
  make_option(c("-Y","--young"),type="character",help="Young interaction data (.csv)"),
  make_option(c("-O","--old"),type="character",help="Old interaction data (.csv)"),
  make_option(c("-C","--CR"),type="character",help="CR interaction data (.csv)"),
  make_option(c("-o","--out"),type="character",help="outpath"),
  make_option(c("-p","--prefix"),type="character",help="sample name"),
  make_option(c("-v","--verbose"),action="store_true",default=TRUE,help="Print extra output")
)

opt <- parse_args(OptionParser(option_list=option_list,usage = "export interaction data to Cytoscape network"))
###########################
library(data.table)
library(xlsx)
###########################

old.df=read.csv(opt$old) 
young.df=read.csv(opt$young)
cr.df=read.csv(opt$cr)
sample.name=opt$prefix
out=opt$out

old.df = "G:\\result\\result\\09_sc_new\\interaction\\data\\BAT\\interaction_count.O.csv"
cr.df = "G:\\result\\result\\09_sc_new\\interaction\\data\\BAT\\interaction_count.CR.csv"
young.df ="G:\\result\\result\\09_sc_new\\interaction\\data\\BAT\\interaction_count.Y.csv"
sample.name = "BAT"
out = "G:\\result\\result\\09_sc_new\\interaction\\data\\BAT\\"


interaction2Cytoscape=function(old.df,young.df,cr.df,sample.name,out){

###########################
old.df=read.csv(old.df) 
young.df=read.csv(young.df)
cr.df=read.csv(cr.df)
  
old.df$X=old.df$celltype.O
young.df$X=young.df$celltype.Y
cr.df$X=cr.df$celltype.CR

########################### Old vs Young  ( O Y )
OY.df=merge(data.frame(old.df),data.frame(young.df),by=c('X'))

OY.out=data.frame(celltype=OY.df$celltype.O,Count=OY.df$number.O.pvalues - OY.df$number.Y.pvalues)
OY.out$Ligand_celltype= unlist(strsplit(as.character(OY.out$celltype) ,split='_') )[seq(1,by=2,length.out = nrow(OY.out))] 
OY.out$Receptor_celltype=unlist(strsplit(as.character(OY.out$celltype) ,split='_'))[seq(2,by=2,length.out = nrow(OY.out))]

OY.out$Regulated='none'
OY.out[OY.out$Count>0,]$Regulated='up'
OY.out[OY.out$Count<0,]$Regulated='down'
OY.out$Count=abs(OY.out$Count)
OY.out=OY.out[,c('Ligand_celltype','Receptor_celltype','Count','Regulated')]

tapply.df1=data.frame(celltype=names(tapply(abs( OY.out$Count),OY.out$Ligand_celltype, sum)),L_Total_Count=unname(tapply(abs(OY.out$Count) ,OY.out$Ligand_celltype, sum)))
tapply.df2=data.frame(celltype=names(tapply(abs( OY.out$Count),OY.out$Receptor_celltype, sum)),R_Total_Count=unname(tapply(abs(OY.out$Count) ,OY.out$Receptor_celltype, sum)))
stat.df=merge(tapply.df1,tapply.df2,by='celltype',all=T)
stat.df$Total_Count=rowSums(stat.df[,2:3],na.rm = T)

write.xlsx(OY.out,paste(out,sample.name,'_OY_network.xlsx',sep = ''),row.names = F)
write.xlsx(stat.df,paste(out,sample.name,'_OY_stat.xlsx',sep = ''),row.names = F)

#############################  CR vs Old  ( CO )

CO.df=merge(data.frame(cr.df),data.frame(old.df),by=c('X'))

CO.out=data.frame(celltype=CO.df$celltype.CR,Count=CO.df$number.CR.pvalues - CO.df$number.O.pvalues)
CO.out$Ligand_celltype= unlist(strsplit(as.character(CO.out$celltype) ,split='_') )[seq(1,by=2,length.out = nrow(CO.out))] 
CO.out$Receptor_celltype=unlist(strsplit(as.character(CO.out$celltype) ,split='_'))[seq(2,by=2,length.out = nrow(CO.out))]

CO.out$Regulated='none'
CO.out[CO.out$Count>0,]$Regulated='up'
CO.out[CO.out$Count<0,]$Regulated='down'
CO.out$Count=abs(CO.out$Count)
CO.out=CO.out[,c('Ligand_celltype','Receptor_celltype','Count','Regulated')]

tapply.df1=data.frame(celltype=names(tapply(abs( CO.out$Count),CO.out$Ligand_celltype, sum)),L_Total_Count=unname(tapply(abs(CO.out$Count) ,CO.out$Ligand_celltype, sum)))
tapply.df2=data.frame(celltype=names(tapply(abs( CO.out$Count),CO.out$Receptor_celltype, sum)),R_Total_Count=unname(tapply(abs(CO.out$Count) ,CO.out$Receptor_celltype, sum)))
stat.df=merge(tapply.df1,tapply.df2,by='celltype',all=T)
stat.df$Total_Count=rowSums(stat.df[,2:3],na.rm = T)

write.xlsx(CO.out,paste(out,sample.name,'_co_network.xlsx',sep = ''),row.names = F)
write.xlsx(stat.df,paste(out,sample.name,'_co_stat.xlsx',sep = ''),row.names = F)

}

interaction2Cytoscape(old.df = "G:\\result\\result\\09_sc_new\\interaction\\data\\Skin\\interaction_count.O.csv",
                      cr.df = "G:\\result\\result\\09_sc_new\\interaction\\data\\Skin\\interaction_count.CR.csv",
                      young.df ="G:\\result\\result\\09_sc_new\\interaction\\data\\Skin\\interaction_count.Y.csv",
                      sample.name = "Skin",
                      out = "G:\\result\\result\\09_sc_new\\interaction\\data\\Skin\\")
