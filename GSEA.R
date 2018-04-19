#############################################################
#   2017.12.22    Kunbang Liu
#
#  Function: 将探针与 Gene symbol对应起来
#
#    result：1.输出 GSEA的输入文件
#            
#    
##############################################################


Probe2Symbol <- function(all_File, gseaPath, sampleInfoFile){

  # 打印分析开始时间
  cat("Start the CEL data Analysis.\n")
  sttime <- print( Sys.time() )
  cat(paste("StartTime:",sttime,sep=" "),"\n\n")


  # 加载数据库和包
  library(AnnotationDbi)
  library(hgu133plus2.db)


  # 读取数据
  all_info <- read.delim(all_File)
  rownames(all_info) <- all_info[,1]
  
  sampleInfo <- read.delim(sampleInfoFile)
  g.cls <- paste(
             paste(nrow(sampleInfo), length(unique(sampleInfo$groupName)), 1),
		     paste('#', paste(unique(sampleInfo$groupName), collapse=" ")),
             paste(paste(sampleInfo$groupName, collapse=" "), "\n", sep=""),
           sep="\n")


  # 过滤质控探针
  all_info <- all_info[grep("^[^A-Z]",rownames(all_info)),]


  # 找探针对应的 gene symbol
  Probe2Genesymbol <- select(hgu133plus2.db, keys=rownames(all_info), 
                      columns=c('SYMBOL', 'GENENAME'), keytype='PROBEID')
  

  # 合并输出数据
  result <- merge(Probe2Genesymbol, all_info, by.x=1 , by.y=1)
  result <- result[!duplicated(result[,2]),c(2:(ncol(result)-6))]
  colnames(result)[1:2] <- c("NAME", "Description")
  result[,2] <- "na"
  
  write.table(result, paste(gseaPath, "/data.gct", sep=""),
              sep="\t", quote=F, row.names=F)

  write.table(g.cls, paste(gseaPath, "/group.cls", sep=""),
              row.names=F, col.names=F, quote=F)

  # 收尾
  cat("\nCongratulations\n")
  endtime<- print(Sys.time())
  cat(paste("EndTime:",endtime,sep=" "),"\n")
  print(difftime(endtime,sttime))
  rm(list=ls())

}
