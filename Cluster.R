####################################################
#  2017.12.23   Kunbang Liu
#
#  Function: 对差异探针和样本进行聚类
#
#    result：1.输出聚类热图
#    
##############################################################


Cluster <- function(ClusterPath){


  ### 分析开始时间
  cat("Start the CEL data Analysis.\n")
  sttime <- print( Sys.time() )
  cat(paste("StartTime:",sttime,sep=" "),"\n\n")


  ### 加载所需软件包
  library(gplots)


  ### 读取输入文件
  datas <- as.matrix(read.table(paste(ClusterPath,"/cluster.txt",sep=""), 
                                header = T, sep ="\t", row.names = 1) )


  ### 调用 pdf 函数，将图形保存到 Cluster_heatmapPlot.pdf 中
  pdf(paste(ClusterPath, "/Cluster_heatmapPlot.pdf", sep=""),
      height=6, width=7
     )


  ### 使用 heatmap 函数画图
  heatmap(datas, 
          labRow=F, 
          col=bluered(300), 
          margins=c(max(nchar(colnames(datas)))+2,2)
         )

  
  ### 设置图片外边界
  box(which = "inner",col = "gray21")

  ### 最后关闭图形函数
  dev.off()
  
  
  ### 收尾
  cat("\nCongratulations\n")
  endtime<- print(Sys.time())
  cat(paste("EndTime:",endtime,sep=" "),"\n")
  print(difftime(endtime,sttime))
  rm( list=ls() )
  
}
 
