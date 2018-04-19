#############################################################
#   2017.12.22    Kunbang Liu
#
#  Function: 对差异基因分类
#
#    result：1.输出上下调统计文件
#            2.输出差异基因柱状图
#    
##############################################################



###### 子函数：用于设定y轴的最大轴
maxx <- function(x){ 

  # 差异个数小于10，设y轴最大值为10
  if( x<11 ) ret <- 10;
 
  # 差异个数大于10，取大于该值前两位数被5整除的最小整数
  if( x>10 ){
    ret <- signif(x,1)
    if( (x-ret) >0 ) ret <- ret+0.05*10**nchar(as.character(x))
   }
   ret
}


###### 绘图主函数
Plot <- function( inputdegFile, inputallFile, PlotPath, sampleInfoFile ){


  cat("Start the Histogram.\n")
  sttime <- print( Sys.time() )
  cat(paste("StartTime:",sttime,sep=" "),"\n\n")
  
  # 设置文件名
  barplot_name <- paste(PlotPath,"/barPlot.pdf",sep="")
  volplot_name <- paste(PlotPath,"/volPlot.pdf",sep="")
  scaplot_name <- paste(PlotPath,"/scaPlot.pdf",sep="")
  file_name    <- paste(PlotPath,"/deg_summary.txt",sep="")
  mergePlot_name <- paste(PlotPath,"/merge_vol_sca_Plots.pdf",sep="")
  
  
  # 读入数据
  datas      <- read.delim(inputdegFile, row.names=1)
  all.data   <- read.delim(inputallFile, row.names=1)
  sampleInfo <- read.delim(sampleInfoFile)
  
  
  # 整理数据
  groups     <- unique(sampleInfo$groupName)
  matchs     <- paste(groups[1], "vs", groups[2])
  s1 <- as.character(sampleInfo[sampleInfo[,"groupName"] == groups[1],"sampleName"])
  s2 <- as.character(sampleInfo[sampleInfo[,"groupName"] == groups[2],"sampleName"])
  x <- apply(all.data[,s1,drop=F], 1, function(x) mean(x, na.rm = T))
  y <- apply(all.data[,s2,drop=F], 1, function(x) mean(x, na.rm = T))
  x.y <- cbind(x,y) 
  
  
  # 获得差异基因上下调的个数
  up_num   <- nrow(datas[datas$logFC >= 1,])
  down_num <- nrow(datas[datas$logFC <= -1,])
  summ     <- cbind(up_num, down_num)
  rownames(summ) <- matchs
  colnames(summ) <- c("up", "down")
  
  
  # 设置画图颜色、图形边界长度、坐标标签相关参数
  barcols <- c("royalblue")
  volcolor <- rep(1, nrow(all.data))
  volcolor[all.data[,"logFC"] >= 1 & all.data[,"P.Value"] <= 0.05] <- 2
  volcolor[all.data[,"logFC"] <= -1 & all.data[,"P.Value"] <= 0.05] <- 3
  scacolor <- rep(1, nrow(x.y))
  scacolor[all.data[,"logFC"] >= 1 & all.data[,"P.Value"] <= 0.05]<- 2           #上调
  scacolor[all.data[,"logFC"] <= -1 & all.data[,"P.Value"] <= 0.05] <- 3          #下调
  ylim <- c(0,maxx(max(summ)))
  xlim <- c(0,length(rownames(summ))*2+0.5)


###################### 柱状图
  # 设置图形参数
  pdf(barplot_name, width = 6, height = 4)
  par(mai=c(0.5,0.7,0.5,0.5))


  # 画图
  barplot(summ, 
                ylim=ylim,
                xlim=xlim, 
                col=barcols, 
                width = 0.5,
                axes=F, 
                border=0,
                space = 1)


  # 绘纵坐标轴
  axis(2,c(seq(0, ylim[2], ylim[2]/5)), cex.axis=0.8)


  # 绘等距横线
  abline(a=0, h=c(seq(0,ylim[2],ylim[2]/5)), 0, col = "azure4")
  

  # 设置图片外边界
  box(which = "inner",col = "gray21")

  
  # 将等距横线放置到柱形图后面的方法，比较笨的方法!!!!
  barplot(summ, 
                ylim=ylim,
                xlim=xlim, 
                col=barcols, 
                width = 0.5,
                axes=F, 
                border=0,
                space = 1, 
                add=T)

  dev.off()


#####################  火山图  
  pdf(volplot_name)
  par(mai=c(1,1.2,0.5,0.5))
  
  plot(all.data[,"logFC"],
             -log10(all.data[,"P.Value"]),
	         xlab="log2 (Fold Change)", 
             ylab="-log10 (p.value)", 
	         main="Volcano Plot", 
	         cex.lab=1.5, 
	         cex.main=1.5, 
	         col=volcolor,
	         pch=20, 
	         cex=0.4)
	   
  box(which = "inner",col = "gray21")

  # 绘差异基因分割线
  abline( v=c(log2(0.5),-log2(0.5)), col=8)
  abline( h=-log10(0.05), 
          text(ceiling(max(all.data[,"logFC"]))/2+2,
		       -log10(0.05+0.015), 
			   cex=0.8,
		       labels="p.value = 0.05"), 
		  cex=1, 
		  col=8)

  dev.off()

  
  
#####################  散点图  
  pdf(scaplot_name)
  par(mai=c(1,1.2,0.5,0.5))

  plot(x.y[,2], 
             x.y[,1],
             col=scacolor,
	         pch=20, 
	         main="Scatter Plot", 
	         axes=T, 
	         cex= 0.4, 
             cex.axis=1, 
	         cex.lab=1.5, 
	         cex.main=1.5, 
	         xlab=groups[2], 
	         ylab=groups[1])

  box(which = "inner",col = "gray21")
  dev.off()
  
  
  
################# 合并 火山图
  pdf(mergePlot_name, width=14, height=6)
  par(mfrow=c(1,2), mai=c(1,1.2,0.5,0.5))
  
  plot(all.data[,"logFC"],
             -log10(all.data[,"P.Value"]),
	         xlab="log2 (Fold Change)", 
             ylab="-log10 (p.value)", 
	         main="Volcano Plot", 
	         cex.lab=1.5, 
	         cex.main=1.5, 
	         col=volcolor,
	         pch=20, 
	         cex=0.4)
	   
  box(which = "inner",col = "gray21")

  abline( v=c(log2(0.5),-log2(0.5)), col=8)
  abline( h=-log10(0.05), 
          text(ceiling(max(all.data[,"logFC"]))/2+2,
		       -log10(0.05+0.015), 
			   cex=0.8,
		       labels="p.value = 0.05"), 
		  cex=1, 
		  col=8)

#################### 合并 散点图
  plot(x.y[,2], 
             x.y[,1],
             col=scacolor,
	         pch=20, 
	         main="Scatter Plot", 
	         axes=T, 
	         cex= 0.4, 
             cex.axis=1, 
	         cex.lab=1.5, 
	         cex.main=1.5, 
	         xlab=groups[2], 
	         ylab=groups[1])

  box(which = "inner",col = "gray21")

  dev.off()


  # 输出差异基因分类表
  summ <- cbind(0,summ)
  summ[1,1] <- rownames(summ)
  colnames(summ)[1] <- "groups"
  write.table(summ, file=file_name, sep="\t", quote=F, row.names=F)
  

  # 收尾
  cat("\nCongratulations\n")
  endtime<- print(Sys.time())
  cat(paste("EndTime:",endtime,sep=" "),"\n")
  print(difftime(endtime,sttime))
  rm( list=ls() )


}
