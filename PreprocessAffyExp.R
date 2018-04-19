#############################################################
#   2017.12.22    Kunbang Liu
#
#  Function: 对 Affy 芯片数据（CEL文件）进行预处理
#
#    result：1.输出质控图（箱线图，降解图，密度曲线图）
#            2.输出芯片表达原始数据表
#            3.输出有表达值的探针的统计数据表
#            4.输出显著差异探针数据表
#            5.输出聚类数据输入表
#    
##############################################################



PreprocessAffyExp <- function(sampleInfoFile, CelDataPath, MidfilePath){

  # 打印分析开始时间
  cat("Start the CEL data Analysis.\n")
  sttime <- print( Sys.time() )
  cat(paste("StartTime:",sttime,sep=" "),"\n\n")

  # 加载 R包
  library(gcrma)
  library(limma)
  library(affy)

  sampleInfo  <- read.delim(sampleInfoFile, header=T, sep="\t", dec=".")
  cols <- sampleInfo$group+1
  

  # 对芯片数据提取
  data.cel <- ReadAffy(filenames=sampleInfo$FileName, 
                       sampleNames=sampleInfo$sampleName,
                       celfile.path=CelDataPath
                       )
  raw.data <- log2(exprs(data.cel))

  
  # 数据归一化
  eset <- gcrma(data.cel)
  norm.data <- exprs(eset)

  
  # 绘制 RNA 降解曲线
  RNAdeg <- AffyRNAdeg(data.cel)
  pdf(paste(MidfilePath, "/QC/RNA_degradation_plot.pdf", sep=""),
      width=8, height=7)
  plotAffyRNAdeg(RNAdeg, col=cols)
  legend("topleft", legend=sampleNames(data.cel), lty=1, col=cols, box.col="transparent", xpd=TRUE)
  
  # 设置图片外边界
  box(which = "inner",col = "gray21")

  # 最后关闭图形函数
  dev.off()
  
  
  rm(RNAdeg)
  
  
  # 对归一化前后的数据绘制箱线图
  pdf(paste(MidfilePath, "/QC/Boxplot.pdf", sep=""), width=20, height=8)
  layout( matrix(1:2, 1, 2) )
  boxplot(raw.data, 
          main="Raw Intensity", 
		  ylab="Log2 (intensity)",
          col=cols)
		  
  boxplot(norm.data, 
          main="Norm Intensity", 
		  ylab="Log2 (intensity)",
          col=cols)
  
  # 设置图片外边界
  box(which = "inner",col = "gray21")

  # 最后关闭图形函数
  dev.off()


  # 绘制密度曲线
  pdf(paste(MidfilePath, "/QC/DensityPlot.pdf", sep=""))
  layout( matrix(1:2, 2, 1) )  
  plotDensity(raw.data, lty=1:3, col=cols, main="Raw Intensity", xlab="Log2 (intensity)")
  legend("topright", legend=sampleNames(data.cel), lty=seq(1,length(cols[cols==2])),
          cex=0.7, col=cols, box.col="transparent", xpd=TRUE)

  plotDensity(norm.data, lty=1:3, col=cols, main="Norm Intensity", xlab="Log2 (intensity)")
  legend("topright", legend=sampleNames(data.cel), lty=seq(1,length(cols[cols==3])),
          cex=0.7, col=cols, box.col="transparent", xpd=TRUE)
  
  # 设置图片外边界
  box(which = "inner",col = "gray21")

  # 最后关闭图形函数
  dev.off()
  

  # 筛选出至少一个芯片有表达值的探针(选取“表达”基因)
  eset.mas5calls <- exprs( mas5calls(data.cel) )
  PMA <- apply(eset.mas5calls, 1, function(x)any(x=="P"))
  present.probes <- names(PMA[PMA])
  
  rm(eset.mas5calls,PMA)
  
   
  # 分组比对
  groups <- unique(sampleInfo$groupName)
  G1vsG2 <- paste(groups[1], "-", groups[2], sep="")
  design <- model.matrix(~ 0+factor(sampleInfo$group))             #构建实验设计信息
  colnames(design) <- groups
  fit <- lmFit(eset, design)                                       #线性模型拟合
  cont.matrix <- makeContrasts(contrasts=G1vsG2, levels=design )   #构建对比模型，比较两实验条件下表达数据
  fit2 <- contrasts.fit(fit, cont.matrix)                          #根据比对模型进行差值运算
  fit2 <- eBayes(fit2)                                             #贝叶斯检验

  
  # 整理输出数据
  # raw.chip.exp <- cbind(0, raw.data)
  # raw.chip.exp[,1] <- rownames(raw.data)
  # colnames(raw.chip.exp)[1] <- "ProbeSetID"
  
  stat.inf <- topTable(fit2, coef=1, number=Inf)                   #生成所有基因检验的结果报表
  stat.inf <- stat.inf[rownames(norm.data),]
  all.inf <- cbind(0,norm.data[present.probes,], stat.inf[present.probes,])
  # all.inf <- cbind(0,norm.data, stat.inf)
  all.inf[,1] <- rownames(all.inf)
  colnames(all.inf)[1] <- "ProbeSetID"


  # 筛选显著差异探针( logFC大于等于 1或小于等于 -1，且 P值小于等于0.05)
  deg <- all.inf[all.inf$P.Value<=0.05 & (all.inf$logFC <= -1 | all.inf$logFC >= 1),]


  # 过滤质控探针
  deg <- deg[grep("^[^A-Z]",rownames(deg)),]


  # 聚类文件
  cluster <- deg[,0:ncol(norm.data)+1]


  # 输出结果
  # write.table(raw.chip.exp,
        # paste(MidfilePath, "/Preprocess/raw_chip_exp.txt", sep=""),
        # append=F, sep="\t", quote=F, dec=".", row.names=F, col.names=T
             # )        
  write.table(all.inf,
        paste(MidfilePath, "/Preprocess/all_expression_info.txt", sep=""),
        append=F, sep="\t", quote=F, dec=".", row.names=F, col.names=T
             )        
  write.table(deg,
        paste(MidfilePath, "/Preprocess/diff_expressed_genes.txt", sep=""),
        append=F, sep="\t", quote=F, dec=".", row.names=F, col.names=T
             )             
  write.table(cluster,
        paste(MidfilePath, "/Cluster/cluster.txt", sep=""),
        append=F, sep="\t", quote=F, dec=".", row.names=F, col.names=T
             )


  # 收尾
  cat("\nCongratulations\n")
  endtime<- print(Sys.time())
  cat(paste("EndTime:",endtime,sep=" "),"\n")
  print(difftime(endtime,sttime))
  rm(list=ls())

}
