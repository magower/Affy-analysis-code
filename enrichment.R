#############################################################
#   2017.12.22    Kunbang Liu
#
#  Function: 获得 GO、KEGG列表
#
#    result：1.输出差异富集的统计文件
#            2.富集柱状图
#            3.富集气泡图
#    
##############################################################


Enrichment <- function(deg_File, richPath){


  ### 分析开始时间
  cat("Start the CEL data Analysis.\n")
  sttime <- print( Sys.time() )
  cat(paste("StartTime:",sttime,sep=" "),"\n\n")


  # 加载数据库和包
  library(DOSE)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(clusterProfiler)
  
  
  # 读取数据
  deg_info <- read.delim(deg_File)
  rownames(deg_info) <- deg_info[,1]
  

  # 找探针对应的 entrez_gene_id
  Probe2symbol <- select(hgu133plus2.db, keys=rownames(deg_info), 
                  columns=c('SYMBOL', 'ENTREZID'), keytype='PROBEID')

  gene <- unique(Probe2symbol[,3])
  rm(Probe2symbol, deg_info)

  
  # 富集啦
  ego_cc <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     readable = TRUE)

  # 画图:条形图 ，按照P值大小排列；点图，按照富集数目大小排列
  pdf(paste(richPath, "/EnrichmentGO_CC.pdf", sep=""), height=4, width=10)
  print( barplot(ego_cc, showCategory=10, title="EnrichmentGO_CC") )
  dev.off()  
  
  pdf(paste(richPath, "/EnrichmentGO_CC_dot.pdf", sep=""), height=4, width=10)
  print( dotplot(ego_cc, title="EnrichmentGO_CC_dot") )
  dev.off()
  
  # 输出结果
  write.table(as.data.frame(ego_cc), 
              file=paste(richPath, "/Enrichment_CC.txt", sep=""),
              row.names=F, quote=F, sep="\t")

  rm(ego_cc)
  


  # 富集啦
  ego_bp <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     readable = TRUE)
					 
  pdf(paste(richPath, "/EnrichmentGO_BP.pdf", sep=""), height=4, width=10)
  print( barplot(ego_bp, showCategory=10, title="EnrichmentGO_BP") )
  dev.off()

  pdf(paste(richPath, "/EnrichmentGO_BP_dot.pdf", sep=""), height=4, width=10)
  print( dotplot(ego_bp, title="EnrichmentGO_BP_dot") )
  dev.off()

  write.table(as.data.frame(ego_bp), 
              file=paste(richPath, "/Enrichment_BP.txt", sep=""),
              row.names=F, quote=F, sep="\t")

  rm(ego_bp)

  
  # 富集啦
  ego_mf <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     readable = TRUE)
					 
  pdf(paste(richPath, "/EnrichmentGO_MF.pdf", sep=""), height=4, width=10)
  print( barplot(ego_mf, showCategory=10, title="EnrichmentGO_MF") )
  dev.off()

  pdf(paste(richPath, "/EnrichmentGO_MF_dot.pdf", sep=""), height=4, width=10)
  print( dotplot(ego_mf, title="EnrichmentGO_MF_dot") )
  dev.off()

  write.table(as.data.frame(ego_mf), 
              file=paste(richPath, "/Enrichment_MF.txt", sep=""),
              row.names=F, quote=F, sep="\t")

  rm(ego_mf)
  
  
			  
  # 富集啦
  kk <- enrichKEGG(gene = gene,
                   organism = "hsa",
                   pvalueCutoff = 0.01,
                   minGSSize = 1)
				   

  pdf(paste(richPath, "/Enrichment_KEGG.pdf", sep=""), height=4, width=10)
  print( barplot(kk, showCategory=10, title="Enrichment_KEGG") )
  dev.off()

  pdf(paste(richPath, "/Enrichment_KEGG_dot.pdf", sep=""), height=4, width=10)
  print( dotplot(kk, title="Enrichment_KEGG_dot") )
  dev.off()

  write.table(as.data.frame(kk@result), 
              file=paste(richPath, "/Enrichment_KEGG.txt", sep=""),
              row.names=F, quote=F, sep="\t")



  ### 收尾
  cat("\nCongratulations\n")
  endtime<- print(Sys.time())
  cat(paste("EndTime:",endtime,sep=" "),"\n")
  print(difftime(endtime,sttime))
  rm( list=ls() )
  
}



