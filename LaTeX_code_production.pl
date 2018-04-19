#!/usr/bin/env perl

############################################################
#    2017.12.26    Kunbang Liu
#
#  Function: 创建 .tex 后缀的文件，
#            再调用 latex生成结果报告
#
############################################################


use Cwd;
use strict;
use warnings;
use Getopt::Long;
use File::Copy::Recursive qw(fcopy);
use FindBin qw( $Bin );
use lib "$Bin";


my $PWD = getcwd;
my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
print "The results report...$hour:$min:$sec\n\n";


### 获得输入参数
my ($outdir, $barPath, $QCPath, 
    $sampleInfoFile, $ClusterPath, $gseaPath, $richPath);

GetOptions(
  "od:s"          => \$outdir,
  "plotPath:s"     => \$barPath,
  "QCPath:s"      => \$QCPath,
  "sample:s"      => \$sampleInfoFile,
  "ClusterPath:s" => \$ClusterPath,
  "GSEA:s"        => \$gseaPath,
  "rich:s"        => \$richPath,
);


#----------- 获得所需数据和文件路径 ---------------------------#
my $allPlot    = "$barPath/merge_vol_sca_Plots.pdf";
my $boxPlot    = "$QCPath/Boxplot.pdf";
my $deg_summ   = "$barPath/deg_summary.txt";
my $RNAdegPDF  = "$QCPath/RNA_degradation_plot.pdf";
my $DensityPDF = "$QCPath/DensityPlot.pdf";
my $ClusterPDF = "$ClusterPath/Cluster_heatmapPlot.pdf";
my $richbar    = "$richPath/Enrichment_KEGG.pdf";
my $richdot    = "$richPath/Enrichment_KEGG_dot.pdf";


open TXT, "$gseaPath/Summary_KEGG.txt" or warn "KEGG.txt:$!";
my @info  = <TXT>;
my $KEGG_info = join "", @info;
close TXT;

open TXT, "$sampleInfoFile" or warn "sampleInfo.txt:$!";
my @lines = <TXT>;
close TXT;


#----------- 将 LaTeX代码写入文件 ---------------------------#
my $essay = &Essay($#lines, $boxPlot, $RNAdegPDF, $DensityPDF, $allPlot,
             $ClusterPDF, $deg_summ, $richbar, $richdot, $KEGG_info);

open HC ,"> $outdir/Result_Report.tex" or die "Result_Report.tex: $!";
print HC $essay,"\n";
close HC;


#----------- 项目 PDF 时所需的字体文件拷贝 -----------------#
fcopy("$Bin/msyhl.latexfont.ttc","$PWD") unless -f "$PWD/msyhl.latexfont.ttc";


######## 运行时间超过1分钟则自动退出,防止文件缺失时自动等待重新输入文件的意外
eval {
  local $SIG{ALRM} = sub { die "timeout\n" };
  alarm(60);  # 设置超时时间
	

#------------ 生成pdf -----------------------------------------#
  system ("xelatex -output-directory=$outdir $outdir/Result_Report.tex > $outdir/Result_Report.out 2>&1");

  alarm(0);
};


if($@)
{
  die $@ unless $@ eq "timeout\n";
  die "Timed out!\n";     #超时，终止程序
}


print "\nConguatulation\n\n";

($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime();
print "Final result arrangement END:\t$hour:$min:$sec\n";



###################### 子函数 #################################
#---------------- 整体 latex 代码模块生产  -------------------#
sub Essay{
    my ($samplenum, $boxPDF, $RNAPDF, $denPDF, $allPDF, $ClusterPDF,
        $degsumm, $Richbar, $Richdot, $KEGG) = @_;
   
   # 开头是"\"的行均为 LaTeX 代码
    my $tex_head = '\documentclass[a4paper,11pt]{article}
\usepackage{geometry}                        %设置页面式样宏包
\usepackage{fancyhdr}                        %设置页眉
\usepackage{Sweave}
\usepackage{ctex}                            %中文支持
\geometry{left=3.0cm,right=3.0cm,top=3.5cm,bottom=3.5cm} %设置页边距
\setCJKmainfont{msyhl.latexfont.ttc}         %设置字体
\setmainfont{msyhl.latexfont.ttc}            %设置字体
\newsavebox{\headpic}
\sbox{\headpic}{\includegraphics[height=2cm,width=5cm]{/lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code/logo.png}}
\title{Affy 芯片分析结果报告}

\begin{document}

  \headheight 36pt                           %设置页眉高度
  \pagestyle{fancy}        
  \lhead{\usebox{\headpic}}                  %页眉图片
  \rhead{\textbf{朴诚、奋勉、求实、创新}}    %页眉文字
  \linespread{1.3}\selectfont                %页眉线的粗细
  \renewcommand\tablename{表}                %将表名 “Table”换成中文“表”
  \maketitle                                 %打印 title

';

    my $tex_end = '

  \clearpage                                 %防止浮游表格和图片位置错乱
\end{document}
';
    
#  \noindent{尊敬的 XXX 老师，}\\\
    my $part1 = '

  \noindent{本次检测的芯片名称是 Affymetrix GeneChip® Human Genome U133 Plus 2.0 Array（简称：hgu133plus2），该芯片包含了47,000个转录本，代表了38,500个明晰的人类基因。数据库来源于GeneBank、dbEST、RefSeq、UniGene database、Washington University EST trace repository、NCBI human genome assembly。}\\\

  \noindent{我们已经完成 '.$samplenum.' 个样品的 hgu133plus2 芯片的检测。}\\\

  \noindent{结果总结如下：}\\\


  \vspace{4 ex}                              %设置4个字母“e”宽度的垂直空白


  \noindent{一、分析内容}\\\
  \indent{a) 芯片质控}\\\
  \indent{b) 聚类分析}\\\
  \indent{c) 差异分析}\\\
  \indent{d) 差异基因的GO和Pathway分析}\\\
  \indent{e) 基因集富集分析(GSEA)}


  \newpage                                   %换页

  \indent{  }\\\
  \noindent{二、基本结果}\\\
  \indent{a) 芯片质控}\\\
';


    my $box = '  \indent{    箱线图（见图1），不同的实验组以不同的颜色标记，同组的生物学重复样本以相同的颜色标记，展现各个样本归一化前后总体探针的表达水平。}\\\
  \begin{figure}[!hbp]                       %添加图片
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.9\textwidth]{'.$boxPDF.'}
  \caption{箱线图}
  \label{图1}
  \end{figure}
  
  \vspace{3 ex}
';


    my $deg = '
  \indent{    RNA降解曲线（见图2），通常，合格的芯片，其降解曲线呈平滑递增的形态。}\\\
  \begin{figure}[!hbpt]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.7\textwidth]{'.$RNAPDF.'}
  \caption{RNA降解图}
  \label{图2}
  \end{figure}

  \newpage  
';


    my $density = '
  \indent{密度曲线（见图3），查看所有芯片中是否有异常过高的的密度值，如果出现双峰，则表明芯片可能存在“人为空白（spatical artifical）”的质量缺陷。}\\\
  \begin{figure}[!hbp]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.9\textwidth]{'.$denPDF.'}
  \caption{密度曲线：\footnotesize{展现归一化前后探针密度。}}
  \label{图3}
  \end{figure}

  \newpage

';


    my $part2 ='
  \indent{b) 差异基因}\\\
  \indent{两组样本中基因表达差异倍数(FC)>=2或<=0.5且P值<=0.05时,定义为显著差异基因。火山图(见图4)和散点图(见图4)中红色点表示上调基因,绿色点表示下调基因。}
';


    my $degtable = &DiffTableCode($degsumm);
    
    my $bar = '
  \begin{figure}[!hbp]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=1\textwidth]{'.$allPDF.'}
  \caption{火山图（左）和散点图（右）}
  \label{图4}
  \end{figure}

';


    my $desc_Clu = '
  \indent{聚类图(见图5)每一行为一个基因，每一列为一个样本，可直观的反映样本聚类是否符合预期分组。}
';


    my $Clu = '
  \begin{figure}[!hbp]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.60\textwidth]{'.$ClusterPDF.'}
  \caption{聚类图}
  \label{图5}
  \end{figure}

  \newpage
';


    my $part3 ='
  \indent{c) 差异富集分析}\\\
  \indent{  挑选差异富集结果中P值最小（最可信）的前10项富集结果绘制条形图(见图6)；选取富集数目最多的前10项结果绘制气泡图(见图7)。}
';


   my $Bar = '
  \begin{figure}[!hbp]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.75\textwidth]{'.$Richbar.'}
  \caption{富集条形图：\footnotesize{按照P值大小排列。}}
  \label{图6}
  \end{figure}

';


   my $dot = '
  \begin{figure}[!hbp]
  \setlength{\abovecaptionskip}{-0.cm}
  \centering
    \includegraphics[width=0.75\textwidth]{'.$Richdot.'}
  \caption{富集气泡图：\footnotesize{按照富集数目大小排列。}}
  \label{图7}
  \end{figure}

';


    my $part4 ='
  \vspace{3 ex}
  \indent{d) GESA——KEGG中部分结果}\\\
';


  $tex_head.$part1.$box.$deg.$density.$part2.$degtable.
        $bar.$desc_Clu.$Clu.$part3.$Bar.$dot.$part4.$KEGG.$tex_end;
}



#------------------ 差异基因数 表格 latex 代码生产 ------------------#
sub DiffTableCode{

    my ($bartxt) = @_;
    
    # LaTeX 表头代码
    my $tablecode = '  \begin{table}[ht]                  %表格
    \centering                            %居中
    \footnotesize                         %字体小号
    \setlength{\belowcaptionskip}{5pt}    %表名与图片的间距
    \caption{差异基因检出数量}            %表名
    \begin{tabular}{|c|c|c|}
    \hline
    分组 & 上调基因 & 下调基因\\\
    \hline
';

    open DIFF,"< $bartxt" or warn "$bartxt: $!";
    
    # LaTeX表中数据代码
    while(<DIFF>){
        next unless /\d/;
        chomp;
        s/_/\\_/g;
        my ($group, $num1, $num2) = split(/\t/);
        $tablecode .= "    $group \& $num1 \& $num2\\\\\n    \\hline\n";
    }
    
    # LaTeX 表尾代码
    $tablecode .= "    \\end\{tabular\}\n  \\end\{table\}\n";
    
    close DIFF;
    
    #返回表代码和数据行数
    $tablecode;
}



1;
