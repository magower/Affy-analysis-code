#!/usr/bin/env perl

############################################################
#    2017.12.28    Kunbang Liu
#
#  Function: 整理 Result文件夹
#            
#
############################################################


use Cwd;
use strict;
use threads;
use warnings;
use Getopt::Long;
use File::Copy::Recursive qw(fcopy dircopy);


my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
print "Final result arrangement...$hour:$min:$sec\n\n";



### 获得输入参数
my ($Redir, $Pldir, $QCdir, $Cldir, $GSdir, $Endir, $Prdir, $Oddir, $Codir);

GetOptions(
  "Re:s"  => \$Redir,
  "Pl:s"  => \$Pldir,
  "QC:s"  => \$QCdir,
  "Cl:s"  => \$Cldir,
  "GS:s"  => \$GSdir,
  "En:s"  => \$Endir,
  "Pr:s"  => \$Prdir,
  "Od:s"  => \$Oddir,
  "Co:s"  => \$Codir,
);



### 创建文件目录
mkdir( "$Oddir/A.芯片信息" ) or warn $! unless -d "$Oddir/A.芯片信息";
mkdir( "$Oddir/C.差异结果" ) or warn $! unless -d "$Oddir/C.差异结果";
mkdir( "$Oddir/C.差异结果/Enrichment" ) or warn $! 
                             unless -d "$Oddir/C.差异结果/Enrichment";
mkdir( "$Oddir/D.GSEA" ) or warn $! unless -d "$Oddir/D.GSEA";
mkdir( "$Oddir/D.GSEA/KEGG" ) or warn $! unless -d "$Oddir/D.GSEA/KEGG";
mkdir( "$Oddir/D.GSEA/GO_CC" ) or warn $! unless -d "$Oddir/D.GSEA/GO_CC";
mkdir( "$Oddir/D.GSEA/GO_BP" ) or warn $! unless -d "$Oddir/D.GSEA/GO_BP";
mkdir( "$Oddir/D.GSEA/GO_MF" ) or warn $! unless -d "$Oddir/D.GSEA/GO_MF";



### 文件和目录复制
# 1.芯片信息
# fcopy( "$Prdir/raw_chip_exp.txt", "$Oddir/A.芯片信息/Raw_Exp.txt" );
fcopy( "$Prdir/all_expression_info.txt", "$Oddir/A.芯片信息/All_Exp_Info.txt" );

# 2.质控报告
dircopy( $QCdir , "$Oddir/B.质控报告" );

# 3.差异结果
fcopy( "$Prdir/diff_expressed_genes.txt", "$Oddir/C.差异结果/Deg_Exp_info.txt" );
fcopy( "$Cldir/Cluster_heatmapPlot.pdf", "$Oddir/C.差异结果/Cluster.pdf" );
fcopy( "$Pldir/scaPlot.pdf", "$Oddir/C.差异结果/scatter_Plot.pdf" );
fcopy( "$Pldir/volPlot.pdf", "$Oddir/C.差异结果/volcano_Plot.pdf" );
dircopy( $Endir , "$Oddir/C.差异结果/Enrichment" );

# 4.GSEA
my @dir = ( "KEGG", "GO_CC", "GO_BP", "GO_MF" );
my @dir2 = glob "$GSdir/c*";
foreach my $dir( @dir2 ){

  my @files = glob "$dir/*.png $dir/*.html $dir/*.xls";
  print "$dir\t$#files\n";
  
  # 多线程文件复制
  foreach( @files ){
    my $thr1 = threads->new( \&FC, $_, "$Oddir/D.GSEA/KEGG" )
	                                   if $dir =~ /c2.KEGG.Gsea/;
    my $thr2 = threads->new( \&FC, $_, "$Oddir/D.GSEA/GO_CC" )
	                                   if $dir =~ /c5.GO_CC.Gsea/;
    my $thr3 = threads->new( \&FC, $_, "$Oddir/D.GSEA/GO_BP" )
	                                   if $dir =~ /c5.GO_BP.Gsea/;
    my $thr4 = threads->new( \&FC, $_, "$Oddir/D.GSEA/GO_MF" )
	                                   if $dir =~ /c5.GO_MF.Gsea/;
  }

  # 等待线程返回
  foreach my $thread(threads->list()){
    $thread->join();
  }

}

# 线程子函数
sub FC {
  my ($file, $out) = @_;
  fcopy($file, $out) or warn $!;
}
  
# 5.Report  
fcopy( "$Redir/Result_Report.pdf", "$Oddir/E.Report.pdf" ) or warn $!;

# 6.readme
fcopy( "$Codir/readme", "$Oddir/readme" ) or warn $!;


print "\nConguatulation\n\n";

($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime();
print "Final result arrangement END:\t$hour:$min:$sec\n";





1;
