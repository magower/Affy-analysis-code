#!/usr/bin/env perl

########################################################################
#   Author: Kunbang Liu
#
#     Date: 2017.12.28    
#
# Function: Analysis of affy chips, 
#           including differential analysis,cluster analysis,,GSEA, etc
#
#     Tool: 1.OS:CentOS Linux release 7.3.1611 (Core)
#           2.R version 3.4.1
#           3.Perl 5 version 16 subversion 3 (v5.16.3)
#           4.XeTeX 3.14159265-2.6-0.99998 (TeX Live 2017)
#           5.JAVA version "1.8.0_151"
#
#
#########################################################################



use Cwd;
use strict;
use threads;
use warnings;
use FindBin qw( $Bin );
use lib "$Bin";


###### Set Path and Parameters
my $Rpath  = "/usr/bin/R";
my $CodePath = '/lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code';
my $ProjectPath = my $RunPath = getcwd;
$ProjectPath =~ s/\/run//;
my $CelDataPath = $ProjectPath."/CelData";
my $MidfilePath  = $ProjectPath."/Midfile";
my $logPath      = $RunPath."/log";
my $sampleInfoFile = $RunPath."/sampleInfo.txt";
my $projectlogFile = "$logPath/project.log";


###### Create Dir
mkdir $MidfilePath unless -d $MidfilePath;
mkdir $logPath     unless -d $logPath;


###### Print project start time
open LOG, "> $projectlogFile" or warn "log:$!";

my $stime = localtime();
print "\n"."#-" x 20,"#\nProject start: $stime\n\n";
print LOG "\n"."#-" x 20,"#\nProject start: $stime\n\n";


###### Analysis CEL Data: Preprocess
if(1){

  print "#-" x 20,"#\nAffy Preprocess analysis...\n";
  print LOG "#-" x 20,"#\nAffy Preprocess analysis...\n";
  

  # Create Dir
  mkdir "$MidfilePath/QC"         unless -d "$MidfilePath/QC";
  mkdir "$MidfilePath/Preprocess" unless -d "$MidfilePath/Preprocess";
  mkdir "$MidfilePath/Cluster"    unless -d "$MidfilePath/Cluster";
  

  # print Input parameters
  print LOG "\$sampleInfoFile: $sampleInfoFile\n";
  print LOG "\$CelDataPath   : $CelDataPath\n";
  print LOG "\$MidfilePath   : $MidfilePath\n";


  # Use R 
  open RPROG, "|$Rpath --vanilla --slave" or warn "$Rpath:$!";
  print RPROG <<"CODE";
  zz <- file("$logPath/PreprocessAffyExp.Rout", open="wt")
  sink(zz)
  sink(zz, type="message")
  source("$CodePath/PreprocessAffyExp.R")
  PreprocessAffyExp( "$sampleInfoFile", "$CelDataPath", 
                     "$MidfilePath" )
  q()
CODE
  close RPROG;
  
  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "Affy Preprocess analysis END:\t$hour:$min:$sec\n\n";
  print LOG "Affy Preprocess analysis END:\t$hour:$min:$sec\n\n";
}



###### Clustering of differentially genes
if(1){

  print "#-" x 20,"#\nCluster analysis...\n";
  print LOG "#-" x 20,"#\nCluster analysis...\n";
  mkdir "$MidfilePath/Cluster" unless -d "$MidfilePath/Cluster";


  # Set Input parameters
  my $ClusterPath = "$MidfilePath/Cluster";
  print LOG "\$ClusterPath   : $ClusterPath\n";


  # Use R 
  open RPROG, "|$Rpath --vanilla --slave" or warn "$Rpath:$!";
  print RPROG <<"CODE";
  zz <- file("$logPath/Cluster.Rout", open="wt")
  sink(zz)
  sink(zz, type="message")
  source("$CodePath/Cluster.R")
  Cluster( "$ClusterPath" )
  q()
CODE
  close RPROG;


  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "Cluster analysis END:\t$hour:$min:$sec\n\n";
  print LOG "Cluster analysis END:\t$hour:$min:$sec\n\n";

}




###### Drawing histogram for differentially genes
if(1){

  print "#-" x 20,"#\nDraw Plot...\n";
  print LOG "#-" x 20,"#\nDraw Plot...\n";
  mkdir "$MidfilePath/degPlot" unless -d "$MidfilePath/degPlot";
  
  
  # Set Input parameters
  my $PlotPath = "$MidfilePath/degPlot";
  my $inputdegFile = "$MidfilePath/Preprocess/diff_expressed_genes.txt";
  my $inputallFile = "$MidfilePath/Preprocess/all_expression_info.txt";
  print LOG "\$PlotPath       : $PlotPath\n";
  print LOG "\$inputdegFile   : $inputdegFile\n";
  print LOG "\$inputallFile   : $inputallFile\n";
  print LOG "\$sampleInfoFile : $sampleInfoFile\n";
  
  
  # Use R
  open RPROG, "|$Rpath --vanilla --slave" or warn "$Rpath:$!";
  print RPROG <<"CODE";
  zz <- file("$logPath/degPlot.Rout", open="wt")
  sink(zz)
  sink(zz, type="message")
  source("$CodePath/degPlot.R")
  Plot( "$inputdegFile", "$inputallFile", "$PlotPath", "$sampleInfoFile" )
  q()
CODE
  close RPROG;


  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "Draw Plot END:\t$hour:$min:$sec\n\n";
  print LOG "Draw Plot END:\t$hour:$min:$sec\n\n";

}



###### Enrichment analysis of differential genes
if(1){

  print "#-" x 20,"#\nEnrichment analysis...\n";
  print LOG "#-" x 20,"#\nEnrichment analysis...\n";
  mkdir "$MidfilePath/Enrichment" unless -d "$MidfilePath/Enrichment";

  
  # Set Input parameters
  my $deg_File = "$MidfilePath/Preprocess/diff_expressed_genes.txt";
  my $richPath = "$MidfilePath/Enrichment";
  print LOG "\$deg_File : $deg_File\n";
  print LOG "\$richPath : $richPath\n";
  
  
  # Use R
  open RPROG, "|$Rpath --vanilla --slave" or warn "$Rpath:$!";
  print RPROG <<"CODE";
  zz <- file("$logPath/enrichment.Rout", open="wt")
  sink(zz)
  sink(zz, type="message")
  source("$CodePath/enrichment.R")
  Enrichment( "$deg_File", "$richPath" )
  q()
CODE
  close RPROG;


  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "Enrichment analysis END:\t$hour:$min:$sec\n\n";
  print LOG "Enrichment analysis END:\t$hour:$min:$sec\n\n";

}



###### Gene Set Enrichment Analysis
if(1){

  print "#-" x 20,"#\nGSEA...\n";
  print LOG "#-" x 20,"#\nGSEA...\n";
  mkdir "$MidfilePath/GSEA" unless -d "$MidfilePath/GSEA";


  # Set Input parameters
  my $KEGG_info;
  my $all_File = "$MidfilePath/Preprocess/all_expression_info.txt";
  my $gseaPath = "$MidfilePath/GSEA";
  my $gctFile  = "$MidfilePath/GSEA/data.gct";
  print LOG "\$all_File       : $all_File\n";
  print LOG "\$gseaPath       : $gseaPath\n";
  print LOG "\$sampleInfoFile : $sampleInfoFile\n";


  # Use R
  open RPROG, "|$Rpath --vanilla --slave" or warn "$Rpath:$!";
  print RPROG <<"CODE";
  zz <- file("$logPath/GSEA.Rout", open="wt")
  sink(zz)
  sink(zz, type="message")
  source("$CodePath/GSEA.R")
  Probe2Symbol( "$all_File", "$gseaPath", "$sampleInfoFile" )
  q()
CODE
  close RPROG;


  # 整理文件数据
  open GCT, "$gctFile" or warn "data.gct:$!";
  my @gct = <GCT>;
  my $sample_num = ($gct[0] =~ tr/\t/\t/) - 1;
  close GCT;

  open OUT, "> $gctFile" or warn "data.gct:$!";
  print OUT "\#1.2\n$#gct\t$sample_num\n";
  print OUT @gct;
  close OUT;
  

  # 获得比较组
  my $match;
  open TXT, "$gseaPath/group.cls" or warn "group.cls:$!";
  
  foreach ( <TXT> ){
    next unless /^\#/;
    chomp;
    s/^\# //;
    s/ /_versus_/;
    $match = $_;
  }
  close TXT;


  # java命令
  my @java = (
    'java -Djava.util.prefs.userRoot=/mnt/icfs/work/bioinfo/test/kunbang/ -Xmx1024m -cp /lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code/gsea2-2.2.4.jar xtools.gsea.Gsea -res ./data.gct -cls ./group.cls#C_versus_N -out ./ -rpt_label c5.GO_BP -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -gui false',
    'java -Djava.util.prefs.userRoot=/mnt/icfs/work/bioinfo/test/kunbang/ -Xmx1024m -cp /lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code/gsea2-2.2.4.jar xtools.gsea.Gsea -res ./data.gct -cls ./group.cls#C_versus_N -out ./ -rpt_label c5.GO_CC -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.cc.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -gui false',
    'java -Djava.util.prefs.userRoot=/mnt/icfs/work/bioinfo/test/kunbang/ -Xmx1024m -cp /lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code/gsea2-2.2.4.jar xtools.gsea.Gsea -res ./data.gct -cls ./group.cls#C_versus_N -out ./ -rpt_label c5.GO_MF -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.mf.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -gui false',
    'java -Djava.util.prefs.userRoot=/mnt/icfs/work/bioinfo/test/kunbang/ -Xmx1024m -cp /lustre/analysis/microarray/bioinfo-RNA-pro/6.test/kunbang/grad_article/chip_analysis/Analysis_code/gsea2-2.2.4.jar xtools.gsea.Gsea -res ./data.gct -cls ./group.cls#C_versus_N -out ./ -rpt_label c2.KEGG -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.2.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -gui false'
  );

  
  chdir( $gseaPath ) or warn $!;
  
  # 采用多线程进行 GSEA
  # 创建线程
  foreach(@java){
    s/\#\w+? /\#$match /;
    my $thr = threads->new(\&gseawork, $_)
  }


  # 等待线程结束，并返回过程信息
  open INFO, "> $gseaPath/java.log" or die "java.log:$!";
  foreach my $thread(threads->list()){
    my $result = $thread->join();  #等待线程返回值
    print INFO $result;
    print INFO "\#\n" x 10;
  }
  close INFO;

  
  # 线程子函数
  sub gseawork{
    my ($command) = @_;
    my $info = `$command`;
  }

  
  # 记录 GSEA中 KEGG结果概要信息，存放于 Summary_KEGG.txt文件
  my $GSEAkeggPath = glob "c2.KEGG*";
  
  open HTML, "$GSEAkeggPath/index.html" or warn "html:$!";
  open TXT,  "> Summary_KEGG.txt" or warn "Summary_KEGG:$!";
  
  foreach my $line( <HTML> ){

    next unless $line =~ /Enrichment/;
    $line =~ s/(.?)Dataset details.+/$1/;
    $line =~ s$<\/?i>$$g;
    $line =~ s$</?h\d>$\t$g;
    $line =~ s$<\/?[^a><b]+?>$\n$g;
    $line =~ s$<\/?[^><]+?>$$g;
    $line =~ s$\n[SDG].+$$g;
    $line =~ s/^\n+//g;
    $line =~ s/\n+/\\\\\n/g;
    $line =~ s/\n([^\t\n\\]+)/\n\\indent{$1}/g;
    $line =~ s/\t([^\t]+)\t\\\\/\\\\\n\\noindent{$1}\\\\/g;
    $line =~ s/(%)/\\$1/g;
	$line =~ s/\\\\\n\t\\\\\n//;
   
    $KEGG_info .= $line;
  }
  print TXT $KEGG_info;

  close TXT;
  close HTML;
  
  
  # 将路径切回 run 目录
  chdir( $RunPath ) or warn $!;

  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "GSEA END:\t$hour:$min:$sec\n\n";
  print LOG "GSEA END:\t$hour:$min:$sec\n\n";

}



###### Draw results report PDF
if(1){

  print "#-" x 20,"#\nThe results report...\n";
  print LOG "#-" x 20,"#\nThe results report...\n";
  mkdir "$MidfilePath/Report" unless -d "$MidfilePath/Report";
  
  
  # Set Input parameters
  my @options = (
  "-od $MidfilePath/Report",
  "-plotPath $MidfilePath/degPlot",
  "-QCPath $MidfilePath/QC",
  "-sample $sampleInfoFile",
  "-ClusterPath $MidfilePath/Cluster",
  "-GSEA $MidfilePath/GSEA",
  "-rich $MidfilePath/Enrichment",
  );
  my $option = join " ", @options;
  
  print LOG "\$option : $option\n";
  system "perl $CodePath/LaTeX_code_production.pl $option > $logPath/latex.log 2>&1";


  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "The results report END:\t$hour:$min:$sec\n\n";
  print LOG "The results report END:\t$hour:$min:$sec\n\n";

}



###### Finishing the result
if(1){

  print "#-" x 20,"#\nFinal result arrangement...\n";
  print LOG "#-" x 20,"#\nFinal result arrangement...\n";
  mkdir "$ProjectPath/Result" unless -d "$ProjectPath/Result";
  
  
  # Set Input parameters
  my @options = (
    "-Re $MidfilePath/Report",
    "-Pl $MidfilePath/degPlot",
    "-QC $MidfilePath/QC",
    "-Cl $MidfilePath/Cluster",
    "-GS $MidfilePath/GSEA",
    "-En $MidfilePath/Enrichment",
    "-Pr $MidfilePath/Preprocess",
    "-Od $ProjectPath/Result",
    "-Co $CodePath",
  );
  my $option = join " ", @options;
  
  print LOG "\$option : $option\n";
  system "perl $CodePath/result2disk.pl $option > $logPath/result2disk.out 2>&1";


  my ( $sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst ) = localtime();
  print "Final result arrangement END:\t$hour:$min:$sec\n\n";
  print LOG "Final result arrangement END:\t$hour:$min:$sec\n\n";

}



# Print project end time
my $etime = localtime();
print "Project end: $etime\n","#-" x 20,"#\n\n";
print LOG "Project end: $etime\n","#-" x 20,"#\n\n";

close LOG;


1;
