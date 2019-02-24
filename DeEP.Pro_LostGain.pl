#!/usr/bin/perl
use warnings;
use diagnostics;
use strict;
use YAML::Tiny;
use Getopt::Long;


my $cfg_file = '';
my @files=();
my @files_reg=();
my @Myfiles=();
my @Myfiles_reg=();
my @Myfiles_gc=();
my $gc_gem='';



GetOptions(
    'cfg|c=s'   => \$cfg_file,
	);

 if(!-s$cfg_file){
  warn "No configuration file. Can not continue.\n";
 &usage();
 exit(1);
 }

sub usage{
 print STDERR <<EOF;
 
*******************************************************************************************
								Gained & Lost Genes											  
This script is used to identify the genomic regions correlated with traits of interest.
It calculates the parameters for the markov chain model using Maximum likehood, and generates 
a table with the gaines and lost genes, according with this model. 
*******************************************************************************************

 USAGE: DeEP_HR.pl -c <cfg_file.yml>
 
 Options 
	--cfg <cfg_file.yml>
	   -c <cfg_file.yml>      
					path of the configuration file to use. 
		                	For example /home/me/cfg-file.yml
EOF

}

  # Create a YAML file
    my $yaml = YAML::Tiny->new();


  # Open the config	
  $yaml= YAML::Tiny->read($cfg_file) or die "Couldn't read file: Error($!) : Errstr =", YAML::Tiny->errstr;
	
    # Reading properties 
   my $dir_out= $yaml->[0]->{A_General}->{dir_out} or die "Fill dir_out in the configuration file ($cfg_file)";
   #my $Rscript_path =  $yaml->[0]->{N_Rscript}->{Rscript_path} or die "Fill Rscript_path in the configuration file ($cfg_file)";
   my $clonaltree = $yaml->[0]->{K_ClonalTree_output}->{final_tree} or die "Fill final_tree in the configuration file ($cfg_file)";
   #my $tree =  $yaml->[0]->{A_General}->{Phylogeny_tree} or die "Fill Phylogeny_tree in the configuration file ($cfg_file)";
   #my $traits= $yaml->[0]->{A_General}->{Traits} or die "Fill Traits in the configuration file ($cfg_file)";
   my $Rcomb= $yaml->[0]->{L_Result}->{Summary_File} or die "Fill Summary_File in the configuration file ($cfg_file)";
    my $xmfa = $yaml->[0]->{C_MAUVE_output}->{xmfa} or die "Fill xmfa  in the configuration file ($cfg_file)";
    
open IN, $xmfa;
  my @xmfa= <IN>;
   close IN;

my @namesxmfa=();
my $allnames=();
my @sequences_num=();
my %hash_tree;
#Sequence1File

foreach my $xmfaline(@xmfa){
  if($xmfaline=~m/#Sequence(\d+)File.+\/(.+)$/){
  my $name=$2;
  my $sequence_num=$1;
  push (@sequences_num,$sequence_num);
  push(@namesxmfa,$name);
  $allnames=join(',',@namesxmfa);
  $hash_tree{$sequence_num}=$name;
 }
}


my $l_names=scalar(@namesxmfa);
my $max='';
my @max=();

   open IN, $clonaltree;
  my @lines_clonaltree= <IN>;
chomp(@lines_clonaltree);
 close IN;
my $tree_concensus='';

  foreach my $line_clonaltree(@lines_clonaltree){
print "Tree: $line_clonaltree\n";
   while($line_clonaltree=~m/(\d+)\:/g){
   my $seq=$1;
my $name_seq='';
push(@max,$seq);
if(exists $hash_tree{$seq}){
    $name_seq=$hash_tree{$seq}."\:";
 $line_clonaltree=~s/$seq\:/$name_seq/g;}   
}
$tree_concensus=$line_clonaltree;
}
print "Tree with species names\: $tree_concensus\n";

@max=sort(@max);
$max=$max[-1];
#print "este es el max $max\n";
my %hash_nodes;

my $count=$max+1;
 my $count1=0;
my @part=();

my $s=$tree_concensus;
my @code_shape=split('',$s);

foreach my $code (@code_shape){
if($code=~m/\(/){
$count--;
my $letter="\*$count\*";
$count1=$count-1;
push (@part,$letter);
 }
 elsif($code=~m/\)/){
$count1++; 
my $letter1="\*$count1\*";
push (@part,$letter1);
 }   
else{push (@part,$code);}
  }
 
print $count;
my $tree_code =join('',@part);
print $tree_code."\n";
my $species='';

for (my $i=$max;$i>$l_names;$i--){
if($tree_code=~m/\*$i\*(.+)\*$i\*/){
$species=$1;
$species=~s/\*\d+\*//g;
$species=~s/\:\d+\.\d+//g;

} 
print "the node $i have this species $species\n";
$hash_tree{$i}=$species;
}

my $node_num='';
my $block_Rcomb='';

open IN, $Rcomb;
 my @Rcomb= <IN>;
   close IN;
foreach my $Rcomb_line(@Rcomb){
if($Rcomb_line=~m/(\d+)\t\d+\t(\d+)\t/g){
$block_Rcomb=$1;
$node_num=$2;
print "$block_Rcomb\t$node_num\t$hash_tree{$node_num}\n";
}
}
