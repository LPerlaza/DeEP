#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use YAML::Tiny;
use Getopt::Long;

my $cfg_file = '';
my @files=();
my @Myfiles=();
my @Myfiles_ready=();
my $Myfiles='';


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
                                    DeEP.Pro: MultiAlignment											  
This script is used to make a multiple genome alignment, strip out and split the 
collinear blocks on the alignment. This script runs progressiveMAUVE, stripSubsetLCBs, 
and blocksplit.pl. Additionally, it generates a binary matrix with the present and absent 
of conserved regions on the genomes.
*******************************************************************************************

 USAGE: DeEP_MultiAlign.pl -c <cfg_file.yml>
 
 Options 
	--cfg <cfg_file.yml>
	   -c <cfg_file.yml>     path of the configuration file to use. For example /home/me/cfg-file.yml.
EOF
}

# Create a YAML file
my $yaml = YAML::Tiny->new();


# Open the config	
$yaml= YAML::Tiny->read($cfg_file) or die "Couldn't read file: Error($!) : Errstr =", YAML::Tiny->errstr;
	
# Reading properties
my $dir_files = $yaml->[0]->{A_General}->{dir_files} or die "Fill dir_files in the configuration file ($cfg_file)";
my $name = $yaml->[0]->{A_General}->{job_name} or die "Fill job_name in the configuration file ($cfg_file)";
my $path_mauve = $yaml ->[0]->{B_MAUVE}->{path_mauve} or die "Fill path_mauve in the configuration file ($cfg_file)"; 
my $length = $yaml ->[0]->{D_ssLCBs}->{length} or die "Fill length in the configuration file ($cfg_file)";
my $path_ssLCBs = $yaml ->[0]->{D_ssLCBs}->{path_ssLCBs} or die "Fill path_ssLCBs in the configuration file ($cfg_file)";
my $path_blocksplit = $yaml ->[0]->{F_Blocksplit}->{path_blocksplit} or die "Fill path_blocksplit in the configuration file ($cfg_file)";
my $cut=$yaml->[0]->{D_ssLCBs}->{length} or die "Fill length in the configuration file ($cfg_file)";


my $name_xmfa=$name.".xmfa";
my $name_backbone=$name.".backbone";
my $name_tree=$name.".tree";

#creating the output folder
system("mkdir $dir_files\_output");
print "The output folder: $dir_files\_output has been created\n";
$yaml->[0]->{A_General}->{dir_out} = "$dir_files\_output";

my $dir_out= $yaml->[0]->{A_General}->{dir_out} or die "Fill dir_out in the configuration file ($cfg_file)";


#########################################################################################
#progressiveMAUVE
#########################################################################################

#opening the dir_file and processing de files
opendir(DIR, $dir_files) or die "Error in opening dir $dir_files\n";

while (my $file= readdir(DIR)){
 next if $file =~/\.\.?$/;
 next if $file !~/(.fna|.fa|.fasta)$/;
 
 my $fullpath_file=$dir_files."/".$file;
 push(@Myfiles,$fullpath_file); 
}

if(@Myfiles<=0){die "No files to do MAUVE: Error($!)"}
else{
  print "\tRunning progressiveMauve with files: @Myfiles\n";
  my $cmd_Mauve= $path_mauve."/progressiveMauve  --output=".$dir_out."/".$name_xmfa." --backbone-output=".$dir_out."/".$name_backbone." --output-guide-tree=".$dir_out."/".$name_tree." @Myfiles";
  print "\tYou are asking to run this command: $cmd_Mauve\n";
  system($cmd_Mauve);
}

closedir(DIR);

#adding the results to the configuration file
$yaml->[0]->{C_MAUVE_output}->{xmfa} = $dir_out."/".$name_xmfa;
print "$dir_out\/$name_xmfa was created\n";

$yaml->[0]->{C_MAUVE_output}->{tree} = $dir_out."/".$name_tree;
print "$dir_out\/$name_tree was created\n";

$yaml->[0]->{C_MAUVE_output}->{backbone} = $dir_out."/".$name_backbone;
print "$dir_out\/$name_backbone was created\n";

$yaml->write( $cfg_file );

#Adding the names of the genomes to the configuration file
my $xmfa = $yaml->[0]->{C_MAUVE_output}->{xmfa} or die "Fill xmfa  in the configuration file ($cfg_file)";

open IN, $xmfa;
  my @xmfa= <IN>;
   close IN;

my @namesxmfa=();
my $allnames=();

foreach my $xmfaline(@xmfa){
  if($xmfaline=~m/Sequence.+\/(.+)$/){
  my $name=$1;
  push(@namesxmfa,$name);
  $allnames=join(',',@namesxmfa);
 }
}

$yaml->[0]->{A_General}->{num_genomes}=scalar(@namesxmfa);
$yaml->[0]->{A_General}->{names}= $allnames;
$yaml->write( $cfg_file ); 

###########################################################################################
#stripSubsetLCBs
###########################################################################################	

# Reading properties
my $name_xmfa1 = $yaml->[0]->{C_MAUVE_output}->{xmfa} or die "Fill xmfa in the configuration file ($cfg_file)";
my $name_bbcols='';
my $name_core='';

if(!-s $name_xmfa1){
   warn "$name_xmfa is empty. There was a problem running Mauve!\n";
}
else{
   print "\tRunning stripSubsetLCBs on $name_xmfa\n";
   $name_bbcols= $name_xmfa1.".bbcols";
   $name_core= $name.".core_alignment.xmfa";
   my $cmd_stripSubsetLCBs= $path_ssLCBs."/stripSubsetLCBs $name_xmfa1 ".$name_bbcols." ".$dir_out."/".$name_core." $length";
   print "\t You are asking to run $cmd_stripSubsetLCBs\n";
   system($cmd_stripSubsetLCBs);
}

$yaml->[0]->{E_ssLCBs_output}->{bbcols} = $name_bbcols;
$yaml->[0]->{E_ssLCBs_output}->{core} = $dir_out."/".$name_core;
$yaml->write( $cfg_file );

#stripSubsetLCBs: http://gel.ahabs.wisc.edu/mauve/snapshots/2012/2012-06-07/linux-x64/. 64bit only.

############################################################################################
#Blocksplit.pl
############################################################################################
#This script is found http://clonalorigin.googlecode.com/svn/trunk/warg/scripts/
#this script is in clonalorigin/warg/scripts/

# Reading properties
   my $name_xmfa2 = $yaml->[0]->{E_ssLCBs_output}->{core} or die "Fill core in the configuration file ($cfg_file)";
   
if(!-s$name_xmfa2){
   warn "$name_xmfa is empty. There was a problem running Mauve!\n";
  }
  else{
   print "\tRunning blocksplit on $name_xmfa\n";
   my $cmd_blocksplit= $path_blocksplit."/blocksplit.pl $name_xmfa2";
 print "\t You are asking to run $cmd_blocksplit\n";
   system($cmd_blocksplit);
}

my @blocks=();

opendir(DIR, $dir_out) or die "Error in opening dir $dir_out\n";

while (my $file1= readdir(DIR)){
 next if $file1 =~/\.\.?$/;
 next if $file1 !~/(core_alignment\.xmfa\.\d+)$/;

 my $fullpath=$dir_out."/".$file1;
 push(@blocks,$fullpath);
}

my $num_blocks=scalar(@blocks);

 $yaml->[0]->{F_Blocksplit}->{num_blocks} = $num_blocks;
 $yaml->write( $cfg_file );

##########################################################################################
#BINMatrixMAUVE
##########################################################################################

# Reading properties
my $backbone=$yaml->[0]->{C_MAUVE_output}->{backbone} or die "Fill dir_out in the configuration file ($cfg_file)";
my $names=$yaml->[0]->{A_General}->{names} or die "Fill names in the configuration file ($cfg_file)";
 
print "You are asking to create a binary matrix with present and absent of conserved regions based on MAUVE results.\n";
my @block_name=();
my @coordinates=();
my @coordt1=();
my @coordt2=();
my @name_genome=();
my %hash_all;
my %hash;

#blocksplit names of the files
foreach my $block(@blocks){

#reading the core_alignment block
open IN, $block;
  my @block_file= <IN>;
   close IN;
   


foreach my $line_block(@block_file){
if($line_block=~m/^>.+:(\d+)\-(\d+).+\/(.+)\n/){
my $coordt1=$1;
my $coordt2=$2;
my $name_genome=$3;
my $coordts="$1\-$2";

push (@coordt1,$coordt1);
push (@coordt2,$coordt2);
push (@coordinates,$coordts);
push (@name_genome,$name_genome);
$block=~m/.+\/(.+core_alignment\.xmfa\.\d+)$/;
my  $justname=$1;
push (@block_name,$justname);
$hash{$coordts}=$justname;
$hash_all {$justname}{$coordts}= $name_genome;
}}}

#Creating the output file
 my $file_out=$backbone.".binary";
open OUT, ">$file_out" or die "Cannot open $file_out for writing\n";
print "$file_out was created\n";
print OUT "Conserved_Region\tCollinear_Block\tGenome_Ref\tStart_BC\tEnd_BC\tStart_CR\tEnd_CR\tLength_CR\t";

#Reading the backbone file
open IN, $backbone;
  my @lines= <IN>;
   close IN;
   
   if(length($lines[0]) > 0) 
 {
 print "$backbone is ok!\n";
  }
if(length($lines[0]) < 1) {die "$backbone is empty\n";}
   
my @genoms=split(',',$names);
my $nam=join("\t",@genoms);
print OUT "$nam\n";

my $num_seq='';
my @value=();
my $value='';
my $start='';
my $end='';
my @p=();
my %index;
my @index=();
my $index='';
my $i='';
my $q='';
my @back_coor1=();
my @back_coor2=();
my @genoms_name=();
my @binary_new=();
my $line_bin='';
my @lines_bin=(); 
my @binary=();
my %BIN;
my $to='';
my @to=();

if ( $lines[0]=~m/seq(\d+)/ ) 
{
   my @num=split('\t', $lines[0]);
   $num_seq=scalar@num; 
}
##################

shift @lines;  
chomp (@lines); 
  
foreach my $lines(@lines) 
{
  $lines=~s/$lines/\t$lines/;
  $lines=~s/\t([1-9])\t/\t0$1\t/g;
  $lines=~s/\t\-([1-9])\t/\t\-0$1\t/g;
     
     while($lines=~m/(\-*\d{2,})\t(\-*\d{2,})/g){
       $start=$1;
       $end=$2;
       $to=abs($end-$start);
       
       @p=split('\t',$lines);
       @index{@p}=(0..$#p);
	   $index=($index{$end})/2;
	   
	if ($to>=$cut ){
    	push (@to,$to);
    	push(@back_coor1,$start);
    	push(@back_coor2,$end);    
    	push(@genoms_name,$genoms[$index]);
	         
       @binary=split('\t',$lines);
       @lines_bin = map {$binary[$_]} grep {$_ & 1} 1..$#binary;
       $line_bin=join("\t",@lines_bin);
	   $line_bin=~s/(\d{2,})/1/g;
	   $line_bin=~s/-1/1/g;	
    	push(@binary_new,$line_bin);
    
       }
       }
       }
       
#Extracting info comparing the two files (backbone, and blocks).       
    
    my $k= scalar(@block_name);
	my $u='';
	my $j=scalar(@back_coor1);
	my $long1='';
	my $long2='';
	my $key='';
	my $block_collinear='';
	my $h=scalar(@to);



if ($j>=$k){$long1=$j;$long2=$k}else{$long2=$j;$long1=$k}

my $count=0;
for($u=0;$u<$long1;$u++){
		for ($i=0; $i<$long2; $i++)
	{ #print "$u\n";
	  #print "$i\n";
	$key="$coordt1[$i]\-$coordt2[$i]";
  	$block_collinear=$hash{$key};
	if( $coordt1[$i] <= $back_coor1[$u] && $coordt2[$i] >= $back_coor2[$u] && $hash_all{$block_collinear}{$key} eq $genoms_name[$u]){		    	 
  	$count++;
  	print OUT "ConserRegion\.$count\t$block_collinear\t$hash_all{$block_collinear}{$key}\t$coordt1[$i]\t$coordt2[$i]\t$back_coor1[$u]\t$back_coor2[$u]\t$to[$u]\t$binary_new[$u]\n";
     }
   }

  }


$yaml->[0]->{G_Result}->{Binary_File} = $file_out;
$yaml->write( $cfg_file );


close OUT;
      	  		
      	  		

