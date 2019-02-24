#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use YAML::Tiny;
use Getopt::Long;

my $cfg_file = '';

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
  				HR (Homologous Recombination)											  
This script is used to identify the regions involved in homologous recombination. This 
script runs ClonalFrame, and GetClonalTree.pl. Additionally, it generates a list of collinear
blocks and their reference sites which probability of recombination is higher than your 
threshold.
*******************************************************************************************

 USAGE: DeEP_HR.pl -c <cfg_file.yml>
 
 Options 
	--cfg <cfg_file.yml>
	   -c <cfg_file.yml>      
							path of the configuration file to use. 
		                	For example /home/me/cfg-file.yml.
EOF
}

# Create a YAML file
my $yaml = YAML::Tiny->new();

# Open the config	
$yaml= YAML::Tiny->read($cfg_file) or die "Couldn't read file: Error($!) : Errstr =", YAML::Tiny->errstr;
	
# Reading properties
my $name  = $yaml->[0]->{A_General}->{job_name}  or die "Fill name in the configuration file ($cfg_file)";
my $parallel = $yaml->[0]->{H_ClonalFrame}->{parallel}  or die "Fill parallel (ClonalFrame) in the configuration file ($cfg_file)";
my $path_clonalframe = $yaml ->[0]->{H_ClonalFrame}->{path_ClonalFrame}  or die "Fill path_ClonalFrame in the configuration file ($cfg_file)";
my $name_core=$yaml->[0]->{E_ssLCBs_output}->{core}  or die "Fill core in the configuration file ($cfg_file)"; 
my $dir_out= $yaml->[0]->{A_General}->{dir_out} or die "Fill dir_out in the configuration file ($cfg_file)";
my $y = $yaml ->[0]->{H_ClonalFrame}->{iterations}  or die "Fill iterations in the configuration file ($cfg_file)";
my $z = $yaml ->[0]->{H_ClonalFrame}->{interval}  or die "Fill interval in the configuration file ($cfg_file)";
my $x = $yaml ->[0]->{H_ClonalFrame}->{burn_in}  or die "Fill burn_in in the configuration file ($cfg_file)";
my $others = $yaml ->[0]->{H_ClonalFrame}->{others}  or die "Fill Others in the configuration file ($cfg_file)";
my $cutoff  = $yaml->[0]->{H_ClonalFrame}->{threshold} or die "Fill threshold in the configuration file ($cfg_file)";

print "Starting main program. ClonalFrame\n";
my @childs=();
my $name_clonal= '';
my $name_core_v='';
my $cf_stdout='';

if(!-s $name_core){
   warn "$name_core is empty. There was a problem running stripSubsetLCBs!\n";
}
else{
 print "\tRunning ClonalFrame on $name_core\n";
 for (my $i=1 ;$i<$parallel+1;$i++){
  my $pid = fork();
  if ($pid) {

   push(@childs, $pid);
  }
  elsif ($pid == 0) {
   sub1($i);
   exit 0;
  }
  else {
   die "couldn't fork: $!\n";
  }
 }

 foreach (@childs) {
  my $tmp = waitpid($_, 0);
 }
 
 print "End of main program. ClonalFrame\n";
} 
 
sub sub1 {
 my $num = shift;
 $name_clonal= $name.".clonalframe.out.".$num;
 $cf_stdout="cf_stdout.".$num;
 my $cmd_ClonalFrame= $path_clonalframe."/ClonalFrame ".$others." -x ".$x." -y ".$y." -z ".$z." ".$name_core." ".$dir_out."/".$name_clonal." > ".$dir_out."/".$cf_stdout;
 print "$cmd_ClonalFrame of  $num\n";
 system($cmd_ClonalFrame);
}

my $tree=$dir_out."/".$name.".clonalframe.out";


##########################################################################################
#Verifying if you run enough iterations
##########################################################################################

#consensus tree?
my $tree1= $tree.".1";
my $tree2= $tree.".2";
my $tree3= $tree.".3";
my $tree4= $tree.".4";
my $tree5= $tree.".5";
my $q='';
my $s='';



if(!-s $tree1){
   die "$tree1 is empty. There was a problem running ClonalFrame!\n";
}
else{
print "Comparing the first 5 results from ClonalFrame.\n"; 
print "First: $tree1\n";
print "Second: $tree2\n";
print "Third: $tree3\n";
print "Fourth: $tree4\n";
print "Fifth: $tree5\n";
}

print "Comparing parameters (theta, rho, delta, nu).\n";

#Extracting data form the first file.
#Theta comparison

open IN, $tree1;
my @theta;
my $theta;
while (<IN>) {
  if (/#theta/../#nu/) {
    next if /#theta/ || /#nu/;
push (@theta, $_);  
}
}
chomp (@theta);
close IN;


open IN, $tree1;
my @nu;
my $nu;
while (<IN>) {
  if (/#nu/../#delta/) {
    next if /#nu/ || /#delta/;
push (@nu, $_);  
}
}
chomp (@nu);
close IN;


open IN, $tree1;
my @delta;
my $delta;
while (<IN>) {
  if (/#delta/../#R/) {
    next if /#delta/ || /#R/;
push (@delta, $_);  }
}
chomp (@delta);
#print "#delta values are: @delta\n";
close IN;

open IN, $tree1;
my @R;
my $R;
while (<IN>) {
  if (/#R/../#end/) {
    next if /#R/ || /#end/;
push (@R, $_);  }
}
chomp (@R);
close IN;

#extrating the tree n.1

open FILE, $tree1;
my @lines1 = <FILE>;
my @contree=();
my $node='';
my $line1=$lines1[1];
while($line1=~m/(\d+)\:/g){
$node=$1;
push(@contree,$node);
}
chomp (@contree);
close FILE;

#####################
#Extracting data from the other files
  my $long= scalar(@contree);
  my $t=scalar(@theta);
  my $n=scalar (@nu);
  my $d=scalar(@delta);
  my $r=scalar(@R);

my @array=($t,$n,$d,$r);
my $max = (sort { $b <=> $a } @array)[0];

#print "this is the max $max\n";
my @list =($tree2,$tree3,$tree4,$tree5);

  foreach my $list_comparison(@list){

open IN, $list_comparison;
my @theta1;
my $theta1;
while (<IN>) {
  if (/#theta/../#nu/) {
    next if /#theta/ || /#nu/;
push (@theta1, $_);  }
}
chomp (@theta1);
close IN;

open IN, $list_comparison;
my @nu1;
my $nu1;
while (<IN>) {
  if (/#nu/../#delta/) {
    next if /#nu/ || /#delta/;
push (@nu1, $_);  }
}
chomp (@nu1);
close IN;

open IN, $list_comparison;
my @delta1;
my $delta1;
while (<IN>) {
  if (/#delta/../#R/) {
    next if /#delta/ || /#R/;
push (@delta1, $_);  }
}
chomp (@delta1);
close IN;

open IN, $list_comparison;
my @R1;
my $R1;
while (<IN>) {
  if (/#R/../#end/) {
    next if /#R/ || /#end/;
push (@R1, $_);  }
}
chomp (@R1);
close IN;

open FILE, $list_comparison;
my @lines_comparison = <FILE>;
my @contree_comparison=();
my $node_comparison='';
my $line_comparison=$lines_comparison[1];
while($line_comparison=~m/(\d+)\:/g){
$node_comparison=$1;
push(@contree_comparison,$node_comparison);
}
close FILE;

print "Comparing parameters $tree1  vs. $list_comparison\n";

    for ($q=0;$q<$max;$q++){
    chomp;
    print "#theta= $theta[$q] should be equal to #theta= $theta1[$q]\n";
    print "#nu= $nu[$q] should be equal to #nu= $nu1[$q]\n";
    print "#delta= $delta[$q] should be equal to #delta= $delta1[$q]\n";
    print "#R= $R[$q] should be equal to #R= $R[$q]\n";
    if( $theta[$q] eq $theta1[$q] && $nu[$q] eq $nu1[$q] && $delta[$q] eq $delta1[$q] && $R[$q] eq $R1[$q])
    {
    print "*** All parameters coincide ***\n";
    }
    else {die "Parameters don't coincide. You will have to make a upper number of interations.\n"}
    }
    
print "Comparing nodes: $tree1  vs. $list_comparison\n";

     for ($s=0;$s<$long;$s++){
    print "$contree[$s] should be equal to $contree_comparison[$s]\n";
	if($contree[$s] eq $contree_comparison[$s])
	{
     print "nodes coincide\n";
	}
     else{ warn "Nodes don't coincide. You may have to make a upper number of iterations (Optional). Trees should be similar but don't have to be identical, if the parameters are roughly the same.\n"}
    }
}
print "A concensus tree was found\n";
$yaml->[0]->{I_ClonalFrame_output}->{clonal_concensus} = $tree.".1";
$yaml->write( $cfg_file );


##################################################################################################################################

  my $clonal_census  = $yaml->[0]->{I_ClonalFrame_output}->{clonal_concensus} or die "Fill clonal_concensus in the configuration file ($cfg_file)";
  my $path_clonaltree = $yaml->[0]->{J_ClonalTree}->{path_clonaltree} or die "Fill path_clonaltree in the configuration file ($cfg_file)";
  my $tree_con=''; 

if(!-s$clonal_census){
   warn "$clonal_census is empty. There was a problem running ClonalFrame!\n";
  }
  else{
   print "\tRunning getClonalTree on $clonal_census\n";
   $tree_con=$dir_out."/".$name.".clonaltree.nwk";
   my $cmd_clonaltree= $path_clonaltree."/getClonalTree $clonal_census $tree_con";
 print "\t You are asking to run $cmd_clonaltree\n";
   system($cmd_clonaltree);
}

$yaml->[0]->{K_ClonalTree_output}->{final_tree} = $tree_con;
$yaml->write( $cfg_file );

#this script is in clonalorigin/warg/scripts/
###################################################################################################################################
#Aqui lo que debo hacer es leer el archivo de salida, las partes: #consevents, #consinfo, #poly
#y generar una tabla de bloques vs recombinación homologa (1/0). Utilizando un umbral que el usuario escoja.

 my $file_out= $dir_out."/".$name.".Rcomb_".$cutoff.".out";
open OUT, ">$file_out" or die "Cannot open $file_out for writing\n";
print "$file_out was created\n";
print OUT "Collinear_Block\tRefsite\tNode\tRecombination_probability\n";


for my $num_clonalout (1..$parallel){
my $file_clonalout=$tree.".".$num_clonalout;

#Recombination probability.
#Extraer la tabla

open IN, $file_clonalout;
my @table_consevents;
my $table_consevents;
while (<IN>) {
  if (/#consevents/../#consinfo/) {
    next if /#consevents/ || /#consinfo/;
push (@table_consevents, $_);  }
}
chomp (@table_consevents);
close IN;


open IN, $file_clonalout;
my @table_poly;
my $table_poly;
while (<IN>) {
  if (/#poly/../#ll/) {
    next if /#poly/ || /#ll/;
push (@table_poly, $_);  }
}
chomp (@table_poly);
close IN;

open IN, $file_clonalout;
my @table_blocks;
my $table_blocks;
while (<IN>) {
  if (/#blocks/../#theta/) {
    next if /#blocks/ || /#theta/;
push (@table_blocks, $_);  }
}
#chomp (@table_blocks);
close IN;
#####
#creating the output file



#extracting the probabilities of recombination for each reference site
my $Rcomb='';
my $Susti='';
my $count=0;
my @Rcomb=();
my @branch=();
my $data=0;

foreach my $LineR_prob(@table_consevents){
$count++;#counting of branches. counting for line.(solo una linea).
while ($LineR_prob=~m/(\d.\d+e[\+|\-]\d+)\s(\d.\d+e[\+|\-]\d+)/g){
$data++;# position of the reference site. counting each time it found a value per line (varios por linea). 
$Rcomb=$2;
$Susti=$1;
#print "$Rcomb\n";
push (@branch,$count);
push (@Rcomb,$Rcomb);
}
}

#extracting the reference sites for each block.
my $refe_site='';
my $value_end='';
my @value_end=();
my $realstart='';
my $h='';
my $rest='';
my $realend='';
my @poly=();
my @j=();
my @seq_refsite=();
my $u='';
my $g='';
my @data_num=();


#aqui me esta generando un error que no entiendo... VERIFICAR
#Argument "" isn't numeric in addition (+) at ./comparison.pl line 362 (#2)
    #(W numeric) The indicated string was fed as an argument to an operator
    #that expected a numeric value instead.  If you're fortunate the message
    #will identify which operator was so unfortunate.

for my $i (0..$#table_blocks) {
  $value_end += $table_blocks[$i];
  push(@value_end,$value_end);
}

my $l=scalar(@table_blocks);
unshift(@value_end,0);
unshift(@table_poly,0);
push (@value_end,0);

####
my @blocks_num=();
my %hash; 
my $data_num=0;
 
 for (my $j=1;$j<=$l;$j++){
$g=$j-1;
$u=$j+1;
$realend=$value_end[$j];
$realstart=($realend-$table_blocks[$g]+1);
@seq_refsite=@table_poly[$realstart..$realend];
foreach my $refsite(@seq_refsite){
$data_num++; #position of the reference site. Solo cuenta una vez cada refence site. El anterior cuanta 9 veces cada refsite
push (@data_num,$data_num);
push (@poly,$refsite);
push (@blocks_num,$j);
my $blockID="$j\.$data_num";
$hash{$refsite}=$blockID; #$j aca es el numero de bloque
}
}
my @Rcomb_new=();
my $IDRH=0;
foreach my $Rcomb_line(@Rcomb){
my $Rcomb_ID=$Rcomb_line."-".$IDRH++;
push (@Rcomb_new,$Rcomb_ID);}


 my $repetitive=scalar(@table_consevents);
 my @branch_all=();
 my @branch_number=();
 
 
my $branch_numb=0;
for (1..$repetitive){
$branch_numb++;
foreach my $line_ProbaRC(@poly){
push (@branch_all,$line_ProbaRC);
push (@branch_number,$branch_numb);
}}

#branch_all (@poly) es la que relaciona a traves de $refsite los dos hashes. 
my %hash_all;
@hash_all{@Rcomb_new}=@branch_all;

my %hash_branch;
@hash_branch{@Rcomb_new}=@branch_number;



foreach my $key (keys %hash_all){
$key=~m/(.+)\-/;
#print "$key\n";
my $qw=$1;
#print "$key\n";
if($qw>=$cutoff){
my $key1=$hash_all{$key};
my $block_real=$hash{$key1};
$block_real=~m/(\d+)\./;
print OUT "$1\t$hash_all{$key}\t$hash_branch{$key}\t$qw\n";
}

}
}

$yaml->[0]->{L_Result}->{Summary_File} = $file_out;
$yaml->write( $cfg_file );


#my $xmfa = $yaml->[0]->{C_MAUVE_output}->{xmfa} or die "Fill xmfa  in the configuration file ($cfg_file)";
my $clonaltree = $yaml->[0]->{K_ClonalTree_output}->{final_tree} or die "Fill final_tree in the configuration file ($cfg_file)";
my $names_species = $yaml->[0]->{A_General}->{name} or die "Fill name in the configuration file ($cfg_file)";

my $file_tree_names="$dir_out\/$name\.clonaltree_tipslabels.nwk";
   open (TREE, '>', $file_tree_names) or die "Cannot open $file_tree_names for writing\n";
       
my @species_names=split(',',$names_species);
my $last=scalar(@species_names)+1;
my @nodes_num_specie=(1..$last);

my %hash_tree;
@hash_tree{@nodes_num_specie}=@species_names;

my $max1='';
my @max1=();
 
  open IN, $clonaltree;
  my @lines_clonaltree= <IN>;
chomp(@lines_clonaltree);
 close IN;
my $tree_concensus='';
    
  foreach my $line_clonaltree(@lines_clonaltree){
print "$line_clonaltree\n";
   while($line_clonaltree=~m/([\(|\)|\,])(\d+)\:/g){
    my $sing="\\".$1;
    my $seq=$2;
my $name_seq=''; 
push(@max1,$seq);
if(exists $hash_tree{$seq}){
    $name_seq=$hash_tree{$seq}."\:";
 $line_clonaltree=~s/$sing$seq\:/$sing$name_seq/g;
$line_clonaltree=~s/\\//;
}
}
$tree_concensus=$line_clonaltree;
}
print TREE "$tree_concensus\n";
 
$yaml->[0]->{I_ClonalFrame_output}->{clonal_concensus} = $file_tree_names;
$yaml->write( $cfg_file );

