#!/usr/bin/perl
use warnings;
use diagnostics;
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
								GC(GC content)											  
This script is used to identify the regions involved in Lateral Gene Transfer (LGT) using 
the GC content. This script runs window-acgt (GLIMMER program), and uses a R-Script to 
make a t-test and a p-value adjustment. Additionally, it generates a binary matrix where 
1 represents the regions involved in LGT and 0 the regions that are not.
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
   my $dir_files = $yaml->[0]->{A_General}->{dir_files} or die "Fill dir_files in the configuration file ($cfg_file)";
   my $path_window = $yaml->[0]->{M_window_acgt}->{path_window} or die "Fill path_window in the configuration file ($cfg_file)";
   my $window_len  = $yaml->[0]->{D_ssLCBs}->{length} or die "Fill length in the configuration file ($cfg_file)";
   my $window_skip = $yaml ->[0]->{M_window_acgt}->{window_skip} or die "Fill window_skip in the configuration file ($cfg_file)"; 
   my $dir_out= $yaml->[0]->{A_General}->{dir_out} or die "Fill dir_out in the configuration file ($cfg_file)";
   my $Rscript_path =  $yaml->[0]->{N_Rscript}->{Rscript_path} or die "Fill dir_out in the configuration file ($cfg_file)";
my $parallel = $yaml->[0]->{N_Rscript}->{parallel}  or die "Fill parallel (N_Rscript) in the configuration file ($cfg_file)";

#contenido de GC de los genomas 
my @Genomes_names=();

opendir(DIR, $dir_files) or die "Error in opening dir $dir_files\n";
while (my $file= readdir(DIR)){
 next if $file =~/\.\.?$/;
 next if $file !~/(.fna|.fa|.fasta)$/;
  my $fullpath_file=$dir_files."/".$file;
 push(@Genomes_names,$file); 
 push(@Myfiles,$fullpath_file); 
}
my $list_genomes=join(',',@Genomes_names);

my $num_genomes=scalar(@Myfiles);

 foreach my $files(@Myfiles){
 my $file_out=$files.".gc";
  print "\tRunning window-acgt on this file: $files\n";
 my $cmd_window= $path_window."window-acgt -p $window_len $window_skip < $files >> $file_out";
print "\t You are asking to run this command: $cmd_window\n";
 system($cmd_window);
}

closedir(DIR);
my @Regions_names=();

#contenido de GC de las regiones
opendir(DIR, $dir_out) or die "Error in opening dir $dir_out\n";
while (my $file_reg= readdir(DIR)){
 next if $file_reg =~/\.\.?$/;
 next if $file_reg !~/(core_alignment.xmfa\.\d+)$/;
  my $fullpath_file_reg=$dir_out."/".$file_reg;
 push(@Myfiles_reg,$fullpath_file_reg); 
 push(@Regions_names,$file_reg);
}

my $list_regions=join(',',@Regions_names);
my $num_region=scalar(@Myfiles_reg);

 foreach my $files_reg(@Myfiles_reg){
 my $file_out_reg=$files_reg.".gc";
  print "\tRunning window-acgt on this file: $files_reg\n";
 my $cmd_window= $path_window."window-acgt -p $window_len $window_skip < $files_reg >> $file_out_reg";
print "\t You are asking to run this command: $cmd_window\n";
 system($cmd_window);

}
closedir(DIR);

########
                   my $gc_file_Reg='';
                  
#Regions
 ##############
   foreach my $Reg_vector(@Myfiles_reg){   
   $gc_file_Reg=$Reg_vector.".gc";
        
open IN, $gc_file_Reg;
 my @lines = <IN>;
 my $file_out_edited=$gc_file_Reg;
open TAB, ">$file_out_edited" or die "Cannot open $file_out_edited for writing\n";
    
         foreach my $line(@lines){
         if($line=~/^\>/)
            { $line= "#end\n$line";
            print TAB $line;
             }
             else{print TAB $line;
  		}
             }
print TAB "#end\n";
close IN;
close TAB;
}
#####
#Reading the outputs to then compare.
#Genomes  

my $genome_name_out='';

   foreach my $gc_genome(@Myfiles){ 
   
   my $gc_file_Gm=$gc_genome.".gc";
	open IN, $gc_file_Gm;
  	my @gc_vector_GM= <IN>;
   	close IN;
   	$gc_genome=~m/.+\/(.+)$/g;
   	$genome_name_out=$1;

    my $file_vectors_out="$dir_out\/vector_gc_$genome_name_out";
     open ($file_vectors_out, '>', $file_vectors_out) or die "Cannot open $file_vectors_out for writing\n";
	system("chmod a+w $file_vectors_out");

     print "$file_vectors_out was created\n";
     print $file_vectors_out "\"$gc_genome\"\t";
      
    foreach my $line_gc_GM(@gc_vector_GM){
	next if $line_gc_GM=~/^>/;
	next if $line_gc_GM=~/^P/; 
       $line_gc_GM=~m/(\d+\.\d+)$/;
      my $vector_GM=$1;
      print "Extracting \%GC vector of:\n$gc_genome\n";
      print $file_vectors_out "$vector_GM\t";
     }
      print $file_vectors_out "\n";
	
		   my $vector_Reg='';
		   my $genome_nameblock='';
		   my $name_vector_reg='';
		   my $gc_file_Reg1='';

#Regions
 ##############
   foreach my $Reg_vector(@Myfiles_reg){ 
   
   $gc_file_Reg1=$Reg_vector.".gc";
   $gc_file_Reg1=~m/.+\/(.+)$/;
   $name_vector_reg=$1;

my @vector;	 
 my $vector; 
open IN, $gc_file_Reg1;
	while (<IN>) {
	if (/.+$genome_name_out/.. /\#end/){
	if(/\s(\d+\.\d+)$/){
	  push (@vector,$1);
	 }
	 }
	}
	chomp (@vector);
	$vector=join("\t",@vector);
	

my $print_file="$dir_out\/vector_gc_$genome_name_out";
print $print_file "\"$name_vector_reg\"\t$vector\n";
print "Extracting \%GC vector of:\n$name_vector_reg\.\.\.\.$genome_name_out\n"; 
close IN;
}
	}


my @Myfiles_vector=();

#contenido de GC de las regiones
opendir(DIR, $dir_out) or die "Error in opening dir $dir_out\n";
while (my $file_vector= readdir(DIR)){
 next if $file_vector =~/\.\.?$/;
 next if $file_vector !~/(vector\_gc)/;
  my $fullpath_file_vector=$dir_out."/".$file_vector;
 push(@Myfiles_vector,$fullpath_file_vector); 
}


my @RSCRIPTS=();


foreach my $gc_analysis(@Myfiles_vector){
##########################################################################
# R.SCRIPT
##########################################################################

my $cmd_Rscript='';
my $FDR_file="$gc_analysis\_FDR.pdf";
my $R_file="$gc_analysis\_GC.R";
my $file_out_analysis= "$gc_analysis\_LT\_statistics.txt";


open R_SCRIPT,">$R_file" or die "Cannot write $R_file script\n";
print "Creating $R_file\n";

        my $R_script= "

library('fdrtool')
line <- readLines('$gc_analysis')
Vectors <- lapply(line,function(x) unlist(strsplit(x,split='\\t')))
l=length(Vectors)
names_vectors=sapply(Vectors, '[[', 1)
names(Vectors)=names_vectors
Genome=as.numeric(Vectors[[1]][-1])
names(Genome)=Vectors[[1]][1]
names_regions=names_vectors[-1]

matrix_gc<-matrix(data=NA,nrow=length(Vectors)-1,ncol=1)
rownames(matrix_gc)=names_regions

for(i in names_regions){ 
b=Genome
a=as.numeric(Vectors[[i]][-1])
names(a)=Vectors[[i]][1]
c=Genome
# Observed difference
diff.observed = mean(b) - mean(a)
number_of_permutations = 1000
diff.random = NULL
for (p in 1 : number_of_permutations) {
  # Sample from the combined dataset
  a.random = sample (c, length(a), TRUE)
  b.random = sample (c, length(b), TRUE)
  # Null (permuated) difference
  diff.random[p] = mean(b.random) - mean(a.random)
}
pvalue_low = sum(diff.random >= diff.observed) / number_of_permutations
pvalue_high =  sum(diff.random <= diff.observed) / number_of_permutations
pvalue <- min(pvalue_high,pvalue_low)
#print (pvalue)
matrix_gc[i,]<- pvalue  
}

y=matrix_gc[!is.na(matrix_gc),]
names_y=names(matrix_gc[!is.na(matrix_gc),])

pdf('$FDR_file')
FDR=fdrtool(y,statistic=c('pvalue'))
dev.off()


matrix_FDR=cbind(FDR\$pval,FDR\$qval)
colnames(matrix_FDR)=c('pvalue','qvalue')
decision=data.frame(LT=rep(NA,length=length(matrix_FDR[,'pvalue'])), row.names=c(names_y))

for(i in rownames(matrix_FDR)){
  if(matrix_FDR[i, 'qvalue' ]<=0.05){
    decision[i,]='LT'
  }else{
    decision[i,]='NO'
  }
}

matrix_results=cbind(matrix_FDR,decision)

write.table(matrix_results, file = '$file_out_analysis')

";


        print R_SCRIPT $R_script;
        close R_SCRIPT;
        system("chmod a+x $R_file"); 

	$cmd_Rscript="$Rscript_path $R_file --save";
	system( $cmd_Rscript);
	push(@RSCRIPTS,$cmd_Rscript);
}


use Parallel::ForkManager;
#### debo agregar aqui el numero de paralelos que escoja el usuario
####Lo logro hacer que esto corra adecuadamente????

#my $manager = new Parallel::ForkManager($parallel);
#foreach my $command (@RSCRIPTS) { 
#$manager->start($command) and next;
#print "You are asking to run this command: $command\n";
 #system( $command );
 #$manager->finish; }; 
#print "Waiting for GC.R script to finish...\n";
 # $pm->wait_all_children;
  #print "All GC.R scripts done!\n";
