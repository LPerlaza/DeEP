#!/usr/bin/perl
#######################################################################
## Main script to run DeEP.Pro
#######################################################################
###ESTE SCRIPT DEBE SER MODIFICADO PARA QUE APLIQUE PARA DEeP.Pro
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

my $task='';
my $DeEP_path='';
my $conf='';
my $module= '';

GetOptions(
    'path|p=s' => \$DeEP_path,
    'task|t=s' => \$task,
    'module|m=s' => \$module, 
    'conf|c=s' => \$conf,
    );

if(!$DeEP_path){
  warn "No path of the scripts folder specified. Can not execute DeEP.\n";
  &usage();
exit(1);
}

if(!$module && !$task){
  warn "No task or module defined. Can not continue.\n";
  &usage();
exit(1);
}

if(!$conf){
  warn "No configuration file specified. Can not continue.\n";
  &usage();
exit(1);
}

sub usage{
 print STDERR <<EOF;
 This is the main script to run DeEP
 
	--path    -p  This is the path where DeEP.Pro's Scripts are.
	--task    -t  This is the task you want to run. Script by Script. Separated by commas. 
                  1 = DeEP.Pro_MultiAlign.pl
                  2 = DeEP.Pro_HR.pl
                  3 = DeEP.Pro_GC.pl
		  4 = DeEP.Pro_CorRegions.pl
		  5 = DeEP.Pro_BackCoor.pl

	--module  -m  This is the module you want to run. Several scripts as module.
 		  MiEvo = Make the whole microevolutionary analysis (equivalent "-t 1,2,3")
		  all = Run the whole pipeline. DeEP.Pro complete.

	--conf    -c  This is the configuration file in yaml.
EOF
}

 my %hash=   (    1 => 'DeEP.Pro_MultiAlign.pl',
                  2 => 'DeEP.Pro_HR.pl',
                  3 => 'DeEP.Pro_GC.pl',
                  4 => 'DeEP.Pro_CorRegions.pl',
                  5 => 'DeEP.Pro_BackCoor.pl',
);

my @alltask=();
push (@alltask,$task);

if ($module =~m/P\/A/g) {$task= "1,2,3,11"; push (@alltask,$task);};
if ($module =~m/HR/g){ $task= "4,5,6,7,8,11";push (@alltask,$task);};
if ($module =~m/LT/g){ $task= "9,11"; push (@alltask,$task);};
if ($module =~m/DP/g){ $task= "10,11"; push (@alltask,$task);};
if ($module =~m/MiEvo/g){ $task= "1,2,3,4,5,6,7,8,9,10,11"; push (@alltask,$task);};
if ($module =~m/AR/g){ $task= "12"; push (@alltask,$task);};
if ($module =~m/PIC/g){ $task= "13"; push (@alltask,$task);};
if ($module =~m/all/g){ $task= "1,2,3,4,5,6,7,8,9,10,11,12,13"; push (@alltask,$task);};

chomp (@alltask);
my $alltask=join(',',@alltask);
$alltask=~s/\s/,/;

my $Script='';
my @task=();
my $eachtask='';
my $cmd='';

if($alltask=~m/./) {
@task=split(',',$alltask);
foreach $eachtask(@task){
$Script=$hash{$eachtask};
$cmd= $DeEP_path."/".$Script." -c $conf\n";
print "You are asking to run $cmd\n";
system($cmd);
}
}


