#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readGOTree{
    my $pathin=$_[0];
    my $gonames=$_[1];
    my $godef=$_[2];
    my $gospace=$_[3];
    my $connections=$_[4];

    open(my $input, $pathin);

    my $line=<$input>;

    my $currentterm="NA";

    while($line){
	chomp $line;

	my @s=split(": ", $line);
	
	if($s[0] eq "id"){
	    $currentterm=$s[1];
	}

	if($s[0] eq "is_a"){
	    my @t=split(" ! ", $s[1]);
	    my $othergo=$t[0];

	    if(exists $connections->{$currentterm}){
		$connections->{$currentterm}{$othergo}=1;
	    } else{
		$connections->{$currentterm}={$othergo=>1};
	    }
	}

	if($s[0] eq "name"){
	    $gonames->{$currentterm}=$s[1];
	}

	if($s[0] eq "def"){
	    my @t=split("\"", $s[1]);
	    $godef->{$currentterm}=$t[1];
	}

	if($s[0] eq "namespace"){
	    $gospace->{$currentterm}=$s[1];
	}
	
	
	$line=<$input>;
    }

    close($input);
}

#######################################################################

sub readGeneAnnotations{
    my $pathin=$_[0];
    my $annot=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;

	my $prefix=substr $line, 0, 1;

	if($prefix ne "!"){
	    my @s=split("\t", $line);
	    
	    my $gene=$s[2];
	    my $cat=$s[4];
	    
	    if(exists $annot->{$gene}){
		$annot->{$gene}{$cat}=1;
	    } else{
		$annot->{$gene}={$cat=>1};
	    }
	}
	
	$line=<$input>;
    }

    close($input);
}

#######################################################################

sub propagateGOAnnotations{
    my $annot=$_[0];
    my $connections=$_[1];

    my $nbdone=0;
    
    foreach my $gene (keys %{$annot}){
	my @origcat=keys %{$annot->{$gene}};
	
	foreach my $gocat (@origcat){
	    addGOAnnotation($gene, $annot->{$gene}, $gocat, $connections);
	}

	$nbdone++;

	print $nbdone." genes done.\n";
    }
}

#######################################################################

sub addGOAnnotation{
    my $gene=$_[0];
    my $thisannot=$_[1];
    my $thisgo=$_[2];
    my $connections=$_[3];
   
    if(exists $connections->{$thisgo}){
	#print "adding connections for ".$thisgo." for ".$gene."\n";
	
	foreach my $othergo (keys %{$connections->{$thisgo}}){
	 #print "looking at connection ".$othergo." for ".$thisgo."\n";
	    
	    $thisannot->{$othergo}=1;
	    addGOAnnotation($gene, $thisannot, $othergo, $connections);
	}
    }
}

#######################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script propagates GO annotations. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

########################################################################
########################################################################

my %parameters;

$parameters{"pathGOTree"}="NA";
$parameters{"pathInputGeneAnnotation"}="NA";
$parameters{"pathOutputGOCategories"}="NA";
$parameters{"pathOutputGeneAnnotation"}="NA";

my %defaultvalues;
my @defaultpars=("pathGOTree", "pathInputGeneAnnotation", "pathOutputGOCategories", "pathOutputGeneAnnotation");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

#####################################################################################
#####################################################################################

print "Reading gene ontology tree...\n";

my %gonames;
my %godef;
my %gospace;
my %connections;

readGOTree($parameters{"pathGOTree"}, \%gonames, \%godef, \%gospace, \%connections);

my $nbc=keys %connections;

print "Found connections for ".$nbc." GO categories.\n";

print "Done.\n";

#####################################################################################

print "Writing simplified output for GO categories...\n";

open(my $outputgo, ">".$parameters{"pathOutputGOCategories"});

print $outputgo "ID\tGOSpace\tName\n";

foreach my $id (keys %gonames){
    if($gospace{$id} ne "external"){
	print $outputgo $id."\t".$gospace{$id}."\t".$gonames{$id}."\n";
    }
}

close($outputgo);

print "Done.\n";

#####################################################################################

print "Reading gene annotations...\n";

my %geneannot;

readGeneAnnotations($parameters{"pathInputGeneAnnotation"}, \%geneannot);

my $nbg=keys %geneannot;

print "Found annotations for ".$nbg." genes.\n";

print "Done.\n";

#####################################################################################

print "Propagating gene annotations...\n";

propagateGOAnnotations(\%geneannot, \%connections);

print "Done.\n";

#####################################################################################

print "Writing output for gene annotations...\n";

open(my $output, ">".$parameters{"pathOutputGeneAnnotation"});

print $output "GeneName\tGOID\n";

foreach my $gene (keys %geneannot){
    foreach my $id (keys %{$geneannot{$gene}}){
	print $output $gene."\t".$id."\n";
    }
}

close($output);

print "Done.\n";

#####################################################################################
