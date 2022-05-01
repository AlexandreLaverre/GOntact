#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }

    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;

	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";

	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;

	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}


#########################################################################################

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
    }
}

#########################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script extracts enhancer sequences. \n";
    print "\n";
    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################################
##########################################################################################

my %parameters;

$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathCoordinates"}="NA";
$parameters{"coordConvention"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathGenomeSequence", "pathCoordinates", "coordConvention", "pathOutput");

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

print "Reading genome sequence...\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);

print "Done.\n";

#####################################################################################

print "Reading input coordinates and writing output for GC content...\n";

my $coordtype=$parameters{"coordConvention"};

my $offsetstart=0;
my $offsetend=0;

if($coordtype eq "0_open_end"){
    $offsetstart=1;
    $offsetend=0;
} else{
    if($coordtype eq "1_closed_end"){
	$offsetstart=0;
	$offsetend=0;
    } else{
	print "Weird! unknown coordinate convention: ".$coordtype."\n";
	exit(1);
    }
}

print "Adding offset ".$offsetstart." for start coordinates to make them one-based, closed-end.\n";
print "Adding offset ".$offsetend." for end coordinates to make them one-based, closed-end.\n";

open(my $output, ">".$parameters{"pathOutput"});
open(my $input, $parameters{"pathCoordinates"});

my $line=<$input>; ## no header

while($line){
    chomp $line;
    my @s=split("\t", $line);

    my $chr=$s[0];
    my $start=$s[1]+$offsetstart;
    my $end=$s[2]+$offsetend;

    my $id=$chr.":".$start."-".$end;

    my $seq=substr $genome{$chr}, ($start-1), ($end-$start+1);

    $seq=uc $seq;

    writeSequence($seq, $id, $output);

    $line=<$input>;
}

close($output);
close($input);

print "Done.\n";

#####################################################################################
