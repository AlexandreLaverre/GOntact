#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readGOCategories{
    my $pathin=$_[0];
    my $spacego=$_[1];
    my $gospace=$_[2];
    my $godesc=$_[3];

    open(my $input, $pathin);

    my $line=<$input>;
    chomp $line;
    my %header;
    my @s=split("\t", $line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $id=$s[$header{"ID"}];
	my $space=$s[$header{"GOSpace"}];
	my $name=$s[$header{"Name"}];

	if(exists $spacego->{$space}){
	    push(@{$spacego->{$space}}, $id);
	} else{
	    $spacego->{$space}=[$id];
	}

	if(exists $gospace->{$id}){
	    print "Weird! already saw ".$id." for ".$gospace->{$id}." and ".$space."\n";
	    exit(1);
	} else{
	    $gospace->{$id}=$space;
	}

	$godesc->{$id}=$name;
	
	$line=<$input>;
    }

    close($input);
}

##########################################################################

sub readGOAnnotations{
    my $pathin=$_[0];
    my $gospace=$_[1];
    my $space=$_[2];
    my $genego=$_[3];
    my $gogene=$_[4];

    open(my $input, $pathin);

    my $line=<$input>;
    chomp $line;
    my %header;
    my @s=split("\t", $line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $gene=$s[$header{"GeneName"}];
	my $go=$s[$header{"GOID"}];

	if(exists $gospace->{$go}){
	    if($gospace->{$go} eq $space){
		
		if(exists $genego->{$gene}){
		    push(@{$genego->{$gene}}, $go);
		} else{
		    $genego->{$gene}=[$go];
		}
		
		if(exists $gogene->{$go}){
		    push(@{$gogene->{$go}}, $gene);
		} else{
		    $gogene->{$go}=[$gene];
		}
	    }
	} else{
	    print "Weird! ".$go." doesn't belong to any GO space.\n";
	    exit(1);
	}

	$line=<$input>;
    }

    close($input);
}

#######################################################################

sub identifyRedundantCategories{
    my $gogene=$_[0];
    my $connected=$_[1];

    my %nbgenes;

    foreach my $go (keys %{$gogene}){
	my $nbg=@{$gogene->{$go}};

	if(exists $nbgenes{$nbg}){
	    push(@{$nbgenes{$nbg}}, $go); 
	} else{
	    $nbgenes{$nbg}=[$go];
	}
    }

    my %connected;

    foreach my $nbg (keys %nbgenes){
	my @goids=@{$nbgenes{$nbg}};

	my $n=@goids;

	if($n>1){
	    for(my $i=0; $i<($n-1); $i++){
		my $go1=$goids[$i];

		for(my $j=($i+1); $j<$n; $j++){
		    my $go2=$goids[$j];

		    my $allfound=1;
		    foreach my $gene1 (@{$gogene->{$go1}}){
			my $found=0;

			foreach my $gene2 (@{$gogene->{$go2}}){
			    if($gene1 eq $gene2){
				$found=1;
				last;
			    }
			}

			if($found==0){
			    $allfound=0;
			    last;
			}
		    }

		    if($allfound==1){
			## same genes

			if(exists $connected->{$go1}){
			    $connected->{$go1}{$go2}=1;
			} else{
			    $connected->{$go1}={$go2=>1, $go1=>1};
			}
			
			if(exists $connected->{$go2}){
			    $connected->{$go2}{$go1}=1;
			} else{
			    $connected->{$go2}={$go1=>1, $go2=>1};
			}
		    }
		}
	    }
	}
    }
}

#######################################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];

    my $nbconnected=keys %{$refconnected};

    my $round=0;

    while($nbconnected>0){
	
	foreach my $key (keys %{$refconnected}){
	    addToCluster($refconnected,$refclusters,$refclustid,$key);
	}
	
	$round++;
	
	$nbconnected=keys %{$refconnected};
    }
}

##############################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
	
	## find the cluster that contains this key
	
	my $indexcluster="NA";
	
	if(exists $refclustid->{$key}){
	    $indexcluster=$refclustid->{$key};
	}
	
	## if there isn't any
	
	if($indexcluster eq "NA"){
	    my $nbclusters=keys %{$refclusters};
	    $indexcluster=$nbclusters+1;
	    $refclusters->{$indexcluster}={$key=>1};
	    $refclustid->{$key}=$indexcluster;
	}
		
	foreach my $connection (keys %{$refconnected->{$key}}){

	    ## check if this island is already in the cluster
	    
	    if(!(exists $refclusters->{$indexcluster}{$connection})){
		$refclusters->{$indexcluster}{$connection}=1;
		$refclustid->{$connection}=$indexcluster;
		addToCluster($refconnected,$refclusters,$refclustid,$connection);
	    }
	}

	## after we've checked all of its connections, remove it from the connected islands

	delete $refconnected->{$key};
    }
    
}

#######################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script identifies redundant GO annotations.\n";
    print "\n";

    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }

    print "\n";
}

####################################################################################
####################################################################################

my %parameters;

$parameters{"pathGOCategories"}="NA";
$parameters{"pathGOAnnotations"}="NA";
$parameters{"GOSpace"}="NA";
$parameters{"pathOutputSimplifiedAnnotations"}="NA";
$parameters{"pathOutputGOClusters"}="NA";

my @defaultpars=("pathGOCategories", "pathGOAnnotations", "GOSpace", "pathOutputSimplifiedAnnotations", "pathOutputGOClusters");

my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;

    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];

    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
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

####################################################################################
####################################################################################

print "Reading GO categories...\n";

my %spacego;
my %gospace;
my %godesc;

readGOCategories($parameters{"pathGOCategories"}, \%spacego, \%gospace, \%godesc);

my $nbspace=keys %spacego;
my $nbgo=keys %gospace;

print "There are ".$nbgo." categories in ".$nbspace." GO spaces.\n";

print "Done.\n";

####################################################################################

print "Reading GO annotations...\n";

my %genego;
my %gogene;

my $space=$parameters{"GOSpace"};

print "GO space: ".$space."\n";

readGOAnnotations($parameters{"pathGOAnnotations"}, \%gospace, $space, \%genego, \%gogene);

my $nbgenes=keys %genego;
my $nbgo=keys %gogene;

print "There are ".$nbgenes." genes and ".$nbgo." GO annotations.\n";

print "Done.\n";

####################################################################################

print "Identifying redundant GO categories...\n";

my %connected;

identifyRedundantCategories(\%gogene, \%connected);

my $nbcon=keys %connected;

print "There are ".$nbcon." GO categories with potential redundancies.\n";

print "Done.\n";

####################################################################################

print "Extracting clusters...\n";

my %clusters;

extractClusters(\%connected, \%clusters);

my $nbclust=keys %clusters;

print "Done.\n";

####################################################################################

print "Writing output...\n";

open(my $outclust, ">".$parameters{"pathOutputGOClusters"});
open(my $outannot, ">".$parameters{"pathOutputSimplifiedAnnotations"});

print $outannot "GeneName\tGOID\n";

my %inclusters;

my %geneclust;

foreach my $idclust (keys %clusters){
    my $written=0;
    
    foreach my $go (keys %{$clusters{$idclust}}){
	print $outclust "GOCluster".$idclust."\t".$go."\t".$godesc{$go}."\n";

	$inclusters{$go}=$idclust;

	if($written==0){
	    foreach my $gene (@{$gogene{$go}}){
		print $outannot $gene."\tGOCluster".$idclust."\n";
	    }
	    
	    $written=1;
	}
    }
}

foreach my $gene (keys %genego){
    foreach my $go (@{$genego{$gene}}){
	if(!exists $inclusters{$go}){
	    print $outannot $gene."\t".$go."\n";
	}
    }
}

close($outclust);
close($outannot);

print "Done.\n";

####################################################################################

