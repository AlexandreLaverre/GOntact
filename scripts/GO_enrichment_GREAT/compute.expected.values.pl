use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readRegulatoryRegions{
    my $pathin=$_[0];
    my $regions=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;
    chomp $line;
    my %header;
    my @s=split("\t", $line);
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    my %duplicated;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $gene=$s[$header{"gene_name"}];
	my $chr=$s[$header{"chr"}];
	my $start=$s[$header{"region_start"}];
	my $end=$s[$header{"region_end"}];

	if(exists $regions->{$gene}){
	    $duplicated{$gene}=1;
	} else{
	    $regions->{$gene}={"chr"=>$chr, "start"=>$start, "end"=>$end};
	}
	
	$line=<$input>;
    }

    close($input);

    ## we remove duplicated gene names

    my $nbd=keys %duplicated;
    
    print "Found ".$nbd." duplicated gene names, we remove them.\n";

    foreach my $d (keys %duplicated){
	delete $regions->{$d};
    }
}

####################################################################################

sub readGOCategories{
    my $pathin=$_[0];
    my $gocat=$_[1];

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

	if(exists $gocat->{$space}){
	    push(@{$gocat->{$space}}, $id);
	} else{
	    $gocat->{$space}=[$id];
	}
	
	$line=<$input>;
    }

    close($input);
}

####################################################################################

sub readGOAnnotations{
    my $pathin=$_[0];
    my $genego=$_[1];
    my $gogene=$_[2];

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
	
	$line=<$input>;
    }
    
    close($input);

}

####################################################################################

sub computeExpectedValues{
    my $gocat=$_[0]; ## GO categories by space
    my $gogene=$_[1]; ## GO annotations - gene
    my $regions=$_[2];
    my $totalsizes=$_[3];

    my %goregions;

    foreach my $space (keys %{$gocat}){
	foreach my $cat (@{$gocat->{$space}}){
	    if(exists $gogene->{$cat}){
		foreach my $gene (@{$gogene->{$cat}}){
		    if(exists $regions->{$gene}){
			my $chr=$regions->{$gene}{"chr"};
			my $start=$regions->{$gene}{"start"};
			my $end=$regions->{$gene}{"end"};

			if(exists $goregions{$cat}){
			    push(@{$goregions{$cat}{"chr"}}, $chr);
			    push(@{$goregions{$cat}{"start"}}, $start);
			    push(@{$goregions{$cat}{"end"}}, $end);
			} else{
			    $goregions{$cat}={"chr"=>[$chr], "start"=>[$start], "end"=>[$end]};
			}
		    }
		}
	    }
	}
    }

    foreach my $cat (keys %goregions){
	my %hashpos;

	my $nbr=@{$goregions{$cat}{"chr"}};

	for(my $i=0; $i<$nbr; $i++){
	    my $chr=${$goregions{$cat}{"chr"}}[$i];
	    my $start=${$goregions{$cat}{"start"}}[$i];
	    my $end=${$goregions{$cat}{"end"}}[$i];

	    for(my $pos=$start; $pos<=$end; $pos++){
		my $key=$chr.":".$pos;
		$hashpos{$key}=1;
	    }
	}

	my $nb=keys %hashpos;

	$totalsizes->{$cat}=$nb;
    }    
}

####################################################################################
####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts promoters and enhancers found in neighboring regulatory regions.\n";
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
$parameters{"pathRegulatoryRegions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGOCategories", "pathGOAnnotations", "pathRegulatoryRegions", "pathOutput");

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

my %gocat;
readGOCategories($parameters{"pathGOCategories"}, \%gocat);

print "Done.\n";

####################################################################################

print "Reading GO annotations...\n";

my %genego;
my %gogene;

readGOAnnotations($parameters{"pathGOAnnotations"}, \%genego, \%gogene);

my $nbgene=keys %genego;
my $nbgo=keys %gogene;

print "There are ".$nbgene." genes and ".$nbgo." GO categories.\n";

####################################################################################

print "Reading regulatory regions...\n";

my %regions;
readRegulatoryRegions($parameters{"pathRegulatoryRegions"}, \%regions);

my $nbgr=keys %regions;

print "There are ".$nbgr." genes with regulatory regions.\n";

print "Done.\n";

####################################################################################

print "Computing expected values...\n";

my %totalsizes;

computeExpectedValues(\%gocat, \%gogene, \%regions, \%totalsizes);

print "Done.\n";

####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "ID\tGOSpace\tNbBases\n";

foreach my $space (keys %gocat){
    foreach my $go (@{$gocat{$space}}){
	if(exists $totalsizes{$go}){
	    print $output $go."\t".$space."\t".$totalsizes{$go}."\n";
	}
    }
}

close($output);

print "Done.\n";

####################################################################################
