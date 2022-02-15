use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
####################################################################################

sub readGOEnhancers{
    my $pathin=$_[0];
    my $goenh=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;
    my $prefix=substr $line, 0, 1;

    while($prefix eq "#"){
	$line=<$input>;
	$prefix=substr $line, 0, 1;
    }
    
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

	my $term=$s[$header{"GOTerm"}];
	my @enh=split(",", $s[$header{"ForegroundEnhancers"}]);

	$goenh->{$term}={};

	foreach my $e (@enh){
	    my $prefix=substr $e, 0, 3;
	    if($prefix eq "chr"){
		$e=substr $e, 3;
	    }

	    my @t=split(":", $e);

	    if(@t==3){
		$e=$t[0].":".$t[1]."-".$t[2];
	    }
	    
	    $goenh->{$term}{$e}=1;
	}
	
	$line=<$input>;
    }

    close($input);

}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script compares gene-enhancer associations between two methods.\n";
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

$parameters{"pathGOEnhancers1"}="NA";
$parameters{"pathGOEnhancers2"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGOEnhancers1", "pathGOEnhancers2", "pathOutput");

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

print "Reading GO-enhancer associations...\n";

my %goenh1;
readGOEnhancers($parameters{"pathGOEnhancers1"}, \%goenh1);

my $nbgo1=keys %goenh1;

print "There are ".$nbgo1." GO categories in the first file.\n";

my %goenh2;
readGOEnhancers($parameters{"pathGOEnhancers2"}, \%goenh2);

my $nbgo2=keys %goenh2;

print "There are ".$nbgo2." GO categories in the first file.\n";

print "Done.\n";

####################################################################################

print "Comparing associations and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GOTerm\tNbEnhancersMethod1\tNbEnhancersMethod2\tNbInCommon\n";

foreach my $go (keys %goenh1){
    my $nb1=keys %{$goenh1{$go}};

    if(exists $goenh2{$go}){
	my $nb2=keys %{$goenh2{$go}};

	my $nbc=0;

	foreach my $e (keys %{$goenh1{$go}}){
	    if(exists $goenh2{$go}{$e}){
		$nbc++;
	    }
	}

	print $output $go."\t".$nb1."\t".$nb2."\t".$nbc."\n";
    } else{
	print $output $go."\t".$nb1."\t0\t0\n";
    }
}

foreach my $go (keys %goenh2){
    my $nb2=keys %{$goenh2{$go}};

    if(!exists $goenh1{$go}){
	print $output $go."\t0\t".$nb2."\t0\n";
    }
}

close($output);

print "Done.\n";

####################################################################################
