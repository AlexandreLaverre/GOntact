use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
use strict;

####################################################################################
####################################################################################

sub readElementCoordinates{
    my $pathin=$_[0];
    my $okchromo=$_[1];
    my $elements=$_[2];

    open(my $input, $pathin);

    my $line=<$input>; ## no header

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[0];
	my $prefix=substr $chr, 0, 3;

	if($prefix eq "chr"){
	    $chr=substr $chr, 3;
	}

	if(exists $okchromo->{$chr}){
	    my $start=$s[1]+0;
	    my $end=$s[2]+0;
	    my $id=$s[3];

	    my $midpos=floor (($start+$end)/2); ## lower median position

	    $elements->{$id}={"chr"=>$chr, "start"=>[$midpos], "end"=>[$midpos]};
	}

	$line=<$input>;
    }

    close($input);
}

####################################################################################

sub readRegulatoryRegions{
    my $pathin=$_[0];
    my $regions=$_[1];
    my $okchromo=$_[2];

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

	my $prefix=substr $chr, 0, 3;

	if($prefix eq "chr"){
	    $chr=substr $chr, 3;
	}

	my $start=$s[$header{"region_start"}];
	my $end=$s[$header{"region_end"}];

	my $excludestart=$s[$header{"exclude_start"}];
	my $excludeend=$s[$header{"exclude_end"}];

	if(exists $regions->{$gene}){
	    $duplicated{$gene}=1;
	} else{

	    if($excludestart eq "NA"){
		$regions->{$gene}={"chr"=>$chr, "start"=>[$start], "end"=>[$end]};
	    } else{
		$excludestart=$s[$header{"exclude_start"}]+0;
		$excludeend=$s[$header{"exclude_end"}]+0;

		if($excludestart>$start){
		    if(exists $regions->{$gene}){
			push(@{$regions->{$gene}{"start"}}, $start);
			push(@{$regions->{$gene}{"end"}}, ($excludestart-1));
		    } else{
			$regions->{$gene}={"chr"=>$chr, "start"=>[$start], "end"=>[$excludestart-1]};
		    }
		}

		if($excludeend<$end){
		    if(exists $regions->{$gene}){
			push(@{$regions->{$gene}{"start"}}, ($excludeend+1));
			push(@{$regions->{$gene}{"end"}}, $end);
		    } else{
			$regions->{$gene}={"chr"=>$chr, "start"=>[$excludeend+1], "end"=>[$end]};
		    }
		}
	    }
	}

	$okchromo->{$chr}=1;

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

sub orderCoordinates{
    my $unordered=$_[0];
    my $ordered=$_[1];

    my %hashcoords;

    foreach my $id (keys %{$unordered}){
	my $chr=$unordered->{$id}{"chr"};

	my $nb=@{$unordered->{$id}{"start"}};

	for(my $i=0; $i<$nb; $i++){
	    my $start=${$unordered->{$id}{"start"}}[$i];
	    my $end=${$unordered->{$id}{"end"}}[$i];

	    if(exists $hashcoords{$chr}){
		if(exists $hashcoords{$chr}{$start}){
		    push(@{$hashcoords{$chr}{$start}{"end"}}, $end);
		    push(@{$hashcoords{$chr}{$start}{"id"}}, $id);
		} else{
		    $hashcoords{$chr}{$start}={"end"=>[$end], "id"=>[$id]};
		}
	    } else{
		$hashcoords{$chr}={$start=>{"end"=>[$end], "id"=>[$id]}};
	    }
	}
    }

    foreach my $chr (keys %hashcoords){
	$ordered->{$chr}={"start"=>[], "end"=>[], "id"=>[]};

	my @startpos=keys %{$hashcoords{$chr}};
	my @orderedstart=sort {$a<=>$b} @startpos;

	foreach my $start (@orderedstart){
	    my $nbpos=@{$hashcoords{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nbpos; $i++){
		my $end=${$hashcoords{$chr}{$start}{"end"}}[$i];
		my $id=${$hashcoords{$chr}{$start}{"id"}}[$i];

		push(@{$ordered->{$chr}{"start"}}, $start);
		push(@{$ordered->{$chr}{"end"}}, $end);
		push(@{$ordered->{$chr}{"id"}}, $id);
	    }
	}
    }
}

####################################################################################

sub overlapCoordinates{
    my $coords1=$_[0];
    my $coords2=$_[1];
    my $overlap=$_[2];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nb1=@{$coords1->{$chr}{"start"}};
	    my $nb2=@{$coords2->{$chr}{"start"}};

	    my $firstj=0;

	    for(my $i=0; $i<$nb1; $i++){
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];

		my $j=$firstj;

		while($j<$nb2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}

		$firstj=$j;

		while($j<$nb2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];

		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);

		    if($M<=$m){
			if(exists $overlap->{$id1}){
			    $overlap->{$id1}{$id2}=1;
			} else{
			    $overlap->{$id1}={$id2=>1};
			}
		    }

		    $j++;
		}
	    }
	}
    }
}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script computes gene ontology stats for a set of input genomic elements.\n";
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

$parameters{"pathInputElements"}="NA";
$parameters{"pathBackgroundElements"}="NA";
$parameters{"pathGOCategories"}="NA";
$parameters{"pathGOAnnotations"}="NA";
$parameters{"pathRegulatoryRegions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathInputElements", "pathBackgroundElements", "pathGOCategories", "pathGOAnnotations", "pathRegulatoryRegions", "pathOutput");

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

print "Reading regulatory regions...\n";

my %regions;
my %okchromo;
readRegulatoryRegions($parameters{"pathRegulatoryRegions"}, \%regions, \%okchromo);

my $nbgr=keys %regions;

print "There are ".$nbgr." genes with regulatory regions.\n";

my %orderedregions;

orderCoordinates(\%regions, \%orderedregions);

print "Done.\n";

my $nbchr=keys %okchromo;

print "We keep elements on ".$nbchr." standard chromosomes.\n";

####################################################################################

print "Reading genomic coordinates for input elements...\n";

my %elements;

readElementCoordinates($parameters{"pathInputElements"}, \%okchromo, \%elements);

my $nbel=keys %elements;

print "Found ".$nbel." elements.\n";

my %orderedelements;

orderCoordinates(\%elements, \%orderedelements);

print "Done.\n";

####################################################################################

print "Computing overlap between elements and regions...\n";

my %overlap;

overlapCoordinates(\%orderedelements, \%orderedregions, \%overlap);

my $nbov=keys %overlap;

print "There are ".$nbov." elements that overlap with regions.\n";

print "Done.\n";

####################################################################################

my $usingbg=0;
my $nbelbg=0;

my %orderedbgelements;
my %bgelements;

if($parameters{"pathBackgroundElements"} ne "NA"){
    $usingbg=1;
    
    print "Reading genomic coordinates for background elements...\n";
    
    
    readElementCoordinates($parameters{"pathBackgroundElements"}, \%okchromo, \%bgelements);
    
    $nbelbg=keys %bgelements;
    
    print "Found ".$nbelbg." background elements.\n";
    
    orderCoordinates(\%bgelements, \%orderedbgelements);
    
    print "Done.\n";
}

####################################################################################

print "Computing overlap between elements and regions...\n";

my %overlap;

overlapCoordinates(\%orderedelements, \%orderedregions, \%overlap);

my $nbov=keys %overlap;

print "There are ".$nbov." input elements that overlap with regions.\n";

my $nbovbg=0;
my %overlapbg;

if($usingbg==1){
    overlapCoordinates(\%orderedbgelements, \%orderedregions, \%overlapbg);
    
    $nbovbg=keys %overlapbg;
    
    print "There are ".$nbovbg." background elements that overlap with regions.\n";
}

print "Done.\n";

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

print "Computing observed values...\n";

my %elgo;

foreach my $idel (keys %overlap){
    foreach my $idgene (keys %{$overlap{$idel}}){
	if(exists $genego{$idgene}){
	    foreach my $go (@{$genego{$idgene}}){
		if(exists $elgo{$go}){
		    $elgo{$go}{$idel}=1;
		} else{
		    $elgo{$go}={$idel=>1};
		}
	    }
	}
    }
}

my %bgelgo;

if($usingbg==1){
    foreach my $idel (keys %overlapbg){
	foreach my $idgene (keys %{$overlapbg{$idel}}){
	    if(exists $genego{$idgene}){
		foreach my $go (@{$genego{$idgene}}){
		    if(exists $bgelgo{$go}){
			$bgelgo{$go}{$idel}=1;
		    } else{
			$bgelgo{$go}={$idel=>1};
		    }
		}
	    }
	}
    }
}

print "Done.\n";

####################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "#NbTotalInputElements\t".$nbel."\n";
print $output "#NbInputElementsInRegions\t".$nbov."\n";

if($usingbg==1){
    print $output "#NbTotalBackgroundElements\t".$nbelbg."\n";
    print $output "#NbBackgroundElementsInRegions\t".$nbovbg."\n";

    print $output "ID\tGOSpace\tNbAssociatedInputElements\tNbAssociatedBackgroundElements\n";
} else{
   print $output "ID\tGOSpace\tNbAssociatedInputElements\n";
}

foreach my $space (keys %gocat){
    foreach my $go (@{$gocat{$space}}){
	my $nbinputel=0;
	
	if(exists $elgo{$go}){
	    $nbinputel=keys %{$elgo{$go}};
	}

	if($usingbg==1){
	    
	    my $nbbgel=0;
	    
	    if(exists $bgelgo{$go}){
		$nbbgel=keys %{$bgelgo{$go}};
	    }
	    
	    print $output $go."\t".$space."\t".$nbinputel."\t".$nbbgel."\n";
	} else{
	    print $output $go."\t".$space."\t".$nbinputel."\n";
	}
    }
}

close($output);

print "Done.\n";

####################################################################################
