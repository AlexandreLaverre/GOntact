use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

####################################################################################
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
	my $start=$s[$header{"region_start"}];
	my $end=$s[$header{"region_end"}];

	if(exists $regions->{$gene}){
	    $duplicated{$gene}=1;
	} else{
	    $regions->{$gene}={"chr"=>$chr, "start"=>$start, "end"=>$end};
	    $okchromo->{$chr}=1;
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

    print "assigning regions to GO\n";
    
    foreach my $space (keys %{$gocat}){
	foreach my $cat (@{$gocat->{$space}}){
	    if(exists $gogene->{$cat}){
		foreach my $gene (@{$gogene->{$cat}}){
		    if(exists $regions->{$gene}){
			my $chr=$regions->{$gene}{"chr"};

			## we are only adding ungapped regions
			
			my $nbug=@{$regions->{$gene}{"ungappedstart"}};

			for(my $i=0; $i<$nbug; $i++){
			    my $start=${$regions->{$gene}{"ungappedstart"}}[$i];
			    my $end=${$regions->{$gene}{"ungappedend"}}[$i];
			    
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
    }

    print "computing total size\n";
    my $nbdone=0;
    
    foreach my $cat (keys %goregions){
	my %blocks;
	makeBlocks($goregions{$cat}, \%blocks);

	my $totalsize;
	my $nbr=@{$blocks{"chr"}};

	for(my $i=0; $i<$nbr; $i++){
	    my $start=${$blocks{"start"}}[$i];
	    my $end=${$blocks{"end"}}[$i];

	    $totalsize+=($end-$start+1);
	}


	$totalsizes->{$cat}=$totalsize;
    }    
}

####################################################################################

sub makeBlocks{
    my $coords=$_[0];
    my $blocks=$_[1];

    $blocks->{"chr"}=[];
    $blocks->{"start"}=[];
    $blocks->{"end"}=[];
    
    my %hashpos;

    my $nb=@{$coords->{"chr"}};

    for(my $i=0; $i<$nb; $i++){
	my $chr=${$coords->{"chr"}}[$i];
	my $start=${$coords->{"start"}}[$i];
	my $end=${$coords->{"end"}}[$i];
	
	if(exists $hashpos{$chr}){
	    if(exists $hashpos{$chr}{$start}){
		if($end>$hashpos{$chr}{$start}){
		    $hashpos{$chr}{$start}=$end;
		}
	    } else{
		$hashpos{$chr}{$start}=$end;
	    }
	} else{
	    $hashpos{$chr}={$start=>$end};
	}
    }

    foreach my $chr (keys %hashpos){
	my @uniquestart = keys %{$hashpos{$chr}};
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	my $nb=@sortedstart;
	
	my $currentstart=$sortedstart[0];
	my $currentend=$hashpos{$chr}{$currentstart};

	for(my $i=1; $i<$nb; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$hashpos{$chr}{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){  
		
		## we only change the end if it's larger than the current position
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$blocks->{"chr"}},$chr);
		push(@{$blocks->{"start"}},$currentstart);
		push(@{$blocks->{"end"}},$currentend);
		
		$currentstart=$thisstart;
		$currentend=$thisend;
		
	    }
	}

	## last block
	push(@{$blocks->{"chr"}},$chr);
	push(@{$blocks->{"start"}},$currentstart);
	push(@{$blocks->{"end"}},$currentend);
	
    }
}

#####################################################################################

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

####################################################################################

sub removeNs{
    my $genome=$_[0];
    my $regions=$_[1];
    my $okchromo=$_[2];
    my $stats=$_[3];

    my $totalN=0;
    my $totalL=0;
    
    foreach my $chr (keys %{$okchromo}){
	print "removing Ns for ".$chr."\n";
	
	if(!exists $genome->{$chr}){
	    print "Weird! cannot find genome sequence for ".$chr."\n";
	    exit(1);
	}
	
	my $seq=$genome->{$chr};
	my $l=length $seq;

	$totalL+=$l;
	
	my %hashN;

	for(my $i=0; $i<$l; $i++){
	    my $base=uc (substr $seq, $i, 1);

	    if($base eq "N"){
		$hashN{$i+1}=1;
	    }
	}

	my $nbN=keys %hashN;

	print "there are ".$nbN." N bases out of ".$l." on chr ".$chr."\n";

	$totalN+=$nbN;
	
	foreach my $gene (keys %{$regions}){
	    if($regions->{$gene}{"chr"} eq $chr){
		$regions->{$gene}{"ungappedstart"}=[];
		$regions->{$gene}{"ungappedend"}=[];

		my $rstart=$regions->{$gene}{"start"};
		my $rend=$regions->{$gene}{"end"};

		my $currentstart="NA";
		my $currentend="NA";

		for(my $i=$rstart; $i<=$rend; $i++){
		    if(!exists $hashN{$i}){
			if($currentstart eq "NA"){
			    $currentstart=$i;
			    $currentend=$i;
			} else{
			    if($i==($currentend+1)){
				$currentend=$i;
			    } else{
				push(@{$regions->{$gene}{"ungappedstart"}}, $currentstart);
				push(@{$regions->{$gene}{"ungappedend"}}, $currentend);

				$currentstart=$i;
				$currentend=$i;
			    }
			}
		    }
		}

		## last ungapped block
		
		if($currentstart ne "NA"){
		    push(@{$regions->{$gene}{"ungappedstart"}}, $currentstart);
		    push(@{$regions->{$gene}{"ungappedend"}}, $currentend);
		}
	    }
	}
    }

    $stats->{"totalL"}=$totalL;
    $stats->{"totalN"}=$totalN;
}

####################################################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts gene ontology stats for regulatory regions.\n";
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
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathRegulatoryRegions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGOCategories", "pathGOAnnotations", "pathGenomeSequence", "pathRegulatoryRegions", "pathOutput");

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
my %okchromo;

readRegulatoryRegions($parameters{"pathRegulatoryRegions"}, \%regions, \%okchromo);

my $nbgr=keys %regions;
my $nbchr=keys %okchromo;

print "There are ".$nbgr." genes with regulatory regions on ".$nbchr." chromosomes.\n";

print "Done.\n";

####################################################################################

print "Reading genome sequence...\n";

my %genome;

readFasta($parameters{"pathGenomeSequence"}, \%genome);

print "Done.\n";

####################################################################################

print "Removing Ns from regulatory regions...\n";

my %stats;
removeNs(\%genome, \%regions, \%okchromo, \%stats);

my $l=$stats{"totalL"};
my $n=$stats{"totalN"};

print "There are ".$l." bases in total on these ".$nbchr." chromosomes. There were ".$n." N bases.\n";

print "Done.\n";

####################################################################################

print "Computing expected values...\n";

my %totalsizes;

computeExpectedValues(\%gocat, \%gogene, \%regions, \%totalsizes);

print "Done.\n";

####################################################################################

if($nbgr>0){
    print "Writing output...\n";
    
    open(my $output, ">".$parameters{"pathOutput"});
    print $output "#TotalChromosomeLength\t".$l."\n";
    print $output "#NbNBases\t".$n."\n";
    
    print $output "ID\tGOSpace\tNbNonNBases\n";
    
    foreach my $space (keys %gocat){
	foreach my $go (@{$gocat{$space}}){
	    if(exists $totalsizes{$go}){
		print $output $go."\t".$space."\t".$totalsizes{$go}."\n";
	    }
	}
    }
    
    close($output);
    
    print "Done.\n";
}

####################################################################################
