#!/usr/bin/perl -w

# miR_homologs.pl -- given a chromosomal-sorted .bam file of query miRNAs aligned against a reference, examine for conformity to Meyers et al. criteria for MIRNAs, assuming the aligned query is the mature miRNA
# Type perldoc miR_homologs.pl for full description; see the License and full documentation after the __END__ symbol in this script, or see the README

use Getopt::Long;
use strict;

my $version = "0.2";

# usage statement
my $usage = "
$0 $version
USAGE: $0 [options] [alignments\.bam] [genome\.fasta]
OPTIONS:
--outdir [string] : Name of directory to be created to store results\.  Defaults to miR_homologs_[time] if not provided by user
--foldwindow [integer] : Size of window to fold\.  Window is centered on the alignment\.  Default: 300
--maxhpsep [integer] : maximum allowed separation of a base pair to span during hairpin prediction \(option -L for RNALfold\; must be between 50 and 2000\;default: 300\)
--minfracpaired [float] : minimum fraction of paired nts allowable in a hairpin structure\; default: 0.67
--minntspaired [integer] : minimum number of base-pairs in an accetable hairpin structure\; default: 30
";

# initial option defintions
my $outdir = '';
my $foldwindow = 300;
my $maxhpsep = 300;
my $minfracpaired = 0.67;
my $minntspaired = 30;

# ensure presence of samtools
# check for required installation of samtools and get version, or quit and complain
my ($samtools_version,$full_samtools_test_text) = get_samtools_version();
if($samtools_version eq "not found") {
    die "samtools not found\n$full_samtools_test_text\n\n$usage\n";
}

# get options from the command line
GetOptions ('outdir=s' => \$outdir,
	    'foldwindow=i' => \$foldwindow,
	    'maxhpsep=i' => \$maxhpsep,
	    'minfracpaired=f' => \$minfracpaired,
	    'minntspaired=i' => \$minntspaired);

# ensure presence and readability of the genome file and bam file pop'd off the end of the command
my $genome = pop @ARGV;
unless(-r $genome) {
    die "genome file $genome is not readable\n\nFATAL\n$usage\n";
}

# ensure bam file is readable
my $bamfile = pop @ARGV;
unless(-r $bamfile) {
    die "bamfile $bamfile not readable\.\n$usage\n";
}

# check that the options make some sense
unless(($foldwindow > 50) and ($foldwindow < 1000)) {
    die "FATAL: Option --foldwindow must be more than 50 and less than 1000\n$usage\n";
}
unless(($maxhpsep >= 50) and ($maxhpsep <= 2000)) {
    die "Option --maxhpsep must be a number between 50 and 2000\n\n$usage\n";
}
unless(($minfracpaired > 0) and ($minfracpaired <= 1)) {
    die "Option --minfracpaired must be a number greater than zero and less than or equal to one\n\n$usage\n";
}
unless(($minntspaired > 0) and ($minntspaired <= $maxhpsep)) {
    die "Option --minntspaired must be number greater than zero and less than or equal to option --maxhpsep\n\n$usage\n";
}

# Make the output directory, unless it already exists
unless($outdir) {
    my $time = time;
    $outdir = "miR_homologs" . "_$time";
}
if(-d $outdir) {
    die "FATAL: Output directory $outdir already exists\!  Pick another name\n$usage\n";
} else {
    system "mkdir $outdir";
}


# Then, check for presence of the .bai index file
my $expected_index = "$bamfile" . "\.bai";
unless(-e $expected_index) {
    die "FATAL: The provided \.bam file does not appear to have been indexed -- I could not find the expected index file $expected_index\nUse samtools to sort the data by chromosomal location, and then index the file, before using this program\n\n$usage\n";
}

# And also for the .fai index of the genome file
my $expected_faidx = "$genome" . "\.fai";
unless(-e $expected_faidx) {
    print STDERR "Expected genome index $expected_faidx for genome file $genome not found\.  Creating it using samtools faidx";
    system "samtools faidx $genome";
    print STDERR " done\n\n";
}

# get yur coordinates
my @clusters = get_miR_homo_clusters(\$expected_faidx,\$bamfile,\$foldwindow); 
# these entries are in the format of Chr:start-stop:strand:query_name:qleft-qright .. where strand is either "+" or "-"
# qleft is the left-end mapping coordinate, qright is the right-end coordinate

# get yur folding regions
my %to_fold = get_folding_regions_mh(\$genome,\@clusters);
# key = @cluster entr ..

# fold em

my %qualifying_hairpins = folder_mh(\%to_fold,\$maxhpsep,\$minntspaired,\$minfracpaired);  ## each entry has locus as key and an array of entries.  Each entry is tab-delim w/ brax, local-start, structure_details, deltaG, and strand
# key = locus, values in an anonymous array with each entry tab-delimited ... start, structure_details,deltaG

# Process Results
my %true_hps = hairpin_coords_mh(\%to_fold,\%qualifying_hairpins);  ## structure same as %qualifying hairpins

my %true_hps_trimmed = remove_redundant_hps_mh(\%true_hps); ## structure same as %qualifying hairpins
my %true_hps_3 = remove_nonoverlapped_hps_mh(\%true_hps_trimmed);  ## structure same as %qualifying hairpins

my %hp_clusters = get_hp_clusters_mh(\%true_hps_3,\$genome);


# Final analysis -- candidate mature miRNAs and their putative star seqences. Output if pass
my @ok = mir_output (\%hp_clusters,\$genome,\$bamfile,\$outdir);

# Overlapping loci .. if any overlap and on same strand
my %overlap = check_overlap(\@ok);

my $n_large = scalar @ok;
my $n_redundant; ## all those listed as values in the overlap hash are redundant
my $loc;

my %redundant = ();
while(($loc) = each %overlap) {
    my @reds = @{$overlap{$loc}};
    foreach my $r (@reds) {
	++$n_redundant;
	$redundant{$r} = 1;
    }
}

# open a Results.txt file
my $res_file = "$outdir" . "\/" . "Results\.txt";
open(RES, ">$res_file");
# print a header
print RES "\#Locus\tName\tHP\tStrand\tOther_names\n";

print STDERR "\nFound $n_large miRNAs\n";
if($n_redundant) {
    print STDERR "Some are overlapping, and likely represent the same MIRNA gene aligning to more than one qualifying query\n";
    my $true_n = $n_large - $n_redundant;
    print STDERR "After accounting for overlapping loci, there are a total of $true_n MIRNA loci found\n\n";
}
foreach $loc (@ok) {
    if(exists($redundant{$loc})) {
	next;
    }
    if($loc =~ /^(\S+):(\d+)-(\d+):(\S):(\S+):/) {
	my $shortloc = "$1" . ":" . "$2" . "-" . "$3";
	print RES "$shortloc\t$5\tMIRNA\t$4\t";
    } else {
	die "FATAL: in main program failed to parse locus name $loc\n";
    }
    if(exists($overlap{$loc})) {
	my @reds = @{$overlap{$loc}};
	foreach my $r (@reds) {
	    if($r =~ /^(\S+):(\d+)-(\d+):(\S):(\S+):/) {
		print RES "$5";
		print RES ",";
	    } else {
		die "FATAL: in main program failed to parse locus name $r\n";
	    }
	}
	print RES "\n";
    } else {
	print RES "\.\n";
    }
}
close RES;


##### here be sub-routines
sub get_samtools_version {
    my $version = "not found";
    my $full_text;
    open(MESSAGE, "samtools 2>&1 |");
    while (<MESSAGE>) {
	$full_text .= $_;
	if($_ =~ /^Version/) {
	    $version = $_;
	    last;
	}
    }
    close MESSAGE;
    return ($version,$full_text);
}

sub get_miR_homo_clusters {
    my($faidx,$bamfile,$winsize) = @_;  ## references, all scalars
    my @clusters = ();
    my $last_one = "null";
    my @fields = ();
    my $query_length;
    my $right_end;
    my $center;
    my $start;
    my $stop;
    my @faidx_fields = ();
    my $chr_length;
    my $strand;
    open(SAM, "samtools view $$bamfile |");
    while (<SAM>) {
	chomp;
	# ignore headers
	if($_ =~ /^@/) {
	    next;
	}
	@fields = split ("\t", $_);
	# ignore unmapped queries
	if($fields[1] & 4) {
	    next;
	}
	# get the length of the query from the CIGAR string
	$query_length = parse_cigar ($fields[5]);
	
	# get the right end
	$right_end = $fields[3] + $query_length - 1;
	
	# get the strand
	if($fields[1] & 16) {
	    $strand = "-";
	} else {
	    $strand = "+";
	}
	
	# get the nominal window, centered on the middle of the alignment
	$center = $fields[3] + (int (0.5 * $query_length));
	$start = $center - (int(0.5*$$winsize));
	$stop = $center + (int(0.5*$$winsize));
	
	# ensure start and stop are still within the chromosome
	if($start < 1) {
	    $start = 1;
	}
	open(FAIDX,"$$faidx");
	while(<FAIDX>) {
	    chomp;
	    @faidx_fields = split ("\t", $_);
	    if($faidx_fields[0] eq $fields[2]) {
		$chr_length = $faidx_fields[1];
		last;
	    }
	}
	close FAIDX;
	if($stop > $chr_length) {
	    $stop = $chr_length;
	}
	
	# build entry
	$last_one = "$fields[2]" . ":" . "$start" . "-" . "$stop" . ":" . "$strand" . ":" . "$fields[0]" . ":" . "$fields[3]" . "-" . "$right_end";
	push(@clusters,$last_one);
    }
    close SAM;
    return @clusters;
}
	
	
sub parse_cigar {
    my($cigar) = @_;
    if($cigar eq "\*") {
	die "FATAL in sub-routine \'parse_cigar\' : The CIGAR string is missing from at least one of the mappings, so I can't caluclate read length\!\n";
    }
    # per the SAM specification, "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ."
    my $read_length = 0;
    while($cigar =~ /(\d+)M/g) {
	$read_length += $1;
    }
    while($cigar =~ /(\d+)I/g) {
	$read_length += $1;
    }
    while($cigar =~ /(\d+)S/g) {
	$read_length += $1;
    }
    while($cigar =~ /(\d+)\=/g) {
	$read_length += $1;
    }
    while($cigar =~ /(\d+)X/g) {
	$read_length += $1;
    }
    return $read_length;
}

sub get_folding_regions_mh {
    my($gen_file,$info) = @_;  ## passed as references .. scalar, array
    my $locus;

    my %to_fold = ();
    my $chr;
    my $strand;
    my $start;
    my $stop;
    my $locus_size;

    my $get_size;

    my $middle;
    my $get_start;
    my $get_stop;
    my $gen_line;
    my $chr_seq;
    my $ok;
    my $subseq;

    my $reported_seq;
    # for progress tracking
    my $n_to_get = scalar (@$info);
    my $five_percent = int (0.05 * $n_to_get);
    my $x;
    print STDERR "\tProgress in sub-routine \"get_folding_regions_mh\" \(dot = five percent\): ";
    
    foreach $locus (@$info) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	# parse locus name
	if($locus =~ /^(\S+):(\d+)-(\d+):(\S):/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	    $strand = $4;
	} else {
	    die "Fatal in sub-routine get_folding_regions: could not parse locus name $locus\n";
	}
	
	# call samtools faidx to get the FASTA formatted version of the whole thing
	$subseq = '';
	my $sam_coords = "$chr" . ":" . "$start" . "-" . "$stop";
	open(FAIDX, "samtools faidx $$gen_file $sam_coords |");
	# parse the FASTA output to get the sequence as a single string
	while (<FAIDX>) {
	    chomp;
	    if($_ =~ /^>/) {
		next;
	    }
	    $_ =~ s/\s//g;
	    # ensure upper case
	    my $fa_line = uc ($_);
	    $subseq .= $fa_line;
	}
	close FAIDX;

	# convert it to RNA
	$subseq =~ s/T/U/g;
	
	# if strand is -, reverse comp it
	if($strand eq "-") {
	    $reported_seq = reverse $subseq;
	    $reported_seq =~ tr/AUCG/UAGC/;
	} else {
	    $reported_seq = $subseq;
	}
	
	## add to the to_fold hash
	$to_fold{$locus} = $reported_seq;
	
    }
    # finish progress tracking
    print STDERR " Done\n";
    # return the hash
    return %to_fold;
}

sub folder_mh {
    my($to_fold,$L,$min_paired,$min_frac_paired) = @_;  ## passed by reference, hash and scalar, respectively
    my $locus;
    my $fold_seq;
    my $brax;
    #my $delta_G;  # deprecated as of v. 0.2
    my %output = ();  
    my $structure_details;
    my $start;
    my $entry;  ## tab-delimited .. brax, start, helix_info (e.g. 123-150,180-200), strand "Watson" or "Crick"
    my $revcomp;
    
    my $adj_st;
    my %regions = ();
    my $brax_section;
    my $left;
    my $right;

    # first, get some information to enable a crude progress bar
    my $n_loci_to_fold = scalar ( keys %$to_fold);
    my $x = 0;
    my $five_percent = int(0.05 * $n_loci_to_fold);
    print STDERR "\n\tProgress in sub-routine \"folder_mh\" \(dot = 5 percent\): ";
    
    while(($locus) = each %$to_fold) {
	# progress tracking
	++$x;
	if($x == $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	
	open(RNALFOLD, "echo $$to_fold{$locus} | RNALfold -d 2 -noLP -L $$L |");
	while (<RNALFOLD>) {
	    chomp;
	    if($_ =~ /^[^\.\(\)]/) {
		## this is a line that does not have a structure .. most likely the last line of output which just has the input sequence .. ignore it
		next;
	    }

	    
	    if($_ =~ /^\S+/) {
		$brax = $&;
	    } else {
		die "FATAL in sub-routine \'folder\' : failed to parse brackets from RNALfold output line:\n$_\n";
	    }
	    
	 # delta G not tracked as of v 0.2
	    #if($_ =~ /\s+\((.*\d+.*)\)/) {
	#	$delta_G = $1;
	#	$delta_G =~ s/\s//g;
	#    } else {
	#	die "FATAL in sub-routine \'folder\' : failed to parse deltaG from RNALfold output line: $_\n";
	#    }
	    
	    if($_ =~ /(\d+)\s*$/) {
		$start = $1;
	    } else {
		die "FATAL in sub-routine \'folder\' : failed to parse start position from RNALfold output line $_\n";
	    }

	    %regions = get_distinct_hps($brax);
	    while(($left, $right) = each %regions) {
		$brax_section = substr($brax,($left - 1),($right-$left+1));
		
		# send the structure to the general structure evaluation sub-routine, which returns zero to reject, one to keep
		$structure_details = evaluate_structure_general($brax_section,$$min_paired,$$min_frac_paired);
		
		# if it is OK, adjust the coordinates, and add to output
		unless($structure_details eq "bogus") {
		    # structure details are 123-150,180-200 .. e.g. start and stop positions of the helix of interest, one-based coordinates, relative to the brackets themselves.
		    $adj_st = $left + $start - 1;
		    $entry = "$brax_section\t$adj_st\t$structure_details";
		    push(@{$output{$locus}}, $entry);
		}
	    }
	}
	close RNALFOLD;
    }
    # close progress
    print STDERR " Done\n";
    return %output;
}

sub evaluate_structure_general {
    ## this is the primary evaluation for putative hairpin secondary structures
    ## input is a scalar, dot-bracket string.  Other two inputs also scalar
    my($brax,$min_paired,$min_paired_frac) = @_;
    my @chars = split ('', $brax);
    my $char;
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    my $i;
    my $left_true_stop;
    my $right_true_start;
    my $pairs;
    my $frac_paired;
    my $to_return;  ## left_start-left_true_stop,right_true_start-right_stop
    
    # make a hash where key is 1-based positions of ( and value is one-based position of corresponding )
    my %left_right = get_left_right($brax);
    
    # find the first run of "(" paired residues
    $i = 0;
    foreach $char (@chars) {
	++$i;
	if($char eq "\(") {
	    unless($left_start) {
		$left_start = $i;
	    }
	    $left_stop = $i;  ## simply updates each time
	} elsif ($char eq "\)") {
	    last;  ## thus, only the first 'left' run will be tracked
	}
    }
    
    # now, find the last run of "(" paired residues
    $i = 0;
    foreach $char (@chars) {
	++$i;
	if($char eq "\(") {
	    $right_start = '';
	    $right_stop = '';
	} elsif ($char eq "\)") {
	    unless($right_start) {
		$right_start = $i;
	    }
	    $right_stop = $i;
	}
    }
    
    # find the helical region shared by the first '(' block and the last ')' block.
    # at the same time, track the number of paired positions in the bottom helix
    
    # first, ensure there actually is a helical region at all
    if(($i > 0) and ($left_stop)) {
    
	for($i = $left_start; $i <= $left_stop; ++$i) {
	    if(exists($left_right{$i})) {
		if(($left_right{$i} >= $right_start) and
		   ($left_right{$i} <= $right_stop)) {
		    ++$pairs;
		    $left_true_stop = $i;
		    $right_true_start = $left_right{$i};
		}
	    }
	}
	# Note that structures can have multiple, indpendent helices... these will have uninitialized values
	#  for $pairs, $right_true_start and $left_true_stop ... if this is the case cease analysis of the structure
	if(($pairs) and ($right_true_start) and ($left_true_stop)) {
	    
	    # calculate the fraction of nts in the bottom stem that are paired
	    $frac_paired = (2 * $pairs) / (($right_stop - $right_true_start + 1) + ($left_true_stop - $left_start + 1));
	    
	    # For reference, here are the stats for Arabidopsis thaliana miRNAs (miRBase 18).
	    # Note that 13 MIRNAs were excluded from analysis because of short side-branches at the ends. . untrimmed hps in miRBase
	    # Stat       Number of pairs     fraction paired
	    # mean        62.37769784         0.748874515
	    # stdev       37.36382455         0.137359546
	    # 10th p-tile 32                  0.600252207
	    # 25th p-tile 39.25               0.714555596
	    # 50th p-tile 50                  0.77694859
	    # 75th p-tile 70.75               0.829021927
	    # 90th p-tile 114.3               0.879151515
	    
	    # min       22                  0.171990172
	    # max       213                 1
	    
	    # n         278
	    
	    if(($frac_paired >= $min_paired_frac) and
	       ($pairs >= $min_paired)) {
		$to_return = "$left_start" . "-$left_true_stop" . "," . "$right_true_start" . "-" . "$right_stop";
	    } else {
		$to_return = "bogus";
	    }
	} else {
	    $to_return = "bogus";
	}
    } else {
	$to_return = "bogus";
    }
    return $to_return;
}

sub get_left_right {
    my($brax) = @_;
    my @chars = split('', $brax);
    my $char;
    my @left = ();
    my %hash = ();
    my $i = 0;
    my $key;

    foreach $char (@chars) {
	++$i;
	if($char eq "\(") {
	    push(@left,$i);
	} elsif ($char eq "\)") {
	    $key = pop @left;
	    $hash{$key} = $i;
	}
    }
    
    return %hash;
}

sub hairpin_coords_mh {
    my($to_fold,$input_hps) = @_;  ## passed by reference, both hashes
    my %true_hps;
    my @inputs = ();
    my @fields;
    my $locus;
    my $in_helix_coords;
    my $local_offset;
    my $brax;
    my $inp;
    #my $delta_G;
    ## coordinates within the initial query region.
    my $left_1_start;
    my $left_1_stop;
    my $right_1_start;
    my $right_1_stop;
    my $start_1;
    my $stop_1;
    ##    
    my $start_true;  ## start_true and stop_true are chromosomal.  start is the lower number for Watson-strand RNAs, higher for Crick - strand RNAs
    my $stop_true;
    my $left_true_start;  ## the 5' arm of the helix.  start < stop for Watson, stop > start for Crick
    my $left_true_stop;
    my $right_true_start;
    my $right_true_stop;  ## the 3' arm of the helix.  start < stop for Watson, stop > start for Crick
    ##
    my $brax_true;  ##  format: $start_true-$stop_true .. start < stop for Watson, start > stop for Crick
    my $helix_true;  ## format: $left_true_start-$left_true_stop,$right_true_start-$right_true_stop
    my $true_entry;  ## format: $brax\t$brax_true\thelix_true\t$delta_G
    
    my $strand_folded;
    my $folded_start;
    my $folded_stop;
    while(($locus) = each %$input_hps) {
	# parse the locus name to see which strand was folded
	if($locus =~ /^(\S+):(\d+)-(\d+):(\S):/) {
	    if($4 eq "+") {
		$strand_folded = "Watson";
	    } elsif ($4 eq "-") {
		$strand_folded = "Crick";
	    }
	    $folded_start = $2;
	    $folded_stop = $3;
	} 
	
	@inputs = @{$$input_hps{$locus}};
	foreach $inp (@inputs) {
	    @fields = split ("\t", $inp);
	    #$delta_G = pop @fields;
	    
	    ## HERE
	    $in_helix_coords = pop @fields;
	    $local_offset = pop @fields;
	    $brax = pop @fields;
	    
	    if($in_helix_coords =~ /^(\d+)-(\d+),(\d+)-(\d+)$/) {
		$left_1_start = $local_offset + $1 - 1;
		$left_1_stop = $local_offset + $2 - 1;
		$right_1_start = $local_offset + $3 -1;
		$right_1_stop = $local_offset + $4 - 1;
		$start_1 = $local_offset;
		$stop_1 = $start_1 + (length $brax) - 1;
	    } else {
		die "Fatal error in sub-routine hairpin_coords: failed to parse local hairpin coordinates $in_helix_coords\n";
	    }
	    
	    ## true coordinates depend on the polarity which was folded
	    if($strand_folded eq "Watson") {
		$start_true = $folded_start + $start_1 - 1;
		$stop_true = $folded_start + $stop_1 - 1;
		$left_true_start = $folded_start + $left_1_start - 1;
		$left_true_stop = $folded_start + $left_1_stop - 1;
		$right_true_start = $folded_start + $right_1_start - 1;
		$right_true_stop = $folded_start + $right_1_stop - 1;
	    } elsif ($strand_folded eq "Crick") {
		$start_true = $folded_stop - $start_1 + 1;
		$stop_true = $folded_stop - $stop_1 + 1;
		$left_true_start = $folded_stop - $left_1_start + 1;
		$left_true_stop = $folded_stop - $left_1_stop + 1;
		$right_true_start = $folded_stop - $right_1_start +1;
		$right_true_stop = $folded_stop - $right_1_stop + 1;
	    }
	    
	    $brax_true = "$start_true" . "-" . "$stop_true";
	    $helix_true = "$left_true_start" . "-" . "$left_true_stop" . "," . "$right_true_start" . "-" . "$right_true_stop";
	    #$true_entry = "$brax\t$brax_true\t$helix_true\t$delta_G\t$strand_folded";
	    $true_entry = "$brax\t$brax_true\t$helix_true\t$strand_folded";
	    
	    ## TEST
#	    print STDERR "\tOUTPUT: $true_entry\n";
	    ##
	    
	    push(@{$true_hps{$locus}}, $true_entry);
	    
	}
    }
    return %true_hps;
}
sub remove_redundant_hps_mh {
    my($input_hps) = @_;  ## the %true_hps hash, passed by reference
    my $locus;
    my %output_hps;  ## the cleaned version to be returned
    my @entries = ();
    my $entry;
    my $n = 0;
    my @fields = ();
    my $helix_coords;
    my $hp_length;
    my @tmpW = ();
    my @tmpC = ();
    my $tmp_entry;
    my %to_deleteW = ();
    my %to_deleteC = ();
    my $x;

    my $delta_G_per_nt;
    
    my @q_fields = ();
    my @q_ranges = ();
    my @q_5p_range = ();
    my @q_3p_range = ();
    
    my $y;
    my $sub_entry;

    my @s_fields = ();
    my @s_ranges = ();
    my @s_5p_range = ();
    my @s_3p_range = ();
    
    my $five_p_overlap;
    my $three_p_overlap;
    
    my $strand;
    
    my @entriesW = ();
    my @entriesC = ();

    my $q_size;
    my $s_size;

    while(($locus) = each %$input_hps) {
	if($locus =~ /^(\S+):(\d+)-(\d+):(\S):/) {
	    if($4 eq "+") {
		$strand = "Watson";
	    } elsif ($4 eq "-") {
		$strand = "Crick";
	    }
	} else {
	    die "FATAL in sub-routine remove_redundant_hps_mh : could not parse locus $locus to find the strand\n";
	}

	@tmpW = ();  ## reset at each locus
	%to_deleteW = ();  ## reset at each locus
	@tmpC = ();  ## reset at each locus
	%to_deleteC = ();  ## reset at each locus
	
	@entriesW = ();
	@entriesC = ();
	
	
	@entries = @{$$input_hps{$locus}};
	foreach $entry (@entries) {
	    
	    # we must discriminate between Watson and Crick derived hairpins
	    
	    @fields = split ("\t", $entry);
	    $helix_coords = $fields[2];
	    $hp_length = length $fields[0];  ## the brackets
	    #$delta_G_per_nt = $fields[3] / $hp_length;
	    $tmp_entry = "$helix_coords";
	 #   $strand = $fields[4];
	    if($strand eq "Watson") {
		push(@tmpW, $tmp_entry);
		push(@entriesW, $entry);
	    } else {
		push(@tmpC, $tmp_entry);
		push(@entriesC, $entry);
	    }
	}
	# now analyze the Watson strand ones
	for($x = 0; $x < (scalar @tmpW); ++$x) {
	    if(exists($to_deleteW{$x})) {
		next;
	    }
	    $tmp_entry = $tmpW[$x];
	    #@q_fields = split ("\t", $tmp_entry);

	    @q_ranges = split (",", $tmp_entry);
	    ## q_ranges[0] is the 5p, [1] is the 3p
	    @q_5p_range = split ("-", $q_ranges[0]);
	    @q_3p_range = split ("-", $q_ranges[1]);
	    # in the above, [0] is the start, [1] is the stop.
	    $q_size = abs($q_5p_range[0] - $q_3p_range[1]);
	    
	    # go through all pairwise combinations
	    for($y = 0; $y < (scalar @tmpW); ++$y) {
		if($x == $y) {
		    next;  ## don't compare the same entry!
		}
		if(exists($to_deleteW{$y})) {
		    next;  ## don't both with entries already on the delete list
		}
		$sub_entry = $tmpW[$y];
		#@s_fields = split ("\t", $sub_entry);

		@s_ranges = split (",", $sub_entry);
		@s_5p_range = split ("-", $s_ranges[0]);
		@s_3p_range = split ("-", $s_ranges[1]);
		$s_size = abs($s_5p_range[0] - $s_3p_range[1]);
		
		# is there overlap in both the 5p and 3 ranges?

		$five_p_overlap = range_overlap(\@q_5p_range,\@s_5p_range);
		$three_p_overlap = range_overlap(\@q_3p_range,\@s_3p_range);
		## zero returned for no overlap, 1 for overlap
		if(($five_p_overlap) and ($three_p_overlap)) {
		    # we will delete the shorter one
		    if($q_size <= $s_size) {
			$to_deleteW{$x} = 1;
		    } else {
			$to_deleteW{$y} = 1;
		    }
		}

	    }
	}
	## add the surviving hairpins to the output hash
	for($x = 0; $x < (scalar @tmpW); ++$x) {
	    unless(exists($to_deleteW{$x})) {
		push (@{$output_hps{$locus}}, $entriesW[$x]);
	    }
	}
	
	## Now for the Crick-strand hairpins

	for($x = 0; $x < (scalar @tmpC); ++$x) {
	    if(exists($to_deleteC{$x})) {
		next;
	    }
	    $tmp_entry = $tmpC[$x];
	    @q_fields = split ("\t", $tmp_entry);
	    ## q_fields[1] is the query delta_g per nt
	    @q_ranges = split (",", $q_fields[0]);
	    ## q_ranges[0] is the 5p, [1] is the 3p
	    @q_5p_range = split ("-", $q_ranges[0]);
	    @q_3p_range = split ("-", $q_ranges[1]);
	    # in the above, [0] is the start, [1] is the stop.
	    $q_size = abs($q_5p_range[0] - $q_3p_range[1]);
    
	    # go through all pairwise combinations
	    for($y = 0; $y < (scalar @tmpC); ++$y) {
		if($x == $y) {
		    next;  ## don't compare the same entry!
		}
		if(exists($to_deleteC{$y})) {
		    next;  ## don't both with entries already on the delete list
		}
		$sub_entry = $tmpC[$y];
		@s_fields = split ("\t", $sub_entry);
		## s_fields[1] is the subject deltaG per nt
		@s_ranges = split (",", $s_fields[0]);
		@s_5p_range = split ("-", $s_ranges[0]);
		@s_3p_range = split ("-", $s_ranges[1]);
		$s_size = abs($s_5p_range[0] - $s_3p_range[1]);
		
		# is there overlap in both the 5p and 3 ranges?

		$five_p_overlap = range_overlap(\@q_5p_range,\@s_5p_range);
		$three_p_overlap = range_overlap(\@q_3p_range,\@s_3p_range);
		## zero returned for no overlap, 1 for overlap
		if(($five_p_overlap) and ($three_p_overlap)) {
		    # delete the shorter of the two
		    if($q_size <= $s_size) {
			$to_deleteC{$x} = 1;
		    } else {
			$to_deleteC{$y} = 1;
		    }
		}
	    }
	}
	## add the surviving hairpins to the output hash
	for($x = 0; $x < (scalar @tmpC); ++$x) {
	    unless(exists($to_deleteC{$x})) {
		push (@{$output_hps{$locus}}, $entriesC[$x]);
	    }
	}
    }
    return %output_hps;
}

sub remove_nonoverlapped_hps_mh {
    my($input_hps) = @_;  ## passed by reference
    my $locus;
    my %output_hps;  ## the cleaned version to be returned
    my $strand;
    my @entries = ();
    my $entry;
    my $qleft;
    my $qright;
    my @fields = ();
    my $helix_coords;
    my $ok;
    my $lowfive;
    my $highfive;
    my $lowthree;
    my $highthree;
    while(($locus) = each %$input_hps) {
	# parse the strand
	if($locus =~ /^(\S+):(\d+)-(\d+):(\S):(\S+):(\d+)-(\d+)$/) {
	    if($4 eq "+") {
		$strand = "Watson";
	    } elsif($4 eq "-") {
		$strand = "Crick";
	    }
	    $qleft = $6;
	    $qright = $7;
	}
	@entries = @{$$input_hps{$locus}};
	foreach $entry (@entries) {
	    $ok = 0;
	    @fields = split ("\t", $entry);
	    $ok = 0;
	    if($fields[2] =~ /(\d+)-(\d+)\,(\d+)-(\d+)/) {  ## (\d+)-(\d+),(\d+)-(\d+)
		if($2 > $1) {
		    $lowfive = $1;
		    $highfive = $2;
		    $lowthree = $3;
		    $highthree = $4;
		} else {
		    $lowfive = $2;
		    $highfive = $1;
		    $lowthree = $4;
		    $highthree = $3;
		}
		if(($qleft >= $lowfive) and
		   ($qleft <= $highfive)) {
		    $ok = 1;
		}
		if(($qleft >= $lowthree) and
		   ($qleft <= $highthree)) {
		    $ok = 1;
		}
		if(($qright >= $lowfive) and
		   ($qright <= $highfive)) {
		    $ok = 1;
		}
		if(($qright >= $lowthree) and
		   ($qright <= $highthree)) {
		    $ok = 1;
		}
	    }
	    if($ok) {
		push(@{$output_hps{$locus}}, $entry);
	    }
	}
    }
    return %output_hps;
}

sub range_overlap {
    my($q_range,$s_range) = @_;  ## passed by reference .. ranges are arrays, number-number
    my $answer = 0;
    my @q_sorted = sort {$a <=> $b} @$q_range;
    my @s_sorted = sort {$a <=> $b} @$s_range;  ## sorting is necessary to deal with minus-strand hairpins
    my $x;
    for($x = $q_sorted[0]; $x <= $q_sorted[1]; ++$x) {
	if(($x >= $s_sorted[0]) and
	   ($x <= $s_sorted[1])) {
	    $answer = 1;
	}
    }
    return $answer;
}
    

sub range_overlap_count {
    my($q_range,$s_range) = @_;  ## passed by reference .. ranges are arrays, number-number
    my $answer = 0;
    my @q_sorted = sort {$a <=> $b} @$q_range;
    my @s_sorted = sort {$a <=> $b} @$s_range;  ## sorting is necessary to deal with minus-strand hairpins
    my $x;
    for($x = $q_sorted[0]; $x <= $q_sorted[1]; ++$x) {
	if(($x >= $s_sorted[0]) and
	   ($x <= $s_sorted[1])) {
	    ++$answer;
	}
    }
    return $answer;
}

sub mir_output {
    my ($hp_hash,$genome,$bamfile,$out_dir) = @_;  ## passed by reference.  hash, scalar, scalar, scalar, scalar, hash, hash
    my $outfile;
    my @output = ();
    my $locus;  ## the padded locus -- keys from the hp_hash
    my $locus_data;
    my @l_data_fields = ();
    
    my $brax;
    my @brax_coords;
    my $strand;
    
    my %hp_seq_hash = ();
    my %brax_hash = ();
    my @hp_letters = ();
    my @brax_chars = ();

    my $i;
    my $j;
    
    my $chr;
    my $subseq;
    my $loc_start;
    my $loc_stop;
    my $displayseq;

    my %sense = ();
    my %antisense = ();
    my @samfields = ();
    my $read_stop;
    my $read_length;
    my $mapping_polarity;
    my $map_coord;
    my $key;
    my @sense_names = ();
    my @antisense_names = ();
    my $total_mappings = 0;
    
    my @candidates = ();
    
    my @mc_limits = ();
    
    my $candidate_brax;
    my $candidate_brax_chopped;
    
    my @candidates_two = ();
    
    my $fail;
    my $n_unpaired_cand;
    my $n_left_paired_cand;
    my $n_right_paired_cand;
    
    my $star_coord;
    
    my @star_limits = ();
    my $star_brax;
    my $star_brax_chopped;
    my $n_unpaired_star;
    my $n_left_paired_star;
    my $n_right_paired_star;
    
    my %miRNA = ();
    my %miRNA_star = ();
    
    my @successes = ();  # the query names, reset for each locus
    
    ## for progress bar
    my $n_to_process = scalar ( keys %$hp_hash);
    my $five_percent = int ($n_to_process / 20);
    my $progress = 0;
    print STDERR "\tProgress in sub-routine \"mir_output\" \(dot = five percent\): ";
    
    
    my $query_name;
    my %qnames = ();
    my $qstart;
    my $qstop;
    while(($locus,$locus_data) = each %$hp_hash) {
	# progress tracking
	++$progress;
	if($progress >= $five_percent) {
	    print STDERR "\.";
	    $progress = 0;
	}
        # parse the data
	@l_data_fields = split ("\t", $locus_data);
	$brax = $l_data_fields[0];
	@brax_coords = split ("-",$l_data_fields[1]);
	# sort so that the lower coordinate is always in the [0]
	@brax_coords = sort {$a <=> $b} @brax_coords;
	#$strand = $l_data_fields[4];  ## either 'Watson' or 'Crick'
	
	# retrieve the genome sequence corresponding to the padded hp / the whole cluster
	if($locus =~ /^(\S+):(\d+)-(\d+):(\S):(\S+):(\d+)-(\d+)/) {
	    $chr = $1;
	    $loc_start = $2;
	    $loc_stop = $3;
	    if($4 eq "+") {
		$strand = "Watson";
	    } elsif ($4 eq "-") {
		$strand = "Crick";
	    }
	    $query_name = $5;
	    $qstart = $6;
	    $qstop = $7;
	} else {
	    die "FATAL: in sub-routine \"hp_output\" could parse chr name from $locus\n";
	}
	# call samtools faidx
	$subseq = '';  ## reset
	my $sam_coords = "$chr" . ":" . "$loc_start" . "-" . "$loc_stop";
	open(FAIDX, "samtools faidx $$genome $sam_coords |");
	while (<FAIDX>) {
	    chomp;
	    if($_ =~ /^>/) {
		next;
	    }
	    $_ =~ s/\s//g;
	    # ensure upper case
	    my $fa_line = uc ($_);
	    $subseq .= $fa_line;
	}
	close FAIDX;
	
	$subseq =~ s/T/U/g;
	if($strand eq "Crick") {
	    $displayseq = reverse $subseq;
	    $displayseq =~ tr/ACUG/UGAC/;
	} else {
	    $displayseq = $subseq;
	}
	
	
	# code the padded hp sequence and the brackets into hashes, keyed by coordinates, for easy retrieval
	@hp_letters = split ('', $displayseq);
	%hp_seq_hash = ();
	$j = 0;
	if($strand eq "Watson") {
	    for($i = $loc_start; $i <= $loc_stop; ++$i) {
		$hp_seq_hash{$i} = $hp_letters[$j];
		++$j;
	    }
	} elsif ($strand eq "Crick") {
	    for($i = $loc_stop; $i >= $loc_start; --$i) {
		$hp_seq_hash{$i} = $hp_letters[$j];
		++$j;
	    }
	}
	
	@brax_chars = split ('', $brax);
	%brax_hash = ();
	$j = 0;
	if($strand eq "Watson") {
	    for($i = $brax_coords[0]; $i <= $brax_coords[1]; ++$i) {
		$brax_hash{$i} = $brax_chars[$j];
		++$j;
	    }
	} elsif ($strand eq "Crick") {
	    for($i = $brax_coords[1]; $i >= $brax_coords[0]; --$i) {
		$brax_hash{$i} = $brax_chars[$j];
		++$j;
	    }
	}
	# check whether candidates span a loop and if not whether they have four or fewer un-paired residues (ignoring the last two nts)
	
	# our candidate is noted in the locus name itself and captured already in $qstart and $qstop
	if($strand eq "Watson") {
	    $map_coord = "$qstart" . "-" . "$qstop";
	    $key = $qstart;
	} elsif ($strand eq "Crick") {
	    $map_coord = "$qstop" . "-" . "$qstart";
	    $key = $qstop;
	}

	$candidate_brax = '';
	$fail = 0;
	$n_unpaired_cand = 0;
	$n_left_paired_cand = 0;
	$n_right_paired_cand = 0;
	@mc_limits = split ("-", $map_coord);

	if($mc_limits[0] < $mc_limits[1]) {
	    for($i = $mc_limits[0]; $i <= $mc_limits[1]; ++$i) {
		if(exists($brax_hash{$i})) {
		    $candidate_brax .= $brax_hash{$i};
		} else {
		    # if the candidate is partially within the padded region, consider all parts of the padded region unpaired
		    $candidate_brax .= "\.";
		}
	    }
	} elsif ($mc_limits[0] > $mc_limits[1]) {
	    for($i = $mc_limits[0]; $i >= $mc_limits[1]; --$i) {
		if(exists($brax_hash{$i})) {
		    $candidate_brax .= "$brax_hash{$i}";
		} else {
		    $candidate_brax .= "\.";
		}
	    }
	}
		
	$candidate_brax_chopped = substr($candidate_brax,0,((length $candidate_brax) - 2));
	while($candidate_brax_chopped =~ /\./g) { ## mismatches at the last two nts are ignored -- they are not part of the miR/miR* duplex
	    ++$n_unpaired_cand;
	}
	
	
	while($candidate_brax =~ /\(/g) {
	    ++$n_left_paired_cand;
	}
	while($candidate_brax =~ /\)/g) {
	    ++$n_right_paired_cand;
	}
	if($n_unpaired_cand > 4) {
	    $fail = 1;
	}
	
	if(($n_left_paired_cand > 0 ) and
	   ($n_right_paired_cand > 0)) {
	    $fail = 1;
	}
	
	# if the candidate is still viable, find the map_coordinates for the star sequence
	unless($fail) {
	    $star_coord = get_star_coord(\$map_coord,\%brax_hash,\$candidate_brax);
	    if($star_coord eq "fail") {
		$fail = 1;
	    }
	}
	
	# first check whether the star sequence has four or fewer mismatches, omitting the last two nts
	unless($fail) {
	    @star_limits = split ("-", $star_coord);
	    $star_brax = '';
	    $n_unpaired_star = 0;
	    $n_left_paired_star = 0;
	    $n_right_paired_star = 0;
	    
	    if($star_limits[0] < $star_limits[1]) {
		for($i = $star_limits[0]; $i <= $star_limits[1]; ++$i) {
		    if(exists($brax_hash{$i})) {
			$star_brax .= $brax_hash{$i};
		    } else {
			# if the candidate is partially within the padded region, consider all parts of the padded region unpaired
			$star_brax .= "\.";
		    }
		}
	    } elsif ($star_limits[0] > $star_limits[1]) {
		for($i = $star_limits[0]; $i >= $star_limits[1]; --$i) {
		    if(exists($brax_hash{$i})) {
			$star_brax .= "$brax_hash{$i}";
		    } else {
			$star_brax .= "\.";
		    }
		}
	    }
	    $star_brax_chopped = substr($star_brax,0,((length $star_brax) - 2));
	    while($star_brax_chopped =~ /\./g) { ## mismatches at the last two nts are ignored -- they are not part of the miR/miR* duplex
		++$n_unpaired_star;
	    }
	    while($star_brax =~ /\(/g) {
		++$n_left_paired_star;
	    }
	    while($star_brax =~ /\)/g) {
		++$n_right_paired_star;
	    }
	    if($n_unpaired_star > 4) {
		$fail = 1;
	    }
	    if(($n_left_paired_star > 0 ) and
	       ($n_right_paired_star > 0)) {
		$fail = 1;
	    }
	}
	
	# begin output if it the locus is a winner
	unless($fail) {
	    push(@output,$locus);
	    my $clean_locus = $locus;
	    $clean_locus =~ s/:/_/g;
	    $outfile = "$$out_dir" . "\/" . "$clean_locus" . "\.txt";
	    open(OUT, ">$outfile");
	    print OUT "$locus\n";
	    print OUT "$displayseq\n";
	    if($strand eq "Watson") {
		for($i = $loc_start; $i <= $loc_stop; ++$i) {
		    if(exists($brax_hash{$i})) {
			print OUT "$brax_hash{$i}";
		    } else {
			print OUT " ";
		    }
		}
		print OUT "\n";
		for($i = $loc_start; $i <= $loc_stop; ++$i) {
		    if($i == $key) {
			for($j = $loc_start; $j < $i; ++$j) {
			    print OUT "\.";
			}
			if($map_coord =~ /(\d+)-(\d+)/) {
			    # get length here as well
			    $read_length = $2 - $1 + 1;
			    for($j = $1; $j <= $2; ++$j) {
				print OUT "$hp_seq_hash{$j}";
			    }
			} else {
			    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
			}
			for($j = ($2 + 1); $j <= $loc_stop; ++$j) {
			    print OUT "\.";
			}
			print OUT "\t";
			print OUT "$query_name\n";
		    }
		}
	    } elsif ($strand eq "Crick") {
		for($i = $loc_stop; $i >= $loc_start; --$i) {
		    if(exists($brax_hash{$i})) {
			print OUT "$brax_hash{$i}";
		    } else {
			print OUT " ";
		    }
		}
		print OUT "\n";
		for($i = $loc_stop; $i >= $loc_start; --$i) {
		    if($i == $key) {
			for($j = $loc_stop; $j > $i; --$j) {
			    print OUT "\.";
			}
			if($map_coord =~ /(\d+)-(\d+)/) {
			    # get length here as well
			    $read_length = $1 - $2 + 1;
			    for($j = $1; $j >= $2; --$j) {
				print OUT "$hp_seq_hash{$j}";
			    }
			} else {
			    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
			}
			for($j = ($2 - 1); $j >= $loc_start; --$j) {
			    print OUT "\.";
			}
			print OUT "\t";
			print OUT "$query_name\n";
		    }
		}
	    }
	}
    }
    print STDERR " Done\n";  ## closes progress
    return @output;
}

sub get_star_coord {
    my ($mir,$brax_hash,$mir_brax) = @_; ## passed by reference
    if ($$mir_brax =~ /^\.+$/) {
	return "fail";
    }
    my $strand;
    # determine strand of the bracket
    my @mir_coords = split ("-",$$mir);
    if($mir_coords[0] < $mir_coords[1]) {
	$strand = "+";
    } else {
	$strand = "-";
    }
    
    # hash the brackets: keys: ( , values: )
    my %leftpairs = ();
    my %rightpairs = ();
    my @left = ();
    my @b_pos = keys (%$brax_hash);
    if($strand eq "+") {
	@b_pos = sort {$a <=> $b} @b_pos;
    } elsif ($strand eq "-") {
	@b_pos = sort {$b <=> $a} @b_pos;
    }
    my $c;
    foreach my $b (@b_pos) {
	if($$brax_hash{$b} eq "\(") {
	    push(@left,$b);
	} elsif($$brax_hash{$b} eq "\)") {
	    $c = pop @left;
	    $leftpairs{$c} = $b;
	    $rightpairs{$b} = $c;
	}
    }
    
    my $i;
    my $j;
    my $star_start;
    my $star_stop;
    # lookup the correct arm
    if($$mir_brax =~ /^[\(\.]+$/) {
	# candidate miRNA on 5p arm, use leftpairs hash

	if($strand eq "+") {
	    $j = 0;
	    $i = ($mir_coords[1]) - 2;
	    until(exists($leftpairs{$i})) {
		--$i;
		++$j;
	    }
	    $star_start = $leftpairs{$i} - $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($leftpairs{$i})) {
		++$i;
		++$j;
	    }
	    $star_stop = $leftpairs{$i} + 2 + $j;
	    
	} elsif($strand eq "-") {
	    $j = 0;
	    $i = ($mir_coords[1]) + 2;
	    until(exists($leftpairs{$i})) {
		++$i;
		++$j;
	    }
	    $star_start = $leftpairs{$i} + $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($leftpairs{$i})) {
		--$i;
		++$j;
	    }
	    $star_stop = $leftpairs{$i} - 2 - $j;
	}
    } elsif ($$mir_brax =~ /^[\)\.]+$/) {
	# candidate miRNA on 3p arm, use rightpairs hash
	if($strand eq "+") {
	    $j = 0;
	    $i = $mir_coords[1] - 2;
	    until(exists($rightpairs{$i})) {
		--$i;
		++$j;
	    }
	    $star_start = $rightpairs{$i} - $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($rightpairs{$i})) {
		++$i;
		++$j;
	    }
	    $star_stop = $rightpairs{$i} + 2 + $j;
	} elsif ($strand eq "-") {
	    $j = 0;
	    $i = $mir_coords[1] + 2;
	    until(exists($rightpairs{$i})) {
		++$i;
		++$j;
	    }
	    $star_start = $rightpairs{$i} + $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($rightpairs{$i})) {
		--$i;
		++$j;
	    }
	    $star_stop = $rightpairs{$i} - 2 - $j;
	}
    } else {
	return "fail";
    }
    
    my $answer = "$star_start" . "-" . "$star_stop";
    return $answer;
}

sub get_hp_clusters_mh {
    my ($input_hash,$genome) = @_;  ## passed by reference
    my %output = ();
    my @entries = ();
    my $orig;
    my $chr;
    my $entry;
    my @en_fields = ();
    my @hp_coords = ();
    my $padded_start;
    my $padded_stop;
    my $chr_len;
    my $new_hp_cluster;
    my $fai_file = "$$genome" . "\.fai";
    my @orig_fields = ();
    while(($orig) = each %$input_hash) {
	if($orig =~ /^(\S+):(\d+)-(\d+):/) {
	    $chr = $1;
	} else {
	    die "FATAL in sub-routine \"get_hp_clusters\" : could not parse chromosome name from $orig\n";
	}
	@orig_fields = split (":",$orig);
	@entries = @{$$input_hash{$orig}};
	foreach $entry (@entries) {
	    @en_fields = split ("\t", $entry);
	    @hp_coords = split ("-", $en_fields[1]);
	    if($hp_coords[0] < $hp_coords[1]) {
		$padded_start = $hp_coords[0] - 15;
		$padded_stop = $hp_coords[1] + 15;
	    } else {
		$padded_start = $hp_coords[1] - 15;
		$padded_stop = $hp_coords[0] + 15;
	    }
	    
	    # get the length of the chromosome
	    unless (-r $fai_file) {
		die "FATAL in sub-routine \"get_hp_clusters\" : could not open fai file $fai_file\n";
	    }
	    open(FAI, "$fai_file");
	    $chr_len = 0;
	    while (<FAI>) {
		if($_ =~ /^$chr\t(\d+)\t/) {
		    $chr_len += $1;
		}
	    }
	    close FAI;
	    unless($chr_len) {
		die "Fatal in sub-routine get_folding_regions : failed to get chromosome length for $chr\n";
	    }
	    
	    if($padded_stop > $chr_len) {
		$padded_stop = $chr_len;
	    }
	    if($padded_start < 1) {
		$padded_start = 1;
	    }
	    
	    $new_hp_cluster = "$chr" . ":" . "$padded_start" . "-" . "$padded_stop" . ":" . "$orig_fields[2]" . ":" . "$orig_fields[3]" . ":" . "$orig_fields[4]";
	    $output{$new_hp_cluster} = $entry;
	}
    }
    return %output;
}

sub check_overlap {
    my($ok) = @_; ## array reference
    my %output = ();
    my %in_another = ();
    my @okdup = @ok;
    my $q_chr;
    my $q_start;
    my $q_stop;
    my $q_strand;
    my $s_chr;
    my $s_start;
    my $s_stop;
    my $s_strand;
    my $overlap;
    foreach my $okloc (@ok) {
	if($okloc =~ /^(\S+):(\d+)-(\d+):(\S)/) {
	    $q_chr = $1;
	    $q_start = $2;
	    $q_stop = $3;
	    $q_strand = $4;
	    foreach my $oklocdup (@okdup) {
		if($oklocdup eq $okloc) {
		    next;
		} elsif ($oklocdup =~ /^(\S+):(\d+)-(\d+):(\S)/) {
		    $s_chr = $1;
		    $s_start = $2;
		    $s_stop = $3;
		    $s_strand = $4;
		    
		    $overlap = 0;
		    if(($q_chr eq $s_chr) and
		       ($q_strand eq $s_strand)) {
			
			if(($q_start >= $s_start) and
			   ($q_start <= $s_stop)) {
			    $overlap = 1;
			}
			if(($q_stop >= $s_start) and
			   ($q_stop <= $s_stop)) {
			    $overlap = 1;
			}
		    }
		    
		    if(exists($output{$oklocdup})) {
			$overlap = 0;  ##
		    }
		    if(exists($in_another{$okloc})) {
			$overlap = 0;
		    }
		    
		    if($overlap) {
			push(@{$output{$okloc}}, $oklocdup);
			$in_another{$oklocdup} = 1;
		    }
		}
	    }
	}
    }
    return %output;
}

sub get_distinct_hps {
    # given a brax, return a hash with keys as starts and values as stops for all distinct hairpins in the input
    my($brax) = @_;  ## simply passed as scalar
    my @chars = split('', $brax);
    my $i = 0;
    my %pairs = get_left_right($brax);
    my %rl_pairs = ();
    my $r;
    my $l;
    while(($l,$r) = each %pairs) {
	$rl_pairs{$r} = $l;
    }
    my $last_right;  # right means ")"
    my $first_left;  # left means "("
    my %regions = (); ## left(key), right(value)
    foreach my $ch (@chars) {
	++$i;  ## one-based position
	if($ch eq "\(") {
	    if($last_right) {
		$regions{"$rl_pairs{$last_right}"} = $last_right;
		$last_right = '';
	    }
	} elsif ($ch eq "\)") {
	    $last_right = $i;
	}
    }
    # the last one
    if($last_right) {
	if(exists($rl_pairs{$last_right})) {
	    $regions{"$rl_pairs{$last_right}"} = $last_right;
	}
    }
    return %regions;
}

			
__END__
=head1 LICENSE

miR_homologs.pl

Copyright (C) 2012 Michael J. Axtell                                                             
                                                                                                 
This program is free software: you can redistribute it and/or modify                             
it under the terms of the GNU General Public License as published by                             
the Free Software Foundation, either version 3 of the License, or                                
(at your option) any later version.                                                              
                                                                                                 
This program is distributed in the hope that it will be useful,                                  
    but WITHOUT ANY WARRANTY; without even the implied warranty of                                   
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                    
GNU General Public License for more details.                                                     
                                                                                                 
You should have received a copy of the GNU General Public License                                
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 SYNOPSIS

Annotation of putative MIRNA loci based on alignments of known mature miRNAs, and analysis of predicted secondary structure.

=head1 VERSIONS

0.2 : THIS VERSION.  Released June 12, 2012 in ShortStack v 0.1.3 package.  Fixes major bug in hairpin structure parsing that was causing inappropriate rejection of valid hairpins.  Mirrors the same bug fix in ShortStack v 0.1.3

0.1 : Initial release. May 4, 2012

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 CITATION

If you use miR_homologs in your work, please cite 

Axtell MJ. (2012) ShortStack: Comprehensive annotation and quantification of small RNA genes.  In prep.

A manuscript describing the ShortStack package will be submitted sometime in the Spring/Summer of 2012, so check Pubmed first or look for an update!


=head1 INSTALL

install samtools from <http://samtools.sourceforge.net/> and ensure that samtools is in your PATH

install the ViennaRNA package <http://www.tbi.univie.ac.at/~ivo/RNA/> and ensure that RNALfold is in your PATH

ensure the script is executable                                                                  
                                                                                                 
    chmod +x miR_homologs.pl                                                         
                                                                                                 
ensure the script is in your PATH (examples):                                                    
                                                                                                 
    sudo cp miR_homologs.pl /usr/bin/                                                
                                                                                                 
OR just for one session assuming script is in your working directory:                            
                                                                                                 
    PATH=$PATH:.                                                                                 
                                                                                                 
ensure 'perl' is located in /usr/bin/ .. if not, edit line 1 of script accordingly                 

=head1 USAGE
                                                                                                                             
miR_homologs.pl [options] [alignments.bam] [genome.fasta]
 

=head1 QUICK START 

1. Install miR_homologs.pl and required third-party tools per above instructions

2. Gather a set of mature miRNAs to use as queries against your reference genome.  It's better to use a non-redundant set (in terms of sequences).

3. Ensure the chromosome names of the reference genome are short and sweet, containing no whitespace or metacharacters (see below)

4. Align your reads to the reference genome, and output the results in sam/bam format.  Suggested aligner is bowtie 1 (0.12.7) but method that outputs in sam/bam format is fine. 

5. Ensure your alignments have intact CIGAR strings, not just "*" placeholders (see SAM specification).  Bowtie 1's SAM output is fine; check for other aligners. 

6. Convert the alignment to a .bam file sorted by chromosomal position.   See samtools documentation for methods.

7. To run with default parameters, call "miR_homologs.pl [in.bam] [genome.fasta]".  See OPTIONS below for other options and run modes.

=head1 TEST

Some Arabidopsis test data can be found at http://axtelldata.bio.psu.edu/data/ShortStack_TestData/

1.  Athaliana_genome.tgz : The "TAIR10" Arabidopsis thaliana (ecotype-Col-0) genome assembly including the plastid and mitochrondria, and it's .fai index.  Retrieved from Phytozome.  This is the assembly to which the .bam files in this directory were mapped.

2.  col_leaf_ok.bam[.bai] : Sorted and indexed small RNA-seq alignments in BAM format.  Derived from wild-type rosette leaves -- Liu et al. (2012) Plant Physiology PMID: 22474216.  This alignment contains 26,523,213 mapped reads, 14,351,052 of which were "uniquely" mapped (just one alignment), and a total of 104,980,568 alignments.  The small RNA sizes range from 15-27nts.  To create this alignment, the raw .csfasta and .QV.qual files were combined to make a colorspace-fastq formatted file, adapters were trimmed along with the corresponding quality values (including the hybrid 3' color and Q value), and mapped using bowtie 0.12.7.  The bowtie settings were -C -v 1 --best --strata -k 50 --col-keepends -S, which allow zero or one mismatch, keeping only the best scoring 'stratum', and retaining only the first 50 alignments observed, and outputting in sam format.  SAM lines corresponding to unmapped reads were filtered out.  The SAM file was then processed with Prep_bam.pl (included in ShortStack package) to add the NH:i: tags to each alignment, and to output a chromosomal-sorted alignment in the BAM format.

3. ath_mb18_ShortStack_loci.txt : Coordinates for Arabidopsis thaliana MIRNA hairpin sequences, as determined by taking the top-scoring hit from a blastn search using miRBase 18 ath- hairpins as queries against the reference genome.  This file is useful as input for a ShortStack run in --count mode.

4.  ath_mature_nr.bam[.bai] : Sorted and indexed alignments of all non-redundant mature Arabidopsis thaliana miRNAs from miRBase18 against the A. thaliana reference genome.  Mapped with bowtie 0.12.7 using settings -f -v 0 -m 20 ... perfect matches only, no alignments reported if more than 20 were observed.  These alignments are useful for testing the miR_homologs.pl helper script.

Some Tests:

A) full de-novo annotation run:

    ./ShortStack.pl col_leaf_ok.bam Athaliana_167.fa

B) count mode run to quantify and annotate known miRBase MIRNA loci:

./ShortStack.pl --count ath_mb18_ShortStack_loci.txt col_leaf_ok.bam Athaliana_167.fa

C) Analyze annotated miRBase mature miRNAs for acceptable structure with miR_homologs.pl:

    ./miR_homologs.pl ath_mature_nr.bam Athaliana_167.fa


=head1 OPTIONS

--outdir [string] : Name of directory to be created to receive results of the run.  Deafults to "miR_homologs_[time]", where time is the number of non-leap seconds since Jan 1, 1970 UCT, if not provided   
                                      
--foldwindow [integer] : Size of window to fold.  Window is centered on the alignment.  Default: 300

--maxhpsep [integer] : Maximum allowed span for a base-pair during hairpin search (Option -L for RNALfold).  Default = 300.  Must be between 50 and 2000.

--minfracpaired [float] : Minimum fraction of paired nucleotides required within a valid hairpin structure.  Default = 0.67.  Allowed values are greater than 0 and less than or equal to 1.

--minntspaired [integer] : Minimum absolute number of paired nucleotides required within a valid hairpin structure.  Default = 30.  Allowed values are greater than zero and less than or equal to --maxhpsep


=head1 KEY FORMATTING REQUIREMENTS AND ASSUMPTIONS

=head2 Input .bam file

The aligned queries in the input .bam file must be sorted by chromosomal location, and indexed using the samtools index command -- specifically, miR_homologs will look for the [prefix].bam.bai file in the same directory as the input [prefix].bam file.   

Each mapped read must have the CIGAR string set (column 6 in the SAM specification) -- miR_homologs.pl determines the small RNA lengths by parsing the CIGAR string .. if any mappings have "*" entered instead of a valid CIGAR string miR_homologs.pl will exit and complain.

Finally, it is critical that the chromosome names referenced in the .bam file correspond exactly to those present in the genome.fasta file.

=head2 Input genome.fasta file

It is critical that this be the precise genome to which the reads in the input .bam file were mapped.

Additionally, the chromosome names in the FASTA headers must be kept SIMPLE.  Specifically, miR_homologs.pl at several points parses locus information by the regex /^(\S+):(\d+)-(\d+)$/ or some variant thereof, where the first pattern is the chromosome name.  Therefore, the chromosome names must match (\S+) .. e.g. a single string of one or more non-white-space characters, with no metacharacters.  So, ">Chr1" in your reference genome is good, but ">Chr1 | XM00023 | this is a bunch of annotation blah blah blah" is bad.  This same concern applies to the input .bam file, so your chromosome names should be shortened BEFORE mapping your reads, so that they are short and they are exactly reflected in the .bam file.

If not already present, a .fai index file for the genome will be created using samtools faidx at the beginning of the run.

=head1 SUGGESTED WORKFLOW

=head2 Queries

Typcially queries should be a set of non-redundant mature miRNA sequences, all upper-case, with U's turned into T's prior to alignment, in FASTA format.  The FASTA headers should be short and purged of whitespace and metacharacters.  For instance, ">ath-miR172B* MIMAT00101" should not be used because it has both whitespace and a metacharacter ('*').  Convert to ">ath-miR172bstar".

=head2 Reference (genome)

FASTA headers should be short and devoid of whitespace, as described above.  Repeat-masking not necessary if strategy below is used.  Index the genome file by typing 

    samtools index [genome]

Since the suggested workflow uses bowtie 0.12.7 as the aligner, you also need to build a bowtie index for the genome.

    ./bowtie-build [genome] [prefix]

=head2 Alignment

We currently use bowtie version 1 (0.12.7) for alignments.  Settings -f -v 0 -a -m 20 -S (fasta formatted input, retain only exact matches, retain all valid alignments subject to a cap of 20, queries with 20 or more alignments are suppressed, output in SAM format).  We also find it helpful to scrub the results on the fly with a awk command to omit .sam lines corresponding to unmapped queries.  An example aligner call might look like:

    ./bowtie -f -v 0 -a -m 20 -S [index] [queries.fasta] | awk '$3!="*"' > [output.sam]

=head2 Processing alignments

After the above, we are left with a read-sorted .sam file, and we need to get a chrosomal-sorted .bam file.  Use samtools to manipulate the alignment as follows:

    samtools view -b -S [output.sam] > [output_unsorted.bam]

    samtools sort [output_unsorted.bam] [output_sorted_prefix]
    
    samtools index [output_sorted_prefix.bam]

The resulting [output_sorted_prefix.bam] file is ready for use by miR_homologs.

=head2 Run

    ./miR_homologs.pl [options] [output_sorted_prefix.bam] [genome.fasta]

=head1 OUTPUT

Each MIRNA locus that passes the analysis will have its own simple text file to display the details of the locus.  These text files all show A) the Name and genomic coordinates of the locus, B) the sequence, in RNA form, C) the identified hairpin structure, in dot-bracket notation, and D) the aligned mature miRNA

A brief description of the results, including overlapping loci, is also provided on STDERR.

Finally, a file called "Resuts.txt" is created in the outdir.  This is a tab-delimited file where describing each non-redudndant locus.  The file is compatiable with input into ShortStack under --count mode.

=head1 METHODS

miR_homologs.pl uses code borrowed from ShortStack in order to determine whether an alignment could correspond to a qualifyin MIRNA locus.  Loci that are positively identified conform the criteria outlined in Meyers et al. (2008) Plant Cell 20:3186-3190. PMID: 19074682, EXCEPT that actual expression is not taken into account.  In order to be annotated as a potential MIRNA homolog by miR_homologs.pl, identified locus must meet the following criteria:

1. The putative hairpin must have a valid structure per options --minfracpaired and --minntspaired.

2. The aligned query (the putative mature miRNA) must contain no more than 4 unpaired nts (excepting the 2nts on the 3' end), and they must not span a loop (i.e., no base-pairs to themselves).

3. The predicted miRNA* of the candidate mature miRNA must contain no more than 4 unpaired nts (excepting the 2nts on the 3' end), and it must not span a loop (i.e., no base-pairs to themselves).  Predicted miRNA*'s are based on identifying the small RNA that would form a miR/miR* duplex with a canonical 2nt, 3' overhang.


=cut


