#!/usr/bin/perl -w
# See below the __END__ mark for docmentation, license, citation, etc., or just see the README distributed with this script

use Getopt::Long;
use strict;

###############MAIN PROGRAM BLOCK
##### VERSION
my $version = "0.2.2";
#####

##### get options and validate them

# usage statement and option gathering
my $usage = usage($0,$version);

# Initial option definition, inlcuding default settings
my $outdir = '';
my $reads = '';
my $mindepth = 20;
my $pad = 100;
my $dicermin = 20;
my $dicermax = 24;
my $maxhpsep = 300;
my $minfracpaired = 0.67;
my $minntspaired = 30;
my $minfrachpdepth = 0.67;  ## raised from 0.5 as of version 0.2.1
my $minstrandfrac = 0.8;
my $mindicerfrac = 0.85;
my $count = '';
my $nohp = 0;
my $raw = '';
my $phasesize = 21;
my $phaseFDR = 0.05;
my $inv_file = '';

# get user options from command line
GetOptions ('outdir=s' => \$outdir,
	    'reads=i' => \$reads,
	    'mindepth=i' => \$mindepth,
	    'pad=i' => \$pad,
	    'dicermin=i' => \$dicermin,
	    'dicermax=i' => \$dicermax,
	    'maxhpsep=i' => \$maxhpsep,
	    'minfracpaired=f' => \$minfracpaired,
	    'minntspaired=i' => \$minntspaired,
	    'minfrachpdepth=f' => \$minfrachpdepth,
	    'minstrandfrac=f' => \$minstrandfrac,
	    'mindicerfrac=f' => \$mindicerfrac,
	    'count=s' => \$count,
	    'nohp' => \$nohp,
	    'raw' => \$raw,
	    'phasesize=s' => \$phasesize,
	    'inv_file=s' => \$inv_file,
	    'phaseFDR=f' => \$phaseFDR);

# Validate the options
# default output directory is ShortStack_[time], where time is the number of non-leap seconds since 00:00:00 UCT Jan 1, 1970
unless($outdir) {
    my $time = time;
    $outdir = "ShortStack" . "_$time";
}

# ensure the output directory does not already exist.  If it does, quit and complain
# then create it or die tryin'
if(-d $outdir) {
    die "Fatal: Output directory $outdir already exists\n$usage\n";
} else {
    system "mkdir $outdir";
}

unless(($mindepth > 0) and
       ($mindepth < 1000000)) {
    die "FATAL: mindepth must be more than zero and less than one million\n\n$usage\n";
}

# ensure that option --pad is present and logical -- 0 >= x < 100,000
unless(($pad >= 0) and ($pad <= 100000)) {
    die "Option --pad must be a number greater than or equal to zero and less than one hundred thousand\n\n$usage\n";
}

# ensure that option --dicermin is present and logical -- 15 >= x <= 35, and less than or equal to dicermax
unless(($dicermin >= 15) and ($dicermin <= 35) and ($dicermin <= $dicermax)) {
    die "Option --dicermin must be a number between 15 and 35, and less than or equal to option --dicermax\n\n$usage\n";
}

# ensure that option --dicermax is present and logical -- 15 >= x <= 35, and more than or equal to dicermin
unless(($dicermax >= 15) and ($dicermax <= 35) and ($dicermin <= $dicermax)) {
    die "Option --dicermax must be a number between 15 and 35, and more than or equal to option --dicermin\n\n$usage\n";
}

# ensure that option --maxhpsep is present and logical -- 50 >= x <= 2000
unless(($maxhpsep >= 50) and ($maxhpsep <= 2000)) {
    die "Option --maxhpsep must be a number between 50 and 2000\n\n$usage\n";
}

# ensure that option --minfracpaired is present and logical -- 0 > x <= 1
unless(($minfracpaired > 0) and ($minfracpaired <= 1)) {
    die "Option --minfracpaired must be a number greater than zero and less than or equal to one\n\n$usage\n";
}

# ensure that option --minntspaired is present and logical -- 0 > x <= --maxhpsep
unless(($minntspaired > 0) and ($minntspaired <= $maxhpsep)) {
    die "Option --minntspaired must be number greater than zero and less than or equal to option --maxhpsep\n\n$usage\n";
}

# ensure that option --minfrachpdepth is present and logical -- 0 >= x <= 1
unless(($minfrachpdepth >= 0) and ($minfrachpdepth <= 1)) {
    die "Option --minfrachpdepth must be a number greater than or equal to zero and less than or equal to one\n\n$usage\n";
}

# ensure that option --minstrandfrac is present and logical -- 0.5 >= x <= 1
unless(($minstrandfrac >= 0.5) and ($minstrandfrac <= 1)) {
    die "Option --minstrandfrac must be number greater than or equal to 0.5 and less than or equal to one\n\n$usage\n";
}

# ensure that option --mindicerfrac is present and logical -- 0 >= x <= 1
unless(($mindicerfrac >= 0) and ($mindicerfrac <= 1)) {
    die "Option -mindicerfrac must be number greater than or equal to 0 and less than or equal to one\n$usage\n";
}

# if the --count option was provided, make sure the indicated file can be read
if($count) {
    unless(-r $count) {
	die "FATAL: Option --count : File $count could not be read\n\n$usage\n";
    }
}

# check the --phasesize option to see if it is a number between --dicermin and --dicermax, OR 'all', OR 'none'
if($phasesize =~ /^\d+$/) {
    unless(($phasesize >= $dicermin) and
	   ($phasesize <= $dicermax)) {
	die "FATAL: Option --phasesize must be an integer within the --dicermin to --dicermax size range OR \'all\'\n\n$usage\n";
    }
} else {
    unless(($phasesize eq "all") or
	   ($phasesize eq "none")) {
	die "FATAL: Option --phasesize must be an integer within the --dicermin to --dicermax size range OR \'all\' OR \'none\'\n\n$usage\n";
    }
}

# check the --phaseFDR value to ensure it is between 0 and 1
unless(($phaseFDR >= 0) and
       ($phaseFDR <= 1)) {
    die "FATAL: Option --phaseFDR must be a number between 0 and 1\n\n$usage\n";
}


###############
# check the two required, ordered files from the end of the command line array
my $genome = pop @ARGV;
unless(-r $genome) {
    die "genome file $genome is not readable\n\nFATAL\n$usage\n";
}

# ensure bam file is readable
my $bamfile = pop @ARGV;
unless(-r $bamfile) {
    die "bamfile $bamfile not readable\.\n$usage\n";
}
#############

# check for required installation of samtools and get version, or quit and complain
my ($samtools_version,$full_samtools_test_text) = get_samtools_version();
if($samtools_version eq "not found") {
    die "samtools not found\n$full_samtools_test_text\n\n$usage\n";
}

###########
# Run some checks on the provided bamfile.  
# First, look at just the first non-header line to make sure that the NH:i: tag is present
my($nh) = check_nh($bamfile);
unless($nh) {
    die "FATAL: The provided \.bam file does not appear to contain the NH:i tag, indicating number of mappings for each read\. $0 $version requires the NH:i: tag\n\n$usage\n";
}

# Then, check for presence of the .bai index file
my $expected_index = "$bamfile" . "\.bai";
unless(-e $expected_index) {
    die "FATAL: The provided \.bam file does not appear to have been indexed -- I could not find the expected index file $expected_index\nUse samtools to sort the data by chromosomal location, and then index the file, before using this program\n\n$usage\n";
}


#########
# Check for the existence of the einverted .inv file
# If not provided, issue a warning and proceed
# If a string is provided but no file is readble, quit and complain
if($inv_file) {
    unless(-r $inv_file) {
	die "FATAL: The provided \.inv file $inv_file from user option --inv_file could not be read\n$usage\n";
    }
    # if the run is in 'nohp' mode, than this file is irrelevant, so warn the user about that
    if($nohp) {
	print STDERR "\nWARNING: An einverted \.inv file was provided \($inv_file\) but this run is in nohp mode, so the einverted file will be ignored\n\n";
    }
} else {
    unless($nohp) {
	print STDERR "\nWARNING: No einverted \.inv file was provided \(option --inv_file\)\.  This may decrease the annotation accuracy, especially for loci deriving from large inverted repeats\n\n";
    }
}


#open a log file
my $logfile = "$outdir" . "\/" . "Log\.txt";
open(LOG, ">$logfile");

##### Report to user on initialization of the run
print STDERR "\n$0 $version\n";
print STDERR `date`;
print STDERR "samtools $samtools_version";
print STDERR "Mapped small RNAs: $bamfile\n";
print STDERR "Genome: $genome\n";
print STDERR "Output Directory: $outdir\n";
unless($nohp) {
    print STDERR "einverted \.inv file of Inverted Repeats:";
    if(-r $inv_file) {
	print STDERR " $inv_file\n";
    } else {
	print STDERR " NOT PROVIDED, hairpins from RNALfold only\n";
    }
}
print STDERR "Number of mapped reads in $bamfile:";
if($reads) {
    print STDERR " User Provided: $reads\n";
} else {
    print STDERR " Not provided: Forcing run in --raw mode\n";
    $raw = 1;
}
print STDERR "Clusters:";
if($count) {
    print STDERR " Running in \"count\" mode\. User provided clusters from file $count\n";
} else {
    print STDERR " To be calculated de novo\n";
    print STDERR "Island threshold: $mindepth mappings\n";
    print STDERR "Padding around initial islands: $pad nts\n";
}
print STDERR "Dicer size range: $dicermin to $dicermax\n";
if($nohp) {
    print STDERR "Running in \"nohp\" mode: Hairpins and MIRNAs will not be inferred\n";
} else {
    print STDERR "Maximum allowed separation of a base pair to span during hairpin prediction: $maxhpsep\n";
    print STDERR "Minimum fraction of paired nts allowable in a valid hairpin structure: $minfracpaired\n";
    print STDERR "Minimum number of base pairs within a valid hairpin structure: $minntspaired\n";
    print STDERR "Minimum fraction of coverage in hairpin helix to keep hairpin: $minfrachpdepth\n";
}
print STDERR "Minimum fraction of mappings to assign a polarity to a non-hairpin cluster: $minstrandfrac\n";
print STDERR "Minimum fraction of mappings within the Dicer size range to annotate a locus as Dicer-derived: $mindicerfrac\n";
print STDERR "Cluster type to analyze for phasing: $phasesize\n";
print STDERR "False Discovery Rate for Significantly Phased Clusters: $phaseFDR\n";
print STDERR "Reporting units: ";
if($raw) {
    print STDERR "Raw mappings\n\n";
} else {
    print STDERR "Mappings per million mapped\n\n";
}

### Duplicate the above to the log
##### Report to user on initialization of the run
print LOG "\n$0 $version\n";
print LOG `date`;
print LOG "samtools $samtools_version";
print LOG "Mapped small RNAs: $bamfile\n";
print LOG "Genome: $genome\n";
print LOG "Output Directory: $outdir\n";
unless($nohp) {
    print LOG "einverted \.inv file of Inverted Repeats:";
    if(-r $inv_file) {
	print LOG " $inv_file\n";
    } else {
	print LOG " NOT PROVIDED, hairpins from RNALfold only\n";
    }
}
print LOG "Number of mapped reads in $bamfile:";
if($reads) {
    print LOG " User Provided: $reads\n";
} else {
    print LOG " Not provided: Forcing run in --raw mode\n";
}
print LOG "Clusters:";
if($count) {
    print LOG " Running in \"count\" mode\. User provided clusters from file $count\n";
} else {
    print LOG " To be calculated de novo\n";
    print LOG "Island threshold: $mindepth mappings\n";
    print LOG "Max padding around initial islands: $pad nts\n";
}
print LOG "Dicer size range: $dicermin to $dicermax\n";
if($nohp) {
    print LOG "Running in \"nohp\" mode: Hairpins and MIRNAs will not be inferred\n";
} else {
    print LOG "Maximum allowed separation of a base pair to span during hairpin prediction: $maxhpsep\n";
    print LOG "Minimum fraction of paired nts allowable in a valid hairpin structure: $minfracpaired\n";
    print LOG "Minimum number of base pairs within a valid hairpin structure: $minntspaired\n";
    print LOG "Minimum fraction of coverage in hairpin helix to keep hairpin: $minfrachpdepth\n";
}
print LOG "Minimum fraction of mappings to assign a polarity to a non-hairpin cluster: $minstrandfrac\n";
print LOG "Minimum fraction of mappings within the Dicer size range to annotate a locus as Dicer-derived: $mindicerfrac\n";
print LOG "Cluster type to analyze for phasing: $phasesize\n";
print LOG "False Discovery Rate for Significantly Phased Clusters: $phaseFDR\n";
print LOG "Reporting units: ";
if($raw) {
    print LOG "Raw mappings\n\n";
} else {
    print LOG "Mappings per million mapped\n\n";
}

# Next, check for the presence of .fai index file corresponding to the genome.  If not found, create one with samtools
my $expected_faidx = "$genome" . "\.fai";
unless(-e $expected_faidx) {
    print STDERR "Expected genome index $expected_faidx for genome file $genome not found\.  Creating it using samtools faidx";
    system "samtools faidx $genome";
    print STDERR " done\n\n";
}

###############################################

##### Phase One: Identify Clusters
my %names = ();  ## keys = locusIDs (e.g. Chr1:1-1000) , values = name designations

print STDERR `date`;
print STDERR "Phase One: Identifying Clusters";
my @clusters = ();
if($count) {
    print STDERR " a priori from file $count\n";
    print LOG "Clusters determined a priori from file $count\n";
    @clusters = get_clusters_a_priori($count);
    %names = get_names_countmode($count);  ## if names not present in count file, arbitrary names assigned
} else {
    print STDERR " de novo\n";
    print LOG "Cluster determined de novo\n";
    my @first_clusters = get_clusters($bamfile,$mindepth);
    @clusters = merge_clusters(\@first_clusters,\$pad,\$genome);
}

unless((scalar @clusters) > 0) {
    die "Sorry, no clusters meeting your criteria were found in these data\.\n";
}


##### Phase Two: Gather genomic sequences for input into RNA folding
my @final_clusters = ();
my %hp_clusters = ();
my %miRNAs = ();
my %put_mir = ();  ## only really used when option --putativemir is called

if($nohp) {
    print STDERR "Running in nohp mode -- skipping phases 2, 3, and 4\n\n";
    @final_clusters = @clusters;
    unless($count) {
	# if clusters were de novo, generate arbitrary names for them now if the hairpin routines are being skipped
	%names = get_names_simple(\@final_clusters);
    }
} else {
    print STDERR "\n";
    print STDERR `date`;
    print STDERR "Phase Two: Gathering genome sequence for RNA secondary structure prediction\n";
    
    # if running in --count mode, only the exact interval defined in the a priori file will be folded.  Else, the de novo method takes a larger window
    my %to_fold = ();
    if($count) {
	%to_fold = get_folding_regions_countmode(\$genome,\@clusters);
    } else {
	%to_fold = get_folding_regions(\$genome,\@clusters,\$pad);  ## hash structure: 'sequence', 'start', 'stop'  strand is always Watson
    }
    print STDERR "\n\n";

    #### Phase three: Assess correlation of clusters with qualifying inverted repeats

    print STDERR `date`;
    print STDERR "Phase Three: Identification of hairpin-associated small RNA clusters\n";

    #fold each query and retain only qualifying structures

    print STDERR "\tFolding sequences and identifying qualifying hairpins\.\.\.";
    my %qualifying_hairpins = folder(\%to_fold,\$maxhpsep,\$minntspaired,\$minfracpaired);  ## each entry has locus as key and an array of entries.  Each entry is tab-delim w/ brax, local-start, structure_details, and strand
    
    #convert the coordinates of the qualifying hairpins to the true chromosomal coordinates
    print STDERR "\tConverting hairpin coordinates to true genomic coordinates\.\.\.";
    my %true_hps = hairpin_coords(\%to_fold,\%qualifying_hairpins);  ## structure same as %qualifying hairpins
    print STDERR "Done\n";
    
    ## HERE einverted
    if(-r $inv_file) {
	print STDERR "\tParsing inverted repeats from file $inv_file\.\.\.\n";
	my @initial_irs = parse_inv(\$inv_file,\$minntspaired,\$minfracpaired);  ## each entry is tab-delimited .. [0] pseudo-brax, [1] helix-coords, [2] helix_details, [3] strand, [4] Chr (to be removed during merge)
	print STDERR "\tdone\n";
	my $initial_ir_tally = scalar @initial_irs;
	print STDERR "\t\tFound $initial_ir_tally potential qualifying IRs\n";
	print STDERR "\tMerging einverted-derived inverted repeats with clusters and with RNAL-fold derived hairpins\.\.\.";
	my ($in_irs,$out_irs) = merge_inv(\@initial_irs,\%true_hps,\@clusters);
	print STDERR " done\n";
	print STDERR "\t\tFrom the $in_irs potentially qualifying IRs, $out_irs had overlap with one or more clusters and were merged with RNALfold hairpins\n";
    }

    print STDERR "\tRemoving redundant hairpins on a per-locus basis\.\.\.";
    my %true_hps_trimmed = remove_redundant_hps(\%true_hps); ## structure same as %qualifying hairpins
    print STDERR "Done\n";
    
    
    print STDERR "\tRemoving hairpins that do not overlap clusters\.\.\.";
    my %true_hps_3 = remove_nonoverlapped_hps(\%true_hps_trimmed);  ## structure same as %qualifying hairpins
    print STDERR "Done\n";
    
    print STDERR "\tFiltering hairpins based on expression evidence\.\.\.\n";
    my %filtered_hps = hp_expression(\%true_hps_3,\$bamfile,\$minfrachpdepth);

    print STDERR "\n\n";
    print STDERR `date`;
    print STDERR "Phase 4: Annotation of miRNA and other hpRNA loci\n";
    print STDERR "\tExamining candidate hairpins for miRNA and other hpRNA annotation\.\.\.\n";
    
    my %trash = hp_output (\%filtered_hps,\$genome,\$bamfile,\$reads,\$minfrachpdepth,\$minstrandfrac);  ## the important result of this process is to modify the %filtered_hps hash
    %trash = ();
    print STDERR "\n";
    
    print STDERR "\tFinishing parsing of hairpin-associated clusters\.\.\.";
    if($count) {
	## In count mode, clusters are not modified in position
	## In addition, if there happens to be more than one qualified hairpin within a locus provided a priori, only the first one listed will be used!
	@final_clusters = @clusters;
	%hp_clusters = get_hp_clusters_countmode(\%filtered_hps);
    } else {
	%hp_clusters = get_hp_clusters(\%filtered_hps);
	@final_clusters = get_final_clusters(\%filtered_hps,\@clusters,\%hp_clusters);
	## not in count mode, so clusters have been de novo and haven't been named yet
	%names = get_names_simple(\@final_clusters);
    }
    print STDERR "done\n";
    ## Write the HP and MIRNA files
    print STDERR "\tWriting miRNA and other hpRNA details to files\.\.\.";
    my($n_miRNAs,$n_hpRNAs) = write_files(\%hp_clusters, \$outdir, \%names, \@final_clusters);
    print STDERR " done\n";
    ## Report
    
    my $clus_count = scalar @final_clusters;
    print STDERR "\tTotal Clusters: $clus_count\n\tmiRNA loci: $n_miRNAs\n\tother hpRNA loci: $n_hpRNAs\n\n";

}

## Phase 5: Annotate and Quantify all clusters
print STDERR `date`;
print STDERR "Phase 5: Quantify all clusters\n";
my %quant_master = quant(\@final_clusters,\$bamfile,\$dicermin,\$dicermax,\$minstrandfrac,\$mindicerfrac,\$phasesize,\$phaseFDR,\%names);

# modify the notations for hairpins
my %quant_2 = ();
if($nohp) {
    foreach my $f_clus (@final_clusters) {
	my @qm_fields = split ("\t",$quant_master{$f_clus});
	my @new_fields = @qm_fields[2..((scalar @qm_fields) -1)];
	unshift(@new_fields,"ND");
	unshift(@new_fields,$qm_fields[1]);  ## name
	unshift(@new_fields,$qm_fields[0]);  ## loc-coords
	my $new_entry = join("\t", @new_fields);
	$quant_2{$f_clus} = $new_entry;
    }
} else {
    %quant_2 = mod_quant_hp(\%quant_master,\%hp_clusters);
}

my %quant_3 = ();
# unless --raw was chosen, convert the raw values to mappings per million mapped
if($raw) {
    %quant_3 = %quant_2;
} else {
    %quant_3 = convert_2_mpmm(\%quant_2,\$reads);
}
print STDERR "\n\n";

## Phase 6 .. Create output files and final reports
print STDERR `date`;
print STDERR "Phase 6: Outputting Results\n";

# output a .bed file
# need difft. methods if hairpins/MIRNAs are included
my $bed_text_for_log;
if($nohp) {
    $bed_text_for_log = write_bed_nohp(\@final_clusters,\%quant_3,\$outdir,\$dicermin,\$dicermax,\%names);
} else {
    $bed_text_for_log = write_bed_with_hp(\@final_clusters,\%quant_3,\$outdir,\$dicermin,\$dicermax,\%names,\%hp_clusters);
}

# write rgb information to log
print LOG "\n$bed_text_for_log\n";

# output the master file
my $big_table = "$outdir" . "\/" . "Results\.txt";
open(BIG, ">$big_table");
print BIG "\#Locus\tName\tHP\tStrand\tFrac_Watson\tTotal\tUniques\tRep-total\tDicerCall\tPhaseOffset\tPhase_pval\tPhase_FDR_Call\tShort\tLong";
for(my $d = $dicermin; $d <= $dicermax; ++$d) {
    print BIG "\t$d";
} 
print BIG "\n";
    
foreach my $f_clus (@final_clusters) {
    print BIG "$quant_3{$f_clus}\n";
}
close BIG;
    
# Provide a summary report to the user on STDERR and in the log
print STDERR "\nSummary of Results:\n";
print LOG "\nSummary of Results:\n";
my %final_summary = ();
if($nohp) {
    %final_summary = final_summary_nohp(\%quant_3);
    print STDERR "Dominant Size\tNumber of Clusters\n";
    print LOG "Dominant Size\tNumber of Clusters\n";
    if(exists($final_summary{'N'})) {
	print STDERR "Non-Dicer\t$final_summary{'N'}\n";
	print LOG "Non-Dicer\t$final_summary{'N'}\n";
    } else {
	print STDERR "Non-Dicer\t0\n";
	print LOG "Non-Dicer\t0\n";
    }
    for(my $dsize = $dicermin; $dsize <= $dicermax; ++$dsize) {
	print STDERR "$dsize\t";
	print LOG "$dsize\t";
	if(exists($final_summary{$dsize})) {
	    print STDERR "$final_summary{$dsize}\n";
	    print LOG "$final_summary{$dsize}\n";
	} else {
	    print STDERR "0\n";
	    print LOG "0\n";
	}
    }
    print STDERR "\n";
    print STDERR "Significantly phased clusters \(size = $phasesize FDR = $phaseFDR\):\n";
    print LOG "\n";
    print LOG "Significantly phased clusters \(size = $phasesize FDR = $phaseFDR\):\n";
    if(exists($final_summary{'phased'})) {
	print STDERR "$final_summary{'phased'}\n";
	print LOG "$final_summary{'phased'}\n";
    } else {
	print STDERR "0\n";
	print LOG "0\n";
    }
} else {
    %final_summary = final_summary(\%quant_3);
    print STDERR "Dominant Size\tNOT-HAIRPINS\tHAIRPINS\tMIRNAs\n";
    print LOG "Dominant Size\tNOT-HAIRPINS\tHAIRPINS\tMIRNAs\n";
    if(exists($final_summary{'N'}{'X'})) {
	print STDERR "Non-Dicer\t$final_summary{'N'}{'X'}\t";
	print LOG "Non-Dicer\t$final_summary{'N'}{'X'}\t";
    } else {
	print STDERR "Non-Dicer\t0\t";
	print LOG "Non-Dicer\t0\t";
    }
    if(exists($final_summary{'N'}{'HP'})) {
	print STDERR "$final_summary{'N'}{'HP'}\t";
	print LOG "$final_summary{'N'}{'HP'}\t";
    } else {
	print STDERR "0\t";
	print LOG "0\t";
    }
    if(exists($final_summary{'N'}{'MIRNA'})) {
	print STDERR "$final_summary{'N'}{'MIRNA'}\n";
	print LOG "$final_summary{'N'}{'MIRNA'}\n";
    } else {
	print STDERR "0\n";
	print LOG "0\n";
    }
    
    
    for(my $dsize = $dicermin; $dsize <= $dicermax; ++$dsize) {
	print STDERR "$dsize\t";
	print LOG "$dsize\t";
	if(exists($final_summary{$dsize}{'X'})) {
	    print STDERR "$final_summary{$dsize}{'X'}\t";
	    print LOG "$final_summary{$dsize}{'X'}\t";
	} else {
	    print STDERR "0\t";
	    print LOG "0\t";
	}
	if(exists($final_summary{$dsize}{'HP'})) {
	    print STDERR "$final_summary{$dsize}{'HP'}\t";
	    print LOG "$final_summary{$dsize}{'HP'}\t";
	} else {
	    print STDERR "0\t";
	    print LOG "0\t";
	}
	if(exists($final_summary{$dsize}{'MIRNA'})) {
	    print STDERR "$final_summary{$dsize}{'MIRNA'}\n";
	    print LOG "$final_summary{$dsize}{'MIRNA'}\n";
	} else {
	    print STDERR "0\n";
	    print LOG "0\n";
	}
    }
    print STDERR "\n";
    print STDERR "Significantly phased clusters \(size = $phasesize FDR = $phaseFDR\):\n";
    print LOG "\n";
    print LOG "Significantly phased clusters \(size = $phasesize FDR = $phaseFDR\):\n";
    if(exists($final_summary{'phased'})) {
	print STDERR "$final_summary{'phased'}\n";
	print LOG "$final_summary{'phased'}\n";
    } else {
	print STDERR "0\n";
	print LOG "0\n";
    }
}
print STDERR "Completed\n";
print STDERR `date`;
print LOG "Completed\n";
print LOG `date`;
close LOG;


############### SUB ROUTINE BLOCKS
sub check_nh {
    my($bamfile) = @_;
    my $result = 0;
    open(SAM, "samtools view $bamfile |");
    while (<SAM>) {
	if($_ =~ /^@/) {
	    next;
	} else {
	    if($_ =~ /\tNH:i:\d+/) {
		$result = 1;
	    }
	    last;
	}
    }
    close SAM;
    return $result;
}    
	
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

sub usage {
    my($name,$version) = @_;
    my $usage = "\n$name $version
USAGE: $name \[options\] bam_file genome_file
OPTIONS:
--outdir [string] : name for output directory that ShortStack will create\.  Defaults to \"ShortStack_\[unix-time\]\" if not provided\.
--reads [float] : number of mapped reads in bam_file\.  If not provided, ShortStack will force the analysis into \'raw\' mode, and no reads-per-million adjustments will be made\.
--inv_file [string] : path to .inv file from einverted analysis of the genome_file in question\.  If not provided, warning is issued but analysis proceeds\.
--mindepth [integer] : threshold for calling islands \(in reads per million\; must be more than 0 and less than 1 million\; default: 20\)
--pad [integer] : Number of nucleotides upstream and downstrem to extend initial islands during cluster definition\; default: 100
--dicermin [integer] : smallest size of the Dicer-derived small RNAs \(must be between 15-35 and less than or equal to option --dicermax\; default: 20\)
--dicermax [integer] : largest size of the Dicer-derived small RNAs \(must be between 15-35 and more than or equal to option --dicermin\; default: 24\)
--maxhpsep [integer] : maximum allowed separation of a base pair to span during hairpin prediction \(option -L for RNALfold\; must be between 50 and 2000\;default: 300\)
--minfracpaired [float] : minimum fraction of paired nts allowable in a hairpin structure\; default: 0.67
--minntspaired [integer] : minimum number of base-pairs in an accetable hairpin structure\; default: 30
--minfrachpdepth [float] : minimum fraction of coverage in hairpin helix to keep hairpin\; default: 0.67
--minstrandfrac [float] : minimum fraction of mappings to assign a polarity to a non-hairpin cluster\; default: 0.8
--mindicerfrac [float] : minimum fraction of mappings within the dicer size range to annotate a locus as Dicer-derived\; default: 0.85
--count [string] : File containing a-priori defined clusters\.  Presence of --count triggers \"count\" mode, and de-novo cluster discovery is skipped\. Default: off
--phasesize [string or integer] : Cluster type to examine for significant phasing\.  Must be within the --dicermin to --dicermax range, or \'all\' to look at all dicer clusters, or \'none\' to skip all phasing analyses\. Default: 21
--phaseFDR [float] : Benjamini-Hochberg false discovery rate for calling significantly phased clusters\.  Default: 0.05
--nohp : If --nohp appears on the command line, program executes in \"no_hp\" mode, omitting hairpin and MIRNA inferences\.  Default: off
--raw : If --raw appears on the command line, program reports quantifications in raw mappings, instead of normalizing for library size to mappings per million mapped\. Default: Off

";
    return $usage;
}
      
sub get_clusters {
    my($bamfile,$actual_threshold) = @_;
    my @clusters = ();
    my @fields = ();
    my $last_chr = "null";
    my $last_good_coordinate;
    my $current_start_coordinate;
    my $open = 0;
    my $chr;
    my $position;
    my $last_position = 0;
    my $cluster;
    my $end;
    my $read_length;
    my %depth = ();
    my $i;
    my @occupied = ();
    my $coordinate;
    my @set_of_depths = ();
    my $delta;
    my $analysis_position;
    my $analysis_depth;
    
    # Sanity check and warning / failure if threshold is too low
    # Warn if >= 2 <= 5 mappings
    # fail if less than 2
    if($actual_threshold < 2) {
	die "\nFATAL: With --mindepth set at $actual_threshold , islands will be defined by having as few as $actual_threshold reads\nThat is too low\!  Try again after increasing --mindepth\n\n";
    } elsif ($actual_threshold <= 5) {
	print STDERR "\n\*\*\* WARNING \*\*\*\n";
	print STDERR "At --mindepth setting $actual_threshold, islands are being defined with as few as $actual_threshold reads \.\.\. which seems very low\nConsider aborting this run and increasing --mindepth\n";
    }
    
    
    # for a crude progress counter
    my $nlines = 0; # number of lines
    my $onehunk = 0; # number of 100ks
    my $onemil = 0;
    print STDERR "\tProgress in sub-routine \'get_clusters\': ";
    ##
    
    # initialize depth hash
    for($i = 0; $i < 100; ++$i) {
	$depth{$i} = 0;
    }
    
    #  read the samfile
    open(SAM, "samtools view $bamfile |");
    while (<SAM>) {
	chomp;
	# skip headers
	if($_ =~ /^@/) {
	    next;
	}
	# progress
	++$nlines;
	if($nlines >= 100000) {
	    ++$onehunk;
	    $nlines = 0;
	    print STDERR "\.";
	}
	if($onehunk >= 10) {
	    ++$onemil;
	    $onehunk = 0;
	    print STDERR " $onemil Million mappings examined\n\t";
	}
	##

	# get fields of interest for the current line
	@fields = split ("\t", $_);
	$chr = $fields[2];
	$position = $fields[3];
	# ensure the read is mapped.  If not, go to the next line
	if($fields[1] & 4 ) {
	    next;
	}
	$read_length = parse_cigar($fields[5]);
	
	
	
	# check if you've passed into the next chromosome.  If so, process accordingly
	if(($last_chr ne $chr) and
	   ($last_chr ne "null")) {
	    
	    ## close out last chr
	    $delta = 99;  ## ensures all possible positions at the end of the last chr are analyzed
	    # analyze read depths at the last position up through all positions skipped before the current position	    
	    
	    for($i = 0; $i < $delta; ++$i) {
		# if the delta has gone beyond the 100nt window, break the loop
		if($i > 99) {
		    last;
		} else {
		    $analysis_position = $i + $last_position;
		    $analysis_depth = $depth{$i};
		    if($open) {
			# currently within a valid cluster
			# if the current position is beyond , must close the currently open cluster and then check to see if a new one should be opened
			if($analysis_position > ($last_good_coordinate + 1)) {
			    # close
			    $open = 0;
			    $end = $last_good_coordinate;
			    $cluster = "$last_chr" . ":" . "$current_start_coordinate" . "-" . "$end";
			    push(@clusters, $cluster);
			    
			    # TEST
			    #print STDERR "$cluster\n";
			    # END TEST
			    
			    $last_good_coordinate = '';
			    $current_start_coordinate = '';
			    
			    # open new one if warranted
			    if($analysis_depth >= $actual_threshold) { # 
				$open = 1;
				$current_start_coordinate = $analysis_position;
				$last_good_coordinate = $analysis_position;
			    }
			} else {
			    # or, the position under analysis is the next consecutive one.  Check depth -- if depth OK, modify the last good coordinate
			    if($analysis_depth >= $actual_threshold) {
				$last_good_coordinate = $analysis_position;
			    }
			}
		    } else {
			# no cluster is currently open.  Initiate a new cluster if warranted by the depth
			if($analysis_depth >= $actual_threshold) {
			    $open = 1;
			    $current_start_coordinate = $analysis_position;
			    $last_good_coordinate = $analysis_position;
			}
		    }
		}
	    }
	    ## close cluster if one is still open
	    if($open) {
		$open = 0;
		$end = $last_good_coordinate;
		$cluster = "$last_chr" . ":" . "$current_start_coordinate" . "-" . "$end";
		push(@clusters, $cluster);
		# TEST
		#print STDERR "$cluster\n";
		# END TEST
		$last_good_coordinate = '';
		$current_start_coordinate = '';
	    }
	    
	    # clear out depth hash
	    for($i = 0; $i < 100; ++$i) {
		$depth{$i} = 0;
	    }
	}
	
	# check if you've passed to a new left-most position on the same chr?  If so, analyze all positions that were just passed by.
	if(($last_position != $position) and
	   ($last_chr eq $chr)) {
	    
	    $delta = $position - $last_position;  # should be a positive number of 1 or more
	    # analyze read depths at the last position up through all positions skipped before the current position	    
	    
	    for($i = 0; $i < $delta; ++$i) {
		# if the delta has gone beyond the 100nt window, break the loop
		if($i > 99) {
		    last;
		} else {
		    $analysis_position = $i + $last_position;
		    $analysis_depth = $depth{$i};
		    if($open) {
			# currently within a valid cluster
			# if the current position is not the next one, must close the currently open cluster and then check to see if a new one should be opened
			if($analysis_position > ($last_good_coordinate + 1)) {
			    # close
			    $open = 0;
			    $end = $last_good_coordinate;
			    $cluster = "$chr" . ":" . "$current_start_coordinate" . "-" . "$end";
			    push(@clusters, $cluster);
			    # TEST
			    #print STDERR "$cluster\n";
			    # END TEST
			    $last_good_coordinate = '';
			    $current_start_coordinate = '';
			    
			    # open new one if warranted
			    if($analysis_depth >= $actual_threshold) { 
				$open = 1;
				$current_start_coordinate = $analysis_position;
				$last_good_coordinate = $analysis_position;
			    }
			} else {
			    # or, the position under analysis is within the allowable gap distance.  Check depth -- if depth OK, modify the last good coordinate
			    if($analysis_depth >= $actual_threshold) {
				$last_good_coordinate = $analysis_position;
			    }
			}
		    } else {
			# no cluster is currently open.  Initiate a new cluster if warranted by the depth
			if($analysis_depth >= $actual_threshold) {
			    $open = 1;
			    $current_start_coordinate = $analysis_position;
			    $last_good_coordinate = $analysis_position;
			}
		    }
		}
	    }
	    
	    # move values over to reflect position
	    # reset the array first
	    @set_of_depths = ();
	    for ($i = 0; $i < 100; ++$i) {
		push(@set_of_depths,$depth{$i});
	    }
	    
	    for ($i = 1; $i <= $delta; ++$i) {
		my $junk = shift @set_of_depths;
		$junk = '';
		push(@set_of_depths, 0);
	    }
	    for($i = 0; $i < 100; ++$i) {
		$depth{$i} = $set_of_depths[$i];
	    }
	}
	
	# Add counts correspoinding to the current mapped read, all the way through the end of the read
	for($i = 0; $i < $read_length; ++$i) {
	    ++$depth{$i};
	}
	# update last chr and last position before turning to the next line
	$last_chr = $chr;
	$last_position = $position;
    }
    
    ## close out the end of the final chr
    ## close out last chr
    $delta = 99;  ## ensures all possible positions at the end of the last chr are analyzed
    # analyze read depths at the last position up through all positions skipped before the current position	    
    
    for($i = 0; $i < $delta; ++$i) {
	# if the delta has gone beyond the 100nt window, break the loop
	if($i > 99) {
	    last;
	} else {
	    $analysis_position = $i + $last_position;
	    $analysis_depth = $depth{$i};
	    if($open) {
		# currently within a valid cluster
		# if the current position is beyond the max gap length, must close the currently open cluster and then check to see if a new one should be opened
		if($analysis_position > ($last_good_coordinate + 1)) {
		    # close
		    $open = 0;
		    $end = $last_good_coordinate;
		    $cluster = "$last_chr" . ":" . "$current_start_coordinate" . "-" . "$end";
		    push(@clusters, $cluster);
		    # TEST
		    #print STDERR "$cluster\n";
		    # END TEST

		    $last_good_coordinate = '';
		    $current_start_coordinate = '';
		    
		    # open new one if warranted
		    if($analysis_depth >= $actual_threshold) { # 
			$open = 1;
			$current_start_coordinate = $analysis_position;
			$last_good_coordinate = $analysis_position;
		    }
		} else {
		    # or, the position under analysis is within the allowable gap distance.  Check depth -- if depth OK, modify the last good coordinate
		    if($analysis_depth >= $actual_threshold) {
			$last_good_coordinate = $analysis_position;
		    }
		}
	    } else {
		# no cluster is currently open.  Initiate a new cluster if warranted by the depth
		if($analysis_depth >= $actual_threshold) {
		    $open = 1;
		    $current_start_coordinate = $analysis_position;
		    $last_good_coordinate = $analysis_position;
		}
	    }
	}
    }
    ## close cluster if one is still open
    if($open) {
	$open = 0;
	$end = $last_good_coordinate;
	$cluster = "$last_chr" . ":" . "$current_start_coordinate" . "-" . "$end";
	push(@clusters, $cluster);
	# TEST
	#print STDERR "$cluster\n";
	# END TEST
	$last_good_coordinate = '';
	$current_start_coordinate = '';
    }
    print STDERR " Done\n";
    return @clusters;
}

sub merge_clusters {
    my($input,$pad,$genome) = @_; ## passed my reference .. array and scalar
    my @output = ();

    my $last_start;
    my $last_stop;
    my $this_padded_start;
    my $this_padded_stop;
    my $last_chr = "null";
    my %chr_sizes = ();
    ## grab the chrom sizes, which you need to ensure that you don't pad off the end of the chroms
    # chrom sizes in column 1 from the fai file
    my $fai_file = "$$genome" . "\.fai";
    unless(-e $fai_file) {
	die "Fatal in sub-routine get_folding_regions : expected fai file $fai_file does not exist\n";
    }
    my @fai_fields = ();
    open(FAI, "$fai_file");

    while (<FAI>) {
	@fai_fields = split ("\t", $_);
	$chr_sizes{$fai_fields[0]} = $fai_fields[1];
    }
    close FAI;
    my $this_chr;
    my $entry;
    foreach my $in_clus (@$input) {
	if($in_clus =~ /^(\S+):(\d+)-(\d+)$/) {
	    $this_chr = $1;
	    $this_padded_start = $2 - $$pad;
	    if($this_padded_start < 1) {
		$this_padded_start = 1;
	    }
	    $this_padded_stop = $3 + $$pad;
	    if($this_padded_stop > $chr_sizes{$this_chr}) {
		$this_padded_stop = $chr_sizes{$this_chr};
	    }
	    
	    # special first case
	    if($last_chr eq "null") {
		$last_start = $this_padded_start;
		$last_stop = $this_padded_stop;
	    } elsif ($this_chr ne $last_chr) {
		$entry = "$last_chr" . ":" . "$last_start" . "-" . "$last_stop";
		push(@output,$entry);
		$last_start = $this_padded_start;
		$last_stop = $this_padded_stop;
	    } else {
		if($this_padded_start > $last_stop) {
		    $entry = "$this_chr" . ":" . "$last_start" . "-" . "$last_stop";
		    push(@output,$entry);
		    $last_start = $this_padded_start;
		    $last_stop = $this_padded_stop;
		} else {
		    if($this_padded_start < $last_start) {
			$last_start = $this_padded_start;
		    }
		    if($this_padded_stop > $last_stop) {
			$last_stop = $this_padded_stop;
		    }
		}
	    }
	    $last_chr = $this_chr;
	} else {
	    die "FATAL: in sub-routine \'merge_clusters\' : failed to parse initial locus $in_clus\n";
	}
    }
    return @output;
}
    

sub get_clusters_a_priori {
    my($flatfile) = @_;
    my @clusters = ();
    my @fields = ();
    my %tracker = ();
    open(FILE, "$flatfile");
    while (<FILE>) {
	chomp;
	# ignore comment lines
	if($_ =~ /^\#/) {
	    next;
	}
	if($_ =~ /\t/) {
	    @fields = split ("\t", $_);
	    unless($fields[0] =~ /^\S+:\d+-\d+$/) {
		die "FATAL in sub-routine get_clusters_a_priori : cluster name $fields[0] is not understandable\n";
	    }
	    push(@clusters,$fields[0]);
	    if(exists($tracker{$fields[0]})) {
		die "FATAL in sub-routine get_clusters_a_priori : cluster $fields[0] appears more than once\!  Duplicate loci must be purged prior to analysis\n";
	    } else {
		$tracker{$fields[0]} = 1;
	    }
	} else {
	    unless($_ =~ /^\S+:\d+-\d+$/ ) {
		die "FATAL in sub-routine get_clusters_a_priori : cluster name $_ is not understandable\n";
	    }
	    if(exists($tracker{$_})) {
		die "FATAL in sub-routine get_clusters_a_priori : cluster $fields[0] appears more than once\!  Duplicate loci must be purged prior to analysis\n";
	    } else {
		$tracker{$_} = 1;
	    }
	    push(@clusters,$_);
	}
    }
    close FILE;
    return @clusters;
}

sub get_names_countmode {
    my($flatfile) = @_;
    my %names = ();
    my @fields = ();
    my $n = 0;
    my $entry;
    open(FILE, "$flatfile");
    while (<FILE>) {
	chomp;
	# ignore comment lines
	if($_ =~ /^\#/) {
	    next;
	}
	if($_ =~ /\t/) {
	    @fields = split ("\t", $_);
	    unless($fields[0] =~ /^\S+:\d+-\d+$/) {
		die "FATAL in sub-routine get_names_countmode : cluster name $fields[0] is not understandable\n";
	    }
	    if($fields[1]) {
		$names{$fields[0]} = $fields[1];
	    } else {
		die "FATAL in sub-routine get_names_countmode : file is tab-delmited but no name entry in second column\n";
	    }
	} else {
	    unless($_ =~ /^\S+:\d+-\d+$/ ) {
		die "FATAL in sub-routine get_names_countmode : cluster name $_ is not understandable\n";
	    }
	    ++$n;
	    $entry = "Cluster_" . "$n";
	    $names{$_} = $entry;
	}
    }
    close FILE;
    return %names;
}
		

sub get_folding_regions {
    my($gen_file,$info,$pad) = @_;  ## passed as references .. scalar, array, scalar, and scalar, respectively
    my $locus;

    my %to_fold = ();
    my $chr;
    my $start;
    my $stop;
    my $locus_size;
    my $unp_locus_size_3x;
    my $get_size;
    my $middle;
    my $get_start;
    my $get_stop;
    my $gen_line;

    my $ok;
    my $subseq;

    my $chr_length;
    my $fai_file = "$$gen_file" . "\.fai";
    
    # for progress tracking
    my $n_to_get = scalar (@$info);
    my $five_percent = int (0.05 * $n_to_get);
    my $x;
    print STDERR "\tProgress in sub-routine \"get_folding_regions\" \(dot = five percent\): ";
    
    foreach $locus (@$info) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	# parse locus name
	if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    die "Fatal in sub-routine get_folding_regions: could not parse locus name $locus\n";
	}

	$locus_size = $stop - $start + 1;
	$unp_locus_size_3x = int (3 * ($locus_size - (2 * $$pad)));  ## the unpadded region is the locus size - (2 * $pad).

	if($locus_size >= 1000) {
	    $get_size = $locus_size;
	} elsif($unp_locus_size_3x > 1000) {
	    $get_size = 1000;
	} elsif($unp_locus_size_3x < 250) {
	    $get_size = 250;
	} else {
	    $get_size = $unp_locus_size_3x;
	}
	
	# get the coordinates .. centered on the middle of the cluster
	# first get the middle of the cluster
	$middle = int(($locus_size / 2) + $start);
	$get_start = int ($middle - ($get_size / 2));
	$get_stop = int ($middle + ($get_size / 2));
	
	# ensure that the get_start is not less than 1
	if($get_start < 1) {
	    $get_start = 1;
	}
	
	# ensure that the get_stop does not exceed the chromosome length
	# the chromosome lengths are in the second column of the .fai file
	unless(-e $fai_file) {
	    die "Fatal in sub-routine get_folding_regions : expected fai file $fai_file does not exist\n";
	}
	open(FAI, "$fai_file");
	$chr_length = 0;
	while (<FAI>) {
	    if($_ =~ /^$chr\t(\d+)\t/) {
		$chr_length += $1;
	    }
	}
	close FAI;
	unless($chr_length) {
	    die "Fatal in sub-routine get_folding_regions : failed to get chromosome length for $chr\n";
	}
	if($get_stop > $chr_length) {
	    $get_stop = $chr_length;
	}
	# get the sequence of interest
	# first, build the query
	my $query = "$chr" . ":" . "$get_start" . "-" . "$get_stop";
	# call samtools faidx
	$subseq = '';  ## reset
	open(FAIDX, "samtools faidx $$gen_file $query |");
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
	
	## add to the to_fold hash
	$to_fold{$locus}{'sequence'} = $subseq;
	$to_fold{$locus}{'start'} = $get_start;
	$to_fold{$locus}{'stop'} = $get_stop;
    }
    # finish progress tracking
    print STDERR " Done\n";
    # return the hash
    return %to_fold;
}


sub get_folding_regions_countmode {
    my($gen_file,$info) = @_;  ## passed as references .. scalar, array
    my $locus;

    my %to_fold = ();
    my $chr;
    #my $chr_name;
    my $start;
    my $stop;
    my $locus_size;
    my $locus_size_3x;
    my $get_size;
    #my @val_fields = ();
    #my $get_strand;
    my $middle;
    my $get_start;
    my $get_stop;
    my $gen_line;
    my $chr_seq;
    my $ok;
    my $subseq;
    #my $got_seq;
    
    # for progress tracking
    my $n_to_get = scalar (@$info);
    my $five_percent = int (0.05 * $n_to_get);
    my $x;
    print STDERR "\tProgress in sub-routine \"get_folding_regions_countmode\" \(dot = five percent\): ";
    
    foreach $locus (@$info) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	# parse locus name
	if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    die "Fatal in sub-routine get_folding_regions: could not parse locus name $locus\n";
	}
	
	
	# call samtools faidx to get the FASTA formatted version of the whole thing
	$subseq = '';
	open(FAIDX, "samtools faidx $$gen_file $locus |");
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
	
	## add to the to_fold hash
	$to_fold{$locus}{'sequence'} = $subseq;
	$to_fold{$locus}{'start'} = $start;
	$to_fold{$locus}{'stop'} = $stop;
    }
    # finish progress tracking
    print STDERR " Done\n";
    # return the hash
    return %to_fold;
}

sub folder {
    my($to_fold,$L,$min_paired,$min_frac_paired) = @_;  ## passed by reference, hash and scalar, respectively
    my $locus;
    my $fold_seq;
    my $brax;
    my $delta_G;
    my %output = ();  
    my $structure_details;
    my $start;
    my $entry;  ## tab-delimited .. brax, start, helix_info (e.g. 123-150,180-200), strand "Watson" or "Crick"
    my $revcomp;
    my %regions = ();
    my $brax_section;
    my $left;
    my $right;

    my $adj_st;
    
    # first, get some information to enable a crude progress bar
    my $n_loci_to_fold = scalar ( keys %$to_fold);
    my $x = 0;
    my $five_percent = int(0.05 * $n_loci_to_fold);
    print STDERR "\n\tProgress in sub-routine \"folder\" \(dot = 5 percent\): ";
    
    while(($locus) = each %$to_fold) {
	# progress tracking
	++$x;
	if($x == $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	
	# First, the Watson Strand is folded
	open(RNALFOLD, "echo $$to_fold{$locus}{'sequence'} | RNALfold -d 2 -noLP -L $$L |");
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
	    
	    ## Note, as of v.0.1.3, deltaG is no longer used. 
	    ##if($_ =~ /\s+\((.*\d+.*)\)/) {
	    #	$delta_G = $1;
	    #	$delta_G =~ s/\s//g;
	    #} else {
	    #	die "FATAL in sub-routine \'folder\' : failed to parse deltaG from RNALfold output line: $_\n";
	    #}
	    
	    if($_ =~ /(\d+)\s*$/) {
		$start = $1;
	    } else {
		die "FATAL in sub-routine \'folder\' : failed to parse start position from RNALfold output line $_\n";
	    }
	    
	    # new v0.1.3 .. split the RNALfold-provided structure into distinct hairpins .. this fixes a major problem in earlier versions with false negatives on hairpin discovery!
	    %regions = get_distinct_hps($brax);
	    while(($left, $right) = each %regions) {
		$brax_section = substr($brax,($left - 1),($right-$left+1));
	    
		# send the structure to the general structure evaluation sub-routine, which returns zero to reject, one to keep
		$structure_details = evaluate_structure_general($brax_section,$$min_paired,$$min_frac_paired);

		# if it is OK, adjust the coordinates, and add to output
		unless($structure_details eq "bogus") {
		    # structure details are 123-150,180-200 .. e.g. start and stop positions of the helix of interest, one-based coordinates, relative to the brackets themselves.
		    $adj_st = $left + $start - 1;
		    $entry = "$brax_section\t$adj_st\t$structure_details\t" . "Watson";
		    push(@{$output{$locus}}, $entry);
		}
	    }
	}
	close RNALFOLD;
	
	# Now, examine the Crick strand .. revcomp the sequence
	$revcomp = reverse $$to_fold{$locus}{'sequence'};
	$revcomp =~ tr/ACUG/UGAC/;

	open(RNALFOLD, "echo $revcomp | RNALfold -d 2 -noLP -L $$L |");
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
	    
	    ## Note, as of v.0.1.3, deltaG is no longer used. 
	    #if($_ =~ /\s+\((.*\d+.*)\)/) {
	    #	$delta_G = $1;
	    #	$delta_G =~ s/\s//g;
	    #} else {
	    #	die "FATAL in sub-routine \'folder\' : failed to parse deltaG from RNALfold output line: $_\n";
	    #}
	    
	    if($_ =~ /(\d+)\s*$/) {
		$start = $1;
	    } else {
		die "FATAL in sub-routine \'folder\' : failed to parse start position from RNALfold output line $_\n";
	    }
	    
            # new v0.1.3 .. split the RNALfold-provided structure into distinct hairpins .. this fixes a major problem in earlier versions with false negatives on hairpin discovery!
	    %regions = get_distinct_hps($brax);
	    while(($left, $right) = each %regions) {
		$brax_section = substr($brax,($left - 1),($right-$left+1));
		
		# send the structure to the general structure evaluation sub-routine, which returns zero to reject, one to keep
		$structure_details = evaluate_structure_general($brax_section,$$min_paired,$$min_frac_paired);
		
		# if it is OK, adjust the coordinates, and add to output
		unless($structure_details eq "bogus") {
		    
		    # structure details are 123-150,180-200 .. e.g. start and stop positions of the helix of interest, one-based coordinates, relative to the brackets themselves.
		    $adj_st = $left + $start - 1;
		    $entry = "$brax_section\t$adj_st\t$structure_details\t" . "Crick";
		    push(@{$output{$locus}}, $entry);
		}
	    }
	}
	close RNALFOLD;
	
	

## TEST .. output to flat file
#	if(@{$output{$locus}}) {
#	    print "$locus\n";
#	    foreach $entry (@{$output{$locus}}) {
#		print"\t$entry\n";
#	    }
#	}
	## end TEST
    }
    print STDERR " Done\n";
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

sub hairpin_coords {
    my($to_fold,$input_hps) = @_;  ## passed by reference, both hashes
    my %true_hps;
    my @inputs = ();
    my @fields;
    my $locus;
    my $in_helix_coords;
    my $local_offset;
    my $brax;
    my $inp;
    #my $delta_G;  deprecated v0.1.3
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
    
    while(($locus) = each %$input_hps) {
	
	## TEST
#	print STDERR "Locus: $locus\n";
#	print STDERR "FoldedRegion: $$to_fold{$locus}{'start'} to $$to_fold{$locus}{'stop'} strand: $$to_fold{$locus}{'strand'}\n";
	##
	
	@inputs = @{$$input_hps{$locus}};
	foreach $inp (@inputs) {
	    @fields = split ("\t", $inp);
	    $strand_folded = pop @fields;
	    #$delta_G = pop @fields;
	    $in_helix_coords = pop @fields;
	    $local_offset = pop @fields;
	    $brax = pop @fields;
	    
	    ## TEST
#	    print STDERR "\tINPUT: $inp\n";
	    ##
	    
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
		$start_true = $$to_fold{$locus}{'start'} + $start_1 - 1;
		$stop_true = $$to_fold{$locus}{'start'} + $stop_1 - 1;
		$left_true_start = $$to_fold{$locus}{'start'} + $left_1_start - 1;
		$left_true_stop = $$to_fold{$locus}{'start'} + $left_1_stop - 1;
		$right_true_start = $$to_fold{$locus}{'start'} + $right_1_start - 1;
		$right_true_stop = $$to_fold{$locus}{'start'} + $right_1_stop - 1;
	    } elsif ($strand_folded eq "Crick") {
		$start_true = $$to_fold{$locus}{'stop'} - $start_1 + 1;
		$stop_true = $$to_fold{$locus}{'stop'} - $stop_1 + 1;
		$left_true_start = $$to_fold{$locus}{'stop'} - $left_1_start + 1;
		$left_true_stop = $$to_fold{$locus}{'stop'} - $left_1_stop + 1;
		$right_true_start = $$to_fold{$locus}{'stop'} - $right_1_start +1;
		$right_true_stop = $$to_fold{$locus}{'stop'} - $right_1_stop + 1;
	    }
	    
	    $brax_true = "$start_true" . "-" . "$stop_true";
	    $helix_true = "$left_true_start" . "-" . "$left_true_stop" . "," . "$right_true_start" . "-" . "$right_true_stop";
	    $true_entry = "$brax\t$brax_true\t$helix_true\t$strand_folded";
	    
	    ## TEST
#	    print STDERR "\tOUTPUT: $true_entry\n";
	    ##
	    
	    push(@{$true_hps{$locus}}, $true_entry);
	    
	}
    }
    return %true_hps;
}

sub remove_redundant_hps {
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
    
    #my $delta_G_per_nt;
    
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
	## TEST
	#print STDERR "LOCUS: $locus\n";
	#print STDERR "INPUTS:\n";
	##
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
	    $strand = $fields[3];
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
	    #@q_fields = split ("\t", $tmp_entry);

	    @q_ranges = split (",", $tmp_entry);
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

sub range_overlap {
    my($q_range,$s_range) = @_;  ## passed by reference .. ranges are arrays, number-number
    my $answer = 0;
    my @q_sorted = sort {$a <=> $b} @$q_range;
    my @s_sorted = sort {$a <=> $b} @$s_range;  ## sorting is necessary to deal with minus-strand hairpins
    my $x;
    # first check .. don't do the time-consuming enumeration if you can easily tell there is no overlap at all.
    if(($q_sorted[1] < $s_sorted[0]) or
       ($s_sorted[1] < $q_sorted[0])) {
	$answer = 0;
    } else {
	for($x = $q_sorted[0]; $x <= $q_sorted[1]; ++$x) {
	    if(($x >= $s_sorted[0]) and
	       ($x <= $s_sorted[1])) {
		$answer = 1;
		last;
	    }
	}
    }
    return $answer;
}
    
sub remove_nonoverlapped_hps {
    my($input_hash) = @_; ## passed by reference
    my $locus;
    my %scrubbed = ();
    my @loc_range = ();
    my @entries = ();
    my $entry;
    my @en_fields = ();
    my @hp_ranges = ();
    my @range1 = ();
    my @range2 = ();
    my $overlap1;
    my $overlap2;
    
    while(($locus) = each %$input_hash) {
	## TEST
	#print STDERR "LOCUS: $locus\n";
	##
	# get the range of the actual locus
	@loc_range = ();  ## clear at each iteration
	if($locus =~ /\S+:(\d+)-(\d+)$/) {
	    push(@loc_range,$1);
	    push(@loc_range,$2);
	} else {
	    die "Fatal in sub-routine remove_nonoverlapped_hps -- reg ex failure to larse locus name $locus\n";
	}
	@entries = @{$$input_hash{$locus}};
	## TEST
	#print STDERR "INPUT:\n";
	##
	foreach $entry (@entries) {
	    ## TEST
	    #print STDERR "\t$entry\n";
	    ##
	    @en_fields = split ("\t", $entry);
	    @hp_ranges = split (",", $en_fields[2]);
	    @range1 = split ("-", $hp_ranges[0]);
	    @range2 = split ("-", $hp_ranges[1]);
	    $overlap1 = range_overlap_count(\@loc_range,\@range1);
	    $overlap2 = range_overlap_count(\@loc_range,\@range2);

	    if(($overlap1 >= 20) or
	       ($overlap2 >= 20)) {
		
		## criteria for overlap is at least 20nts of overlap between at least one of the two helical arms and the cluster
		
		push(@{$scrubbed{$locus}}, $entry);
		
		## TEST
		#print STDERR "\t\tKEPT\n";
		##
		
	    }
	}
    }
    return %scrubbed;
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

sub hp_expression {
    my($input_hps,$bamfile,$min_frac) = @_; # passed by reference, hash and three scalars
    
    my $cluster;
    my @hps = ();
    my $hp_entry;
    my @hp_fields = ();
    my @hp_arms = ();
    my @fivep = ();
    my @threep = ();
    
    my $strand;
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    my $left_flank_start;
    my $right_flank_stop;
    
    my $chr;
    my $sam_start;
    my $sam_query;
    
    my @fields = ();
    my $position;
    my $read_length;
    my $read_strand;
    
    my %total_coverage = ();
    my %left_arm_coverage = ();
    my %right_arm_coverage = ();
    my $i;
    
    my $ok_strand_left_cov_sum;
    my $ok_strand_right_cov_sum;
    
    my $sum;
    my $qualifying_sum;
    my $ratio;
    
    my %output_hash = ();
    
    my $sum_left;  # 0.2.1
    my $sum_right; # 0.2.1
    
    my %original_locus = (); ## 0.2.1
    my $orig_coverage = 0;
    my $orig_start; ## 0.2.1
    my $orig_stop; ## 0.2.1
    my $sum_hp_orig; ## 0.2.1
    my $qual_sum_hp_orig;  ## 0.2.1
    my $final_ratio; ## 0.2.1

    # for progress tracking
    my $n_to_examine = scalar (keys %$input_hps);  ## this is a bit crude, as loci will have varying numbers of hairpins within them
    my $five_percent = int (0.05 * $n_to_examine);
    my $x = 0;
    print STDERR "\n\tProgress within sub-routine \"hp_expression\" \(dots = five percent\): ";
    
    while(($cluster) = each %$input_hps) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    print STDERR ".";
	    $x = 0;
	}
	# first, get the coverage total across the entire, original locus
	%original_locus = ();
	$orig_coverage = 0;
	$qual_sum_hp_orig = 0;
	if($cluster =~ /^\S+:(\d+)-(\d+)$/) {
	    $orig_start = $1;
	    $orig_stop = $2;
	} else {
	    die "FATAL: Failed to parse cluster name $cluster in sub-routine hp_expression\n";
	}
	for($i = $orig_start; $i <= $orig_stop; ++$i) {
	    $original_locus{$i} = 0;
	}
	open(SAM, "samtools view $$bamfile $cluster |");
	while (<SAM>) {
	    chomp;
	    # skip headers, just in case
	    if($_ =~ /^@/) {
		next;
	    }
	    
	    # get fields of interest for the current line
	    @fields = split ("\t", $_);
	    $position = $fields[3];
	    # ensure the read is mapped.  If not, go to the next line
	    if($fields[1] & 4 ) {
		next;
	    }
	    $read_length = parse_cigar ($fields[5]);
	    for($i = $position; $i < ($position + $read_length); ++$i) {
		if(exists($original_locus{$i})) {
		    ++$orig_coverage;
		    ++$original_locus{$i};
		}
	    }
	}
	close SAM;
	
	@hps = @{$$input_hps{$cluster}};
	foreach $hp_entry (@hps) {
	    @hp_fields = split ("\t", $hp_entry);
	    @hp_arms = split (",", $hp_fields[2]);
	    @fivep = split ("-", $hp_arms[0]);
	    @threep = split ("-", $hp_arms[1]);
	    
	    if($fivep[0] < $fivep[1]) { ## sense 
		$strand = "+";
		$left_start = $fivep[0];
		$left_stop = $fivep[1];
		$right_start = $threep[0];
		$right_stop = $threep[1];
	    } else {
		## antisense
		$strand = "-";
		$left_start = $threep[1];
		$left_stop = $threep[0];
		$right_start = $fivep[1];
		$right_stop = $fivep[0];
	    }
	    
	    $left_flank_start = $left_start - ($left_stop - $left_start);
	    $right_flank_stop = $right_stop + ($right_stop - $right_start);
	    
	    # build query region for samtools view
	    if($cluster =~ /^(\S+):/) {
		$chr = $1;
	    } else {
		die "FATAL: in hp_expression sub-routine failure to parse chr from $cluster\n";
	    }
	    

	    $sam_start = $left_flank_start;
	    if($sam_start < 1) {
		$sam_start = 1;  ## in case you are at the edge of the chr
	    }
	    $sam_query = "$chr" . ":" . "$sam_start" . "-" . "$right_flank_stop";
	    
	    # initialize tracking hashes
	    %total_coverage = ();
	    %left_arm_coverage = ();
	    %right_arm_coverage = ();
	    for($i = $left_flank_start; $i <= $right_flank_stop; ++$i) {
		$total_coverage{$i} = 0;
		if(($i >= $left_start) and
		   ($i <= $left_stop)) {
		    $left_arm_coverage{$i} = 0;
		}
		if(($i >= $right_start) and
		   ($i <= $right_stop)) {
		    $right_arm_coverage{$i} = 0;
		}
	    }
	    
            # get the alignments
	    $ok_strand_left_cov_sum = 0;
	    $ok_strand_right_cov_sum = 0;
	    
	    open(SAM, "samtools view $$bamfile $sam_query |");
	    while (<SAM>) {
		chomp;
		# skip headers, just in case
		if($_ =~ /^@/) {
		    next;
		}
	
		# get fields of interest for the current line
		@fields = split ("\t", $_);
		$position = $fields[3];
		# ensure the read is mapped.  If not, go to the next line
		if($fields[1] & 4 ) {
		    next;
		}
		$read_length = parse_cigar ($fields[5]);
		
		if($fields[1] & 16) {
		    $read_strand = "-";
		} else {
		    $read_strand = "+";
		}

		# track
		for($i = $position; $i <($position + $read_length); ++$i) {
		    if(exists($total_coverage{$i})) {
			++$total_coverage{$i};
		    }
		    if (exists($left_arm_coverage{$i})) {
			++$left_arm_coverage{$i};	
			if($read_strand eq $strand) {
			    ++$ok_strand_left_cov_sum;
			}
		    }
		    if (exists($right_arm_coverage{$i})) {
			++$right_arm_coverage{$i};
			if($read_strand eq $strand) {
			    ++$ok_strand_right_cov_sum;
			}
		    }
		}
	    }
	    close SAM;
	    
	    # calculate the result
	    $sum = 0;
	    $qualifying_sum = 0;
	    $sum_left = 0;
	    $sum_right = 0;
	    $sum_hp_orig = 0;  ## total coverage withing the HP arms only in coordinates within the original locus
	    for($i = $left_flank_start; $i <= $right_flank_stop; ++$i) {
		$sum += $total_coverage{$i};
		if(exists($left_arm_coverage{$i})) {
		    $sum_left += $left_arm_coverage{$i};
		    if(exists($original_locus{$i})) {
			$sum_hp_orig += $original_locus{$i};
		    }
		}
		if(exists($right_arm_coverage{$i})) {
		    $sum_right += $right_arm_coverage{$i};
		    if(exists($original_locus{$i})) {
			$sum_hp_orig += $original_locus{$i};
		    }
		}
	    }
	    $qualifying_sum = $sum_left + $sum_right;
	    # protect against undefined errors
	    if($sum == 0) {
		$ratio = 0;
	    } elsif (($ok_strand_left_cov_sum == 0) or  ## new in 0.2.1 .... has to be at least some reads on both arms of the correct strand to make it past this point
		     ($ok_strand_right_cov_sum == 0)) {
		$ratio = 0;
	    } elsif ((($ok_strand_left_cov_sum + $ok_strand_right_cov_sum) / $qualifying_sum) <= 0.5) {
		$ratio = 0;  ## a simple majority of coverage in the hp_arms has to be on the same strand as the hairpin itself
	    } else {
		$ratio = $qualifying_sum / $sum;
	    }
	    
	    # is the observed ratio greater than the user's minimum fraction?
	    if($ratio >= $$min_frac) {
		# OK to keep the entry for now
		$hp_entry .= "\t$ratio\t$sum_hp_orig\t$orig_coverage";
		push(@{$output_hash{$cluster}}, $hp_entry);
		$qual_sum_hp_orig += $sum_hp_orig;
	    }
	}
	# final check .. does the qual_sum_hp_orig (which is from the surviving hairpins to this point) justify splitting the original locus?
	## protect against undefined errors
	if(exists($output_hash{$cluster})) {
	    if($orig_coverage == 0) {
		$final_ratio = 0;
	    } else {
		$final_ratio = $qual_sum_hp_orig / $orig_coverage;
	    }
	    if($final_ratio < $$min_frac) {
		delete $output_hash{$cluster};
	    }
	}
    }
    
    # close progress tracking
    print STDERR "  Done\n";
    return %output_hash;
}

sub get_hp_clusters {
    my ($input_hash) = @_;  ## passed by reference
    my %output = ();
    my @entries = ();
    my $orig;
    my $chr;
    my $entry;
    my @en_fields = ();
    my @hp_coords = ();
    my $new_hp_cluster;
    my $start;
    my $stop;
    
    while(($orig) = each %$input_hash) {
	if($orig =~ /^(\S+):/) {
	    $chr = $1;
	} else {
	    die "FATAL in sub-routine \"get_hp_clusters\" : could not parse chromosome name from $orig\n";
	}
	@entries = @{$$input_hash{$orig}};
	foreach $entry (@entries) {
	    @en_fields = split ("\t", $entry);
	    @hp_coords = split ("-", $en_fields[1]);
	    if($hp_coords[0] < $hp_coords[1]) {
		$start = $hp_coords[0];
		$stop = $hp_coords[1];
	    } else {
		$start = $hp_coords[1];
		$stop = $hp_coords[0];
	    }
	    
	    $new_hp_cluster = "$chr" . ":" . "$start" . "-" . "$stop";
	    $output{$new_hp_cluster} = $entry;
	}
    }
    return %output;
}

sub get_hp_clusters_countmode {
    my ($input_hash) = @_;  ## passed by reference
    my %output = ();
    my @entries = ();
    my $orig;
    my $orig_size;
    my $hp_size;
    my @en_fields = ();
    my @hp_coords = ();
    
    while(($orig) = each %$input_hash) {
	@entries = @{$$input_hash{$orig}};
	if((scalar @entries) == 1) {
	    if($orig =~ /^\S+:(\d+)-(\d+)/) {
		$orig_size = $2 - $1 + 1;
	    } else {
		# failsafe
		$orig_size = 1000000000000000;
	    }
	    @en_fields = split ("\t", $entries[0]);
	    @hp_coords = split ("-", $en_fields[1]);
	    $hp_size = abs ($hp_coords[1] - $hp_coords[0]);
	    if(($hp_size / $orig_size) > 0.9) {  ## 90% rule .. in count mode, HP or MIRNA is kept only if there is just one within the original cluster that accounts for more than 90% of the original cluster's length
		$output{$orig} = $entries[0];  ## no padding, only the first entry is kept if there was more than one hp
	    }
	}
    }
    return %output;
}

sub get_final_clusters {
    my($input_hp_hash,$orig_clus,$hp_clus_hash) = @_; # passed by reference ..hash-hash-array-hash
    my $orig;
    my @output = ();
    my $chr;
    my @entries = ();
    my @en_fields = ();
    my @hp_coords = ();
    my $start;
    my $stop;
    my $chr_len;
    my $new_hp_cluster;
    my $entry;
    my %failsafe = ();

    
    foreach $orig (@$orig_clus) {
	if(exists($$input_hp_hash{$orig})) {
	    if($orig =~ /^(\S+):/) {
		$chr = $1;
	    } else {
		die "FATAL in sub-routine \"get_hp_clusters\" : could not parse chromosome name from $orig\n";
	    }
	    @entries = @{$$input_hp_hash{$orig}};
	    foreach $entry (@entries) {
		@en_fields = split ("\t", $entry);
		@hp_coords = split ("-", $en_fields[1]);
		if($hp_coords[0] < $hp_coords[1]) {
		    $start = $hp_coords[0];
		    $stop = $hp_coords[1];
		} else {
		    $start = $hp_coords[1];
		    $stop = $hp_coords[0];
		}
		
		$new_hp_cluster = "$chr" . ":" . "$start" . "-" . "$stop";
		
		# check to see if you've got the name right
		unless(exists($$hp_clus_hash{$new_hp_cluster})) {
		    die "FATAL in sub-routine \"get_final_clusters\" : incorrect name for $new_hp_cluster\n";
		}
		
		# add to array, after checking failsafe
		unless(exists($failsafe{$new_hp_cluster})) {
		    push(@output,$new_hp_cluster);
		}
		$failsafe{$new_hp_cluster} = 1;
	    }
	} else {
	    unless(exists($failsafe{$orig})) {
		push(@output,$orig);
	    }
	    $failsafe{$orig} = 1;
	}
    }
    return @output;
}

sub hp_output {  ## modified as of 0.2.1
    my ($hp_hash,$genome,$bamfile,$n_mapped_reads,$minhpfrac,$minstrandfrac) = @_;  ## passed by reference.  hash, scalar, scalar, scalar, scalar, scalar
    my $outfile;
    my %output = ();
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
    
    my $mmmr;
    
    my @entries = (); ## 0.2.1
    my $ratio; ## 0.2.1
    my $sam_query; ## 0.2.1
    my @kept_entries = (); ## 0.2.1
    my $hp_call; ## 0.2.1
    my $output_string; ## 0.2.1
    my $orig_cov; ## 0.2.1
    my $hp_cov_sum;  ## 0.2.1
    my $ke; ## 0.2.1
    
    my $corrected_strand_frac; ## 0.2.1
    
    my @olines = (); ## 0.2.1
    
    my $earlier_fail;  ## 0.2.1
    
    ## for progress bar
    my $n_to_process = scalar ( keys %$hp_hash);
    my $five_percent = int ($n_to_process / 20);
    my $progress = 0;
    
    print STDERR "\tProgress in sub-routine \"hp_output\" \(dot = five percent\): ";
    
    while(($locus) = each %$hp_hash) {
	# progress tracking
	++$progress;
	if($progress >= $five_percent) {
	    print STDERR "\.";
	    $progress = 0;
	}

	@entries = @{$$hp_hash{$locus}};
	@kept_entries = ();
	foreach $locus_data (@entries) {
	    
	    # parse the data
	    @l_data_fields = split ("\t", $locus_data);
	    $brax = $l_data_fields[0];
	    @brax_coords = split ("-",$l_data_fields[1]);
	    # sort so that the lower coordinate is always in the [0]
	    @brax_coords = sort {$a <=> $b} @brax_coords;
	    $strand = $l_data_fields[3];  ## either 'Watson' or 'Crick'
	    $ratio = $l_data_fields[4];
	    ### hp_cov = $l_data_fields[5]; Not needed yet, but true
	    ### orig_cov = $l_data_fields[6]; Not needed yet, but true
	    
	    # retrieve coordinates
	    if($locus =~ /^(\S+):/) {
		$chr = $1;
		$loc_start = $brax_coords[0];
		$loc_stop = $brax_coords[1];
		$sam_query = "$chr" . ":" . "$loc_start" . "-" . "$loc_stop";
	    } else {
		die "FATAL: in sub-routine \"hp_output\" could not parse chr name from $locus\n";
	    }
	    # call samtools faidx
	    $subseq = '';  ## reset
	    open(FAIDX, "samtools faidx $$genome $sam_query |");
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
	    
	    
	    # code the hp sequence and the brackets into hashes, keyed by coordinates, for easy retrieval
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
	    
	    # get all mappings that have any overlap with the interval, and track in a data structure
	    #  that keeps 'em sorted by left-most position and distinguishes sense from antisense
	    
	    %sense = ();  ## NOTE, sense and antisense are relative the hairpin direction, not (necessarily) the genome
	    %antisense = ();
	    @sense_names = ();
	    @antisense_names = ();
	    $total_mappings = 0;
	    %miRNA = (); 
	    %miRNA_star = ();
	    
	    open(SAM, "samtools view $$bamfile $sam_query |");
	    while (<SAM>) {
		chomp;
		# ignore headers
		if($_ =~ /^@/) {
		    next;
		}
		# get fields
		@samfields = split ("\t", $_);
		# ignore unmapped reads
		if($samfields[1] & 4) {
		    next;
		}
		
		$read_length = parse_cigar ($samfields[5]);
		
		$read_stop = $read_length + $samfields[3] - 1;
		# for purposes of hairpin output and MIRNA search, consider only those mappings that both start and stop within the locus
		if(($read_stop > $loc_stop) or
		   ($samfields[3] < $loc_start)) {
		    next;
		}
		
		++$total_mappings;
		
		
		
		# decide whether the read is sense or antisense wrt to the hairpin direction
		if($samfields[1] & 16) {
		    # mapping is antisense relative to genome
		    if($strand eq "Watson") {
			# mapping is antisense relative to hairpin
			$mapping_polarity = "antisense";
		    } elsif ($strand eq "Crick") {
			# mapping is sense relative to hairpin
			$mapping_polarity = "sense";
		    }
		} else {
		    # mapping is sense relative to genome
		    if($strand eq "Watson") {
			# mapping is sense relative to hairpin
			$mapping_polarity = "sense";
		    } elsif($strand eq "Crick") {
			$mapping_polarity = "antisense";
		    }
		}
		
		# the naming convention depends solely upon whether the hairpin is Watson or Crick.
		# All mappings for a Watson hairpin are start-stop, while all mappings for a Crick hairpin are
		# keyed as stop-start.  Doesn't matter whether the mapping itself is sense or antisense wrt to the hairpin
		
		if($strand eq "Watson") {
		    $map_coord = "$samfields[3]" . "-" . "$read_stop";
		    $key = $samfields[3];
		} elsif ($strand eq "Crick") {
		    $map_coord = "$read_stop" . "-" . "$samfields[3]";
		    $key = $read_stop;
		}
		
		if($mapping_polarity eq "sense") {
		    unless(exists($sense{$key}{$map_coord})) {
			push(@sense_names,$map_coord);
		    }
		    ++$sense{$key}{$map_coord};
		} elsif ($mapping_polarity eq "antisense") {
		    unless(exists($antisense{$key}{$map_coord})) {
			push(@antisense_names, $map_coord);
		    }
		    ++$antisense{$key}{$map_coord};
		}
	    }
	    close SAM;
	    
	    # Begin assessment of whether the hairpin can be annotated as a miRNA
	    @candidates = ();
	    # First, gather all sense mappings that account for 20% or more of ALL mappings at the locus
	    foreach $map_coord (@sense_names) {
		$key = $map_coord;
		$key =~ s/-\d+$//g;
		if(($sense{$key}{$map_coord} / $total_mappings) > 0.2) {
		    push (@candidates,$map_coord);
		}
	    }
	    # tracking
	    if(@candidates) {
		my $n_cands = scalar @candidates;
		${$output{$sam_query}}[0] = 1;
		${$output{$sam_query}}[1] = $ratio;
		${$output{$sam_query}}[2] = $n_cands;
		${$output{$sam_query}}[3] = 0;
		${$output{$sam_query}}[4] = 0;
		${$output{$sam_query}}[5] = 0;
		${$output{$sam_query}}[6] = 0;
		${$output{$sam_query}}[7] = 0;
		${$output{$sam_query}}[8] = 0;
		${$output{$sam_query}}[9] = 0;
	    } else {
		${$output{$sam_query}}[0] = 1;
		${$output{$sam_query}}[1] = $ratio;
		${$output{$sam_query}}[2] = 0;
		${$output{$sam_query}}[3] = 0;
		${$output{$sam_query}}[4] = 0;
		${$output{$sam_query}}[5] = 0;
		${$output{$sam_query}}[6] = 0;
		${$output{$sam_query}}[7] = 0;
		${$output{$sam_query}}[8] = 0;
		${$output{$sam_query}}[9] = 0;
	    }
	    # check whether candidates span a loop and if not whether they have four or fewer un-paired residues (ignoring the last two nts)
	    if(@candidates) {
		foreach $map_coord (@candidates) {
		    $candidate_brax = '';
		    $fail = 0;
		    $n_unpaired_cand = 0;
		    $n_left_paired_cand = 0;
		    $n_right_paired_cand = 0;
		    @mc_limits = split ("-", $map_coord);
		    #@mc_limits = sort {$a <=> $b} @mc_limits;
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
		    } else {
			++${$output{$sam_query}}[3];
		    }
		    
		    if(($n_left_paired_cand > 0 ) and
		       ($n_right_paired_cand > 0)) {
			$fail = 1;
		    } else {
			++${$output{$sam_query}}[4];
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
			} else {
			    ++${$output{$sam_query}}[5];
			}
			if(($n_left_paired_star > 0 ) and
			   ($n_right_paired_star > 0)) {
			    $fail = 1;
			} else {
			    ++${$output{$sam_query}}[6];
			}
		    }
		
		    # if the star has survived this far, see if it was mapped, and if so, whether the sum of candidate + star is >= 25% of all mappings
		    unless($fail) {
			if(exists($sense{$star_limits[0]}{$star_coord})) {
			    ++${$output{$sam_query}}[7];
			    if( (($sense{$star_limits[0]}{$star_coord} + $sense{$mc_limits[0]}{$map_coord}) / $total_mappings ) <= 0.25) {
				$fail = 1;
			    } else {
				++${$output{$sam_query}}[8];
			    }
			} else {
			    $fail = 1;
			}
		    }
		    
		    unless($fail) {
			# this duplex is good.  But might be redundant.  Check.
			if((exists($miRNA{$map_coord})) or
			   (exists($miRNA{$star_coord})) or
			   (exists($miRNA_star{$map_coord})) or
			   (exists($miRNA_star{$star_coord}))) {
			    $fail = 1;
			} else {
			    ++${$output{$sam_query}}[9];
			}
		    }
		    
		    unless($fail) {
			# you've got a winner.  Check abundances to decide which one should really be the miRNA, and which should be the star
			if($sense{$mc_limits[0]}{$map_coord} >= $sense{$star_limits[0]}{$star_coord}) {
			    $miRNA{$map_coord} = $star_coord;
			    $miRNA_star{$star_coord} = $map_coord;
			} else {
			    $miRNA{$star_coord} = $map_coord;
			    $miRNA_star{$map_coord} = $star_coord;
			}
		    }
		    
		}
	    } else {  ## closes "if candidates" 
		$fail = 1;
	    }
	    
	    # so, if it failed as a MIRNA, we will only keep it IF the ratio was >= $minhpfrac + (0.5 * (1 - $minhpfrac)).
	    if(${$output{$sam_query}}[9] >= 1) {
		$hp_call = "MIRNA";
	    } elsif ($ratio >= ($$minhpfrac + (0.5*(1-$$minhpfrac)))) {
		$hp_call = "HP";
	    } else {
		$hp_call = "FAIL";
		next;
	    }
	    
	    # begin output
	    $output_string = '';
	    $output_string .= "$sam_query $strand\n";
	    $output_string .= "$displayseq\n";
	    if($strand eq "Watson") {
		for($i = $loc_start; $i <= $loc_stop; ++$i) {
		    if(exists($brax_hash{$i})) {
			$output_string .= "$brax_hash{$i}";
		    } else {
			$output_string .= " ";
		    }
		}
		$output_string .= "\n";
		for($i = $loc_start; $i <= $loc_stop; ++$i) {
		    if(exists($sense{$i})) {
			foreach $map_coord (@sense_names) {
			    if(exists($sense{$i}{$map_coord})) {
				for($j = $loc_start; $j < $i; ++$j) {
				    if(exists($miRNA{$map_coord})) {
					$output_string .= "m";
				    } elsif(exists($miRNA_star{$map_coord})) {
					$output_string .= "\*";
				    } else {
					$output_string .= "\.";
				    }
				}
				if($map_coord =~ /(\d+)-(\d+)/) {
				    # get length here as well
				    $read_length = $2 - $1 + 1;
				    for($j = $1; $j <= $2; ++$j) {
					$output_string .= "$hp_seq_hash{$j}";
				    }
				} else {
				    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
				}
				for($j = ($2 + 1); $j <= $loc_stop; ++$j) {
				    if(exists($miRNA{$map_coord})) {
					$output_string .= "m";
				    } elsif(exists($miRNA_star{$map_coord})) {
					$output_string .= "\*";
				    } else {
					$output_string .= "\.";
				    }
				}
				
				$output_string .= " ";
				# print length and abundance, and if appropriate, designation as a miRNA or miRNA-star
				$output_string .= "l=$read_length ";
				$output_string .= "m=$sense{$i}{$map_coord}";
				if($$n_mapped_reads) {
				    $mmmr = sprintf("%.4f", (($sense{$i}{$map_coord} / $$n_mapped_reads) * 1000000));
				    $output_string .=  " mmmr=$mmmr";
				}
				if(exists($miRNA{$map_coord})) {
				    $output_string .= " miRNA\n";
				} elsif (exists($miRNA_star{$map_coord})) {
				    $output_string .= " miRNA-star\n";
				} else {
				    $output_string .= "\n";
				}
			    }
			}
		    }
		}
		
		# now check for any antisense-mapped reads
		for($i = $loc_start; $i <= $loc_stop; ++$i) {
		    if(exists($antisense{$i})) {
			foreach $map_coord (@antisense_names) {
			    if(exists($antisense{$i}{$map_coord})) {
				for($j = $loc_start; $j < $i; ++$j) {
				    $output_string .= "<";
				}
				if($map_coord =~ /(\d+)-(\d+)/) {
				    # get length here as well
				    $read_length = $2 - $1 + 1;
				    for($j = $1; $j <= $2; ++$j) {
					my $letter = $hp_seq_hash{$j};
					$letter =~ tr/ACUG/UGAC/;
					$output_string .= "$letter";
				    }
				} else {
				    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
				}
				for($j = ($2 + 1); $j <= $loc_stop; ++$j) {
				    $output_string .= "<";
				}
				
				$output_string .= " ";
				# print length and abundance
				$output_string .= "l=$read_length ";
				$output_string .= "m=$antisense{$i}{$map_coord}";
				if($$n_mapped_reads) {
				    $mmmr = sprintf("%.4f", (($antisense{$i}{$map_coord} / $$n_mapped_reads) * 1000000));
				    $output_string .= " mmmr=$mmmr\n";
				} else {
				    $output_string .= "\n";
				}
			    }
			}
		    }
		}
	    } elsif ($strand eq "Crick") {
		for($i = $loc_stop; $i >= $loc_start; --$i) {
		    if(exists($brax_hash{$i})) {
			$output_string .= "$brax_hash{$i}";
		    } else {
			$output_string .= " ";
		    }
		}
		$output_string .= "\n";
		for($i = $loc_stop; $i >= $loc_start; --$i) {
		    if(exists($sense{$i})) {
			foreach $map_coord (@sense_names) {
			    if(exists($sense{$i}{$map_coord})) {
				for($j = $loc_stop; $j > $i; --$j) {
				    if(exists($miRNA{$map_coord})) {
					$output_string .= "m";
				    } elsif(exists($miRNA_star{$map_coord})) {
					$output_string .= "\*";
				    } else {
					$output_string .= "\.";
				    }
				}
				if($map_coord =~ /(\d+)-(\d+)/) {
				    # get length here as well
				    $read_length = $1 - $2 + 1;
				    for($j = $1; $j >= $2; --$j) {
					$output_string .= "$hp_seq_hash{$j}";
				    }
				} else {
				    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
				}
				for($j = ($2 - 1); $j >= $loc_start; --$j) {
				    if(exists($miRNA{$map_coord})) {
					$output_string .= "m";
				    } elsif(exists($miRNA_star{$map_coord})) {
					$output_string .= "\*";
				    } else {
					$output_string .= "\.";
				    }
				}
				
				$output_string .= " ";
				# print length and abundance, and if appropriate, designation as a miRNA or miRNA-star
				$output_string .= "l=$read_length ";
				$output_string .= "m=$sense{$i}{$map_coord}";
				if($$n_mapped_reads) {
				    $mmmr = sprintf("%.4f", (($sense{$i}{$map_coord} / $$n_mapped_reads) * 1000000));
				    $output_string .= " mmmr=$mmmr";
				}
				if(exists($miRNA{$map_coord})) {
				    $output_string .= " miRNA\n";
				} elsif (exists($miRNA_star{$map_coord})) {
				    $output_string .= " miRNA-star\n";
				} else {
				    $output_string .= "\n";
				}
			    }
			}
		    }
		}
		
		# now check for any antisense-mapped reads
		for($i = $loc_stop; $i >= $loc_start; --$i) {
		    if(exists($antisense{$i})) {
			foreach $map_coord (@antisense_names) {
			    if(exists($antisense{$i}{$map_coord})) {
				for($j = $loc_stop; $j > $i; --$j) {
				    $output_string .= "<";
				}
				if($map_coord =~ /(\d+)-(\d+)/) {
				    # get length here as well
				    $read_length = $1 - $2 + 1;
				    for($j = $1; $j >= $2; --$j) {
					my $letter = $hp_seq_hash{$j};
					$letter =~ tr/ACUG/UGAC/;
					$output_string .= "$letter";
				    }
				} else {
				    die "error in sub-routine hp_outout : failed to parse map_coord $map_coord\n";
				}
				for($j = ($2 - 1); $j >= $loc_start; --$j) {
				    $output_string .= "<";
				}
				
				$output_string .= " ";
				# print length and abundance
				$output_string .= "l=$read_length ";
				$output_string .= "m=$antisense{$i}{$map_coord}";
				if($$n_mapped_reads) {
				    $mmmr = sprintf("%.4f", (($antisense{$i}{$map_coord} / $$n_mapped_reads) * 1000000));
				    $output_string .= " mmmr=$mmmr\n";
				} else {
				    $output_string .= "\n";
				}
			    }
			}
		    }
		}
	    }
	    
	    # recode data in "output" hash and append to the $locus_data entry
	    # Precision
	    $earlier_fail = 0;  ## reset
	    if(${$output{$sam_query}}[2] >= 1) {
		$locus_data .= "\tPASS";
	    } else {
		$locus_data .= "\tFAIL";
		$earlier_fail = 1;
	    }
	    
	    # Duplex structure
	    if((${$output{$sam_query}}[3] >= 1) and
	       (${$output{$sam_query}}[4] >= 1) and
	       (${$output{$sam_query}}[5] >= 1) and
	       (${$output{$sam_query}}[6] >= 1)) {
		$locus_data .= "\tPASS";
	    } elsif($earlier_fail) {
		$locus_data .= "\tND";
	    } else {
		$locus_data .= "\tFAIL";
		$earlier_fail = 1;
	    }
	    
	    # Star
	    if((${$output{$sam_query}}[7] >= 1) and
	       (${$output{$sam_query}}[8] >= 1) and
	       (${$output{$sam_query}}[9] >= 1)) {
		$locus_data .= "\tPASS";
	    } elsif ($earlier_fail) {
		$locus_data .= "\tND";
	    } else {
		$locus_data .= "\tFAIL";
		$earlier_fail = 1;
	    }
	    
	    # parse the alignments with respect to strandedness.  Demand minimum of $$minstrandfrac fraction to keep
	    # However, "dyad" mappings (map twice, once on each end, on opposite strands, because of perfect inverted-repeat structure) are not tallied when calculating this ratio.
	    
	    $corrected_strand_frac = get_corrected_strand_frac($output_string);
	    if($corrected_strand_frac < $$minstrandfrac) {
		$hp_call = "FAIL";
	    }
	    
	    # check the output string, and suppress it if it is too big (> 400nts long or > 400 distinct small RNAs)
	    # make an exception for MIRNAs, for which we always print all of the alignments.
	    # this should reduce the memory footprint
	    
	    @olines = split ("\n", $output_string);
	    
	    unless($hp_call eq "MIRNA") {
		if((($loc_stop - $loc_start + 1) > 400) or
		   ((scalar @olines) > 403)) {  ## 3 lines are header, brackets, and sequence
		    $output_string = "$olines[0]\n$olines[1]\n$olines[2]\n xxx Alignments Not Shown - too complex xxx\n";
		}
	    }
	    
	    # now add the call and the output
	    
	    $locus_data .= "\t$hp_call";
	    $locus_data .= "\t$output_string";
	    
	    unless($hp_call eq "FAIL") {
		push(@kept_entries, $locus_data);
	    }
	}
	
	# delete old
	delete $$hp_hash{$locus};
	
	# final check on whether to split the original locus.
	$hp_cov_sum = 0;
	$orig_cov = 0;
	foreach $ke (@kept_entries) {
	    @l_data_fields = split ("\t", $ke);
	    $hp_cov_sum += $l_data_fields[5];
	    $orig_cov = $l_data_fields[6];  ## same for every entry, so this is the easiest way to get it
	}
	unless ($orig_cov == 0) { ## if there were no kept entries at this point
	    if(($hp_cov_sum / $orig_cov) < $$minhpfrac) {
		@kept_entries = ();
	    }
	}

	if((scalar @kept_entries) > 0) {
	    @{$$hp_hash{$locus}} = @kept_entries;
	}
    }
    print STDERR " Done\n";  ## closes progress
    return %output;  ## irrelevant but too much of a pain to modify.
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
    
sub quant {
    my($clus_array,$bamfile,$dicer_min,$dicer_max,$strand_cutoff,$dicer_cutoff,$phasesize,$phaseFDR,$names) = @_; ## passed by reference .. first one is array, hp_hash and miR_hash are hashes, others scalars
    my %output = ();
    my %internal = ();
    my $total;
    my $watson;
    my $repnorm;
    my $i;
    my @fields = ();
    my $loc_start;
    my $loc_stop;
    my $read_length;
    my $nh;
    my $frac_watson;
    my $frac_crick;
    my $strand;
    my $total_norm;
    my $repnorm_norm;
    my $size_norm;
    
    my %phase_hash = ();
    my $in_dicer_range;
    my $dicer_call;
    my $max_d_size;
    my $max_d_proportion;
    
    my %phase_p_values = ();
    my %phase_offsets = ();
    my $p_val;
    my $offset;    
    
    my $uniques;
    
    my $x;
    
    # for progress tracking
    my $n_2_analyze = scalar @$clus_array;
    my $five_percent = int($n_2_analyze/20);
    my $n_analyzed = 0;
    print STDERR "\tProgress in sub-routine \"quant\" \(dot = five percent\): ";
    
    foreach my $locus (@$clus_array) {
        ## progress
	++$n_analyzed;
	if($n_analyzed >= $five_percent) {
	    ## TEST
	    ##last;
	    ## END TEST
	    print STDERR "\.";
	    $n_analyzed = 0;
	}
	
	%internal = (); ## reset each time
	%phase_hash = (); ## reset each time
	$total = 0;
	$watson = 0;
	$repnorm = 0;
	$uniques = 0;
	$internal{'short'} = 0;
	$internal{'long'} = 0;
	for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
	    $internal{$i} = 0;
	}
	
	if($locus =~ /:(\d+)-(\d+)/) {
	    $loc_start = $1;
	    $loc_stop = $2;
	} else {
	    die "FATAL in sub-routine \"quant\" : could not parse coordinates from locus $locus\n";
	}
	# initialize phase hash
	for($i = $loc_start; $i <= $loc_stop; ++$i) {
	    $phase_hash{$i} = 0;
	}
	
	open(SAM, "samtools view $$bamfile $locus |");
	while (<SAM>) {
	    chomp;
	    # ignore header lines, should they be present
	    if($_ =~ /^@/) {
		next;
	    }
	    # get fields
	    @fields = split ("\t", $_);
	    # ignore unmapped reads, should they appear
	    if($fields[1] & 4) {
		next;
	    }
	    $read_length = parse_cigar($fields[5]);

	    # determine the number of total mappings for the read from the NH:i: tag.
	    $nh = '';
	    for($i = ((scalar @fields) - 1); $i >= 0; --$i) {
		if($fields[$i] =~ /^NH:i:(\d+)$/) {
		    $nh = $1;
		    last;
		}
	    }
	    unless($nh) {
		die "FATAL in sub-routine \"quant\": Could not parse NH:i flag from sam line $_\n";
	    }
	    
	    # tally
	    ++$total;
	    if($nh == 1) {
		++$uniques;
	    }
	    unless($fields[1] & 16) {
		++$watson;
	    }
	    
	    $repnorm += (1 / $nh);
	    if($read_length < $$dicer_min) {
		++$internal{'short'};
	    } elsif ($read_length > $$dicer_max) {
		++$internal{'long'};
	    } else {
		++$internal{$read_length};
		## add to phase hash
		## make sure to inlcude reads whose left end is outside the locus.. for those simply add the read length to get a correct register
		if($fields[1] & 16) {
		    if(($fields[3] + 2) < $loc_start) {
			++$phase_hash{($fields[3] + 2 + $read_length)};
		    } else {
			++$phase_hash{($fields[3] + 2)};
		    }
		} else {
		    if($fields[3] < $loc_start) {
			++$phase_hash{($fields[3] + $read_length)};
		    } else {
			++$phase_hash{$fields[3]};
		    }
		}
	    }
	}
	close SAM;
	if($total > 0) {
	    $frac_watson = sprintf ("%.4f", ($watson / $total));
	} else {
	    $frac_watson = 0;
	}
	
	$frac_crick = 1 - $frac_watson;
	
	# Evaluate whether this will be annotated a Dicer-derived locus or not, based on the cutoff supplied
	$in_dicer_range = 0; # reset
	for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
	    $in_dicer_range += $internal{$i};
	}
	
	if($total > 0) {
	    if(($in_dicer_range / $total) >= $$dicer_cutoff) {
		# is a dicer locus.  find the dicer-read size with max n reads, and calc proportion
		# reset variables
		$max_d_size = "null";
		$max_d_proportion = 0;
		for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
		    if(($internal{$i} / $total) > $max_d_proportion) {
			$max_d_proportion = $internal{$i} / $total;
			$max_d_size = $i;
		    }
		}
		# shorten the proportion to three decimal places
		$max_d_proportion = sprintf("%.3f",$max_d_proportion);
		## as of version 0.1.2, dicer_call no longer has colon-delimted information, and the proportion is omitted
		#$dicer_call = "$max_d_size" . ":" . "$max_d_proportion"; ## OLD
		$dicer_call = $max_d_size;
	    } else {
		# Not a Dicer locus
		$dicer_call = "N";
	    }
	} else {
	    # shouldn't be here unless somehow $total was zero
	    $dicer_call = "N";
	}
	
	# see whether phasing should be examined
	unless(($$phasesize =~ /^none$/) or
	       ($dicer_call =~ /^N$/)) {
	    if($$phasesize =~  /^all$/) {
		if (($loc_stop - $loc_start + 1) > (4 * $max_d_size)) {
		    ($p_val,$offset) = eval_phasing(\%phase_hash,\$dicer_call,\$loc_start,\$loc_stop);
		    $phase_p_values{$locus} = $p_val;
		    $phase_offsets{$locus} = $offset;
		}
	    } elsif (($dicer_call == $$phasesize) and
		     (($loc_stop - $loc_start + 1) > (4 * $$phasesize))) {
		
		($p_val,$offset) = eval_phasing(\%phase_hash,\$dicer_call,\$loc_start,\$loc_stop);
		$phase_p_values{$locus} = $p_val;
		$phase_offsets{$locus} = $offset;
	    }
	}
	
        # begin entry
	# [0] : locus name
	$output{$locus} .= "$locus";
	
	# [1] : name from name hash
	$output{$locus} .= "\t$$names{$locus}";
	
	# [2] : STRAND .. in this sub-routine, based solely upon $$strand_cutoff.
	
	if($frac_watson >= $$strand_cutoff) {
	    $strand = "+";
	} elsif($frac_crick >= $$strand_cutoff) {
	    $strand = "-";
	} else {
	    $strand = "\.";
	}
	$output{$locus} .= "\t$strand";
	
	# [3] FRAC_WAT
	$output{$locus} .= "\t$frac_watson";
	
	# calculate mappings, reporting only in raw in this sub-routine
	# [4] TOTAL
	$output{$locus} .= "\t$total";
	
	# [5] UNIQUE MAPPERS
	$output{$locus} .= "\t$uniques";
	
	# [6] REP-TOTAL
	$repnorm_norm = sprintf("%.4f",$repnorm);
	$output{$locus} .= "\t$repnorm_norm";
	
	# [7] DICER ... either 'N' or a number within the dicer_min to dicer_max range
	# was calculated above
	$output{$locus} .= "\t$dicer_call";
	
	### [8], [9], and [10] PHASING -- added later
		
	# [11] SHORT
	$output{$locus} .= "\t$internal{'short'}";
	
	# [12] LONG
	$output{$locus} .= "\t$internal{'long'}";
	
	# [13] through whenenver .. Dicer 
	for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
	    $output{$locus} .= "\t$internal{$i}";
	}

	##TEST
	##print STDERR "$output{$locus}\n";
	##END TEST
    }
    
    # Generate the ordered set of Benjamini-Hochberg values (i/m) * FDR
    my @bh_vals = ();
    for($i = 1; $i <= (scalar (keys %phase_p_values)); ++$i) {
	push(@bh_vals,(($i/scalar (keys %phase_p_values)) * $$phaseFDR));
    }
    
    # Generate the ordered set of loci where phasing p-values were determined, ordered in the basis of ascending p-values
    my @phase_sorted_loci = sort {$phase_p_values{$a} <=> $phase_p_values{$b}} (keys %phase_p_values);
    
    # Go in order, and track whether p-value is significant or not
    my %phase_results = ();
    my $ok = 1;
    my $column_entry;
    
    ## TEST
    #print STDERR "\n\nPhasing Results\n";
    
    for($i = 0; $i < (scalar @phase_sorted_loci); ++$i) {
	my $locus = $phase_sorted_loci[$i];
	if($ok) {
	    unless($phase_p_values{$locus} <= $bh_vals[$i]) {
		$ok = 0;
	    }
	    if($ok) {
		$column_entry = "$phase_offsets{$locus}" . ":" . sprintf("%.3e",$phase_p_values{$locus}) . ":" . "OK-FDR-" . "$$phaseFDR";
	    } else {
		$column_entry = "$phase_offsets{$locus}" . ":" . sprintf("%.3e",$phase_p_values{$locus}) . ":" . "NS-FDR-" . "$$phaseFDR";
	    }
	} else {
	    $column_entry = "$phase_offsets{$locus}" . ":" . sprintf("%.3e",$phase_p_values{$locus}) . ":" . "NS-FDR-" . "$$phaseFDR";
	}
	$phase_results{$locus} = $column_entry;
	## TEST
	#print STDERR "$i\t$locus\tBH:$bh_vals[$i]\t$column_entry\n";
	##
    }

    my @old = ();
    my $old_string;
    my $new_string;
    while((my $locus,$old_string) = each %output) {
	@old = split ("\t", $old_string);
	$new_string = '';  ## reset each time
	$new_string .= $old[0]; ## locus
	$new_string .= "\t";
	$new_string .= $old[1]; ## name
	$new_string .= "\t";
	$new_string .= $old[2]; ## strand
	$new_string .= "\t";
	$new_string .= $old[3]; ## frac_wat
	$new_string .= "\t";
	$new_string .= $old[4]; ## total
	$new_string .= "\t";
	$new_string .= $old[5]; ## uniques
	$new_string .= "\t";
	$new_string .= $old[6]; ## rep-total
	$new_string .= "\t";
	$new_string .= $old[7]; ## dicer call
	$new_string .= "\t";
	
	if(exists($phase_results{$locus})) {
	    #$new_string .= $phase_results{$locus};
	    my @pres = split (":", $phase_results{$locus});
	    $new_string .= "$pres[0]\t$pres[1]\t$pres[2]";
	} else {
	    $new_string .= "ND\tND\tND";
	}
	
	for($i = 8; $i < (scalar @old); ++$i) {
	    $new_string .= "\t$old[$i]";
	}
	$output{$locus} = $new_string;
	## TEST
	## print "$new_string\n";
	## END TEST
    }
    print STDERR " Done\n";
    return %output;
}

sub eval_phasing {
    my($phase_hash,$dcall,$loc_start,$loc_stop) = @_; ## passed by reference
    
    my $i;
    my $j;
    
    my $loc_size = ($$loc_stop - $$loc_start + 1);
    my $N;
    my $m;
    my $n;
    my $k;
    
    my $x;
    my $y;
    
    my $term1;
    my $term2;
    my $term2a;
    my $term2b;
    my $denominator;
    
    my $p;

    my $offset;
    
    my $sum = 0;
    
    # determine the phasing size to examine for this locus, which is $i
#    if($$dcall =~ /^(\d+):/) {
#	$i = $1;
#    } else {
#	die "FATAL in sub-routine \'eval_phasing\'  could not parse dicer call $$dcall\n";
#    }
    $i = $$dcall;
    
    # If the locus is larger than 20 * $i, concatenate the phase hash information beyond the first 20 cycles
    my $exam_stop;
    if($loc_size > (20 * $i)) {
	$exam_stop = $$loc_start + (20 * $i) - 1;
	for($x = $$loc_start + (20 * $i); $x <= $$loc_stop; ++$x) {
	    $y = $x;
	    until($y < ($$loc_start + (20 * $i))) {
		$y = $y - (20 * $i);
	    }
	    $$phase_hash{$y} += $$phase_hash{$x};
	}
    } else {
	$exam_stop = $$loc_stop;
    }
    
    # Determine N, the number of all possible states
    $N = $exam_stop - $$loc_start + 1;
    
    
    # determine the mean frequency of occupancy within the exam window
    for($j = $$loc_start; $j <= $exam_stop; ++$j) {
	$sum += $$phase_hash{$j};
    }
    my $mean = $sum / $N;
    
    # count all positions that are above the average within the exam window and also greater than one
    $n = 0;  ## reset
    for($x = $$loc_start; $x <= $exam_stop; ++$x) {
	if(($$phase_hash{$x} > 0) and
	   ($$phase_hash{$x} > $mean)) {
	    ++$n;
	}
    }
    
    # determine the register to examine -- the register with the highest total number of reads
    my $reg_max = 0;
    my $reg_max_offset = 0;
    my $this_total = 0;
    for($j = 0; $j < $i; ++$j) {
	$this_total = 0; ## reset each time	
	for($x = $j + $$loc_start; $x <= $exam_stop; $x += $i) {
	    $this_total += $$phase_hash{$x};
	}
	if($this_total > $reg_max) {
	    $reg_max = $this_total;
	    $reg_max_offset = $j;
	}
    }
    
    # determine m and k in the register of interest within the examination window
    # Use fuzzy phasing, counting -1 and +1 as 'successes' as well as the exact coordinate
    $m = 0; ## reset first
    $k = 0; ## reset first
    # calculate offset
    $offset = $$loc_start + $reg_max_offset;
    for($x = $offset; $x <= $exam_stop; $x += $i) {
	++$m;
	# add flanking nts, as long as they are actually within the window being examined
	if((($x - 1) >= $$loc_start) and
	   (($x - 1) <= $exam_stop)) {
	    ++$m;
	}
	if((($x + 1) >= $$loc_start) and
	   (($x + 1) <= $exam_stop)) {
	    ++$m;
	}
	    
	if(($$phase_hash{$x} > 1) and
	   ($$phase_hash{$x} > $mean)) {
	    ++$k;
	}
	if((($x - 1) >= $$loc_start) and
	   (($x - 1) <= $exam_stop)) {
	    if(($$phase_hash{($x - 1)} > 1) and
	       ($$phase_hash{($x - 1)} > $mean)) {
		++$k;
	    }
	}
	if((($x + 1) >= $$loc_start) and
	   (($x + 1) <= $exam_stop)) {
	    if(($$phase_hash{($x + 1)} > 1) and
	       ($$phase_hash{($x + 1)} > $mean)) {
		++$k;
	    }
	}
    }
    # calculate the p-value
    $p = 0;
    for($y = $k; $y <= $m; ++$y) {
	$term1 = binom_coeff($m,$y);
	$term2a = $N - $m;
	$term2b = $n - $y;
	$term2 = binom_coeff($term2a,$term2b);
	$denominator = binom_coeff($N,$n);
	$p += ($term1 * $term2) / $denominator;
    }
    
    ## TEST
    #print "loc_start: $$loc_start loc_stop: $$loc_stop reg_max: $reg_max reg_max_offset: $reg_max_offset offset: $offset N: $N n: $n m: $m k: $k p: $p\n";
    #print "$i\t$N\t$n\t$m\t$k\t$p\n";
    ## END TEST
    
    # return
    return($p,$offset);
}

sub binom_coeff {
    my($n,$k) = @_; ## thanks PERL monks!
    my $r=1;
    $r*=$n/($n-$k),$n--while$n>$k;
    return $r;
}

sub mod_quant_hp {
    my($quant_hash,$hp_hash) = @_;  ## passed by reference
    my $locus;
    my $entry;
    my @fields = ();
    my @new_fields = ();
    my $new_entry;
    my %output = ();
    while(($locus,$entry) = each %$quant_hash) {
	@fields = split ("\t", $entry);
	@new_fields = @fields[3..((scalar @fields) - 1)];
	if(exists($$hp_hash{$locus})) {
	    # polarity of the hairpin supersedes that calcuated by the reads alone
	    if($$hp_hash{$locus} =~ /Watson/) {
		unshift(@new_fields,"+");
	    } elsif ($$hp_hash{$locus} =~ /Crick/) {
		unshift(@new_fields,"-");
	    } else {
		die"FATAL error in sub-routine mod_quant_hp : could find strand for HP locus $locus in entry $$hp_hash{$locus}\n";
	    }
	    
	    # now check if its a MIRNA or not
	    if($$hp_hash{$locus} =~ /MIRNA/) {
		unshift(@new_fields,"MIRNA");
	    } else {
		unshift(@new_fields,"HP");
	    }
	} else {
	    # not an HP, leave the original strand alone, and enter "." for the HP designation
	    unshift(@new_fields,$fields[2]);
	    unshift(@new_fields,"\.");
	}
	
	# add the name and the locus
	unshift(@new_fields,$fields[1]);
	unshift(@new_fields,$locus);
	
	# add the new entry
	$new_entry = join("\t", @new_fields);
	$output{$locus} = $new_entry;
    }
    return %output;
}

sub convert_2_mpmm {
    my($in_hash,$reads) = @_;  ## passed by reference
    my $locus;
    my $input_string;
    my %output = ();
    my @old_fields = ();
    my @new_fields = ();
    my $i;
    while(($locus,$input_string) = each %$in_hash) {
	@new_fields = ();  ## reset each time though
	@old_fields = split ("\t", $input_string);
	# 0: locus, 1: name, 2: HP, 3: strand, 4: frac_wat, 5: total, 6: uniques 7: rep-total, 8: dicercall, 9: phasing_offset, 10: phase p-val, 11; phase_FDRcall,  12: short, 13: long, 14 to end, dicer
	for($i = 0; $i <= 4; ++$i) {
	    push(@new_fields,$old_fields[$i]);
	}
	# total
	push(@new_fields, sprintf("%.3f",(($old_fields[5] / $$reads) * 1000000)));
	# uniques
	push(@new_fields, sprintf("%.3f",(($old_fields[6] / $$reads) * 1000000)));
	# rep-total
	push(@new_fields, sprintf("%.3f",(($old_fields[7] / $$reads) * 1000000)));
	# dicer call and phasing
	push(@new_fields,$old_fields[8]);
	push(@new_fields,$old_fields[9]);
	push(@new_fields,$old_fields[10]);
	push(@new_fields,$old_fields[11]);
	# short and long
	push(@new_fields, sprintf("%.3f",(($old_fields[12] / $$reads) * 1000000)));
	push(@new_fields, sprintf("%.3f",(($old_fields[13] / $$reads) * 1000000)));
	# the rest
	for($i = 14; $i < (scalar @old_fields); ++$i) {
	    push(@new_fields, sprintf("%.3f",(($old_fields[$i] / $$reads) * 1000000)));
	}
	my $result = join("\t",@new_fields);
	$output{$locus} = $result;
	
    }
    return %output;
}

sub write_bed_nohp {
    my($clus_array,$quant_hash,$outdir,$dicermin,$dicermax,$names) = @_; ## passed by reference

    # open file
    my $bedfile = "$$outdir" . "\/" . "ShortStack\.bed";
    open(BED, ">$bedfile");
    print BED "track name=ShortStack itemRgb=\"On\"\n";

    my $i;
    my $j;
    my $text_to_return;
    # non-dicer always gray
    print STDERR "\tColor-Scheme in bed file $bedfile : \n";
    $text_to_return .= "\tColor-Scheme in bed file $bedfile : \n";
    my %colors = ();
    $colors{'N'} = "169,169,169";  ## dark gray
    print STDERR "\tNon-Dicer Clusters: Dark Gray RGB: 169,169,169\n";
    $text_to_return .= "\tNon-Dicer Clusters: Dark Gray RGB: 169,169,169\n";
    
    my $n_colors_to_get = $$dicermax - $$dicermin + 1;
    my %rgb_hash = get_colors($n_colors_to_get);
    my $x = 0;
    
    for ($i = $$dicermin; $i <= $$dicermax; ++$i) {
	++$x;
	$colors{$i} = $rgb_hash{$x};
	print STDERR "\t$i Clusters: RGB: $rgb_hash{$x}\n";
	$text_to_return .= "\t$i Clusters: RGB: $rgb_hash{$x}\n";
    }

    my @fcfields = ();
    my $chrom;
    my $bedstart;
    my $bedstop;
    my $name;
    foreach my $f_clus (@$clus_array) {
	@fcfields = split("\t",$$quant_hash{$f_clus});
	if($fcfields[0] =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chrom = $1;
	    $bedstart = $2 - 1;
	    $bedstop = $3;
	} else {
	    die "FATAL: in sub-routine \"write_bed_nohp\" : Could not parse locus name $fcfields[0]\n";
	}
	print BED "$chrom\t$bedstart\t$bedstop\t";
	
	# name
	print BED "$$names{$f_clus}\t";
	
	# score (ignored -- just zero)
	print BED "0\t";
	
	# strand, which is $fcfields[3]
	print BED "$fcfields[3]\t";
	
	# thickStart and thickEnd, which are always the same as bedstart and bedstop
	print BED "$bedstart\t$bedstop\t";
	
	# Determine the dominant size of the cluster from $fcfields[8] and write correct color
	if($fcfields[8] eq "N") {
	    print BED "$colors{'N'}\t";
	} elsif ($fcfields[8] =~ /^(\d+)$/) {
	    print BED "$colors{$1}\t";
	} else {
	    die "FATAL in sub-routine \"write_bed_nohp\" : Could not find the cluster size from entry $fcfields[8]\n";
	}
	
	# in no-hp mode, block count is always 1
	print BED "1\t";
	
	# the block size is the locus size
	my $block_size = $bedstop - $bedstart;
	print BED "$block_size\t";
	
	# the start of the block is the beginning of the locus
	print BED "0\n";
    }
    close BED;
    return($text_to_return);
}

sub write_bed_with_hp {
    my($clus_array,$quant_hash,$outdir,$dicermin,$dicermax,$names,$hp_hash) = @_; ## passed by reference

    # open file
    my $bedfile = "$$outdir" . "\/" . "ShortStack\.bed";
    open(BED, ">$bedfile");
    print BED "track name=ShortStack itemRgb=\"On\"\n";
    
    # determine color scheme .. ROY G BIV, unless the dicer range exceeds seven colors.  in that case, just use red for all dicer-sized ones
    my %colors = ();
    my $i;
    my $j;
    my $text_to_return;
    # non-dicer always gray
    print STDERR "\n\tColor-Scheme in bed file $bedfile : \n";
    $text_to_return .= "\n\tColor-Scheme in bed file $bedfile : \n";
    $colors{'N'} = "169,169,169";  ## dark gray
    print STDERR "\tNon-Dicer Clusters: Dark Gray RGB: 169,169,169\n";
    $text_to_return .= "\tNon-Dicer Clusters: Dark Gray RGB: 169,169,169\n";

    my $n_colors_to_get = $$dicermax - $$dicermin + 1;
    my %rgb_hash = get_colors($n_colors_to_get);
    my $x = 0;
    
    for ($i = $$dicermin; $i <= $$dicermax; ++$i) {
	++$x;
	$colors{$i} = $rgb_hash{$x};
	print STDERR "\t$i Clusters: RGB: $rgb_hash{$x}\n";
	$text_to_return .= "\t$i Clusters: RGB: $rgb_hash{$x}\n";
    }

    my @fcfields = ();
    my $chrom;
    my $bedstart;
    my $bedstop;
    my $name;
    foreach my $f_clus (@$clus_array) {
	@fcfields = split("\t",$$quant_hash{$f_clus});
	if($fcfields[0] =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chrom = $1;
	    $bedstart = $2 - 1;
	    $bedstop = $3;
	} else {
	    die "FATAL: in sub-routine \"write_bed_nohp\" : Could not parse locus name $fcfields[0]\n";
	}
	print BED "$chrom\t$bedstart\t$bedstop\t";
	
	#name
	print BED "$$names{$f_clus}";
	# Check $fcfields[2]  -- if it is 'HP' or 'MIRNA', then
	if($fcfields[2] eq "HP") {
	    print BED "_HP\t";
	} elsif ($fcfields[2] eq "MIRNA") {
	    print BED "_MIRNA\t";
	} else {
	    print BED "\t";
	}
	
	# score (ignored -- just zero)
	print BED "0\t";
	
	# strand, which is $fcfields[3]
	print BED "$fcfields[3]\t";
	
	# thickStart and thickEnd, which are always the same as bedstart and bedstop
	print BED "$bedstart\t$bedstop\t";
	
	# Determine the dominant size of the cluster from $fcfields[8] and write correct color
	if($fcfields[8] eq "N") {
	    print BED "$colors{'N'}\t";
	} elsif ($fcfields[8] =~ /^(\d+)$/) {
	    print BED "$colors{$1}\t";
	} else {
	    die "FATAL in sub-routine \"write_bed_nohp\" : Could not find the cluster size from entry $fcfields[8]\n";
	}
	
	if(($fcfields[2] eq "HP") or
	   ($fcfields[2] eq "MIRNA")) {
	    
	    # block count will be two, representing the two arms of the hp
	    print BED "2\t";
	    # retrieve information on the arm locations
	    my $left_block_size;
	    my $right_block_size;
	    my $left_begin;
	    my $right_begin;
	    if($$hp_hash{$f_clus} =~ /\t(\d+)-(\d+),(\d+)-(\d+)\t/) {
		if($1 < $2) {
		    $left_block_size = $2 - $1 + 1;
		    $right_block_size = $4 - $3 + 1;
		    $left_begin = $1 - $bedstart;
		    $right_begin = $3 - $bedstart;
		} else {
		    $left_block_size = $3 - $4 + 1;
		    $right_block_size = $1 - $2 + 1;
		    $left_begin = $4 - $bedstart;
		    $right_begin = $2 - $bedstart;
		}
		print BED "$left_block_size";
		print BED ",";
		print BED "$right_block_size\t";
		print BED "$left_begin";
		print BED ",";
		print BED "$right_begin\n";
	    } else {
		die "FATAL: in sub-routine write_bed_with_hp : could not parse hairpin information with regex from $$hp_hash{$f_clus}\n";
	    }
	} else {
	    print BED "1\t";
	    # the block size is the locus size
	    my $block_size = $bedstop - $bedstart;
	    print BED "$block_size\t";
	    
	    # the start of the block is the beginning of the locus
	    print BED "0\n";
	}
    }
    close BED;
    return($text_to_return);
}
	
sub get_names_simple {
    my($clusters) = @_;
    my %names = ();
    my $n = 0;
    foreach my $loc (@$clusters) {
	++$n;
	$names{$loc} = "Cluster_" . "$n";
    }
    return %names;
}


sub final_summary_nohp {
    my ($clus,$phasesize)  = @_;
    my %hash = ();
    my @fields = ();
    my $size;
    foreach my $line (keys %$clus) {
	@fields = split ("\t", $$clus{$line});
	
	# dicer call
	if($fields[8] eq "N") {
	    ++$hash{$fields[8]};
	} else {
	    if($fields[8] =~ /^(\d+)$/) {
		$size = $1;
	    } else {
		die "FATAL in sub-routine final_summary_hp : failed to parse dicer size category of $fields[8]\n";
	    }
	    ++$hash{$size};
	}
	
	# phased?
	if($fields[11] =~ /OK/) {
	    ++$hash{'phased'};
	}
    }
    return %hash;
}

sub final_summary {
    my ($clus,$phasesize)  = @_;
    my %hash = ();
    my @fields = ();
    my $hp_call;
    my $dcall;
    foreach my $line (keys %$clus) {
	@fields = split ("\t", $$clus{$line});
	
	# dicer call
	if($fields[8] eq "N") {
	    $dcall = $fields[8];
	} else {
	    if($fields[8] =~ /^(\d+)$/) {
		$dcall = $1;
	    } else {
		die "FATAL in sub-routine final_summary_hp : failed to parse dicer size category of $fields[8]\n";
	    }
	}
	
	# HP call
	if(($fields[2] eq "MIRNA") or
	   ($fields[2] eq "HP")) {
	    $hp_call = "$fields[2]";
	} else {
	    $hp_call = "X";
	}

	++$hash{$dcall}{$hp_call};
	
	# phased?
	if($fields[11] =~ /OK/) {
	    ++$hash{'phased'};
	}
    }
    return %hash;
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

sub get_colors {
    my($n_to_get) = @_;
    my %rgb_hash = ();
    # initial conditions
    my $r = 255;
    my $g = 0;
    my $b = 0;
    my $x = 1;
    $rgb_hash{$x} = "$r" . "," . "$g" . "," . "$b";
    # increase green
    for(my $i = 1; $i <= 255; ++$i) {
	++$g;
	++$x;
	$rgb_hash{$x} = "$r" . "," . "$g" . "," . "$b";
    }
    # decrease red
    for(my $i = 254; $i >=0; --$i) {
	--$r;
	++$x;
	$rgb_hash{$x} = "$r" . "," . "$g" . "," . "$b";
    }
    # increase blue
    for(my $i = 1; $i <= 255; ++$i) {
	++$b;
	++$x;
	$rgb_hash{$x} = "$r" . "," . "$g" . "," . "$b";
    }
    # decrease green
    for(my $i = 254; $i >=0; --$i) {
	--$g;
	++$x;
	$rgb_hash{$x} = "$r" . "," . "$g" . "," . "$b";
    }
    
    my %answer = ();

    # determine even spacing on the $x number line
    my $denom = $n_to_get - 1;
    if($denom <= 0) {
	$answer{'1'} = $rgb_hash{'1'};
    } else {
	my $spacing = int ($x / $denom);
	my $position = 1;
	$answer{$position} = $rgb_hash{$position};
	for(my $y = 2; $y <= $n_to_get; ++$y) {
	    $position += $spacing;
	    if($position > $x) {
		$position = $x;
	    }
	    $answer{$y} = $rgb_hash{$position};
	}
    }
    
    return %answer;
}
	    

sub write_files {
    my($hp_clusters,$outdir,$names,$final_clusters) = @_;  ## passed by reference, hash, scalar, hash, array
    my $n_miRNAs = 0;
    my $n_hpRNAs = 0;
    
    # make the sub-directories
    system "mkdir $$outdir/HP_details";
    system "mkdir $$outdir/MIRNA_details";
    
    # open table files for each and print their headers
    open(MIR_TABLE, ">$$outdir/miRNA_summary\.txt");
    print MIR_TABLE "\#Locus\tName\tmiRNA\tmiRNA_mappings\tmiRNA-star\tmiRNA-star_mappings\tTotal_mappings\n";
    open(HP_TABLE, ">$$outdir/hpRNA_summary\.txt");
    print HP_TABLE "\#Locus\tName\tPrecision\tDuplex\tStar\n";

    my @fields = ();
    my $entry;
    my $locus;
    my $name;
    my @file_lines = ();
    my $file_line;
    
    my @mirs = ();
    my @mir_counts = ();
    my @stars = ();
    my @star_counts = ();
    my $total;
    my $mappings;
    
    my $mir_out;
    my $mir_count_out;
    my $star_out;
    my $star_count_out;
	
    
    foreach $locus (@$final_clusters) {
	if(exists($$hp_clusters{$locus})) {
	    $name = $$names{$locus};
	    @fields = split ("\t", $$hp_clusters{$locus});
	    if($fields[10] eq "MIRNA") {
		++$n_miRNAs;
		open(OUT, ">$$outdir/MIRNA_details/$name\.txt");
		print OUT "$name $fields[11]\n";
		close OUT;
		
		@file_lines = split ("\n", $fields[11]);
		@mirs = ();
		@mir_counts = ();
		@stars = ();
		@star_counts = ();
		$total = 0;
		
		foreach $file_line (@file_lines) {
		    if($file_line =~ /m=(\d+)/) {
			$total += $1;
			$mappings = $1;
			if($file_line =~ /([AUGC]{15,}).*miRNA$/) {
			    push(@mir_counts, $mappings);
			    push(@mirs, $1);
			}
			if($file_line =~ /([AUGC]{15,}).*miRNA-star$/) {
			    push(@star_counts, $mappings);
			    push(@stars, $1);
			}
		    }
		}
		print MIR_TABLE "$locus\t$name\t";
		if((scalar @mirs) > 1) {
		    $mir_out = join (",", @mirs);
		    $mir_count_out = join (",", @mir_counts);
		    print MIR_TABLE "$mir_out\t$mir_count_out\t";
		} else {
		    print MIR_TABLE "$mirs[0]\t$mir_counts[0]\t";
		}
		
		if((scalar @stars) > 1) {
		    $star_out = join(",", @stars);
		    $star_count_out = join(",", @star_counts);
		    print MIR_TABLE "$star_out\t$star_count_out\t";
		} else {
		    print MIR_TABLE "$stars[0]\t$star_counts[0]\t";
		}
		print MIR_TABLE "$total\n";
		
	    } elsif ($fields[10] eq "HP") {
		++$n_hpRNAs;
		open(OUT, ">$$outdir/HP_details/$name\.txt");
		print OUT "$name $fields[11]\n";
		close OUT;
		print HP_TABLE "$locus\t$name\t$fields[7]\t$fields[8]\t$fields[9]\n";
	    } else {
		die "FATAL in sub-routine write_files: failed to understand following entry as MIRNA or HP  found $fields[10] instead:\n$$hp_clusters{$locus}\n";
	    }
	}
    }
    
    # close the tables
    close MIR_TABLE;
    close HP_TABLE;
    
    # delete everything if there were no miRNAs and hpRNAs
    if(($n_miRNAs + $n_hpRNAs) == 0) {
	system "rm -f -r $$outdir/MIRNA_details";
	system "rm -f -r $$outdir/HP_details";
    }
    
    return($n_miRNAs,$n_hpRNAs);
}

sub parse_inv {
    my($inv_file,$minntspaired,$minfracpaired) = @_;  ## passed by reference, all scalars
    my @output = ();
    my $outline;
    my $chr;
    
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    my $left_string;
    my $right_string;
    
    my $i;
    my $j;
    my $k;
    my %left_keyed = ();
    my %right_keyed = ();
    my @left_st = ();
    my @right_st = ();
    
    my $w_brax;
    my $w_pairs;
    my $w_frac_paired;
    
    my $c_brax;
    my $c_pairs;
    my $c_frac_paired;
    
    my $is_pair;
    my $strand;
    
    my @single = ();
    
    my $top_bases;
    my $bottom_bases;
    
    open(INV, "$$inv_file");
    while (<INV>) {
	chomp;
	if($_ =~ /\S/) {
	    push (@single, $_);
	    
	    if((scalar @single) == 4) {
		
		# analyze
		if($single[0] =~ /^(\S+):/) {
		    $chr = $1;
		} else {
		    print STDERR "\n\tWARNING in sub-routine parse_inv: Failed to parse chr name from line $single[0] and skipped the entry\n\n";
		    @single = ();
		    next;
		}
		if($single[1] =~ /(\d+) (\S+) (\d+)/) {
		    $left_start = $1;
		    $left_string = $2;
		    $left_stop = $3;
		} else {
		    print STDERR "\n\tWARNING in sub-routine parse_inv: Failed to parse expected left arm info from line $single[1] and skipped the entry\n\n";
		    @single = ();
		    next;
		}
		if($single[3] =~ /(\d+) (\S+) (\d+)/) {
		    $right_stop = $1;
		    $right_string = $2;  ## top strand, reading 3' to 5' i.e. descending
		    $right_start = $3;
		}
		
		# sanity check:  einverted seems to have a rare bug
		$top_bases = 0;
		$bottom_bases = 0;
		while($left_string =~ /[atgc]/g) {
		    ++$top_bases;
		}
		while($right_string =~ /[atgc]/g) {
		    ++$bottom_bases;
		}
		if(($top_bases != ($left_stop - $left_start + 1))) {
		    print STDERR "\n\tWARNING in sub-routine parse_inv: Possible einverted bug -- On top strand of entry below there are $top_bases nucleotides but the coordinates $chr $left_start to $left_stop do not concur\n\t\tSkipping entry $chr $left_start to $right_stop\n\n";
		    @single = ();
		    next;
		}
		if(($bottom_bases != ($right_stop - $right_start + 1))) {
		    print STDERR "\n\tWARNING in sub-routine parse_inv: Possible einverted bug -- On bottom strand of entry below there are $bottom_bases nucleotides but the coordinates $chr $right_start to $right_stop do not concur\n\t\tSkipping entry $chr $left_start to $right_stop\n\n";
		    @single = ();
		    next;
		}
		
		# convert t's to u's
		$left_string =~ s/t/u/g;
		$right_string =~ s/t/u/g;
		
		# explode them to arrays
		@left_st = split ('', $left_string);
		@right_st = split ('', $right_string);
		
		# key in the sequences of the helix, including gapped positions
		# reset first
		%left_keyed = ();
		%right_keyed = ();
		$i = $left_start;
		$j = $right_stop;
		for($k = 0; $k < (length $left_string); ++$k) {
		    unless($left_st[$k] eq "-") {
			if($right_st[$k] eq "-") {
			    $left_keyed{$i} = "$left_st[$k]\tgap\tgap";
			} else {
			    $left_keyed{$i} = "$left_st[$k]\t$j\t$right_st[$k]";
			}
		    }
		    unless($right_st[$k] eq "-") {
			if($left_st[$k] eq "-") {
			    $right_keyed{$j} = "$right_st[$k]\tgap\tgap";
			} else {
			    $right_keyed{$j} = "$right_st[$k]\t$i\t$left_st[$k]";
			}
		    }
		    unless($left_st[$k] eq "-") {
			++$i;
		    }
		    unless($right_st[$k] eq "-") {
			--$j;
		    }
		}
		
		# Watson first
		$w_brax = '';  ## reset
		$w_pairs = 0;
		$strand = "Watson";
		for($i = $left_start; $i <= $left_stop; ++$i) {
		    $is_pair = is_it_paired(\$left_keyed{$i},\$strand);
		    if($is_pair) {
			$w_brax .= "\(";
			++$w_pairs;
		    } else {
			$w_brax .= "\.";
			}
		}
		for($i = ($left_stop + 1); $i < $right_start; ++$i) {
		    $w_brax .= "\.";
		}
		for($i = $right_start; $i <= $right_stop; ++$i) {
		    $is_pair = is_it_paired(\$right_keyed{$i},\$strand);
		    if($is_pair) {
			$w_brax .= "\)";
		    } else {
			$w_brax .= "\.";
		    }
		}
		
		# Now Crick
		$c_brax = '';
		$c_pairs = 0;
		$strand = "Crick";
		for($i = $right_stop; $i >= $right_start; --$i) {
		    $is_pair = is_it_paired(\$right_keyed{$i},\$strand);
		    if($is_pair) {
			$c_brax .= "\(";
			++$c_pairs;
		    } else {
			$c_brax .= "\.";
		    }
		}
		for($i = ($right_start - 1); $i > $left_stop; --$i) {
			$c_brax .= "\.";
		}
		for($i = $left_stop; $i >= $left_start; --$i) {
		    $is_pair = is_it_paired(\$left_keyed{$i},\$strand);
		    if($is_pair) {
			$c_brax .= "\)";
		    } else {
			$c_brax .= "\.";
		    }
		}
		
		# assemble entries after testing
		# Watson first
		$outline = ''; ## reset
		$strand = "Watson";
		$w_frac_paired = $w_pairs / (0.5 * (($left_stop - $left_start + 1) + ($right_stop - $right_start + 1)));
		if(($w_frac_paired >= $$minfracpaired) and
		   ($w_pairs >= $$minntspaired) and
		   (($right_start - $left_stop + 1) < (0.5 * (($left_stop- $left_start +1) + ($right_stop - $right_start + 1))))) {
		    $outline = "$w_brax\t$left_start" . "-" . "$right_stop\t$left_start" . "-" . "$left_stop" . "," . "$right_start" . "-" . "$right_stop\t$strand\t$chr";
		    
		    push(@output, $outline);
		}
		
		# now Crick
		$outline = '';
		$strand = "Crick";
		$c_frac_paired = $c_pairs / (0.5 * (($left_stop - $left_start + 1) + ($right_stop - $right_start + 1)));
		if(($c_frac_paired >= $$minfracpaired) and
		   ($c_pairs >= $$minntspaired) and
		   (($right_start - $left_stop + 1) < (0.5 * (($left_stop- $left_start +1) + ($right_stop - $right_start + 1))))) {
		    $outline = "$c_brax\t$right_stop" . "-" . "$left_start\t$right_stop" . "-" . "$right_start" . "," . "$left_stop" . "-" . "$left_start\t$strand\t$chr";

		    push(@output, $outline);
		}
		@single = ();
	    }
	}
    }
    close INV;
    return (@output);
}

sub is_it_paired {
    my($string,$strand) = @_;  ## passed by reference
    my @fields = split ("\t", $$string);
    my $paired = 0;
    my $actual_one = $fields[0];
    my $actual_two = $fields[2];
    
    if($$strand eq "Crick") {
	$actual_one =~ s/augc/uacg/g;
	$actual_two =~ s/uacg/augc/g;
    }
    
    if($actual_one eq "a") {
	if($actual_two eq "u") {
	    $paired = 1;
	}
    } elsif ($actual_one eq "u") {
	if(($actual_two eq "a") or ($actual_two eq "g")) {
	    $paired = 1;
	}
    } elsif ($actual_one eq "g") {
	if(($actual_two eq "c") or ($actual_two eq "u")) {
	    $paired = 1;
	}
    } elsif ($actual_one eq "c") {
	if($actual_two eq "g") {
	    $paired = 1;
	}
    }
    return $paired;
}
    
sub merge_inv {
    my($irs,$true_hps,$clusters) = @_; ## by reference, array, hash, array
    my $c_chr;
    my @c_fields = ();
    my @c_coords = ();
    my $ir_entry;
    my @ir_fields = ();
    my @ir_subfields = ();
    my @ir_lefts = ();
    my @ir_rights = ();
    my $left_overlap;
    my $right_overlap;
    my $junk;
    my $new_entry;
    my $in = scalar @$irs;
    my $merged = 0;
    foreach my $clus (@$clusters) {
	@c_fields = split (":", $clus);
	$c_chr = $c_fields[0];
	@c_coords = split ("-", $c_fields[1]);
	foreach $ir_entry (@$irs) {
	    @ir_fields = split ("\t", $ir_entry);
	    ## only continue if the two are on the same chromosome
	    if($ir_fields[4] eq $c_chr) {
		@ir_subfields = split (",", $ir_fields[2]);
		@ir_lefts = split ("-", $ir_subfields[0]);
		@ir_rights = split ("-", $ir_subfields[1]);
		$left_overlap = range_overlap(\@c_coords,\@ir_lefts);
		$right_overlap = range_overlap(\@c_coords,\@ir_rights);
		if(($left_overlap) or ($right_overlap)) {
		    $junk = pop @ir_fields; ## remove Chr from the entry, to conform with the hash style
		    $new_entry = join ("\t", @ir_fields);
		    ++$merged;
		    push(@{$$true_hps{$clus}}, $new_entry);
		}
	    }
	}
    }
    return($in,$merged);
}
		    
sub get_corrected_strand_frac {
    my($output_string) = @_;
    my @output_lines = split ("\n", $output_string);
    my $sense_mappings = 0;
    my $antisense_mappings = 0;
    my %sense_seqs = ();
    my $real_mappings = 0;
    my $this_seq;
    my $this_n;
    my $revseq;
    foreach my $line (@output_lines) {
	if($line =~ /([AUGC]+).* m=(\d+)/) {
	    $this_seq = $1;
	    $this_n = $2;
	    ++$real_mappings;
	    if($line =~ /</) {
		$revseq = reverse $this_seq;
		if(exists($sense_seqs{$revseq})) {
		    $sense_seqs{$revseq} -= $this_n;
		} else {
		    $antisense_mappings += $this_n;
		}
	    } else {
		$sense_seqs{$this_seq} += $this_n;
	    }
	}
    }
    my $corrected_frac;
    while(($this_seq, $this_n) = each %sense_seqs) {
	$sense_mappings += $this_n;
    }
    my $cor_sum = $sense_mappings + $antisense_mappings;
    if($cor_sum == 0) {
	if($real_mappings > 0) {
	    $corrected_frac = 1;  ## this will occur when all reads are in dyads .. e.g. perfect IR.  we'll keep those!
	} else {
	    $corrected_frac = 0;  ## this is a failsafe in case the sub-routine is passed an entry with no reads at all
	}
    } else {
	$corrected_frac = $sense_mappings / $cor_sum;
    }
    return $corrected_frac;
}
    
    

__END__
=head1 LICENSE

ShortStack.pl

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

Annotation and quantification of small RNA genes based upon reference-aligned small RNA sequences

=head1 CITATION

If you use ShortStack in your work, please cite 

Axtell MJ. (2012) ShortStack: Comprehensive annotation and quantification of small RNA genes.  In prep.

A manuscript describing the ShortStack package has been written and submitted as of late June 2012, so check Pubmed first or look for an update!

=head1 VERSIONS

0.2.2 : September 23, 2012.  THIS VERSION.  Minor changes to HP and MIRNA calling during countmode to increase the consistency of calls between de novo and count mode runs.

0.2.1 : September 14, 2012.  Speed improvements during the merging of einverted-derived hairpins with RNALfold-derived hairpins.

0.2.0 : August 29, 2012.  Major update.  Annotation of non-miRNA hairpin RNAS (hpRNAs) is much improved, with many fewer "borderline" cases.  Multiple internal modificaitons changed the method for dealing with these loci.  In addition, as of this version, ShortStack takes in a file of inverted repeats in the form of a .inv file produced by einverted (from the EMBOSS package).  Using inverted repeats complements the usage of RNALfold to find hairpins.  In particular, the use of einverted-derived inverted repeat annotations enables the correct annotation of very large hairpin-derived small RNA genes, such as Arabidopsis thaliana IR71.

0.1.4 : June 28, 2012.   Prep_bam.pl helper script updated to deal with .sam lines corresponding to unmapped reads without breaking.  This also enables "one-step" mapping and .bam file preparation using bowtie (see 'Quick Start' section below).  ShortStack.pl also updated to prevent breaking upon encountering .sam lines corresponding to unmapped reads.  Coloring of .bed files now dynamically adjusts to make a nice 'rainbow' (from red through green to blue) regardless of the size of the dicer range.  Finally, an additional filter was added to prevent excessive splitting of initial de-novo clusters into hairpin-associated clusters ... initial clusters are now prevented from splitting if they would spawn more than --maxsplit child hairpins (default, 3).  This prevents inappropriate splitting of very large clusters that might harbor a number of sub-optimal hairpins.

0.1.3 : June 12 2012.  Critical bug fix:  Fixed issue in analysis of RNALfold-derived secondary structures that was causing a substantial rate of acceptable hairpin structures to be incorrectly rejected .. thus under-reporting hairpin-and MIRNA-derived loci.  Additionally, removed the option --maxtofoldwindow ; all loci are now folded and analyzed unless running in nohp mode.  Other smaller fixes and changes including fixing the column names on the "Hairpin-MIRNA_summary.txt" file ... previously, the "Name" and "Locus" columns lacked headers.

0.1.2 : May 17, 2012.  Changes to simplify downstream analysis of the Results.txt file, and to phasing analysis.  Details: A) changed format of 'DicerCall' column in result to remove the colon-delimted fractions .. now it simply reports size or 'N'  B) Changed format of phasing analysis results to split the phase offset, p-value, and false-discovery rate call into three separate columns.  C) added option of --phasesize none to suppress analysis of phasing at all clusters.

0.1.1 : May 4, 2012.  Added helper script "miR_homologs.pl" to package.  No change to ShortStack.pl code itself except version change.

0.1.0 : Initial release. April 29, 2012

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 INSTALL

install samtools from <http://samtools.sourceforge.net/> and ensure that samtools is in your PATH

install the ViennaRNA package VERSION 1.7.X or 1.8.X <http://www.tbi.univie.tion run, including inverted repeats file:

ShortStack.pl --inv_file Athaliana_167.inv col_leaf_ok.bam Athaliana_167.fa

B) count mode run to quantify and annotate known miRBase MIRNA loci:

ShortStack.pl --count ath_hp_mb19_SStack_Athal_167.txt col_leaf_ok.bam Athaliana_167.fa

C) Analyze annotated miRBase mature miRNAs for acceptable structure with miR_homologs.pl:

miR_homologs.pl ath_mature_nr_mb19.bam Athaliana_167.fa


=head1 OPTIONS

--outdir [string] : Name of directory to be created to receive results of the run.  Deafults to "ShortStack_[time]", where time is the number of non-leap seconds since Jan 1, 1970 UCT, if not provided   
                                      
--reads [integer] : Number of reads (NOT mappings) in the input .bam file.  No default.  Reads are required to output quantifications in mappings per million mapped, instead of in raw rads.  If not provided, the run will be forced into "--raw" mode, because mappings per million mapped reads cannot be calculated.

--inv_file [string] : PATH to an einverted-produced .inv file of inverted repeats within the genome of interest.  Not required but strongly suggested for more complete annotations of hairpin-derived small RNA genes.  Default = {blank}.  Not needed for runs in "nohp" mode.  A typical eniverted run uses default parameters except "-maxrepeat 10000", in order to capture long IRs.

--mindepth [integer] : Minimum depth of mapping coverage to define an 'island'.  Default = 20.  Must be at least 2, more than 5 preferred.

--pad [integer] : Number of nucleotides upstream and downstream to extend initial islands during cluster definition.  Default = 100

--dicermin [integer] : Smallest size in the Dicer size range (or size range of interest).  Deafult = 20.  Must be between 15 and 35, and less than or equal to --dicermax

--dicermax [integer] : Largest size in the Dicer size range (or size range of interest).  Deafult = 24.  Must be between 15 and 35, and more than or equal to --dicermin

--maxhpsep [integer] : Maximum allowed span for a base-pair during hairpin search with RNALfold.  Default = 300.  Must be between 50 and 2000.

--minfracpaired [float] : Minimum fraction of paired nucleotides required within a valid hairpin structure.  Default = 0.67.  Allowed values are greater than 0 and less than or equal to 1.

--minntspaired [integer] : Minimum absolute number of paired nucleotides required within a valid hairpin structure.  Default = 30.  Allowed values are greater than zero and less than or equal to --maxhpsep

--minfrachpdepth [float] : Minimum fraction of nearby coverage within hairpin arms to keep hairpin for later miRNA analysis.  Default = 0.67.  Allowed values between 0 and 1.  See below for details.

--minstrandfrac [float] : Minimum fraction of mappings to one or the other strand call a polarity for non-hairpin clusters.  Also the minimum fraction of "non-dyad" mappings to the sense strand within potential hairpins/miRNAs to keep the locus annotated as a hp or miRNA.  See below for details.  Default = 0.8.  Allowed values between 0.5 and 1.

--mindicerfrac [float] : Minimum fraction of mappings within Dicer size range to annotate a locus as Dicer-derived.  Default = 0.85.  Allowed values between 0 and 1.

--phasesize [integer] : Examine phasing only for clusters dominated by the indicated size range.  Size must be within the bounds described by --dicermin and --dicermax.  Set to 'all' to examine p-values of each locus within the Dicer range, in its dominant size.  Set to 'none' to suppress all phasing analysis.  Default = 21.  Allowed values between --dicermin and --dicermax.

--phaseFDR [float] : False Discovery Rate for phased cluster analysis.  FDR of p-values set to this value for Benjamini-Hochberg method.  See below for details.  Default = 0.05.  Allowed values greater than 0 up to 1.

--count [string] : Invokes count mode, in which user-provided clusters are annotated and quantified instead of being defined de novo.  When invoked, the file provided with --count is assumed to contain a simple list of clusters.  Formatting details below.  Default : Not invoked.

--nohp : If "--nohp" appears on the command line, it invokes running in "no hairpin" mode.  RNA folding, hairpin annotation, and MIRNA annotation will be skipped (likely saving significant time).  Default: Not invoked.

--raw : If "--raw" appears on the command line, it prevents conversion of abundances into mappings per million mapped reads, and instead all tallies in the results will simply be the raw reads.  --raw mode is forced if the user does not provide the number of reads via the --reads option.  Default: Not invoked, unless --reads is left blank.

=head1 KEY FORMATTING REQUIREMENTS AND ASSUMPTIONS

=head2 Input .bam file

The mapped reads in the input .bam file must be sorted by chromosomal location, and indexed using the samtools index command -- specifically, ShortStack will look for the [prefix].bam.bai file in the same directory as the input [prefix].bam file.  

Additionally, each mapping in the .bam file must have the NH:i: tag, which indicates the total number of mappings for that read.  

Finally, each mapped read must have the CIGAR string set (column 6 in the SAM specification) -- ShortStack.pl determines the small RNA lengths by parsing the CIGAR string .. if any mappings (except unmapped reads, which are ignored) have "*" entered instead of a valid CIGAR string ShortStack.pl will exit and complain.

Preparation of proper .bam files can be achieved with the helper script, "Prep_bam.pl", included as part of the ShortStack package.  Prep_bam.pl takes in a READ-sorted .bam, .sam, or .sam.gz file and calculates the NH:i: tags, and finally outputs a chromosomal-sorted and indexed .bam file suitable for use by ShortStack.pl.  It will also warn you if any CIGAR strings are missing.  It can also be used directly to receive and process a SAM-formated alignment stream output by bowtie.

Finally, it is critical that the chromosome names referenced in the .bam file correspond exactly to those present in the genome.fasta file.

=head2 Input genome.fasta file

It is critical that this be the precise genome to which the reads in the input .bam file were mapped.

Additionally, the chromosome names in the FASTA headers must be kept SIMPLE.  Specifically, ShortStack.pl at several points parses clusters by the regex /^(\S+):(\d+)-(\d+)$/ or some variant thereof, where the first pattern is the chromosome name.  Therefore, the chromosome names must match (\S+) .. e.g. a single string of one or more non-white-space characters, with no metacharacters.  So, ">Chr1" in your reference genome is good, but ">Chr1 | XM00023 | this is a bunch of annotation blah blah blah" is bad.  This same concern applies to the input .bam file, so your chromosome names should be shortened BEFORE mapping your reads, so that they are short and they are exactly reflected in the .bam file.

If not already present, a .fai index file for the genome will be created using samtools faidx at the beginning of the run.

=head2 --count file

If running in --count mode, the user-provided file is expected to be a simple text file containing a list of coordinates in the format : [Chr]:[start]-[stop], where Chr is defined in the genome file AND in the .bam file, and start and stop are one-based, inclusive.  The same requirement for short, non-whitespaced chromosome names as discussed above holds true for input --count files.  Comment lines, that begin with '#', are ignored.  Tab-delimited files are also accepted, provided the first column has the coordinates.  The second column in tab-delimted files is assumed to be the names of the clusters, and will be used accordingly.  Any other columns in a tab-delimited input file are ignored.

Importantly, the 'Results.txt' file produced by a previous ShortStack.pl run can be used directly in subsequent runs in --count mode.  This is useful when comparing identical intervals across multiple samples.

=head2 --inv_file 

Unless you are running in "nohp" mode, providing a inv_file will enhance the accuracy of the hairpin annotations.  RNALfold-based folding of clusters will often miss very large inverted repeats that einverted can capture.  To make the .inv file, download and install the EMBOSS package ( http://emboss.sourceforge.net/ ), then run einverted against your genome of interest.  For the purposes of ShortStack analysis, the fasta file can be ignored / deleted.  The .inv file is used for ShortStack analysis.  Note that the inv_file is not required, but ShortStack will warn you if it is missing (unless you are running in 'nohp' mode).

=head1 SUGGESTIONS FOR ADAPTER-TRIMMING AND ALIGNMENTS

The results from ShortStack.pl are strongly affected by how the reads were processed and aligned.  Alignment parameters in particular need to be carefully documented, especially with regard to how multi-mapped small RNA reads are treated.  For instance, if your alignment protocol demanded a unique match to the genome, and suppressed results for multi-mapped reads, that will strongly influence cluster discovery as well as render the comparisons of 'total' and 'rep-total' and 'unique-mappers' in the results meaningless.  There (probably?) is no single 'best' method (although I have some suggestions below), but it is clear that, for comparison of results from different samples, all pre-ShortStack processing steps should be explicit and identical between the samples.

Trimming of adapters, while seemingly mundane and simple, also will have profound effects on ShortStack results.  For instance, if reads were computationally filtered before alignment to only retain those in the Dicer size range, than all clusters will be annotated by ShortStack.pl as being dominated by Dicer-sized small RNAs.  However, because the data were selectively used, this may lead to false results.

=head2 Suggestions for Adapter Trimming

In general, I suggest retaining the broadest possible size range of adapter-trimmed reads.  In our group, we are typically retaining 15-35nt reads for mapping at present.  Inputting all reads for alignment will allow ShortStack to confidently discern clusters dominated by reads in the Dicer size range, from clusters that are not (which will often be degradation fragments from abundant RNA species).  Of course, the range of small RNA sizes will also be dictated by the library construction method used, but my suggestion is that a broad size range is preferable in order to allow confident discrimination of Dicer-derived clusters from non-Dicer-derived clusters.

=head2 Suggestions for Pre-Filtering

In short, I suggest avoiding any pre-filtering of the small RNAs prior to alignment.  Pre-filtering is often used to remove reads mapped to rRNAs, tRNAs, and other abundant RNA species that frequently generate a lot of small RNAs that are not thought to be Dicer-derived.  Because ShortStack discriminates clusters based on the sizes, clusters formed by non-Dicer processes will be readily apparent.  In addition, there is evidence that some tRNA and snoRNA-derived fragments might be biologically meaningful, instead of just random degradation ... so throwing out those data altogether does not seem like the best idea.

=head2 Suggestions for Multi-mapped Reads

In many species, particularly plants, a great deal of small RNAs correspond to repetitive genomic sequences, so in general it seems imprudent to discard multi-mapped reads.  Because ShortStack.pl will report the total alignments within a cluster, the total alignments from uniquely-mapped small RNAs, and the repeat-normalized total alignments, the 'repetitiveness' of each cluster should be readily apparent.  So in general I advocate for alignment parameters that retain multi-mapped reads in all possible positions (e.g. all mappings for the given read).  However, there is a practical limit where storing huge numbers of alignments for a single read becomes prohibitive in terms of file size.  We generally cap the number of allowed mappings for any one read to 50 .. e.g., report only the 1st 50 alignments for a single read.  (In bowtie 1, this would be -k 50).

=head1 OUTPUT

=head2 Results.txt

This is a simple tab-delimited text file.  The first line begins with a "#" (comment) sign, and then lists column headers.  Each subsequent line describes the key traits of a single cluster.

Column 1: Locus : The genome-browser-friendly coordinates of the clusters.  Coordinates are one-based, inclusive (e.g. Chr1:1-100 refers to a 100 nt interval beginning with nt 1 and ending with nt 100).

Column 2: Name : Name of cluster.  Unless the run was in --count mode and the input file of a priori clusters already had names, the names are arbitrarily designated as "Cluster_1", "Cluster_2", etc.

Column 3: HP : Whether this cluster appears to be hairpin-derived or not.  If not, a "." is present.  If it is a hairpin, but NOT qualified as a MIRNA, "HP" is indicated.  MIRNAs are indicated by "MIRNA".  If the run was in "--nohp" mode, than all entries in the column will be "ND" (meaning 'not determined').

Column 4: Strand : The pre-dominant strand from which the small RNA emanate.  If ".", no strand was called.  HPs and MIRNAs always have a polarity, based on the hairpin's originating strand.  Non-HP clusters have their polarity determined by the --minstrandfrac setting.

Column 5: Frac_Wat : Fraction of mappings to the Watson (e.g. +) strand of the cluster.  1 means all were from Watson Strand (e.g. +), 0 means all were from Crick (e.g. -) strand.

Column 6: Total : Total mappings within the cluster, either in raw mappings (for --raw mode) or in mappings per million mapped.

Column 7: Uniquely Mapped Total : Total mappings derived from uniquely mapped reads .. e.g., those with NH:i:1.  In raw mappings (for --raw mode) or in mappings per million mapped.

Column 8: Rep-Total : Repeat normalized total mappings.  Instead of each mapping counting as "1", each mapping instead counts as "1/NH", where NH is the total number of mappings that read had, according to the NH:i: tag.

Column 9: DicerCall : If "N", the cluster was not annotated as dicer-derived, per options --dicermin, --dicermax, and --mindicerfrac.  Otherwise this is a number, within the --dicermin to --dicermax size range, which indicates the most abundant small RNA size within the mappings at that cluster.

Column 10: PhaseOffset : If "ND", phasing p-value was not calculated for this cluster.  Otherwise, the offset is the one-based genomic position with which the cluster appears to be "in-phase" (based on the 5' nt of a sense-mapped small RNA).  Phasing is always in increments identicial to the Dicer size call in column 9.

Column 11: Phase_pval :  If "ND", phasing p-value was not calculated for this cluster.  Otherwise, the p-value is derived from a modified hypergeometric distribution, as described below. 

Column 12: Phase_FDR_Call : If "ND", phasing p-value was not calculated for this cluster.  Otherwise, the FDR call will begin with either "OK" or "NS" (meaning 'not significant') and then list the alpha used in the Benjamini-Hochberg procedure.  

Column 13: Short : The total mappings from reads with lengths less than --dicermin, either in raw reads (--raw mode), or mappings per million mapped.

Column 14: Long : The total mappings from reads with lengths more than --dicermax, either in raw reads (--raw mode), or mappings per million mapped.

Columns 15 - the end : The total mappings from reads with the indicated lengths.  These are the sizes within the Dicer range.

=head2 Log.txt

This is a simple log file which records the key settings and key results from the run.

=head2 ShortStack.bed

This is a .bed file for viewing the clusters on a genome browser.  It follows the .bed specification given at the UCSC broswer site <http://genome.ucsc.edu/FAQ/FAQformat.html>.  Clusters are color-coded based on the dominant size.  Non-Dicer clusters are always dark gray.  Dicer-clusters are RGB-rainbow colored from red (shortest) to green (middle) to blue (longest).  Note that the bed coordinate system is zero-based, and the 'stop' coordinate is the first nt NOT in the interval.  So, a 100 nt interval beginning at base 1 and ending at base 100 would have a start of 0 and a stop of 100 in the bed file.

Hairpins and MIRNAs are graphically indicated: The helical arms will be shown as thick boxes, and the rest of the cluster will be thin lines.

=head2 miRNA_summary.txt

This is a tab-delimited text file that summarizes key features of the loci annotated as MIRNAs, including mature miRNA sequences, miRNA-star sequences, and the numbers of mappings for each and for the entire locus.

=head2 hpRNA_summary.txt

A teb-delimited text file summarizing the non-miRNA hpRNA loci.  Indicates "precision", "duplex", and "star".  A "PASS" for precision indicates that one or more mapped small RNAs accounted for >= 20% of the abundance at the locus.  A "PASS" for "duplex" means that the possible miRNA/miRNA* duplex structure was aberrant in some way, thus disqualifying the locus for further consideration as a miRNA locus.  A "PASS" "star" indicates that the miRNA-star was actually sequenced and the sum of the miRNA/miRNA* abudnance accounted for at least 25% of the locus total. (which would make the locus a MIRNA).  These three parameters are analyzed sequentially, so if an earlier one is "FAIL", the subsequent parameters are simply not determined, indicated by "ND".

=head2 Hairpin and MIRNA detail files

Unless the run was done in --nohp mode, each annotated hairpin-derived and MIRNA locus will have its own simple text file to display the details of the locus.  These text files all show A) the Name and genomic coordinates of the locus, B) the sequence, in RNA form, C) the identified hairpin structure, in dot-bracket notation, and D) all mappings whose start and stop is within the interval being examined.

Reads mapped to the sense strand (sense relative to the hairpin, not necessarily relative to the genome) have "."s as placeholders and are shown in the 5'-->3' orientation.  Reads mapped to the antisense strand (antisense relative to the hairpin, not necessarily relative to the genome) have "<"s as placeholders, and are written in the 3' --> 5' orientation.  Annotated mature miRNAs have "m"s as placeholders, and annotated miRNA*'s have "*"s as placeholders.

After each read, the read length (l) and the number of mappings (m) is shown.  Unless the program was run in --raw mode, the normalized mappings per million mapped reads (mmmr) is also shown.

Note that alignments for very complex HP loci are suppressed (only the sequence and structure will be displayed in the detail file, along with a note indicating the alignment was suppressed).  This has a major impact on reducing the memory footprint of ShortStack.  In addition, since the alignments are meant for visual inspection, very complex alignments aren't really parseable by eye anyway.  Complex loci are defined as those where the hairpin is longer than 400nts AND/OR has more than 400 distinct small RNA sequences.  All MIRNA loci have the alignments presented, regardless of whether it is complex or not.

=head1 KEY METHODS

=head2 de novo Cluster Discovery

Cluster discovery proceeds in two simple steps:

1. The total depth of small RNA coverage at each occupied nucleotide in the genome is examined, and initial 'islands' of coverage are defined as continuous stretches where the read depth is greater than or equal to the threshold depth specified by option --mindepth.  Note that islands could theoretically be as small as one nucleotide, since they depend on total depth of coverage.  Many islands will often be 20-24nts in length, corresponding to a pile of a single small RNA species.

2. The initial islands are then extended on both sides by the distance specified by option --pad.  Islands that overlap after extension are merged.  After all extensions and resultant mergers are performed, the final result is the initial clusters.  If the run is performed in --nohp mode, these are the final clusters.  If hairpins and MIRNAs are being examined, some of the clusters may be adjusted in position and/or split to fully capture the putative hairpin(s) (see below).

=head2 Hairpin and MIRNA analysis in de novo mode

1.  The genomic window to be subject to RNA folding is first determined.  For loci whose length is greater than or equal to 1,000nts, the window of genomic DNA corresponding to the locus itself is used.  For loci less than 1,000nts in length, a window, centered upon the middle of the locus, with a size of the lesser of 1,000nts OR 3 x unpadded_cluster_length is used. (The unpadded cluster length is the length - (2 * option --pad)). If 3 x unpadded_cluster_length is less than 250nts, a window size of 250nts is applied.

2.  Both the top and bottom genomic strands from the window are then subjected to secondary structure prediction using RNALfold (options -d 2 -noLP (--maxhpsep)), which returns a diverse set of often overlapping predicted structures.

3.  The structures are parsed, retaining only those that satisfy options --minfracpaired and --minntspaired.  minfracpaired refers to the fraction of nts in the 'lowest' helix that are paired, NOT to the fraction of ALL nts in the window that are paired.  Same for --minntspaired .. only the positions within the lowest helix are considered.

4.  If an .inv file was provided, all inverted repeats in that file are parsed, and then filtered to also satisfy options -- --minfracpaired and --minntspaired.  In addition, loop lengths are not allowed to be longer than 50% of the helix length of the putative hairpin.  Putative RNA secondary structures in dot-bracket notation are generated from the .inv alignment, not by actual RNA folding.  All G-U alignments are considered paired, in addition to the standard A-U and G-C pairings.  Both strands are used, subject to passing the --minfracpaired and --minntspaired filters.  Inverted-repeats that survive these filters are then filtered to retain only those with overlap to the original clusters, and the resulting set of eniverted-derived hairpins is merged with the RNALfold-derived set.

5.  Redundant hairpins are then removed.  Redundant hairpins are those whose 5' arms and 3' arms overlap.  In pairwise comparisons of redundant hairpins, the longest hairpin is retained.

6.  Hairpins that don't have overlap with the original cluster are then removed.  Because the folding window is often extended substantially around the cluster, there could be many putative hairpins that are not within the original cluster.  To have overlap, at least one of the hairpin's helical arms must have at least 20nts within the original cluster coordinates.

7. The pattern of small RNA expression relative to the remaining hairpins is then examined.  To examine the coverage pattern, the hairpin window (not the original cluster window) is first, temporarily, expanded to include an upstream flanking region equal in length to the 5' arm, and a downstream flanking region equal in length to the 3' arm.  The per-nucleotide depth of coverage is then determined across the extended region.  The sum of coverage on the sense strand (relative to the hairpin direction) within the two arms is divided by the total sum of coverage on both strands of the extended region.  This ratio must be >= to the fraction specififed in option --minfrachpdepth in order to pass this step.  In addition there must be at least some coverage on both arms of the hairpin to continue considering the hairpin.

8.  The per-nucleotide depth of coverage across the entire, original cluster is calculated.  The per-nucleotide coverage within the two arms of all of the child hairpins is also calculated, and summed for the entire suite of remaning candidate hairpins.  The coverage in the hairpin arms has to be at least minfrachpdepth in order to continue considering any of the hairpins.  This prevents inappropriate "splitting" of large initial clusters into a hairpin or hairpins that only account for a minor percentage of the original cluster.

9.  Each potential hairpin that remains is next analyzed to see if it qualifies as a MIRNA.  MIRNA locus annotation is designed to satisfy the criteria for de novo annotation of plant MIRNAs as described in Meyers et al. (2008) Plant Cell 20:3186-3190. PMID: 19074682.  In fact, ShortStack's criteria is a little stricter than Meyers et al., in that ShortStack has an absolute requirement for sequencing of the exact predicted miRNA* sequence for a candidate mature miRNA.  It is important to note that ShortStack's MIRNA annotation method is designed to reduce false positives at the expense of an increased rate of false negatives.  In other words, there are likely many bona fide MIRNA loci that end up being classified as Hairpins, instead of MIRNAs, because they don't quite meet the strict criteria set forth below.  

- Precision: There must be at least one candidate mature miRNA that comprises at least 20% of the total abundance of small RNAs mapped to the hairpin.  

- Duplex:  Candidate mature miRNAs must contain no more than 4 unpaired nts (excepting the 2nts on the 3' end), and they must not span a loop (i.e., no base-pairs to themselves).  Additionally, the predicted miRNA*s of candidate mature miRNAs must contain no more than 4 unpaired nts (excepting the 2nts on the 3' end), and they must not span a loop (i.e., no base-pairs to themselves).  Predicted miRNA*'s are based on identifying the small RNA that would form a miR/miR* duplex with a canonical 2nt, 3' overhang.

- Star: The exact predicted miRNA*s of candidate miRNAs must have at least one mapped read, and the total abundance of any candidate mature miRNA/miRNA* pair must be at least 25% of the total small RNA abundance at the locus.  Finally, redundant candidate mature miRNAs are removed, as the steps above initially might classify a small RNA as both a miRNA* and mature miRNA.  In such cases, the partner with the higher abundance is called the miRNA, the other the miRNA*.

10.  The strand bias of small RNA accumulation at each potential hairpin / MIRNA is then evaluated.  At least -- minstrandfrac of the abundance must be to the sense strand relative to the hairpin / MIRNA or it will be eliminated from consideration as a hairpin / MIRNA locus.  However, "dyads", formed by reads that map to both helix arms on opposite strands are ignored during this analysis of strandedness.  Perfectly base-paired hairpins will cause this mapping pattern, so this correction allows for the recovery of hairpins / MIRNAs where the helical region is extensively / completely paired.

11.  Hairpin candidates that failed MIRNA analysis are then subject to increased scrutiny by increasing --minfrachpdepth to "half the distance" from the original --minfrachpdepth value to 1.  For instance, under the default --minfrachpdepth setting of 0.67, a non-MIRNA hairpin must be at least (0.5 * (1 - 0.67)) + 0.67 = 0.835.  Recall from step 7 that minfrachpdepth refers to the fraction of coverage residing on the helical region of the putative hairpin relative to the surrounding area and the loop.  Non-MIRNA haiprins that fail this increased stringency are rejected from further consideration as hairpins.

12. Finally, step 8 is repeated again with the remaining set of potential MIRNAs / hairpins, again to prevent inappropriate splitting of large initial clusters into a hairpin(s)/MIRNA(s).  Survivors are then annotated as MIRNAs or HPs, as appropriate.  

=head2 Hairpin analysis in --count mode

In --count mode, hairpin analysis differs in that A) it does NOT fold an extended region around the unpadded input cluster and B) in case more than one valid hairpin is returned, no HP or MIRNA annotation will be reported, and C) if just a single valid HP or MIRNA is found, its length must be at least 90% of that of the full, originally input coordinates.  The upshot of this is that essentially the input cluster coordiantes have to be a single HP or MIRNA to get annotated as such in count mode -- there is no 'splitting' of larger clusters, even if there might be some valid ones embedded in there.

=head2 Quantification of clusters

All mappings with either their left-ends, right-ends, or both within the cluster are tallied as being within the cluster.  Thus, for a cluster located at Chr1:1000-2000, reads mapped to 980-1000, 1100-1123, and 2000-2021 are all counted as being within the cluster during quantification.  Note that it's possible to count the same mapping within non-overlapping clusters.

=head2 Analysis of Phasing

'Phasing' describes the periodic mapping of small RNAs to repeating intervals equal to their size.  It occurs when helical RNA is Diced processively from a defined terminus; often the terminus is defined by a prior small RNA slicing event followed by RDRP activity, although some MIRNA hairpins are also phased.  Nearly all documented examples of phased small RNA production (in plants) occur for 21nt small RNAs in 21nt increments, hence the default settings of ShortStack to examine only 21-dominated clusters.  This can be changed with option --phasesize.

ShortStack's basic method to identify phased small RNAs involves calculation of a p-value based on the hypergeometric distribution -- this approach was inspired by Chen et al. (2007) PNAS 104: 3318-3323 PMID: 17360645.  However, ShortStack's method modifies the Chen et al. approach to make it more robust at detecting phasing in highly expressed clusters with a background of non-phased noise; the method also allows phasing analysis in any register within the dicer size range (controlled by option --phasesize), and analyzes regions of arbitrary length.  Finally, ShortStack's analysis of phasing is "fuzzy" -- that it, exactly phased reads, and those +1 and -1 phase are all counted as "phased".

Phasing analysis proceeds as follows:

1. Clusters to be analyzed must be annotated as Dicer-derived and be dominated by the size class indicated by option --phasesize.  If --phasesize is set to 'all', all clusters within the Dicer size range will be analyzed.  Conversely, phasing analysis is suppressed for all clusters if option --phasesize is set to 'none'.

2. Cluster must also have a length of more than 4 x the phase size in question .. so, more than 84nts under the default --phasesize 21 setting.  Clusters that are too short are never examined.

3. Phasing is only analyzed with respect to the dominant size of the cluster.  So, for a cluster dominated by 21mers, only phasing in 21nt increments will be examined.

4. The 5' positions of all sense-mapped small RNAs are tallied as a function of genomic position.  The 3' positions of all antisense-mapped small RNAs are also tallied, after adding 2nts to account for the 2nt, 3' overhangs left by Dicer processing.  After this process, each genomic position within the cluster has a number reflecting the number of small RNA termini at that position.  If the cluster is longer than 20 times the phase (e.g. 20 x 21 for the default settings), reads mapped beyond the 20 x 21 mark are allocated to the beginning of the cluster, keeping it in phase.  For instance, assuming --phasesize of 21, reads in position 420 are assigned at 420, those at 421 get flipped back to 1, 422 back to 2, and so on.  This is necessary because p-value calculation involved calculation of binomial coefficents, which grow too large to calculate (easily) with inputs of more than 500 or so.

5. The average abundance of termini across the locus is calculated from the above representation of the reads.

6. The total abundance in each of the possible phasing registers (there are 21 registers in the default mode of --phasesize 21) is calculated.  The register with the maximum total abundance is the used in p-value determination.  The offset of this register is also noted; the offset is the 1st genomic position representing the 5'-sense position of a phased small RNA.

7.  The p-value within the chosen register is then calculated using the cumulative distribution function (CDF) for the hypergeometric distribution.  Sorry, hard to show equations in plain-text -- see Wikipedia's Hypergeometric distribution entry, under CDF. N (the population size) is the number of nt positions in the locus. m (the number of success states in the population) is the number of possible positions in the phasing register of interest, INLCUDING POSITIONS +1 AND -1 RELATIVE TO THE REGISTER OF INTEREST.  This means phasing is "fuzzy", which is often seen in the known examples of this phenomenon.  n (the number of draws) is defined as the total number of positions with ABOVE AVERAGE abundance.  k (the number of successes) is the number of phased positions (inlduing the fuzzy +1 and -1 positions) with ABOVE AVERAGE abundance.  The p-value is then calculated per the hypergeometric distribution CDF.  NOTE: The restriction of n and k to only above-average abundance works well to eliminate low-level noise and focus on the dominant small RNA pattern within the locus.

8. After all p-values have been caluclated, ShortStack uses the Benjamini-Hochberg method to control for the false discovery rate at a user-specified alpha (default is 0.05).  Significant clusters are noted 'OK' and non-significant clusters are noted 'NS'; in both cases the p-values are reported.   Results of phasing analysis are in columns 10-12 of the Results.txt file.



=cut


