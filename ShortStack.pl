#!/usr/bin/perl -w
# See below the __END__ mark for docmentation, license, citation, etc., or just see the README distributed with this script

use Getopt::Long;
use strict;

###############MAIN PROGRAM BLOCK
##### VERSION
my $version_num = "1.2.0"; 


##### get options and validate them

# usage statement and option gathering
my $usage = usage($version_num);

# If user has no arguments, give help statement and quit
unless($ARGV[0]) {
    die "\n$usage\n";
}

# Initial option definition, inlcuding default settings
my $outdir = '';
my $mindepth = 20;
my $pad = 100;
my $dicermin = 20;
my $dicermax = 24;
my $maxhpsep = 300;
my $minfracpaired = 0.67;
my $minntspaired = 15; ## lowered to 15 from 30 as of version 0.3.0
my $minfrachpdepth = 0.67;  ## raised from 0.5 as of version 0.2.1
my $minstrandfrac = 0.8;
my $mindicerfrac = 0.8;  ## lowered to 0.8 from 0.85 as of 0.4.0
my $minUI = 0.1;  ## parameter added as of version 0.4.0 .. "minimum Uniqueness Index required to attempt RNA folding / hairpin analysis"
my $count = '';
my $nohp = 0;
my $phasesize = 21;
my $inv_file = '';
my $flag_file = "NULL";
my $maxmiRHPPairs = '';  ## parameter added as of version 0.3.0
my $maxmiRUnpaired = '';  ## parameter added as of version 0.3.0
my $miRType = "plant";  ## parameter added as of version 0.3.0
my $maxdGperStem = -0.5;  ## parameter added as of version 0.3.0
my $maxLoopLength = '';  ## parameter added as of version 0.3.0

## New in 1.0.0
my $untrimmedFA;
my $untrimmedFQ;
my $adapter;
my $trimmedFA;
my $trimmedFQ;
my $help;
my $version;
my $bamfile;
my $align_only;

## New in 1.1.0
my $read_group;

## New in 1.2.0 .. support of colorspace data
my $untrimmedCS;
my $trimmedCS;
my $untrimmedCSQV;
my $trimmedCSQV;

# get user options from command line
GetOptions ('outdir=s' => \$outdir,
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
	    'minUI=f' => \$minUI,
	    'count=s' => \$count,
	    'nohp' => \$nohp,
	    'phasesize=s' => \$phasesize,
	    'inv_file=s' => \$inv_file,
	    'flag_file=s' => \$flag_file,
	    'maxmiRHPPairs=i' => \$maxmiRHPPairs,
	    'maxmiRUnpaired=i' => \$maxmiRUnpaired,
	    'miRType=s' => \$miRType,
	    'maxdGperStem=f' => \$maxdGperStem,
	    'maxLoopLength=i' => \$maxLoopLength,
	    'untrimmedFA=s' => \$untrimmedFA,
	    'untrimmedFQ=s' => \$untrimmedFQ,
	    'untrimmedCS=s' => \$untrimmedCS,
	    'untrimmedCSQV=s' => \$untrimmedCSQV,
	    'adapter=s' => \$adapter,
	    'trimmedFA=s' => \$trimmedFA,
	    'trimmedFQ=s' => \$trimmedFQ,
	    'trimmedCS=s' => \$trimmedCS,
	    'trimmedCSQV=s' => \$trimmedCSQV,
	    'bamfile=s' => \$bamfile,
	    'align_only' => \$align_only,
	    'read_group=s' => \$read_group,
	    'version' => \$version,
	    'help' => \$help);

# If help, print and quit
if($help) {
    my $long_help = full_help($version_num);
    die "\n$long_help\n";
}

# If version, print version and quit
if($version) {
    die "ShortStack.pl version $version_num\n";
}

# default output directory is ShortStack_[time], where time is the number of non-leap seconds since 00:00:00 UCT Jan 1, 1970
unless($outdir) {
    my $time = time;
    $outdir = "ShortStack" . "_$time";
}
# ensure the output directory does not already exist.  If it does, quit and complain
# then create it or die tryin'
if(-d $outdir) {
    die "FATAL: Output directory $outdir already exists\n$usage\n";
} else {
    system "mkdir $outdir";
}

#designate a log file
our $logfile = "$outdir" . "\/" . "Log\.txt";

# log the outdir
log_it($logfile, "ShortStack log file at $logfile\n\nBeginning initial checks");

       
# Determine the mode
my $mode;
# Mode 1: Trim, align, then analyze.
# Expect either untrimmedFA or untrimmed FQ, but not both:
if(($untrimmedFA) or ($untrimmedFQ) or ($untrimmedCS)) {
    # Should be mode 1. Ensure adapter specified and valid .. must be 8nts, ATCG
    unless($adapter) {
	log_it($logfile,"FATAL: Option --adapter must be specified in conjunction with untrimmed reads\n\n$usage\n");
	exit;
    }
    ## Adapters can now be comma-delimited
    my $adapter_check = check_adapter($adapter);
    unless($adapter_check) {
	log_it($logfile,"FATAL: Option --adapter must be a string of at least 8 ATGC characters .. case-insensitive or comma delimited set of such\n\n$usage\n");
	exit;
    }
    # Cant have more than one
    if((($untrimmedFA) and ($untrimmedFQ)) or
       (($untrimmedFA) and ($untrimmedCS)) or
       (($untrimmedFQ) and ($untrimmedCS))) {
	log_it($logfile,"FATAL: Must specifiy either --untrimmedFA OR --untrimmedFQ OR --untrimmedCS, but not more than one type\n\n$usage\n");
	exit;
    }
    $mode = 1;
    
    # If both untrimmedCS and untrimmedCSQV are specified, they must have the same number of files .. and they all must be readable
    if(($untrimmedCS) and ($untrimmedCSQV)) {
	my $csqv_check = check_csqv($untrimmedCS,$untrimmedCSQV);
	unless($csqv_check) {
	    exit;
	}
    }
    
    
    # Warn about irrelevant options
    if($trimmedFA) {
	log_it($logfile,"\nWARNING: in mode 1 option --trimmedFA is irrelvant. It is being ignored\n");
    }
    if($trimmedFQ) {
	log_it($logfile,"\nWARNING: in mode 1 option --trimmedFQ is irrelvant. It is being ignored\n");
    }
    if($trimmedCS) {
	log_it($logfile,"\nWARNING: in mode 1 option --trimmedCS is irrelvant. It is being ignored\n");
    }
    if($bamfile) {
	log_it($logfile,"\nWARNING: in mode 1 option --bamfile is irrelvant. It is being ignored\n");
    }
} elsif (($trimmedFA) or ($trimmedFQ) or ($trimmedCS)) {
    # Should be mode 2. 
    # Cant have both
    if((($trimmedFA) and ($trimmedFQ)) or
       (($trimmedFA) and ($trimmedCS)) or
       (($trimmedFQ) and ($trimmedCS))) {
	log_it($logfile,"FATAL: Must specify either --trimmedFA OR --trimmedFQ OR --trimmedCS for mode 2, but not more than one type\n\n$usage\n");
	exit;
    }
    $mode = 2;

    # Warn about irrelevant options --adapter and --bamfile
    if($adapter) {
	log_it($logfile,"\nWARNING: in mode 2 option --adapter is irrelevant. It is being ignored\n");
    }
    if($bamfile) {
	log_it($logfile,"\nWARNING: in mode 2 option --bamfile is ireelevent. It is being ignored\n");
    }
} elsif ($bamfile) {
    # mode 3. 
    $mode = 3;

    # Warn about irrelevant option --adapter
    if($adapter) {
	log_it($logfile,"\nWARNING: in mode 2 option --adapter is irrelevant. It is being ignored\n");
    }
} else {
    log_it($logfile,"\nFATAL: Could not determine mode, aborting\n");
    exit;
}

## Check dependencies

# check for required installation of samtools and get version, or quit and complain
my ($samtools_version,$full_samtools_test_text) = get_samtools_version();
if($samtools_version eq "not found") {
    log_it($logfile, "\nFATAL: samtools not found\n\n");
    exit;
}

# check for required installation of RNALfold, get version
my $RNALfold_check = check_RNALfold();
unless($RNALfold_check) {
    log_it($logfile, "\nFATAL: RNALfold not found\n\n");
    exit;
}

# check for required installation of RNAeval, get version
my $RNAeval_check = check_RNAeval();
unless($RNAeval_check) {
    log_it($logfile, "\nFATAL: RNAeval not found\n\n");
    exit;
}

# If modes one or two, check for required bowtie version 1.x or 0.12.x and bowtie-build same
my $bowtie_check;
my $bowtie_build_check;
if(($mode == 1) or ($mode == 2)) {
    $bowtie_check = check_bowtie();
    unless(($bowtie_check =~ /^0\.12/) or ($bowtie_check =~ /^1\./)) {
	log_it($logfile, "\nFATAL: bowtie version 1.x or 0.12.x not found\n\n");
	exit;
    }
    $bowtie_build_check = check_bowtie_build();
    unless(($bowtie_build_check =~ /^0\.12/) or ($bowtie_build_check =~ /^1\./)) {
	log_it($logfile, "\nFATAL: bowtie version 1.x or 0.12.x not found\n\n");
	exit;
    }
}

# Validate the options

unless(($mindepth > 0) and
       ($mindepth < 1000000)) {
    log_it($logfile,"FATAL: mindepth must be more than zero and less than one million\n\n$usage\n");
    exit;
}

# ensure that option --pad is present and logical -- 0 >= x < 100,000
unless(($pad >= 0) and ($pad <= 100000)) {
    log_it($logfile,"FATAL: Option --pad must be a number greater than or equal to zero and less than one hundred thousand\n\n$usage\n");
    exit;
}

# ensure that option --dicermin is present and logical -- 15 >= x <= 35, and less than or equal to dicermax
unless(($dicermin >= 15) and ($dicermin <= 35) and ($dicermin <= $dicermax)) {
    log_it($logfile, "FATAL: Option --dicermin must be a number between 15 and 35, and less than or equal to option --dicermax\n\n$usage\n");
    exit;
}

# ensure that option --dicermax is present and logical -- 15 >= x <= 35, and more than or equal to dicermin
unless(($dicermax >= 15) and ($dicermax <= 35) and ($dicermin <= $dicermax)) {
    log_it($logfile,"FATAL: Option --dicermax must be a number between 15 and 35, and more than or equal to option --dicermin\n\n$usage\n");
    exit;
}

# ensure that option --maxhpsep is present and logical -- 50 >= x <= 2000
unless(($maxhpsep >= 50) and ($maxhpsep <= 2000)) {
    log_it($logfile, "FATAL: Option --maxhpsep must be a number between 50 and 2000\n\n$usage\n");
    exit;
}

# ensure that option --minfracpaired is present and logical -- 0 > x <= 1
unless(($minfracpaired > 0) and ($minfracpaired <= 1)) {
    log_it($logfile, "FATAL: Option --minfracpaired must be a number greater than zero and less than or equal to one\n\n$usage\n");
    exit;
}

# ensure that option --minntspaired is present and logical -- 0 > x <= --maxhpsep
unless(($minntspaired > 0) and ($minntspaired <= $maxhpsep)) {
    log_it($logfile,"FATAL: Option --minntspaired must be number greater than zero and less than or equal to option --maxhpsep\n\n$usage\n");
    exit;
}

# ensure that option --minfrachpdepth is present and logical -- 0 >= x <= 1
unless(($minfrachpdepth >= 0) and ($minfrachpdepth <= 1)) {
    log_it($logfile, "FATAL: Option --minfrachpdepth must be a number greater than or equal to zero and less than or equal to one\n\n$usage\n");
    exit;
}

# ensure that option --minstrandfrac is present and logical -- 0.5 >= x <= 1
unless(($minstrandfrac >= 0.5) and ($minstrandfrac <= 1)) {
    log_it($logfile, "FATAL: Option --minstrandfrac must be number greater than or equal to 0.5 and less than or equal to one\n\n$usage\n");
    exit;
}

# ensure that option --mindicerfrac is present and logical -- 0 >= x <= 1
unless(($mindicerfrac >= 0) and ($mindicerfrac <= 1)) {
    log_it($logfile, "FATAL: Option -mindicerfrac must be number greater than or equal to 0 and less than or equal to one\n$usage\n");
    exit;
}

# ensure that option --minUI is present and between 0 and 1
unless(($minUI >= 0) and ($minUI <= 1)) {
    log_it($logfile, "FATAL: Option --minUI must be a number between 0 and 1\n$usage\n");
    exit;
}

# if the --count option was provided, make sure the indicated file can be read
if($count) {
    unless(-r $count) {
	log_it($logfile, "FATAL: Option --count : File $count could not be read\n\n$usage\n");
	exit;
    }
    # As of version 0.3.0, count mode also forces nohp mode
    $nohp = 1;
}

# check the --phasesize option to see if it is a number between --dicermin and --dicermax, OR 'all', OR 'none'
if($phasesize =~ /^\d+$/) {
    unless(($phasesize >= $dicermin) and
	   ($phasesize <= $dicermax)) {
	log_it($logfile, "FATAL: Option --phasesize must be an integer within the --dicermin to --dicermax size range OR \'all\'\n\n$usage\n");
	exit;
    }
} else {
    unless(($phasesize eq "all") or
	   ($phasesize eq "none")) {
	log_it($logfile, "FATAL: Option --phasesize must be an integer within the --dicermin to --dicermax size range OR \'all\' OR \'none\'\n\n$usage\n");
	exit;
    }
}

# parse the mirType and the associated settings
unless(($miRType eq "plant") or ($miRType eq "animal")) {
    log_it($logfile, "FATAL: Option --miRType $miRType is invalid\.  Must be either \'plant\' , \'animal\', or unspecified \(unspecified defaults to \'plant\'\)\n\n$usage\n");
    exit;
}
my $maxmiRHPPairs_override = 0;
my $maxmiRUnpaired_override = 0;
my $maxLoopLength_override = 0;

if($maxmiRHPPairs) {
    log_it($logfile, "\n\tNote: Option --maxmiRHPPairs was explicitly set to $maxmiRHPPairs by user: Overides setting provided by option --miRType of $miRType\n\n");
    $maxmiRHPPairs_override = 1;
} elsif ($miRType eq "plant") {
    $maxmiRHPPairs = 150;
} elsif ($miRType eq "animal") {
    $maxmiRHPPairs = 45;
}

if($maxmiRUnpaired) {
    log_it($logfile, "\n\tNote: Option --maxmiRUnpaired was explicitly set to $maxmiRUnpaired by user: Overrides setting provided by option --miRType of $miRType\n\n");
    $maxmiRUnpaired_override = 1;
} elsif ($miRType eq "plant") {
    $maxmiRUnpaired = 5;
} elsif ($miRType eq "animal") {
    $maxmiRUnpaired = 6;
}

if($maxLoopLength) {
    log_it($logfile,"\n\tNote: Option --maxLoopLength was explicitly set to $maxLoopLength by user: Overrides setting provided by option --miRType of $miRType\n\n");
    $maxLoopLength_override = 1;
} elsif ($miRType eq "plant") {
    $maxLoopLength = 100000;  ## an arbitrarily high number
} elsif ($miRType eq "animal") {
    $maxLoopLength = 15;
}

###############
# check that the genome file is at the end of @ARGV
my $genome = pop @ARGV;
unless(-r $genome) {
    log_it($logfile, "FATAL: genome file $genome is not readable\n\nFATAL\n$usage\n");
    exit;
}

#########
# Check for the existence of the einverted .inv file
# If not provided, issue a warning and proceed
# If a string is provided but no file is readble, quit and complain
if($inv_file) {
    unless(-r $inv_file) {
	log_it($logfile,"FATAL: The provided \.inv file $inv_file from user option --inv_file could not be read\n$usage\n");
	exit;
    }
    # if the run is in 'nohp' mode, than this file is irrelevant, so warn the user about that
    if($nohp) {
	log_it($logfile,"\nWARNING: An einverted \.inv file was provided \($inv_file\) but this run is in nohp mode, so the einverted file will be ignored\n\n");
    }
} else {
    unless($nohp) {
	log_it($logfile,"\nWARNING: No einverted \.inv file was provided \(option --inv_file\)\.  This may decrease the annotation accuracy, especially for loci deriving from large inverted repeats\n\n");
    }
}

#### 
# See if the user provided a flag_file, and if so, make sure it is readable
if($flag_file ne "NULL") {
    unless(-r $flag_file) {
	log_it($logfile, "FATAL: The provided flag file $flag_file from user option --flag_file could not be read\n$usage\n");
	exit;
    }
}
##

##### Report to user on initialization of the run
log_it ($logfile,"\n$0 $version_num\n");
log_it ($logfile,`date`);
log_it ($logfile,"Genome: $genome\n");
log_it ($logfile,"Output Directory: $outdir\n");
log_it ($logfile,"samtools $samtools_version");
log_it ($logfile,"RNALfold $RNALfold_check\n");
log_it ($logfile,"RNAeval $RNAeval_check\n");
if(($mode == 1) or ($mode == 2)) {
    log_it ($logfile,"bowtie version $bowtie_check\n");
    log_it ($logfile,"bowtie-build version $bowtie_build_check\n");
}
log_it($logfile,"Pre-analysis mode: $mode ");
if($mode == 1) {
    log_it($logfile,": Trim reads, then align trimmed reads, \n\tUntrimmed Reads: ");
    if($untrimmedFA) {
	log_it($logfile,"$untrimmedFA\n");
    } elsif($untrimmedFQ) {
	log_it($logfile,"$untrimmedFQ\n");
    } else {
	log_it($logfile,"$untrimmedCS\n");
	if($untrimmedCSQV) {
	    log_it($logfile,"\tUntrimmed Qual values: $untrimmedCSQV\n");
	}
    }
	
    log_it($logfile, "\tAdapter: $adapter\n");
} elsif ($mode == 2) {
    log_it($logfile,"Align pre-trimmed reads\n\tTrimmed Reads: ");
    if($trimmedFA) {
	log_it($logfile,"$trimmedFA\n");
    } elsif ($trimmedFQ) {
	log_it($logfile,"$trimmedFQ\n");
    } else {
	log_it($logfile,"$trimmedCS\n");
	if($trimmedCSQV) {
	    log_it($logfile,"\tTrimmed Qual values: $trimmedCSQV\n");
	}
    }
} elsif ($mode == 3) {
    log_it($logfile,"Analyze pre-existing BAM alignment\n\tAlignments: $bamfile\n");
    if($align_only) {
	log_it($logfile,"\nWARNING: option --align_only is irrelevant if you provide a bamfile. It is being ignored and analysis is proceeding\n");
    }
}

if($align_only) {
    log_it($logfile, "\nAlign-only mode: Trimmed reads will be aligned but no analysis will be performed.\n");
} else {
    log_it($logfile, "Read Group: ");
    if($read_group) {
	log_it($logfile, "$read_group\n");
    } else {
	log_it($logfile, "None .. all reads will be analyzed.\n");
    }
    unless($nohp) {
	log_it ($logfile,"einverted \.inv file of Inverted Repeats:");
	if(-r $inv_file) {
	    log_it ($logfile," $inv_file\n");
	} else {
	    log_it ($logfile," NOT PROVIDED, hairpins from RNALfold only\n");
	}
    }
    log_it ($logfile,"Flag file of loci to report if overlapped: ");
    if($flag_file ne "NULL") {
	log_it ($logfile,"$flag_file\n");
    } else {
	log_it ($logfile,"Not Provided\n");
    }
    
    log_it ($logfile,"Clusters:");
    if($count) {
	log_it ($logfile," Running in \"count\" mode\. User provided clusters from file $count\n");
    } else {
	log_it ($logfile," To be calculated de novo\n");
	log_it ($logfile,"Island threshold: $mindepth mappings\n");
	log_it ($logfile,"Padding around initial islands: $pad nts\n");
    }
    log_it ($logfile,"Dicer size range: $dicermin to $dicermax\n");
    if($nohp) {
	log_it ($logfile,"Running in \"nohp\" mode: Hairpins and MIRNAs will not be inferred\n");
    } else {
	log_it ($logfile,"Minimum uniqueness index required to attempt RNA folding during hairpin preduction: $minUI\n");
	log_it ($logfile,"Maximum allowed separation of a base pair to span during hairpin prediction: $maxhpsep\n");
	log_it ($logfile,"Minimum fraction of paired nts allowable in a valid hairpin structure: $minfracpaired\n");
	log_it ($logfile,"Minimum number of base pairs within a valid hairpin structure: $minntspaired\n");
	log_it ($logfile,"Maximum deltaG per stem length for a valid hairpin structure: $maxdGperStem\n");
	log_it ($logfile,"Minimum fraction of coverage in hairpin helix to keep hairpin: $minfrachpdepth\n");
	log_it ($logfile,"microRNA type: $miRType\n");
	log_it ($logfile,"\tMaximum number of base pairs in hairpin for MIRNAs: $maxmiRHPPairs");
	if($maxmiRHPPairs_override) {
	    log_it ($logfile," Overrides normal setting for type $miRType\n");
	} else {
	    log_it ($logfile,"\n");
	}
	log_it ($logfile,"\tMaximum number of unpaired miRNA nts in miR/miR\* duplex: $maxmiRUnpaired");
	if($maxmiRUnpaired_override) {
	    log_it ($logfile," Overrides normal setting for type $miRType\n");
	} else {
	    log_it ($logfile,"\n");
	}
	log_it ($logfile,"\tMaximum loop size for miRNA or other hpRNA loci: $maxLoopLength");
	if($maxLoopLength_override) {
	    log_it ($logfile," Overrides normal setting for type $miRType\n");
	} else {
	    log_it ($logfile,"\n");
	}
    }
    log_it ($logfile,"Minimum fraction of mappings to assign a polarity to a non-hairpin cluster: $minstrandfrac\n");
    log_it ($logfile,"Minimum fraction of mappings within the Dicer size range to annotate a locus as Dicer-derived: $mindicerfrac\n");
    log_it ($logfile,"Cluster type to analyze for phasing: $phasesize\n");
}

# Next, check for the presence of .fai index file corresponding to the genome.  If not found, create one with samtools
my $expected_faidx = "$genome" . "\.fai";
unless(-e $expected_faidx) {
    log_it ($logfile,"Expected genome index $expected_faidx for genome file $genome not found\.  Creating it using samtools faidx");
    system "samtools faidx $genome";
log_it ($logfile," done\n\n");
}


##############################################
# pre-analysis
my @trimmed_FAs = ();
my @trimmed_FQs = ();
my @trimmed_CSs = ();


if($mode == 1) {
    if($untrimmedFA) {
	@trimmed_FAs = trim_FA($untrimmedFA,$adapter);
    } elsif ($untrimmedFQ) {
	@trimmed_FQs = trim_FQ($untrimmedFQ,$adapter);
    } elsif ($untrimmedCS) {
	@trimmed_CSs = trim_CS($untrimmedCS,$adapter,$untrimmedCSQV);  ## if there were input qual files, the entries in the array are concatenated pairs of trimmedCS and trimmedQV with ":::"
	
    } else {
	log_it($logfile,"FATAL: mode 1 but neither --untrimmedFA nor --untrimmedFQ nor --untrimmedCS found\?\n\n");
	exit;
    }
} elsif ($mode == 2) {
    if($trimmedFA) {
	@trimmed_FAs = split (",",$trimmedFA);
    } elsif ($trimmedFQ) {
	@trimmed_FQs = split (",",$trimmedFQ);
    } elsif ($trimmedCS) {
	if($trimmedCSQV) {
	    my @initial_CSs = split (",", $trimmedCS);
	    my @initial_QVs = split (",", $trimmedCSQV);
	    for(my $x = 0; $x < (scalar @initial_CSs); ++$x) {
		my $concat = "$initial_CSs[$x]" . ":::" . "$initial_QVs[$x]";
		push(@trimmed_CSs, $concat);
	    }
	} else {
	    @trimmed_CSs = split (",",$trimmedCS);
	}
    }
}

if(($mode == 1) or ($mode == 2)) {
    my $filetype;
    my $to_align_file;
    my @to_align_files = ();
    if(@trimmed_FAs) {
	$filetype = "a";
	@to_align_files = @trimmed_FAs;
    } elsif (@trimmed_FQs) {
	$filetype = "q";
	@to_align_files = @trimmed_FQs;
    } elsif (@trimmed_CSs) {
	$filetype = "c";
	@to_align_files = @trimmed_CSs;
    } else {
	log_it($logfile,"\nFATAL: Failed to ID filetype between fastq or fasta or csfasta for some reason\n\n");
	exit;
    }
    my @bamfiles = ();
    foreach $to_align_file (@to_align_files) {
	$bamfile = perform_alignment($to_align_file,$genome,$filetype);
	unless(-r $bamfile) {
	    log_it($logfile,"\nFATAL: Failed to read bamfile $bamfile after sub-routine perform_alignment\n");
	    exit;
	}
	log_it($logfile,`date`);
	log_it($logfile,"Alignment of $to_align_file completed: See file $bamfile\n");
	push(@bamfiles, $bamfile);
    }
    
    # If more than one, merge em
    if((scalar @bamfiles) > 1) {
	$bamfile = merge_em(\@bamfiles,\$outdir);
    }
    
    if($align_only) {
	log_it($logfile,"\nRun was in --align_only mode .. terminating.\n");
	exit;
    }
}

################################################
# bam file validation
my %rgs = validate_bam($bamfile,$expected_faidx);
unless(%rgs) {
    log_it($logfile,"\nBAM file validation failed. Aborting run.\n");
    exit;
}

################################################
# If user wanted to only analyze a specific read-group, ensure it is in the bamfile
# Clear hash if it had a null place-holder
if(exists($rgs{'NULL'})) {
    %rgs = ();
}
if($read_group) {
    unless(exists($rgs{$read_group})) {
	log_it($logfile, "\nFATAL: Read group $read_group specified with option --read_group is NOT present in the bamfile\!\n");
	exit;
    }
}

###############################################

##### Phase One: Identify Clusters
my %names = ();  ## keys = locusIDs (e.g. Chr1:1-1000) , values = name designations

log_it ($logfile,`date`);
log_it ($logfile,"Phase One: Identifying Clusters");
my @clusters = ();
if($count) {
    log_it ($logfile," a priori from file $count\n");
    @clusters = get_clusters_a_priori($count);
    %names = get_names_countmode($count);  ## if names not present in count file, arbitrary names assigned
} else {
    log_it ($logfile," de novo\n");

    my @islands = get_islands($bamfile,$mindepth,$expected_faidx,$read_group);
    
    @clusters = merge_clusters(\@islands,\$pad,\$genome);
    %names = get_names_simple(\@clusters);
}

my $phase_one_n_clusters;
if(@clusters) {
    $phase_one_n_clusters = scalar @clusters;
    log_it ($logfile,"Phase One complete: Found $phase_one_n_clusters Clusters\n\n");
} else {
    log_it ($logfile,"Sorry, no clusters meeting your criteria were found in these data\.\n");
    exit;
}

### NEW PHASE 2 -- INITIAL QUANTIFICATION
## OLD Phase 5: Annotate and Quantify all clusters
log_it ($logfile,`date`);
log_it ($logfile,"Phase Two: Quantify all clusters\n");

my %quant_master = quant(\@clusters,\$bamfile,\$dicermin,\$dicermax,\$minstrandfrac,\$mindicerfrac,\$phasesize,\%names, \$read_group);

# in this initial hash, the fields are
#[0] : locus
#[1] : name
#[2] : strand  (+ or -)
#[3] : frac_Watson
#[4] : total
#[5] : unique
#[6] : uniqueness_index
#[7] : dicer_call
#[8] : phase_offset
#[9] : phase_pval
#[10] : short
#[11] : long
#[12 -- end] : dicer size mappings


# later, we will modify the notations for hairpins

# new fields to add
#[0] : locus
#[1] : name
#[2] : HP_call  <<<<<
#[3] : strand
#[4] : frac_Watson
#[5] : total
#[6] : unique
#[7] : uniqueness_index
#[8] : dicer_call
#[9] : phase_offset
#[10] : phase_pval
######### >>>>>>>> OBSOLETED !#[11] : phase_FDR_call
#[11] : pairs           <<<<<<
#[12] : frac_paired    <<<<<<
#[13] : stem_length    <<<<<<
#[14] : loop_length    <<<<<<
#[15] : dGperStem      <<<<<<
#[16] : frac_hp        <<<<<<
#[17] : hp_size_result <<<<<<
#[18] : prec_result    <<<<<<
#[19] : duplex_result  <<<<<<
#[20] : star_result    <<<<<<
#[21 -- end] : dicer size mappings

my @final_clusters = ();
my %final_hps = ();
my %miRNAs = ();

## Deal with hairpins now
if($nohp) {
    log_it ($logfile,"\n");
    log_it ($logfile,`date`);
    log_it ($logfile,"Running in nohp mode: Skipping phases 3, 4, and 5\n");
    @final_clusters = @clusters;


} else {
    log_it ($logfile,"\n");
    log_it ($logfile,`date`);
    log_it ($logfile,"Phase Three: Preparing for RNA secondary structure prediction\n");
    log_it ($logfile,"\n\tFiltering clusters based on repetitveness and DicerCall\.\.\.");

    # get only those clusters with uniqueness values greater than or equal to --minUI and a DicerCall that is NOT "N"
    my @folding_clusters = get_folding_clusters(\@clusters,\%quant_master,\$minUI);
    
    # report 
    my $n_to_fold = scalar @folding_clusters;
    
    log_it ($logfile, "Done\n");
    log_it ($logfile, "\t\tHairpin and miRNA analysis limited to $n_to_fold out of the $phase_one_n_clusters\n");
    log_it ($logfile, "\t\tThe others are omitted from structural analysis because of high repetitiveness and\/or a DicerCall of \"N\"\n");
    
    # Begin folding processes
    log_it ($logfile,"\n\tGathering genomic sequences for folding\n");
    my %to_fold = get_folding_regions(\$genome,\@folding_clusters,\$maxhpsep);  ## hash structure: 'sequence', 'start', 'stop'  strand is always Watson
    
    ## Phase Four -- Folding
    log_it ($logfile, "\n\n");
    log_it ($logfile,`date`);
    log_it ($logfile,"Phase Four: Identification of qualified hairpins overlapping clusters\n");

    #fold each query and retain only qualifying structures

    log_it ($logfile,"\tFolding sequences and identifying qualifying hairpins\.\.\.\n");

    my %qualifying_hairpins = folder(\%to_fold,\$maxhpsep,\$minntspaired,\$minfracpaired,\$maxdGperStem,\$maxLoopLength);

    ## Hash structure:  key is locus.  value is an array.  Each array entry is tab-delimited string containing:
    # 0 : brax
    # 1 : adj_start
    # 2 : details
    # 3 : strand (Watson or Crick)
    #### 0 - 3 are as before  .. new ones are below
    # 4 : pairs
    # 5 : frac_paired
    # 6 : stem_length
    # 7 : loop_length
    # 8 : dGperStem
    
    #convert the coordinates of the qualifying hairpins to the true chromosomal coordinates
    log_it ($logfile,"\tConverting hairpin coordinates to true genomic coordinates\.\.\.");
    my %true_hps = hairpin_coords(\%to_fold,\%qualifying_hairpins);  ## structure same as %qualifying hairpins
    ### same structure as %qualifying_hairpins, except field[1] of each entry is now a start-stop
    
    log_it ($logfile,"Done\n");
    
    ## take in the einverted file
    if(-r $inv_file) {
	log_it ($logfile,"\tParsing inverted repeats from file $inv_file\.\.\.");
	my @initial_irs = parse_inv(\$inv_file,\$minntspaired,\$minfracpaired,\$maxdGperStem,\$maxLoopLength);
	## each line is tab-delimited, with the same fields as in %true_hps, except addition of an extra field ([9]), which is the chr
	
	log_it ($logfile," Done\n");
	my $initial_ir_tally = scalar @initial_irs;
	log_it ($logfile,"\t\tFound $initial_ir_tally potential qualifying IRs\n");
	log_it ($logfile,"\tMerging einverted-derived inverted repeats with clusters and with RNAL-fold derived hairpins\.\.\.");
	my ($in_irs,$out_irs) = merge_inv(\@initial_irs,\%true_hps,\@folding_clusters);
	log_it ($logfile," Done\n");
	log_it ($logfile,"\t\tFrom the $in_irs potentially qualifying IRs, $out_irs had overlap with one or more clusters and were merged with RNALfold hairpins\n");
    }
    
    log_it ($logfile,"\tRemoving redundant hairpins on a per-locus basis\.\.\.");
    my %true_hps_trimmed = remove_redundant_hps(\%true_hps); ## structure same as %true_hps
    log_it ($logfile," Done\n");
    
    log_it ($logfile,"\tRemoving hairpins that do not overlap clusters\.\.\.");
    my %true_hps_3 = remove_nonoverlapped_hps(\%true_hps_trimmed);  ## structure same as %qualifying hairpins
    log_it ($logfile," Done\n");

    
    ## Phase Five : Annotation of hpRNA and microRNA loci
    log_it($logfile, "\n");
    log_it($logfile, `date`);
    log_it($logfile, "Phase Five: Annotation of hpRNA and microRNA loci\n");
    
    log_it ($logfile,"\tFiltering hairpins based on expression evidence\.\.\.");
    my %filtered_hps = hp_expression(\%true_hps_3,\$bamfile,\$minfrachpdepth,\$read_group);  ## radical changes 
    ## Structure of hash .. key is locus, value is a single tab-delimited string.  Fields are as for %true_hps WITH ADDITION of [9] which is the frac_hp_depth

    ## NEW NEW NEW NEW
    log_it ($logfile,"\tFinal filtering of hairpins, and annotation of miRNA and hpRNA loci\.\.\.");
    %final_hps = final_hp(\%filtered_hps,\$genome,\$bamfile,\$minstrandfrac,\$maxmiRHPPairs,\$maxmiRUnpaired,\$read_group);
    log_it ($logfile," Done\n");
	   
    ## structure like %filtered_hps, but expanded:
    # [10] : hp_size_result
    # [11] : prec_result
    # [12] : duplex_result
    # [13] : star_result
    # [14] : hp_call
    # [15] : output string of alignment to be written
    
    log_it ($logfile,"\tFinishing parsing of hairpin-associated clusters\.\.\.");

    @final_clusters = get_final_clusters(\%final_hps,\@clusters);
    log_it($logfile," Done\n");
    ## get names again
    %names = get_names_simple(\@final_clusters);
    
    ## Write the HP and MIRNA files
    log_it ($logfile,"\tWriting miRNA and other hpRNA details to files\.\.\.");
    my($n_miRNAs,$n_hpRNAs) = write_files(\%final_hps, \$outdir, \%names, \@final_clusters);
    log_it ($logfile," Done\n");
    
    my $clus_count = scalar @final_clusters;
    log_it ($logfile,"\n\tTotal Clusters: $clus_count\n\tmiRNA loci: $n_miRNAs\n\tother hpRNA loci: $n_hpRNAs\n\n");
    
    ## now, we must re-quantify .. remove from the %quant_master hash entries for replaced clusters, and recalculate for the hairpin clusters with new coordinates.  Update all names too
    log_it($logfile,"\tRe-quantifying hpRNA and miRNA loci\.\.\.");
    requant(\%quant_master,\@clusters,\@final_clusters,\$bamfile,\$dicermin,\$dicermax,\$minstrandfrac,\$mindicerfrac,\$phasesize,\%names,\$read_group);
    log_it($logfile," Done\n\n");
}

## Phase 6 - Final output
log_it($logfile,`date`);
log_it($logfile,"Phase Six: Finalizing results and writing files\n");

my %quant_2 = ();
if($nohp) {
    foreach my $f_clus (@final_clusters) {
	my @qm_fields = split ("\t",$quant_master{$f_clus});
	my @new_fields = @qm_fields[2..9];
	unshift(@new_fields,"ND");
	unshift(@new_fields,$qm_fields[1]);  ## name
	unshift(@new_fields,$qm_fields[0]);  ## loc-coords
	
	push(@new_fields, "NA"); ## pairs
	push(@new_fields, "NA"); ## frac_paired
	push(@new_fields, "NA"); ## stem_length
	push(@new_fields, "NA"); ## loop_length
	push(@new_fields, "NA"); ## dGperStem
	push(@new_fields, "NA"); ## frac_hp
	push(@new_fields, "NA"); ## hp_size_result
	push(@new_fields, "NA"); ## prec_result
	push(@new_fields, "NA"); ## duplex_result
	push(@new_fields, "NA"); ## star_result
	
	for(my $ii = 10; $ii < (scalar @qm_fields); ++$ii) {
	    push(@new_fields, $qm_fields[$ii]);
	}
	my $new_entry = join("\t", @new_fields);
	$quant_2{$f_clus} = $new_entry;
    }
} else {
    %quant_2 = mod_quant_hp(\%quant_master,\%final_hps);
}

my %quant_3 = ();
# unless --raw was chosen, convert the raw values to mappings per million mapped
# OBSOLETED everything is "RAW" now.
#if($raw) {
%quant_3 = %quant_2;
#} else {
#    %quant_3 = convert_2_mpmm(\%quant_2,\$reads);
#}
log_it ($logfile,"\n\n");

## Assess overlap with flag_file loci, if any

my %quant_4 = ();
unless($flag_file eq "NULL") {
    log_it ($logfile,"\tAssessing overlap of small RNA loci with loci from flag file $flag_file\n");
}
%quant_4 = flag_overlap(\%quant_3, \$flag_file);

## after this conversion, a new column is inserted.  Now the fields are:
#[0] : locus
#[1] : name
#[2] : overlap      <<<<<<<<<
#[3] : HP_call
#[4] : strand
#[5] : frac_Watson
#[6] : total
#[7] : unique
#[8] : uniqueness_index
#[9] : dicer_call
#[10] : phase_offset
#[11] : phase_pval

#[12] : pairs          
#[13] : frac_paired    
#[14] : stem_length    
#[15] : loop_length    
#[16] : dGperStem      
#[17] : frac_hp        
#[18] : hp_size_result 
#[19] : prec_result    
#[20] : duplex_result  
#[21] : star_result   
#[22 -- end] : dicer size mappings


# output two gff3 files .. one for DCL and the other for NON-DCL loci.
# only do this if a de novo run
unless($count) {
    my ($dcl_file,$nfile) = write_gff3s(\@final_clusters,\%quant_4,\$outdir);
    
# write gff3 file information to log
    log_it ($logfile,"\nGFF3 files are $dcl_file and $nfile\n");
}

# output the master file
my $big_table = "$outdir" . "\/" . "Results\.txt";
open(BIG, ">$big_table");
print BIG "\#Locus\tName\tFlagOverlap\tHP\tStrand\tFrac_Watson\tTotal\tUniques\tUI\tDicerCall\tPhaseOffset\tPhase_pval\t";
print BIG "Pairs\tFracPaired\tStemLength\tLoopLength\tdGperStem\tFracCovHP\tHPSizeResult\tPrecisionResult\tDuplexResult\tStarResult\tShort\tLong";
for(my $d = $dicermin; $d <= $dicermax; ++$d) {
    print BIG "\t$d";
} 
print BIG "\n";
    
foreach my $f_clus (@final_clusters) {
    print BIG "$quant_4{$f_clus}\n";
}
close BIG;


summarize(\%quant_4,\$nohp,\$dicermin,\$dicermax,\$logfile);

##########################
# If the bamfile had > 1 read group AND no read_group was specified on the command line, assume
#  the user wants to finish up by separatlely quantify EACH of the read groups. Call ShortStack recursively, and 
#  quietly in --count mode

if(%rgs) {
    if((scalar(keys %rgs)) > 1) {
	unless($read_group) {
	    my $rg;
	    log_it($logfile,"\nBeginning count-mode analysis of each read-group separately.\n");
	    my $common_options = "--dicermin $dicermin --dicermax $dicermax --phasesize $phasesize --flag_file $flag_file --minstrandfrac $minstrandfrac --mindicerfrac $mindicerfrac";
	    while(($rg) = each %rgs) {
		log_it($logfile,"\tWorking on read group $rg ... ");
		my $rgoutdir = "$outdir" . "\/" . "rg_$rg";
		system "$0 --outdir $rgoutdir --count $big_table --bamfile $bamfile --read_group $rg $common_options $genome 2> /dev/null";
		log_it($logfile,"Done .. see $rgoutdir\n");
	    }
	}
    }
}



log_it ($logfile,"\nCompleted\n");
log_it ($logfile, `date`);


############### SUB ROUTINE BLOCKS
	
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
    my($version) = @_;
    my $usage = "\nShortStack.pl version $version

USAGE: ShortStack.pl \[options\] genome.fasta

MODES:
1. Trim, align, and analyze: Requires --untrimmedFA OR --untrimmedFQ OR --untrimmedCS. Also requires --adapter. To stop after alignments, specify --align_only
2. Align, and analyze: Requires --trimmedFA OR --trimmedFQ OR --trimmedCS. To stop after alignments, specify --align_only
3. Analyze: Requires --bamfile

Type \'ShortStack.pl --help\' for full list of options

DOCUMENTATION: type \'perldoc ShortStack.pl\'
";
    return $usage;
}

sub full_help {
    my($version) = @_;
    my $message = "\nShortStack.pl version $version

USAGE: ShortStack.pl \[options\] genome.fasta

MODES:
1. Trim, align, and analyze: Requires --untrimmedFA OR --untrimmedFQ OR --untrimmedCS. Also requires --adapter. To stop after alignments, specify --align_only
2. Align, and analyze: Requires --trimmedFA OR --trimmedFQ OR --trimmedCS. To stop after alignments, specify --align_only
3. Analyze: Requires --bamfile

DOCUMENTATION: type \'perldoc ShortStack.pl\'

OPTIONS:
--help : Print a help message and then quit.

--version : Print the version number and then quit.

--outdir [string] : Name of directory to be created to receive results of the run.  Deafults to \"ShortStack_[time]\", where time is \"UNIX time\" (the number of non-leap seconds since Jan 1, 1970 UCT), if not provided

--untrimmedFA [string] : Path to untrimmed small RNA-seq data in FASTA format. Multiple datasets can be provided as a comma-delimited list.

--untrimmedFQ [string] : Path to untrimmed small RNA-seq data in FASTQ format.  Multiple datasets can be provided as a comma-delimited list.

--untrimmedCS [string] : Path to untrimmed small RNA-seq data is Colorpsace-FASTA format. Multiple datasets can be provided as a comma-delimited list.

--untrimmedCSQV [string] : Path to untrimmed small RNA-seq quality values from SOLiD datasets, corresponding to colorspace-fasta files given in option --untrimmedCS. Multiple datasets can be provided as a comma-delimited list.

--adapter [string] : Sequence of 3' adapter to search for during adapter trimming. Must be at least 8 nts in length, and all ATGC characters. Required if either --untrimmedFA or --untrimmedFQ are specified. Multiple adapters (for when multiple input untrimmedFA/FQ files are specified) can be provided as a comma-delimited list.

--trimmedFA [string] : Path to trimmed and ready to map small RNA-seq data in FASTA format. Multiple datasets can be provided as a comma-delimited list.

--trimmedFQ [string] : Path to trimmed and ready to map small RNA-seq data in FASTQ format. Multiple datasets can be provided as a comma-delimited list.

--trimmedCS [string] : Path to trimmed and ready to map small RNA-seq data in Colorspace-FASTA format. Multiple datasets can be provided as a comma-delimited list.

--trimmedCSQV [string] : Path to trimmed color-space quality value files, corresponding to trimmed colorspace-fasta files provided in option --trimmedCS. Multiple datasets can be provided as a comma-delimited list.

--align_only : Exits program after completion of small RNA-seq data alignment, creating BAM file.

--bamfile [string] : Path to properly formatted and sorted BAM alignment file of small RNA-seq data.

--read_group [string] : Analyze only the indicated read-group. Read-group must be specified in the bam alignment file header. Default = [not active -- all reads analyzed]

--inv_file [string] : PATH to an einverted-produced .inv file of inverted repeats within the genome of interest.  Not required but strongly suggested for more complete annotations of hairpin-derived small RNA genes.  Default = {blank}.  Not needed for runs in \"nohp\" mode or runs in \"count\" mode (because \"count\" mode forces \"nohp\" mode as well).  A typical eniverted run uses default parameters except \"-maxrepeat 10000\", in order to capture long IRs.

--flag_file [string] : PATH to a simple file of genomic loci of interest.  The ShortStack-analyzed small RNA clusters will be analyzed for overlap with the loci in the flag_file .. if there is any overlap (as little as one nt), it will be reported.  Format for this file is describe below.

--mindepth [integer] : Minimum depth of mapping coverage to define an 'island'.  Default = 20.  Must be at least 2, more than 5 preferred.

--pad [integer] : Number of nucleotides upstream and downstream to extend initial islands during cluster definition.  Default = 100

--dicermin [integer] : Smallest size in the Dicer size range (or size range of interest).  Deafult = 20.  Must be between 15 and 35, and less than or equal to --dicermax

--dicermax [integer] : Largest size in the Dicer size range (or size range of interest).  Deafult = 24.  Must be between 15 and 35, and more than or equal to --dicermin

--minUI [float] : Minimum uniqueness index required to attempt RNA folding.  Must be a value between 0 and 1.  Zero forces all clusters to be folded; default: 0.1

--maxhpsep [integer] : Maximum allowed span for a base-pair during hairpin search with RNALfold; Also serves as the maximum size of genomic query to fold with RNALfold .. loci whose unpadded size is more than --maxhpsep will not be analyzed at all with RNALfold.  Default = 300.  Must be between 50 and 2000.

--minfracpaired [float] : Minimum fraction of paired nucleotides required within a valid hairpin structure.  Default = 0.67.  Allowed values are greater than 0 and less than or equal to 1.

--minntspaired [integer] : Minimum absolute number of paired nucleotides required within a valid hairpin structure.  Default = 15.  Allowed values are greater than zero and less than or equal to --maxhpsep

--maxdGperStem [float] : Maximum deltaG / stem length allowed in a valid hairpin structure.  Stem length is 0.5 * (left_stem_length + right_stem_length).  Default = -0.5

--minfrachpdepth [float] : Minimum fraction of corrected coverage within hairpin arms to keep hairpin for further analysis.  Default = 0.67.  Allowed values between 0 and 1.  See below for details.

--miRType [string] : Either \"plant\" or \"animal\".  Defaults to \"plant\".  This option sets --maxmiRHPPairs, --maxmiRUnpaired, and --maxLoopLength to 150, 5, and 100,000 respectively for type \"plant\".  For type \"animal\", the three are instead set to 45, 6, and 15, respectively.

--maxmiRHPPairs [integer] : Maximum number of base pairs in a valid MIRNA hairpin. default: set by --miRType \"plant\" to 150.  --miRType \"animal\" sets to 45 instead.  When provided, user settings will override miRType settings. 

--maxmiRUnpaired [integer] : Maximum number of unpaired miRNA nts in a miRNA/miRNA* duplex. default: set by --miRType \"plant\" to 5.  --miRType \"animal\" instead sets it to 6.  When provided, user settings will override miRType settings.

--maxLoopLength [integer] : maximum allowed loop length for a valid hairpin. default: set by --miRType \"plant\" be essentially unlimited (100,000).  --miRType \"plant\" sets it to 15.  When provided, user settings will override miRType settings.

--minstrandfrac [float] : Minimum fraction of mappings to one or the other strand call a polarity for non-hairpin clusters.  Also the minimum fraction of \"non-dyad\" mappings to the sense strand within potential hairpins/miRNAs to keep the locus annotated as a hp or miRNA.  See below for details.  Default = 0.8.  Allowed values between 0.5 and 1.

--mindicerfrac [float] : Minimum fraction of mappings within Dicer size range to annotate a locus as Dicer-derived.  Default = 0.85.  Allowed values between 0 and 1.

--phasesize [integer] : Examine phasing only for clusters dominated by the indicated size range.  Size must be within the bounds described by --dicermin and --dicermax.  Set to 'all' to examine p-values of each locus within the Dicer range, in its dominant size.  Set to 'none' to suppress all phasing analysis.  Default = 21.  Allowed values between --dicermin and --dicermax.

--count [string] : Invokes count mode, in which user-provided clusters are annotated and quantified instead of being defined de novo.  When invoked, the file provided with --count is assumed to contain a simple list of clusters.  Count mode also forces nohp mode.  Formatting details below.  Default : Not invoked.

--nohp : If \"--nohp\" appears on the command line, it invokes running in \"no hairpin\" mode.  RNA folding, hairpin annotation, and MIRNA annotation will be skipped (likely saving significant time).  Note that --count mode forces --nohp mode as well.  Default: Not invoked.

";
    return $message;
}

sub merge_clusters {
    my($input,$pad,$genome) = @_; ## passed my reference .. array and scalar
    my @output = ();

    my $this_start;
    my $this_stop;
    my $last_start;
    my $last_stop;
    my $last_padded_start;
    my $last_padded_stop;
    my $this_padded_start;
    my $this_padded_stop;
    my $last_chr = "null";
    my %chr_sizes = ();
    ## grab the chrom sizes, which you need to ensure that you don't pad off the end of the chroms
    # chrom sizes in column 1 from the fai file
    my $fai_file = "$$genome" . "\.fai";
    unless(-e $fai_file) {
	log_it($logfile, "\nFatal in sub-routine get_folding_regions : expected fai file $fai_file does not exist\n");
	exit;
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
	    $this_start = $2;
	    $this_padded_start = $this_start - $$pad;
	    if($this_padded_start < 1) {
		$this_padded_start = 1;
	    }
	    $this_stop = $3;
	    $this_padded_stop = $this_stop + $$pad;
	    if($this_padded_stop > $chr_sizes{$this_chr}) {
		$this_padded_stop = $chr_sizes{$this_chr};
	    }
	    
	    # special first case
	    if($last_chr eq "null") {
		$last_padded_start = $this_padded_start;
		$last_padded_stop = $this_padded_stop;
		$last_start = $this_start;
		$last_stop = $this_stop;
	    } elsif ($this_chr ne $last_chr) {
		$entry = "$last_chr" . ":" . "$last_start" . "-" . "$last_stop";
		push(@output,$entry);
		$last_padded_start = $this_padded_start;
		$last_padded_stop = $this_padded_stop;
		$last_start = $this_start;
		$last_stop = $this_stop;
	    } else {
		if($this_padded_start > $last_padded_stop) {
		    ## no overlap between these padded clusters.  Report the last one, trimming off its dangling pads
		    $entry = "$this_chr" . ":" . "$last_start" . "-" . "$last_stop";
		    push(@output,$entry);
		    $last_padded_start = $this_padded_start;
		    $last_padded_stop = $this_padded_stop;
		    $last_start = $this_start;
		    $last_stop = $this_stop;
		} else {
		    # here, same chr, this_padded_start is <= last_padded_stop, so we are merging
		    if($this_padded_start < $last_padded_start) {
			$last_padded_start = $this_padded_start;
		    }
		    if($this_start < $last_start) {
			$last_start = $this_start;
		    }
		    if($this_padded_stop > $last_padded_stop) {
			$last_padded_stop = $this_padded_stop;
		    }
		    if($this_stop > $last_stop) {
			$last_stop = $this_stop;
		    }
		}
	    }
	    $last_chr = $this_chr;
	} else {
	    log_it($logfile, "\nFATAL: in sub-routine \'merge_clusters\' : failed to parse initial locus $in_clus\n");
	    exit;
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
		log_it($logfile, "\nFATAL in sub-routine get_clusters_a_priori : cluster name $fields[0] is not understandable\n");
		exit;
	    }
	    push(@clusters,$fields[0]);
	    if(exists($tracker{$fields[0]})) {
		log_it($logfile, "\nFATAL in sub-routine get_clusters_a_priori : cluster $fields[0] appears more than once\!  Duplicate loci must be purged prior to analysis\n");
		exit;
	    } else {
		$tracker{$fields[0]} = 1;
	    }
	} else {
	    unless($_ =~ /^\S+:\d+-\d+$/ ) {
		log_it($logfile, "\nFATAL in sub-routine get_clusters_a_priori : cluster name $_ is not understandable\n");
		exit;
	    }
	    if(exists($tracker{$_})) {
		log_it($logfile, "\nFATAL in sub-routine get_clusters_a_priori : cluster $fields[0] appears more than once\!  Duplicate loci must be purged prior to analysis\n");
		exit;
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
		log_it($logfile, "\nFATAL in sub-routine get_names_countmode : cluster name $fields[0] is not understandable\n");
		exit;
	    }
	    if($fields[1]) {
		$names{$fields[0]} = $fields[1];
	    } else {
		log_it($logfile, "\nFATAL in sub-routine get_names_countmode : file is tab-delmited but no name entry in second column\n");
		exit;
	    }
	} else {
	    unless($_ =~ /^\S+:\d+-\d+$/ ) {
		log_it($logfile, "\nFATAL in sub-routine get_names_countmode : cluster name $_ is not understandable\n");
		exit;
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
    my($gen_file,$info,$maxhpsep) = @_;  ## passed as references .. scalar, array,and scalar, respectively
    my $locus;

    my %to_fold = ();
    my $chr;
    my $start;
    my $stop;
    my $locus_size;
#    my $unp_locus_size;
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
    log_it($logfile,"\tProgress in sub-routine \"get_folding_regions\" \(dot = five percent\): ");
    
    foreach $locus (@$info) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    log_it($logfile, ".");
	    $x = 0;
	}
	# parse locus name
	if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    log_it($logfile, "\nFatal in sub-routine get_folding_regions: could not parse locus name $locus\n");
	    exit;
	}

	$locus_size = $stop - $start + 1;
	#$unp_locus_size = int($locus_size - (2 * $$pad));  ## the unpadded region is the locus size - (2 * $pad).

	if($locus_size > $$maxhpsep) {
	    next;  ## NOT folded if locus size is larger than the max hp sep
	} else {
	    $get_size = $$maxhpsep;  ## otherwise, a window the size of the maxhpsep will be folded in all cases.
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
	    log_it($logfile, "\nFatal in sub-routine get_folding_regions : expected fai file $fai_file does not exist\n");
	    exit;
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
	    log_it($logfile, "\nFatal in sub-routine get_folding_regions : failed to get chromosome length for $chr\n");
	    exit;
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
    log_it($logfile," Done\n");
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
    log_it($logfile,"\tProgress in sub-routine \"get_folding_regions_countmode\" \(dot = five percent\): ");
    
    foreach $locus (@$info) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    log_it($logfile,".");
	    $x = 0;
	}
	# parse locus name
	if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    log_it($logfile, "\nFatal in sub-routine get_folding_regions: could not parse locus name $locus\n");
	    exit;
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
    log_it($logfile," Done\n");
    # return the hash
    return %to_fold;
}

sub folder {
    my($to_fold,$L,$min_paired,$min_frac_paired,$maxdGperStem,$maxLoopLength) = @_;  ## passed by reference: first is hash, rest are scalar
    my $locus;
    my $fold_seq;
    my $brax;

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
    
    my $pairs;
    my $frac_paired;
    my $stem_length;
    my $loop_length;
    
    my $seq_section;
    my $left_stem_brax;
    my $right_stem_brax;
    my $left_stem_seq;
    my $right_stem_seq;
    
    my $eval_string;
    my $dG;
    my $dGperStem;
    
    # first, get some information to enable a crude progress bar
    my $n_loci_to_fold = scalar ( keys %$to_fold);
    my $x = 0;
    my $five_percent = int(0.05 * $n_loci_to_fold);
    log_it($logfile,"\n\tProgress in sub-routine \"folder\" \(dot = 5 percent\): ");
    
    while(($locus) = each %$to_fold) {
	# progress tracking
	++$x;
	if($x == $five_percent) {
	    log_it($logfile,".");
	    $x = 0;
	}
	
	# First, the Watson Strand is folded
	## compatible with either RNALfold 1.x or 2.x ... removed -noLP
	
	open(RNALFOLD, "echo $$to_fold{$locus}{'sequence'} | RNALfold -L $$L |");
	while (<RNALFOLD>) {
	    chomp;
	    if($_ =~ /^[^\.\(\)]/) {
		## this is a line that does not have a structure .. most likely the last line of output which just has the input sequence .. ignore it
		next;
	    }

	    
	    if($_ =~ /^\S+/) {
		$brax = $&;
	    } else {
		log_it($logfile,"\nFATAL in sub-routine \'folder\' : failed to parse brackets from RNALfold output line:\n$_\n"); 
		exit;
	    }
	    
	    if($_ =~ /(\d+)\s*$/) {
		$start = $1;
	    } else {
		log_it($logfile, "\nFATAL in sub-routine \'folder\' : failed to parse start position from RNALfold output line $_\n");
		exit;
	    }
	    
	    # new v0.1.3 .. split the RNALfold-provided structure into distinct hairpins .. this fixes a major problem in earlier versions with false negatives on hairpin discovery!
	    %regions = get_distinct_hps($brax);
	    while(($left, $right) = each %regions) {
		$brax_section = substr($brax,($left - 1),($right-$left+1));
		$seq_section = substr($$to_fold{$locus}{'sequence'}, ($start + $left - 2), ($right - $left + 1));
		($structure_details,$pairs,$frac_paired,$stem_length,$loop_length) = evaluate_structure_general($brax_section,$$min_paired,$$min_frac_paired,$$maxLoopLength);

		# if it is OK so far, test via RNAeval of just the stem regions
		unless($structure_details eq "bogus") {
		    # structure details are 123-150,180-200 .. e.g. start and stop positions of the helix of interest, one-based coordinates, relative to the brackets themselves.
		    
		    if($structure_details =~ /(\d+)-(\d+)\,(\d+)-(\d+)/) {
			$left_stem_brax = substr($brax_section, ($1 - 1), ($2 - $1 + 1));
			$right_stem_brax = substr($brax_section, ($3 - 1), ($4 - $3 + 1));
			$left_stem_seq = substr($seq_section, ($1 - 1), ($2 - $1 + 1));
			$right_stem_seq = substr($seq_section, ($3 - 1), ($4 - $3 + 1));
			
			$eval_string = "$left_stem_seq" . "\&" . "$right_stem_seq" . "\n" . "$left_stem_brax" . "\&" . "$right_stem_brax";
			
			$dG = dG_from_RNAeval($eval_string);
			
			unless($dG eq "FAIL") {
			    $dGperStem = $dG / $stem_length;
			    if($dGperStem <= $$maxdGperStem) {
				
				## OK
		    
				$adj_st = $left + $start - 1;
				$entry = "$brax_section\t$adj_st\t$structure_details\t" . "Watson";
				$entry .= "\t$pairs\t$frac_paired\t$stem_length\t$loop_length\t$dGperStem";
				
				push(@{$output{$locus}}, $entry);
			    }
			}
		    }
		}
	    }
	}
	close RNALFOLD;
	
	# Now, examine the Crick strand .. revcomp the sequence
	$revcomp = reverse $$to_fold{$locus}{'sequence'};
	$revcomp =~ tr/ACUG/UGAC/;

	open(RNALFOLD, "echo $revcomp | RNALfold -L $$L |");
	while (<RNALFOLD>) {
	    chomp;
	    if($_ =~ /^[^\.\(\)]/) {
		## this is a line that does not have a structure .. most likely the last line of output which just has the input sequence .. ignore it
		next;
	    }

	    
	    if($_ =~ /^\S+/) {
		$brax = $&;
	    } else {
		log_it($logfile,"\nFATAL in sub-routine \'folder\' : failed to parse brackets from RNALfold output line:\n$_\n");
		exit;
	    }
	    
	    if($_ =~ /(\d+)\s*$/) {
		$start = $1;
	    } else {
		log_it($logfile,"\nFATAL in sub-routine \'folder\' : failed to parse start position from RNALfold output line $_\n");
		exit;
	    }
	    
            # new v0.1.3 .. split the RNALfold-provided structure into distinct hairpins .. this fixes a major problem in earlier versions with false negatives on hairpin discovery!
	    %regions = get_distinct_hps($brax);
	    while(($left, $right) = each %regions) {
		$brax_section = substr($brax,($left - 1),($right-$left+1));
		$seq_section = substr($revcomp, ($start + $left - 2), ($right - $left + 1));


		($structure_details,$pairs,$frac_paired,$stem_length,$loop_length) = evaluate_structure_general($brax_section,$$min_paired,$$min_frac_paired,$$maxLoopLength);
		
		# if it is OK, adjust the coordinates, and add to output
		unless($structure_details eq "bogus") {
		    if($structure_details =~ /(\d+)-(\d+)\,(\d+)-(\d+)/) {
			$left_stem_brax = substr($brax_section, ($1 - 1), ($2 - $1 + 1));
			$right_stem_brax = substr($brax_section, ($3 - 1), ($4 - $3 + 1));
			$left_stem_seq = substr($seq_section, ($1 - 1), ($2 - $1 + 1));
			$right_stem_seq = substr($seq_section, ($3 - 1), ($4 - $3 + 1));
			
			$eval_string = "$left_stem_seq" . "\&" . "$right_stem_seq" . "\n" . "$left_stem_brax" . "\&" . "$right_stem_brax";
			$dG = dG_from_RNAeval($eval_string);
			
			unless($dG eq "FAIL") {
			    $dGperStem = $dG / $stem_length;
			    if($dGperStem <= $$maxdGperStem) {
				
				## OK
				
				$adj_st = $left + $start - 1;
				$entry = "$brax_section\t$adj_st\t$structure_details\t" . "Crick";
				$entry .= "\t$pairs\t$frac_paired\t$stem_length\t$loop_length\t$dGperStem";
				
				push(@{$output{$locus}}, $entry);
			    }
			}
		    }
		}
	    }
	}
	close RNALFOLD;
    }
    log_it($logfile," Done\n");
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
    my($brax,$min_paired,$min_paired_frac,$maxLoopLength) = @_;
    my @chars = split ('', $brax);
    my $char;
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    my $i;
    my $left_true_stop;
    my $right_true_start;
    my $pairs = 0;
    my $frac_paired = 0;
    my $details;  ## left_start-left_true_stop,right_true_start-right_stop
    my $stem_length = 0;
    my $loop_length = 0;
    
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
	    
	    # Determine the stem length, which is defined as 0.5 * (left_arm_length + right_arm_length)
	    $stem_length = 0.5 * (($right_stop - $right_true_start + 1) + ($left_true_stop - $left_start + 1));
	    
	    # Determine the loop length
	    $loop_length = $right_true_start - $left_true_stop - 1;
	    
	    # If everything is OK, send it back.  RNAeval will occur in the parent sub-routine
	    
	    if(($frac_paired >= $min_paired_frac) and
	       ($pairs >= $min_paired) and
	       ($loop_length <= $maxLoopLength)) {
		
		$details = "$left_start" . "-$left_true_stop" . "," . "$right_true_start" . "-" . "$right_stop";
	    } else {
		$details = "bogus";
	    }
	} else {
	    $details = "bogus";
	}
    } else {
	$details = "bogus";
    }
    return ($details,$pairs,$frac_paired,$stem_length,$loop_length);
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
	
	@inputs = @{$$input_hps{$locus}};
	foreach $inp (@inputs) {
	    @fields = split ("\t", $inp);
	    $strand_folded = $fields[3];

	    $in_helix_coords = $fields[2];
	    $local_offset = $fields[1];
	    $brax = $fields[0];
	    
	    if($in_helix_coords =~ /^(\d+)-(\d+),(\d+)-(\d+)$/) {
		$left_1_start = $local_offset + $1 - 1;
		$left_1_stop = $local_offset + $2 - 1;
		$right_1_start = $local_offset + $3 -1;
		$right_1_stop = $local_offset + $4 - 1;
		$start_1 = $local_offset;
		$stop_1 = $start_1 + (length $brax) - 1;
	    } else {
		log_it($logfile,"\nFatal error in sub-routine hairpin_coords: failed to parse local hairpin coordinates $in_helix_coords\n");
		exit;
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
	    $true_entry .= "\t$fields[4]\t$fields[5]\t$fields[6]\t$fields[7]\t$fields[8]";
	    
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

	# get the range of the actual locus
	@loc_range = ();  ## clear at each iteration
	if($locus =~ /\S+:(\d+)-(\d+)$/) {
	    push(@loc_range,$1);
	    push(@loc_range,$2);
	} else {
	    log_it($logfile, "\nFatal in sub-routine remove_nonoverlapped_hps -- reg ex failure to larse locus name $locus\n");
	    exit;
	}
	@entries = @{$$input_hash{$locus}};

	foreach $entry (@entries) {

	    @en_fields = split ("\t", $entry);
	    @hp_ranges = split (",", $en_fields[2]);
	    @range1 = split ("-", $hp_ranges[0]);
	    @range2 = split ("-", $hp_ranges[1]);
	    $overlap1 = range_overlap(\@loc_range,\@range1);
	    $overlap2 = range_overlap(\@loc_range,\@range2);

	    if(($overlap1) or
	       ($overlap2)) {
		
		## criteria for overlap is the bare min of a single overlapping nt with one of the two helical arms
		
		push(@{$scrubbed{$locus}}, $entry);
		

		
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
    my($input_hps,$bamfile,$min_frac,$read_group) = @_; # passed by reference, hash and four scalars
   
    my %original_locus_W = ();
    my %original_locus_C = ();
    my $orig_coverage;
    my $orig_start;
    my $orig_stop;
    my @fields = ();
    my $position;
    my $read_length;
 
    my $cluster;
    my @hps = ();
    my $hp_entry;
    my @hp_fields = ();
    my @hp_arms = ();
    my @fivep = ();
    my @threep = ();
    
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    
    my %output_hash = ();
    
    my $i;
    
    my $hp_coverage;
   
    my $hp_cov_frac;
    
    my @alignments = ();
    my $hp_strand;
    
    my %left_names = ();
    my %right_names = ();
    
    my $stop;
    my $readstrand;
    my $readname;
    
    my %left_read_lengths = ();
    
    my $max_frac = 0;
    my $max_string;
    
    # for progress tracking
    my $n_to_examine = scalar (keys %$input_hps);  ## this is a bit crude, as loci will have varying numbers of hairpins within them
    my $five_percent = int (0.05 * $n_to_examine);
    my $x = 0;
    log_it($logfile,"\n\tProgress within sub-routine \"hp_expression\" \(dots = five percent\): ");
    
    
    
    
    while(($cluster) = each %$input_hps) {
	# progress tracking
	++$x;
	if($x >= $five_percent) {
	    log_it($logfile,".");
	    $x = 0;
	}
	# first, get the coverage total across the entire, original locus
	%original_locus_W = ();
	%original_locus_C = ();
	$orig_coverage = 0;
	if($cluster =~ /^\S+:(\d+)-(\d+)$/) {
	    $orig_start = $1;
	    $orig_stop = $2;
	} else {
	    log_it($logfile, "\nFATAL: Failed to parse cluster name $cluster in sub-routine hp_expression\n");
	    exit;
	}
	for($i = $orig_start; $i <= $orig_stop; ++$i) {
	    $original_locus_W{$i} = 0;
	    $original_locus_C{$i} = 0;
	}
	## save all alignments 
	## reset 
	@alignments = ();
	if($$read_group) {
	    open(SAM, "samtools view -F 0x4 -r $$read_group $$bamfile $cluster |");
	} else {
	    open(SAM, "samtools view -F 0x4 $$bamfile $cluster |");
	}
	while (<SAM>) {
	    push(@alignments, $_);
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
		if(exists($original_locus_W{$i})) {  ## exists will be same for both strands, so just using W is fine
		    ++$orig_coverage;
		    if($fields[1] & 16) {
			++$original_locus_C{$i};
		    } else {
			++$original_locus_W{$i};
		    }
		}
	    }
	}
	close SAM;
	
	@hps = @{$$input_hps{$cluster}};
	foreach $hp_entry (@hps) {
	    
	    # reset
	    $max_frac = 0;
	    $max_string = '';
	    
	    @hp_fields = split ("\t", $hp_entry);
	    @hp_arms = split (",", $hp_fields[2]);
	    @fivep = split ("-", $hp_arms[0]);
	    @threep = split ("-", $hp_arms[1]);
	    
	    if($fivep[0] < $fivep[1]) { ## sense 
		$hp_strand = "W";
		$left_start = $fivep[0];
		$left_stop = $fivep[1];
		$right_start = $threep[0];
		$right_stop = $threep[1];
	    } else {
		## antisense
		$hp_strand = "C";
		$left_start = $threep[1];
		$left_stop = $threep[0];
		$right_start = $fivep[1];
		$right_stop = $fivep[0];
	    }
	    
	    ## check the coverage within the current hairpin
	    ## at this point, not yet corrected for any 'dyads' formed by mapping of the same read to opposite strands and arms as
	    ##  would happen for perfect inverted repeats
	    
	    $hp_coverage = 0;  ## reset each time through the loop
	    for($i = $left_start; $i <= $left_stop; ++$i) {
		if($hp_strand eq "W") {
		    if(exists($original_locus_W{$i})) {
			$hp_coverage += $original_locus_W{$i};
		    }
		} else {
		    if(exists($original_locus_C{$i})) {
			$hp_coverage += $original_locus_C{$i};
		    }
		}
	    }
	    for($i = $right_start; $i <= $right_stop; ++$i) {
		if($hp_strand eq "W") {
		    if(exists($original_locus_W{$i})) {
			$hp_coverage += $original_locus_W{$i};
		    }
		} else {
		    if(exists($original_locus_C{$i})) {
			$hp_coverage += $original_locus_C{$i};
		    }
		}
	    }
	    
	    ## determine the fraction of the original locus's coverage in the current hairpin
	    ## protect against dividing by zero
	    if($orig_coverage > 0) {
		$hp_cov_frac = $hp_coverage / $orig_coverage;
	    } else {
		$hp_cov_frac = 0;
	    }
	    
	    ## At this point, the fraction has not been corrected for dyads
	    ## Consider the limiting case .. 100% of coverage is in the hairpin arms, and all
	    ##  mappings are dyads.  Thus, the fraction would be 0.5.  To save time then, we can quit
	    ##  on this hairpin if this preliminary fraction is < 0.5
	    

	    if ($hp_cov_frac < 0.5) {
		next;
	    } else {
		## determine the hp_cov_frac, corrected for dyads
		##  for any possible dyads.
		
		## first, get the read names for any possible dyads
		# reset
		%left_names = ();
		%right_names = ();
		foreach my $samline (@alignments) {
		    chomp $samline;
		    ## ignore headers
		    if($samline =~ /^@/) {
			next;
		    }
		    ## ignore unmapped reads
		    @fields = split ("\t", $samline);
		    if($fields[1] & 4) {
			next;
		    }
		    ## determine stop position
		    $read_length = parse_cigar ($fields[5]);
		    $stop = $fields[3] + $read_length - 1;
		    
		    ## determine strand
		    if($fields[1] & 16) {
			$readstrand = "C";
		    } else {
			$readstrand = "W";
		    }
		    
		    ## is it fully within the left arm?
		    if(($fields[3] >= $left_start) and
		       ($stop <= $left_stop)) {
			$left_names{$fields[0]} = $readstrand;
			$left_read_lengths{$fields[0]} = $read_length;
		    }
		    
		    ## is it fully within the right arm?
		    if(($fields[3] >= $right_start) and
		       ($stop <= $right_stop)) {
			$right_names{$fields[0]} = $readstrand;
		    }
		}
		
		while(($readname,$readstrand) = each %left_names) {
		    if(exists($right_names{$readname})) {
			if($readstrand ne $right_names{$readname}) {
			    ## this is a dyad
			    ## substract its length from the total coverage
			    ## this is the equivalent of just not counting the 
			    ## opposite polarity mapping in the total
			    $read_length = $left_read_lengths{$readname};
			    $orig_coverage -= $read_length;
			}
		    }
		}
		
		## recalculate the fraction and reevaluate
		if($orig_coverage > 0) {
		    $hp_cov_frac = $hp_coverage / $orig_coverage;
		} else {
		    $hp_cov_frac = 0;
		}
		
		## if it meets the threshold, record it
		if($hp_cov_frac >= $$min_frac) {
		    if($hp_cov_frac > $max_frac) {
			$max_frac = $hp_cov_frac;
			$max_string = "$hp_entry" . "\t" . "$max_frac";   
		    }
		}
	    }
	    ## If theres any that qualify, record it.  It will be the one with max fraction.  In case of tie, arbitrary one selected 
	    if($max_frac >= $$min_frac) {
		$output_hash{$cluster} = $max_string;
	    }
	}
    }
    
    # close progress tracking
    log_it($logfile,"  Done\n");
    return %output_hash;
}

#sub get_hp_clusters {
#    my ($input_hash) = @_;  ## passed by reference
#    my %output = ();

#    my $orig;
#    my $chr;

#    my @en_fields = ();
#    my @hp_coords = ();
#    my $new_hp_cluster;
#    my $start;
#    my $stop;
    
#    while(($orig) = each %$input_hash) {
#	if($orig =~ /^(\S+):/) {
#	    $chr = $1;
#	} else {
#	    die "FATAL in sub-routine \"get_hp_clusters\" : could not parse chromosome name from $orig\n";
#	}
#	@en_fields = split ("\t", $input_hash{$orig});
#	@hp_coords = split ("-", $en_fields[1]);
#	if($hp_coords[0] < $hp_coords[1]) {
#	    $start = $hp_coords[0];
#	    $stop = $hp_coords[1];
#	} else {
#	    $start = $hp_coords[1];
#	    $stop = $hp_coords[0];
#	}
	
#	$new_hp_cluster = "$chr" . ":" . "$start" . "-" . "$stop";
#	$output{$new_hp_cluster} = $entry;
#    }
#    return %output;
#}

#sub get_hp_clusters_countmode {
#    my ($input_hash) = @_;  ## passed by reference
#    my %output = ();
#    my $orig;
#    my $orig_size;
#    my $hp_size;
#    my @en_fields = ();
#    my @hp_coords = ();
#    
#    while(($orig) = each %$input_hash) {
#	if($orig =~ /^\S+:(\d+)-(\d+)/) {
#	    $orig_size = $2 - $1 + 1;
#	} else {
#	    # failsafe
#	    $orig_size = 1000000000000000;
#	}
#	@en_fields = split ("\t", $$input_hash{$orig});
#	@hp_coords = split ("-", $en_fields[1]);
#	$hp_size = abs ($hp_coords[1] - $hp_coords[0]);
#	    if(($hp_size / $orig_size) > 0.9) {  ## 90% rule .. in count mode, HP or MIRNA is kept only if there is just one within the original cluster that accounts for more than 90% of the original cluster's length
#		$output{$orig} = $entries[0];  ## no padding, only the first entry is kept if there was more than one hp
#	    }
#	}
#    }
#    return %output;
#}

sub get_final_clusters {
    my($input_hp_hash,$orig_clus) = @_; # passed by reference ..hash-array
    my $orig;
    my @output = ();
    my $chr;
    
    ## new in 0.4.4
    my @pre_output = ();
    
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
		log_it($logfile,"\nFATAL in sub-routine \"get_hp_clusters\" : could not parse chromosome name from $orig\n");
		exit;
	    }
	    @en_fields = split ("\t", $$input_hp_hash{$orig});
	    @hp_coords = split ("-", $en_fields[1]);
	    if($hp_coords[0] < $hp_coords[1]) {
		$start = $hp_coords[0];
		$stop = $hp_coords[1];
	    } else {
		$start = $hp_coords[1];
		$stop = $hp_coords[0];
	    }
	    
	    $new_hp_cluster = "$chr" . ":" . "$start" . "-" . "$stop";
	    
	    # add to array, after checking failsafe
	    unless(exists($failsafe{$new_hp_cluster})) {
		push(@pre_output,$new_hp_cluster);
	    }
	    $failsafe{$new_hp_cluster} = 1;
	    
	    # add new entry to hash and then delete the old one
	    # unless the coordinates of the original cluster are the same
	    unless($new_hp_cluster eq $orig) {
		$$input_hp_hash{$new_hp_cluster} = $$input_hp_hash{$orig};
		delete $$input_hp_hash{$orig};
	    }
	    
	} else {
	    unless(exists($failsafe{$orig})) {
		push(@pre_output,$orig);
	    }
	    $failsafe{$orig} = 1;
	}
    }
    
    ## new in 0.4.4 .. merge overlapping clusters
    my $last_chr = "NULL";
    my $last_end;
    my $last_start;
    my $pre_chr;
    my $pre_end;
    my $pre_start;
    my $last_locus;
    my $size_last;
    my $size_pre;
    foreach my $pre_out_locus (@pre_output) {
	if($pre_out_locus =~ /(\S+):(\d+)-(\d+)/) {
	    $pre_chr = $1;
	    $pre_end = $3;
	    $pre_start = $2;
	    
	    if($last_chr eq "NULL") {
		## first time through, no comparison to make
		$last_locus = $pre_out_locus;
		$last_chr = $pre_chr;
		$last_start = $pre_start;
		$last_end = $pre_end;
		
	    } elsif (($last_chr ne $pre_chr) or
	       ($pre_start > $last_end) or
	       ($last_start > $pre_end)) {
		
		## no overlap, add the last one to the output
		push(@output, $last_locus);
		
		## reset
		$last_locus = $pre_out_locus;
		$last_chr = $pre_chr;
		$last_start = $pre_start;
		$last_end = $pre_end;
	    } else {
		## there is overlap
		## if only one of the two is a HP, keep it.
		## if both, or none,  are HP, keep the longer one
		if((exists($$input_hp_hash{$last_locus})) and
		   (exists($$input_hp_hash{$pre_out_locus}))) {
		    $size_last = $last_end - $last_start +1;
		    $size_pre = $pre_end - $pre_start + 1;
		    
		    if($size_last < $size_pre) {
			## deleted = $last_locus;
			delete $$input_hp_hash{$last_locus};
			
			## judgement suspended on the current one, current gets set as last
			$last_locus = $pre_out_locus;
			$last_chr = $pre_chr;
			$last_start = $pre_start;
			$last_end = $pre_end;
			
		    } else {
			## deleted = $pre_out_locus;
			delete $$input_hp_hash{$pre_out_locus};
			
			## judgement still suspend on the last one, which remains the last one, so no resetting required here
		    }
		    
		} elsif (exists($$input_hp_hash{$last_locus})) {
		    ## deleted = $pre_out_locus;
			## judgement still suspend on the last one, which remains the last one, so no resetting required here
		    
		    
		} elsif (exists($$input_hp_hash{$pre_out_locus})) {
		    ## deleted = $last_locus;
		    
		    ### simply store the current one as last in prep for next time through the loop
			$last_locus = $pre_out_locus;
			$last_chr = $pre_chr;
			$last_start = $pre_start;
			$last_end = $pre_end;
		    
		} else {
		    $size_last = $last_end - $last_start +1;
		    $size_pre = $pre_end - $pre_start + 1;
		    
		    if($size_last < $size_pre) {
			## deleted = $last_locus;
			## judgement suspended on the current one, current gets set as last
			$last_locus = $pre_out_locus;
			$last_chr = $pre_chr;
			$last_start = $pre_start;
			$last_end = $pre_end;
			
			
		    } else {
			## deleted = $pre_out_locus;
			
			## judgement still suspend on the last one, which remains the last one, so no resetting required here

		    }
		}
	    }
	}
    }
    ## don't forget the last one
    push(@output, $last_locus);
    
    ## end new part from 0.4.4
    
    return @output;
}

sub hp_output {  ## modified as of 0.2.1 and then again 0.3.0
    my ($hp_hash,$genome,$bamfile,$n_mapped_reads,$minhpfrac,$minstrandfrac,$maxmiRHPPairs,$maxmiRUnpaired) = @_;  ## passed by reference.  hash, scalar, scalar, scalar, scalar, scalar, scalar,sclara
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
    
    my $n_of_pairs; ## 0.3.0
    
    ## for progress bar
    my $n_to_process = scalar ( keys %$hp_hash);
    my $five_percent = int ($n_to_process / 20);
    my $progress = 0;
    
    log_it($logfile,"\tProgress in sub-routine \"hp_output\" \(dot = five percent\): ");
    
    while(($locus) = each %$hp_hash) {
	# progress tracking
	++$progress;
	if($progress >= $five_percent) {
	    log_it($logfile,"\.");
	    $progress = 0;
	}

	@entries = @{$$hp_hash{$locus}};
	@kept_entries = ();
	foreach $locus_data (@entries) {
	    
	    # parse the data
	    @l_data_fields = split ("\t", $locus_data);
	    $brax = $l_data_fields[0];
	    # get the total number of paired nts in this structure
	    $n_of_pairs = 0;  ## reset each time
	    while ($brax =~ /\(/g) {
		++$n_of_pairs;
	    }
	    
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
		log_it($logfile,"\nFATAL: in sub-routine \"hp_output\" could not parse chr name from $locus\n");
		exit;
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
	    
	    open(SAM, "samtools view -F 0x4 $$bamfile $sam_query |");
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

	    # check whether the hairpin itself qualifies in terms of maxmiRHPPairs.... if it fails that, dump the @candidates array
	    unless($n_of_pairs <= $$maxmiRHPPairs) {
		@candidates = ();
		${$output{$sam_query}}[0] = 0;
	    }
	    
	    
            # check whether candidates span a loop and if not whether they have maxmiRHPUnpaired or fewer un-paired residues (ignoring the last two nts)
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
		    if($n_unpaired_cand > $$maxmiRUnpaired) {
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
		    
		    
		    
		    # first check whether the star sequence has maxmiRUnpaired or fewer mismatches, omitting the last two nts
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
			if($n_unpaired_star > $$maxmiRUnpaired) {
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
				    log_it($logfile,"\nFATAL error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
				    exit;
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
				    log_it($logfile,"\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
				    exit;
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
				    log_it($logfile,"\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
				    exit;
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
				    log_it($logfile, "\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
				    exit;
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
	    $earlier_fail = 0;  ## reset
	    # <= maxmiRPairs ?
	    if(${$output{$sam_query}}[0] >= 1) {
		$locus_data .= "\tPASS";
	    } else {
		$locus_data .= "\tFAIL";
		$earlier_fail = 1;
	    }


            # Precision
	    if(${$output{$sam_query}}[2] >= 1) {
		$locus_data .= "\tPASS";
	    } elsif ($earlier_fail) {
		$locus_data .= "\tND";
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
    log_it($logfile," Done\n");  ## closes progress
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
    
    my $lkeymax;
    my $lkeymin;
    my $rkeymax;
    my $rkeymin;

    foreach my $b (@b_pos) {
	if($$brax_hash{$b} eq "\(") {
	    push(@left,$b);
	} elsif($$brax_hash{$b} eq "\)") {
	    $c = pop @left;
	    $leftpairs{$c} = $b;

            if($lkeymax) {
                if($c > $lkeymax) {
                    $lkeymax = $c;
                }
            } else {
                $lkeymax = $c;
            }
            if($lkeymin) {
                if($c < $lkeymin) {
                    $lkeymin = $c;
                }
            } else {
                $lkeymin = $c;
            }

	    $rightpairs{$b} = $c;
	    
	    if($rkeymax) {
                if($b > $rkeymax) {
                    $rkeymax = $b;
                }
            } else {
                $rkeymax = $b;
            }
            if($rkeymin) {
                if($b < $rkeymin) {
                    $rkeymin = $b;
                }
            } else {
                $rkeymin = $b;
            }
	    
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
		if($i < ($lkeymin - 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_start = $leftpairs{$i} - $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($leftpairs{$i})) {
		++$i;
		++$j;
		if($i > ($lkeymax + 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_stop = $leftpairs{$i} + 2 + $j;
	    
	} elsif($strand eq "-") {
	    $j = 0;
	    $i = ($mir_coords[1]) + 2;
	    until(exists($leftpairs{$i})) {
		++$i;
		++$j;
		if($i > ($lkeymax + 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_start = $leftpairs{$i} + $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($leftpairs{$i})) {
		--$i;
		++$j;
		if($i < ($lkeymin - 1)) {
                    return "fail";
                    last;
                }
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
		if($i < ($rkeymin - 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_start = $rightpairs{$i} - $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($rightpairs{$i})) {
		++$i;
		++$j;
		if($i > ($rkeymax + 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_stop = $rightpairs{$i} + 2 + $j;
	} elsif ($strand eq "-") {
	    $j = 0;
	    $i = $mir_coords[1] + 2;
	    until(exists($rightpairs{$i})) {
		++$i;
		++$j;
		if($i > ($rkeymax + 1)) {
                    return "fail";
                    last;
                }
	    }
	    $star_start = $rightpairs{$i} + $j;
	    
	    $j = 0;
	    $i = $mir_coords[0];
	    until(exists($rightpairs{$i})) {
		--$i;
		++$j;
		if($i < ($rkeymin - 1)) {
                    return "fail";
                    last;
                }
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
    my($clus_array,$bamfile,$dicer_min,$dicer_max,$strand_cutoff,$dicer_cutoff,$phasesize,$names,$read_group) = @_; ## passed by reference .. first one is array, hp_hash and miR_hash are hashes, others scalars
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
    
    my $ui;
    
    my $x;
    
    # for progress tracking
    my $n_2_analyze = scalar @$clus_array;
    my $five_percent = int($n_2_analyze/20);
    my $n_analyzed = 0;
    log_it($logfile, "\tProgress in sub-routine \"quant\" \(dot = five percent\): ");
    
    foreach my $locus (@$clus_array) {
        ## progress
	++$n_analyzed;
	if($n_analyzed >= $five_percent) {
	    ## TEST
	    ##last;
	    ## END TEST
	    log_it($logfile,"\.");
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
	    log_it($logfile, "\nFATAL in sub-routine \"quant\" : could not parse coordinates from locus $locus\n");
	    exit;
	}
	# initialize phase hash
	for($i = $loc_start; $i <= $loc_stop; ++$i) {
	    $phase_hash{$i} = 0;
	}
	if($$read_group) {
	    open(SAM, "samtools view -F 0x4 -r $$read_group $$bamfile $locus |");
	} else {
	    open(SAM, "samtools view -F 0x4 $$bamfile $locus |");
	}
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

	    # determine the number of total mappings for the read from the XX:i tag
	    $nh = '';
	    if($_ =~ /\tXX:i:(\d+)/) {
		$nh = $1;
	    }
	    unless($nh) {
		log_it($logfile,"\nFATAL in sub-routine \"quant\": Could not parse XX:i flag from sam line $_\n");
		exit;
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
	
	# [6] UNIQUENESS INDEX (USED TO BE REP-TOTAL)
	if($total) {
	    $ui = sprintf("%.4f",($repnorm/$total));
	} else {
	    $ui = "NA";
	}
	$output{$locus} .= "\t$ui";
	
	# [7] DICER ... either 'N' or a number within the dicer_min to dicer_max range
	# was calculated above
	$output{$locus} .= "\t$dicer_call";
	
	# [8] and [9] .. phase offset and phase p-value
	if(exists($phase_p_values{$locus})) {
	    $output{$locus} .= "\t$phase_offsets{$locus}\t$phase_p_values{$locus}";
	} else {
	    $output{$locus} .= "\tND\tND";
	}
		
	# [10] SHORT
	$output{$locus} .= "\t$internal{'short'}";
	
	# [11] LONG
	$output{$locus} .= "\t$internal{'long'}";
	
	# [12] through whenenver .. Dicer 
	for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
	    $output{$locus} .= "\t$internal{$i}";
	}
    }

    log_it($logfile," Done\n");
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
    my @hp_hash_fields = ();
    
    while(($locus,$entry) = each %$quant_hash) {
	@fields = split ("\t", $entry);
	@new_fields = @fields[3..9];
	if(exists($$hp_hash{$locus})) {
	    # polarity of the hairpin supersedes that calcuated by the reads alone
	    if($$hp_hash{$locus} =~ /Watson/) {
		unshift(@new_fields,"+");
	    } elsif ($$hp_hash{$locus} =~ /Crick/) {
		unshift(@new_fields,"-");
	    } else {
		log_it($logfile, "\nFATAL error in sub-routine mod_quant_hp : could find strand for HP locus $locus in entry $$hp_hash{$locus}\n");
		exit;
	    }
	    
	    # now check if its a MIRNA or not
	    if($$hp_hash{$locus} =~ /MIRNA/) {
		unshift(@new_fields,"MIRNA");
	    } else {
		unshift(@new_fields,"HP");
	    }

	    ## add all of the hairpin information
	    @hp_hash_fields = split ("\t", $$hp_hash{$locus});
	    push(@new_fields, $hp_hash_fields[4]);  ## Pairs
	    push(@new_fields, sprintf("%.4f",$hp_hash_fields[5]));  ## frac_paired, rounded to 4 decimal places
	    push(@new_fields, $hp_hash_fields[6]);  ## stem_length
	    push(@new_fields, $hp_hash_fields[7]);  ## loop_length
	    push(@new_fields, sprintf("%.4f",$hp_hash_fields[8]));  ## dGperStem, rounded to 4 decimal places
	    push(@new_fields, sprintf("%.4f",$hp_hash_fields[9]));  ## frac_hp, rounded to 4 decimal places
	    push(@new_fields, $hp_hash_fields[10]);  ## hp_size_result
	    push(@new_fields, $hp_hash_fields[11]);  ## precision_result
	    push(@new_fields, $hp_hash_fields[12]);  ## duplex_result
	    push(@new_fields, $hp_hash_fields[13]);  ## star_result

	    
	} else {
	    # not an HP, leave the original strand alone, and enter "." for the HP designation
	    unshift(@new_fields,$fields[2]);
	    unshift(@new_fields,"\.");
	    
	    # and enter NAs for all of the other HP information
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    push(@new_fields, "NA");
	    
	}
	
	# add the name and the locus
	unshift(@new_fields,$fields[1]);
	unshift(@new_fields,$locus);
	
	# add the rest of the original information (the read counts)
	for(my $ii = 10; $ii < (scalar @fields); ++$ii) {
	    push(@new_fields, $fields[$ii]);
	}
	
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
	# hp information
	push(@new_fields,$old_fields[11]);
	push(@new_fields,$old_fields[12]);
	push(@new_fields,$old_fields[13]);
	push(@new_fields,$old_fields[14]);
	push(@new_fields,$old_fields[15]);
	push(@new_fields,$old_fields[16]);
	push(@new_fields,$old_fields[17]);
	push(@new_fields,$old_fields[18]);
	push(@new_fields,$old_fields[19]);
	push(@new_fields,$old_fields[20]);
	
	# short and long
	push(@new_fields, sprintf("%.3f",(($old_fields[21] / $$reads) * 1000000)));
	push(@new_fields, sprintf("%.3f",(($old_fields[22] / $$reads) * 1000000)));
	# the rest
	for($i = 23; $i < (scalar @old_fields); ++$i) {
	    push(@new_fields, sprintf("%.3f",(($old_fields[$i] / $$reads) * 1000000)));
	}
	my $result = join("\t",@new_fields);
	$output{$locus} = $result;
	
    }
    return %output;
}

sub write_gff3s {
    my($clus_array,$quant_hash,$outdir) = @_; ## passed by reference

    # open files
    my $dclfile = "$$outdir" . "\/" . "$$outdir" . "_DCL.gff3";
    my $nfile = "$$outdir" . "\/" . "$$outdir" . "_N.gff3";
    open(DCL, ">$dclfile");
    open(N, ">$nfile");

    my $i;
    my $j;

    my @fcfields = ();
    my $chrom;
    my $start;
    my $stop;
    my $name;
    foreach my $f_clus (@$clus_array) {
	@fcfields = split("\t",$$quant_hash{$f_clus});
	if($fcfields[0] =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chrom = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    log_it($logfile, "\nFATAL: in sub-routine \"write_gff3s\" : Could not parse locus name $fcfields[0]\n");
	    exit;
	}
	
	my $out = "$chrom\t" . "ShortStack\t";
	if($fcfields[9] eq "N") {
	    $out .= "NOT_DCL\t";
	} else {
	    $out .= "DCL\t";
	}
	$out .= "$start\t$stop\t\.\t";
	$out .= "$fcfields[4]\t"; ## strand
	$out .= "\.\t";
	$out .= "ID=$fcfields[1]\;";
	$out .= "Name=$fcfields[1]\;";
	$out .= "DicerCall=$fcfields[9]\;";
	$out .= "HP=";
	if($fcfields[3] eq ".") {
	    $out .= "NONE\;";
	} else {
	    $out .= "$fcfields[3]\;";
	}
	$out .= "FracWatson=$fcfields[5]\;";
	$out .= "UniquenessIndex=$fcfields[8]\n";
	
	if($fcfields[9] eq "N") {
	    print N "$out";
	} else {
	    print DCL "$out";
	}
    }
    close N;
    close DCL;
    return($dclfile,$nfile);
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

    $text_to_return .= "\n\tColor-Scheme in bed file $bedfile : \n";
    $colors{'N'} = "169,169,169";  ## dark gray

    $text_to_return .= "\tNon-Dicer Clusters: Dark Gray RGB: 169,169,169\n";

    my $n_colors_to_get = $$dicermax - $$dicermin + 1;
    my %rgb_hash = get_colors($n_colors_to_get);
    my $x = 0;
    
    for ($i = $$dicermin; $i <= $$dicermax; ++$i) {
	++$x;
	$colors{$i} = $rgb_hash{$x};

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
	    log_it($logfile, "\nFATAL: in sub-routine \"write_bed_nohp\" : Could not parse locus name $fcfields[0]\n");
	    exit;
	}
	print BED "$chrom\t$bedstart\t$bedstop\t";
	
	#name
	print BED "$$names{$f_clus}";
	# Check $fcfields[3]  -- if it is 'HP' or 'MIRNA', then ... as of version 0.3.0, array element 3 is the HP call (used to be 2)
	if($fcfields[3] eq "HP") {
	    print BED "_HP\t";
	} elsif ($fcfields[3] eq "MIRNA") {
	    print BED "_MIRNA\t";
	} else {
	    print BED "\t";
	}
	
	# score (ignored -- just zero)
	print BED "0\t";
	
	# strand, which is $fcfields[4]  .. array element 4 as of version 0.3.0
	print BED "$fcfields[4]\t";
	
	# thickStart and thickEnd, which are always the same as bedstart and bedstop
	print BED "$bedstart\t$bedstop\t";
	
	# Determine the dominant size of the cluster from $fcfields[9] and write correct color .. array element 9 as of version 0.3.0
	if($fcfields[9] eq "N") {
	    print BED "$colors{'N'}\t";
	} elsif ($fcfields[9] =~ /^(\d+)$/) {
	    print BED "$colors{$1}\t";
	} else {
	    log_it($logfile, "\nFATAL in sub-routine \"write_bed_nohp\" : Could not find the cluster size from entry $fcfields[8]\n");
	    exit;
	}
	
	if(($fcfields[3] eq "HP") or
	   ($fcfields[3] eq "MIRNA")) {
	    
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
		## added in 0.2.3 ... adjust the left and right begin coordinates to zero-based
		--$left_begin;
		--$right_begin;
		print BED "$left_block_size";
		print BED ",";
		print BED "$right_block_size\t";
		print BED "$left_begin";
		print BED ",";
		print BED "$right_begin\n";
	    } else {
		log_it($logfile, "\nFATAL: in sub-routine write_bed_with_hp : could not parse hairpin information with regex from $$hp_hash{$f_clus}\n");
		exit;
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
    my ($clus)  = @_;
    my %hash = ();
    my @fields = ();
    my $size;
    foreach my $line (keys %$clus) {
	@fields = split ("\t", $$clus{$line});
	
	# dicer call
	if($fields[9] eq "N") {  ## array element 9 as of version 0.3.0
	    ++$hash{$fields[9]};
	} else {
	    if($fields[9] =~ /^(\d+)$/) {
		$size = $1;
	    } else {
		log_it($logfile, "\nFATAL in sub-routine final_summary_hp : failed to parse dicer size category of $fields[9]\n");
		exit;
	    }
	    ++$hash{$size};
	}
	
	# phased?
	if($fields[12] =~ /OK/) { ## array element 12 as of version 0.3.0
	    ++$hash{'phased'};
	}
    }
    return %hash;
}

sub final_summary {
    my ($clus)  = @_;
    my %hash = ();
    my @fields = ();
    my $hp_call;
    my $dcall;
    foreach my $line (keys %$clus) {
	@fields = split ("\t", $$clus{$line});
	
	# dicer call
	if($fields[9] eq "N") { ## array element 9 as of version 0.3.0
	    $dcall = $fields[9];
	} else {
	    if($fields[9] =~ /^(\d+)$/) {
		$dcall = $1;
	    } else {
		log_it($logfile, "\nFATAL in sub-routine final_summary_hp : failed to parse dicer size category of $fields[9]\n");
		exit;
	    }
	}
	
	# HP call
	if(($fields[3] eq "MIRNA") or  ## array element 3 as of version 0.3.0
	   ($fields[3] eq "HP")) {
	    $hp_call = "$fields[3]";
	} else {
	    $hp_call = "X";
	}

	++$hash{$dcall}{$hp_call};
	
	# phased?
	if($fields[12] =~ /OK/) {  ## array element 12 as of version 0.3.0
	    ++$hash{'phased'};
	}
    }
    return %hash;
}

sub parse_cigar {
    my($cigar) = @_;
    if($cigar eq "\*") {
	log_it($logfile,"FATAL in sub-routine \'parse_cigar\' : The CIGAR string is missing from at least one of the mappings, so I can't caluclate read length\!\n");
	exit;
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
	    if($fields[14] eq "MIRNA") {
		++$n_miRNAs;
		open(OUT, ">$$outdir/MIRNA_details/$name\.txt");
		print OUT "$name $fields[15]\n";
		close OUT;
		
		@file_lines = split ("\n", $fields[15]);
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
		
	    } elsif ($fields[14] eq "HP") {
		++$n_hpRNAs;
		open(OUT, ">$$outdir/HP_details/$name\.txt");
		print OUT "$name $fields[15]\n";
		close OUT;
	    } else {
		log_it($logfile, "\nFATAL in sub-routine write_files: failed to understand following entry as MIRNA or HP  found $fields[14] instead:\n$$hp_clusters{$locus}\n");
		exit;
	    }
	}
    }
    
    # close the tables
    close MIR_TABLE;
    
    # delete everything if there were no miRNAs and hpRNAs
    if(($n_miRNAs + $n_hpRNAs) == 0) {
	system "rm -f -r $$outdir/MIRNA_details";
	system "rm -f -r $$outdir/HP_details";
    }
    
    return($n_miRNAs,$n_hpRNAs);
}

sub parse_inv {
    my($inv_file,$minntspaired,$minfracpaired,$maxdGperStem,$maxLoopLength) = @_;  ## passed by reference, all scalars
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
    
    my $stem_length;
    my $loop_length;
    my $dG;
    my $dGperStem;
    
    my $left_brax;
    my $right_brax;
    my $left_seq;
    my $right_seq;
    
    my $eval_string;
    
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
		    log_it($logfile,"\n\tWARNING in sub-routine parse_inv: Failed to parse chr name from line $single[0] and skipped the entry\n\n");
		    @single = ();
		    next;
		}
		if($single[1] =~ /(\d+) (\S+) (\d+)/) {
		    $left_start = $1;
		    $left_string = $2;
		    $left_stop = $3;
		} else {
		    log_it($logfile,"\n\tWARNING in sub-routine parse_inv: Failed to parse expected left arm info from line $single[1] and skipped the entry\n\n");
		    @single = ();
		    next;
		}
		if($single[3] =~ /(\d+) (\S+) (\d+)/) {
		    $right_stop = $1;
		    $right_string = $2;  ## top strand, reading 3' to 5' i.e. descending
		    $right_start = $3;
		}
		
		# first check that loop criteria are met for this hairpin.  If not, save time by ceasing any further analysis
		$stem_length = 0.5 * (($left_stop - $left_start + 1) + ($right_stop - $right_start + 1));
		$loop_length = $right_start - $left_stop - 1;
		
		unless(($loop_length <= $$maxLoopLength) and
		       ($loop_length < $stem_length)) {
		    @single = ();
		    next;
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
		    log_it($logfile,"\n\tWARNING in sub-routine parse_inv: Possible einverted bug -- On top strand of entry below there are $top_bases nucleotides but the coordinates $chr $left_start to $left_stop do not concur\n\t\tSkipping entry $chr $left_start to $right_stop\n\n");
		    @single = ();
		    next;
		}
		if(($bottom_bases != ($right_stop - $right_start + 1))) {
		    log_it($logfile,"\n\tWARNING in sub-routine parse_inv: Possible einverted bug -- On bottom strand of entry below there are $bottom_bases nucleotides but the coordinates $chr $right_start to $right_stop do not concur\n\t\tSkipping entry $chr $left_start to $right_stop\n\n");
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
		
		$w_frac_paired = $w_pairs / $stem_length;
		if(($w_frac_paired >= $$minfracpaired) and
		   ($w_pairs >= $$minntspaired)) {
		    
		    ## evaluate the stem dG using RNAeval
		    $left_brax = substr($w_brax,0,($left_stop - $left_start + 1));
		    $right_brax = substr($w_brax,($right_start - $left_start),($right_stop - $right_start + 1));
		    ## left seq for w strand is left_string
		    $left_seq = $left_string;
		    $left_seq =~ s/-//g;  ## remove gap symbols!
		    $right_seq = reverse $right_string;  ##  The right string is 3'-5'
		    $right_seq =~ s/-//g;  
		    
		    $eval_string = "$left_seq" . "\&" . "$right_seq" . "\n" . "$left_brax" . "\&" . "$right_brax";
		    
		    $dG = dG_from_RNAeval($eval_string);
		    
		    unless($dG eq "FAIL") {
			$dGperStem = $dG / $stem_length;
			if($dGperStem <= $$maxdGperStem) {
			    
			    $outline = "$w_brax\t$left_start" . "-" . "$right_stop\t$left_start" . "-" . "$left_stop" . "," . "$right_start" . "-" . "$right_stop\t$strand";
			    $outline .= "\t$w_pairs\t$w_frac_paired\t$stem_length\t$loop_length\t$dGperStem\t$chr";
			    
			    push(@output, $outline);
			}
		    }
		}
		# now Crick
		$outline = '';
		$strand = "Crick";
		
		
		
		$c_frac_paired = $c_pairs / $stem_length;
		if(($c_frac_paired >= $$minfracpaired) and
		   ($c_pairs >= $$minntspaired)) {
		    
		    ## evaluate the stem dG using RNAeval
		    $left_brax = substr($c_brax,0,($right_stop - $right_start + 1));
		    $right_brax = substr($c_brax,($right_stop - $left_stop),($left_stop - $left_start + 1));
		    # left_seq is the complement of the right_string
		    $left_seq = $right_string;
		    $left_seq =~ tr/aucg/uagc/;
		    $left_seq =~ s/-//g;
		    # right_seq is the reverse complement of the left_string
		    $right_seq = reverse $left_string;
		    $right_seq =~ tr/aucg/uagc/;
		    $right_seq =~ s/-//g;
		    
		    $eval_string = "$left_seq" . "\&" . "$right_seq" . "\n" . "$left_brax" . "\&" . "$right_brax";

		    $dG = dG_from_RNAeval($eval_string);
		    unless($dG eq "FAIL") {
			$dGperStem = $dG / $stem_length;
			if($dGperStem <= $$maxdGperStem) {
			    $outline = "$c_brax\t$right_stop" . "-" . "$left_start\t$right_stop" . "-" . "$right_start" . "," . "$left_stop" . "-" . "$left_start\t$strand";
			    $outline .= "\t$c_pairs\t$c_frac_paired\t$stem_length\t$loop_length\t$dGperStem\t$chr";

			    push(@output, $outline);
			}
		    }
		}
		@single = ();
	    }
	    #@single = ();
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
	$actual_one =~ tr/augc/uacg/;
	$actual_two =~ tr/uacg/augc/;
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
	    if($ir_fields[-1] eq $c_chr) {
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
    
sub flag_overlap {
    my($in_hash,$file) = @_; ## passed by reference .. hash and scalar
    my %out_hash = ();
    my @ss_loci = keys %$in_hash;
    my $s_locus;
    my $s_chr;
    my $s_start;
    my $s_stop;
    my @flag_data = ();
    my $f_string;
    my @f_fields = ();
    my $f_locus;
    my $f_chr;
    my $f_start;
    my $f_stop;
    my $overlap;
    
    my $n_s_overlapped = 0;
    my $n_f_overlapped = 0;
 
    my %f_over = ();

    my @hash_fields = ();
    my @new_hash_fields = ();
    my $hash_value;
    my $x;
    my $new_string;
    
    if(-r $$file) {
	open(FLAG,"$$file");
	@flag_data = <FLAG>;
	close FLAG;
	foreach $s_locus (@ss_loci) {
	    if($s_locus =~ /^(\S+):(\d+)-(\d+)/) {
		$s_chr = $1;
		$s_start = $2;
		$s_stop = $3;
	    } else {
		log_it($logfile, "\nFATAL in sub-routine flag_overlap: failed to parse small RNA locus name $s_locus\n");
		exit;
	    }
	    $overlap = '';
	    foreach $f_string (@flag_data) {
		chomp $f_string;
		@f_fields = split ("\t", $f_string);
		if($f_fields[0] =~ /^(\S+):(\d+)-(\d+)/) {
		    $f_chr = $1;
		    $f_start = $2;
		    $f_stop = $3;
		} else {
		    log_it($logfile, "\nFATAL in sub-routine flag_overlap: failed to parse flag locus name $f_fields[0]\n");
		    exit;
		}
		if($f_chr ne $s_chr) {
		    next;
		} elsif (($s_start > $f_stop) or ($f_start > $s_stop)) {
		    next;
		} else {
		    $f_over{$f_fields[0]} = 1;
		    if($overlap) {
			$overlap .= ",$f_fields[1]";
		    } else {
			$overlap = "$f_fields[1]";
		    }
		}
	    }
	    @hash_fields = split ("\t", $$in_hash{$s_locus});
	    @new_hash_fields = ();
	    $x = 0;
	    foreach $hash_value (@hash_fields) {
		++$x;
		push(@new_hash_fields, $hash_value);
		if($x == 2) {
		    if($overlap) {
			++$n_s_overlapped;
			push(@new_hash_fields, $overlap);
		    } else {
			push(@new_hash_fields, "\.");
		    }
		}
	    }
	    $new_string = join("\t", @new_hash_fields);
	    $out_hash{$s_locus} = $new_string;
	}
	$n_f_overlapped = scalar (keys %f_over);
	
	# Report to user
	log_it($logfile,"\t\tA total of $n_s_overlapped small RNA loci overlapped with one or more flag_file locus from file $$file\n");
	log_it($logfile,"\t\tA total of $n_f_overlapped loci from flag file $$file overlapped with one or more small RNA locus\n");
    } else {
	# just add a dot to the output hash for all of them
	foreach $s_locus (@ss_loci) {
	    $overlap = "\.";
	    @hash_fields = split ("\t", $$in_hash{$s_locus});
	    @new_hash_fields = ();
	    $x = 0;
	    foreach $hash_value (@hash_fields) {
		++$x;
		push(@new_hash_fields, $hash_value);
		if($x == 2) {
		    push(@new_hash_fields, $overlap);
		}
	    }
	    $new_string = join("\t", @new_hash_fields);
	    $out_hash{$s_locus} = $new_string;
	}
    }
    return %out_hash;
}

sub final_hp { ## First added 0.3.0
    my($in_hash,$genome,$bamfile,$minstrandfrac,$maxmiRHPPairs,$maxmiRUnpaired,$read_group) = @_;  ## by reference .. first one hash, all others scalar
    
    my $orig_locus;
    my $hp_locus;
    my $hp_entry;
    my @hp_fields = ();
    my @hp_arms = ();
    my @fivep = ();
    my @threep = ();
    my $left_start;
    my $left_stop;
    my $right_start;
    my $right_stop;
    my $strand;
    my $chr;
    
    my @alignments = ();
    my @samfields = ();
    my $readstrand;
    
    my $total_reads;
    my $right_strand_reads;
    
    my $read_length;
    my $read_stop;
    
    my %left_dyad = ();
    my %right_dyad = ();
    my $readname;
    
    my $strandfrac;
    
    my %output_hash = ();
    my $result_string = ();
    
    
    
    while(($orig_locus, $hp_entry) = each %$in_hash) {
	@hp_fields = split ("\t", $hp_entry);
	@hp_arms = split (",", $hp_fields[2]);
	@fivep = split ("-", $hp_arms[0]);
	@threep = split ("-", $hp_arms[1]);
	
	## get chr
	if($orig_locus =~ /^(\S+):/) {
	    $chr = $1;
	} else {
	    log_it($logfile, "\nFATAL: Could not parse chromosome name from entry $orig_locus in sub-routine final_hp\n");
	    exit;
	}
	
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
	
	# query for samtools view
	$hp_locus = "$chr" . ":" . "$left_start" . "-" . "$right_stop";
	
	# store alignments in an array while you go through them the first time
	# If it passes the first tests, the array is used again during the miRNA annotation process
	
	@alignments = ();  ## reset each time through
	$total_reads = 0;  ## reset each time through
	$right_strand_reads = 0; ## reset each time through
	%left_dyad = ();
	%right_dyad = ();
	if($$read_group) {
	    open(SAM, "samtools view -F 0x4 -r $$read_group $$bamfile $hp_locus \|");
	} else {
	    open(SAM, "samtools view -F 0x4 $$bamfile $hp_locus \|");
	}
	while (<SAM>) {
	    chomp;
	    push(@alignments, $_);  ## store for later use
	    
	    ## ignore headers
	    if($_ =~ /^@/) {
		next;
	    }
	    
	    @samfields = split ("\t", $_);
	    
            ## ignore unmapped reads
	    if($samfields[1] & 4) {
		next;
	    }
	    
	    ## increment total reads
	    ++$total_reads;
	    
	    ## determine strand of read
	    if($samfields[1] & 16) {
		$readstrand = "-";
	    } else {
		$readstrand = "+";
	    }
	    
	    ## increment if it is the right strand (will amend later for dyads)
	    if($readstrand eq $strand) {
		++$right_strand_reads;
	    }
	    
	    ## determine 1-based stop position of this read
	    $read_length = parse_cigar($samfields[5]);
	    $read_stop = $samfields[3] + $read_length - 1;
	    
	    ## if overlaps with left arm, record in left dyad hash
	    if(($samfields[3] >= $left_start) and
	       ($read_stop <= $left_stop)) {
		$left_dyad{$samfields[0]} = $readstrand;
	    }
	    
	    ## if overlaps with right arm, record in right dyad hash
	    if(($samfields[3] >= $right_start) and
	       ($read_stop <= $right_stop)) {
		$right_dyad{$samfields[0]} = $readstrand;
	    }
	}
	close SAM;
	
	## adjust for any dyads.  For each dyad, remove one count from the **total** but NOT from the right_strand_reads
	## in other words, dyad mappings to the opposite strand are treated as if they didn't exist, but those to the same strand are counted as right
	
	while(($readname) = each %left_dyad) {
	    if(exists($right_dyad{$readname})) {
		if($left_dyad{$readname} ne $right_dyad{$readname}) { ## e.g., on opposite strands
		    --$total_reads;
		}
	    }
	}
	
	## calculate the strand frac
	## protect against undefined errors
	if($total_reads > 0) {
	    $strandfrac = $right_strand_reads / $total_reads;
	} else {
	    $strandfrac = 0;
	}
	
	## test if the strandfrac meets the minimum
	## If so, this locus is either an hpRNA or miRNA
	## If not, this is not a hairpin-associated locus, no entry to the output hash, and the original cluster coordinates will be retained instead
	
	if($strandfrac >= $$minstrandfrac) {
	    
	    ## miRNA vs. hpRNA analysis
	    $result_string = analyze_for_mir(\$hp_entry, \$$genome, \$$maxmiRHPPairs, \$$maxmiRUnpaired, \@alignments, \$strand, \$hp_locus);
	    
	    ## add to output hash
	    $output_hash{$orig_locus} = $result_string;
	}
    }
    
    ## return the output hash
    return %output_hash;
}

sub analyze_for_mir {
    my($hp_entry,$genome,$maxmiRHPPairs,$maxmiRUnpaired,$alignments,$strand,$hp_locus) = @_; ## by reference, all scalars except the alignments, which is an array
    
    ## parse the hp_entry
    my @l_data_fields = split ("\t", $$hp_entry);
    my $brax = $l_data_fields[0];
    
    # get the total number of paired nts in this structure
    my $n_of_pairs = 0;  ## reset each time
    while ($brax =~ /\(/g) {
	++$n_of_pairs;
    }
    
    my @brax_coords = split ("-",$l_data_fields[1]);
    # sort so that the lower coordinate is always in the [0]
    @brax_coords = sort {$a <=> $b} @brax_coords;
    
    ## retrieve the genomic sequence of interest
    my $subseq;
    open(FAIDX, "samtools faidx $$genome $$hp_locus |");
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
    my $displayseq;
    if($$strand eq "-") {
	$displayseq = reverse $subseq;
	$displayseq =~ tr/ACUG/UGAC/;
    } else {
	$displayseq = $subseq;
    }
    
    # code the hp sequence and the brackets into hashes, keyed by coordinates, for easy retrieval
    my @hp_letters = split ('', $displayseq);
    my %hp_seq_hash = ();
    my $j = 0;
    my $i;
    my $loc_start;
    my $loc_stop;
    if($$hp_locus =~ /^\S+:(\d+)-(\d+)$/) {
	$loc_start = $1;
	$loc_stop = $2;
    } else {
	log_it($logfile, "\nFATAL in sub-routine analyze_for_mir : failed to parse locus name $$hp_locus\n");
	exit;
    }
    if($$strand eq "+") {
	for($i = $loc_start; $i <= $loc_stop; ++$i) {
	    $hp_seq_hash{$i} = $hp_letters[$j];
	    ++$j;
	}
    } elsif ($$strand eq "-") {
	for($i = $loc_stop; $i >= $loc_start; --$i) {
	    $hp_seq_hash{$i} = $hp_letters[$j];
	    ++$j;
	}
    }
    
    my @brax_chars = split ('', $brax);
    my %brax_hash = ();
    $j = 0;
    if($$strand eq "+") {
	for($i = $brax_coords[0]; $i <= $brax_coords[1]; ++$i) {
	    $brax_hash{$i} = $brax_chars[$j];
	    ++$j;
	}
    } elsif ($$strand eq "-") {
	for($i = $brax_coords[1]; $i >= $brax_coords[0]; --$i) {
	    $brax_hash{$i} = $brax_chars[$j];
	    ++$j;
	}
    }
    
    # get all mappings that have any overlap with the interval, and track in a data structure
    #  that keeps 'em sorted by left-most position and distinguishes sense from antisense
    
    my %sense = ();  ## NOTE, sense and antisense are relative the hairpin direction, not (necessarily) the genome
    my %antisense = ();
    my @sense_names = ();
    my @antisense_names = ();
    my $total_mappings = 0;
    my %miRNA = (); 
    my %miRNA_star = ();
    
    my @samfields = ();
    my $read_length;
    my $read_stop;

    my $readstrand;
    my $mapping_polarity;
    my $map_coord;
    my $key;
    foreach my $samline (@$alignments) {
	chomp $samline;
	# ignore headers
	if($samline =~ /^@/) {
	    next;
	}
	@samfields = split ("\t", $samline);
	
	#ignore unmapped reads
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
	
	# get readstrand (relative to genome)
	if($samfields[1] & 16) {
	    $readstrand = "-";
	} else {
	    $readstrand = "+";
	}
	
        # decide whether the read is sense or antisense wrt to the hairpin direction
	if($readstrand eq "-") {
	    # mapping is antisense relative to genome
	    if($$strand eq "+") {
		# mapping is antisense relative to hairpin
		$mapping_polarity = "antisense";
	    } elsif ($$strand eq "-") {
		# mapping is sense relative to hairpin
		$mapping_polarity = "sense";
	    }
	} else {
	    # mapping is sense relative to genome
	    if($$strand eq "+") {
		# mapping is sense relative to hairpin
		$mapping_polarity = "sense";
	    } elsif($$strand eq "-") {
		$mapping_polarity = "antisense";
	    }
	}
    
	
	# the naming convention depends solely upon whether the hairpin is Watson or Crick.
	# All mappings for a Watson hairpin are start-stop, while all mappings for a Crick hairpin are
	# keyed as stop-start.  Doesn't matter whether the mapping itself is sense or antisense wrt to the hairpin
	
	if($$strand eq "+") {
	    $map_coord = "$samfields[3]" . "-" . "$read_stop";
	    $key = $samfields[3];
	} elsif ($$strand eq "-") {
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
    
    # Begin assessment of whether the hairpin can be annotated as a miRNA
    my @candidates = ();
    # First, gather all sense mappings that account for 20% or more of ALL mappings at the locus
    foreach $map_coord (@sense_names) {
	$key = $map_coord;
	$key =~ s/-\d+$//g;
	if(($sense{$key}{$map_coord} / $total_mappings) > 0.2) {
	    push (@candidates,$map_coord);
	}
    }
    
    my $prec_result = 0;
    my $hp_size_result = 0;
    my $duplex_result = 0;
    my $star_result = 0;
    if(@candidates) {
	++$prec_result;
    }
    
    unless($n_of_pairs > $$maxmiRHPPairs) {
	++$hp_size_result;
    }
    
    ## analyzing duplex structure
    
    my $candidate_brax;
    my $n_unpaired_cand;
    my $n_left_paired_cand;
    my $n_right_paired_cand;
    my @mc_limits;
    
    my $candidate_brax_chopped;
    
    my $star_coord;
    
    my @star_limits = ();
    my $star_brax;
    my $n_unpaired_star;
    my $n_left_paired_star;
    my $n_right_paired_star;
    
    my $star_brax_chopped;
    
    my $cand_fail;
    
    if(@candidates) {
	foreach $map_coord (@candidates) {
	    $cand_fail = 0;
	    $candidate_brax = '';
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
	    if($n_unpaired_cand > $$maxmiRUnpaired) {
		$cand_fail = 1;
	    }
	    
	    if(($n_left_paired_cand > 0 ) and
	       ($n_right_paired_cand > 0)) {
		$cand_fail = 1;
	    }
	    
	    # find the map_coordinates for the star sequence
	    $star_coord = get_star_coord(\$map_coord,\%brax_hash,\$candidate_brax);
	    if($star_coord eq "fail") {
		$cand_fail = 1;
	    } else {
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
		if($n_unpaired_star > $$maxmiRUnpaired) {
		    $cand_fail = 1;
		}
		if(($n_left_paired_star > 0 ) and
		   ($n_right_paired_star > 0)) {
		    $cand_fail = 1;
		}
		
		## if $cand_fail is still 0 at this point, then this candidate passes the duplex test
		unless($cand_fail) {
		    ++$duplex_result;
		}
	    }
	    
	    ## unless the candidate already failed, check for existence of star sequence
	    unless($cand_fail) {
		if(exists($sense{$star_limits[0]}{$star_coord})) {
		    if( (($sense{$star_limits[0]}{$star_coord} + $sense{$mc_limits[0]}{$map_coord}) / $total_mappings ) <= 0.25) {
			$cand_fail = 1;
		    }
		} else {
		    $cand_fail = 1;
		}
	    }
	    
	    ## unless the candidate already failed, check for redundancy
	    unless($cand_fail) {
		if((exists($miRNA{$map_coord})) or
		   (exists($miRNA{$star_coord})) or
		   (exists($miRNA_star{$map_coord})) or
		   (exists($miRNA_star{$star_coord}))) {
		    $cand_fail = 1;;
		} 
	    }
	    unless($cand_fail) {
		# you've got a winner.  Check abundances to decide which one should really be the miRNA, and which should be the star
		++$star_result;
		if($sense{$mc_limits[0]}{$map_coord} >= $sense{$star_limits[0]}{$star_coord}) {
		    $miRNA{$map_coord} = $star_coord;
		    $miRNA_star{$star_coord} = $map_coord;
		} else {
		    $miRNA{$star_coord} = $map_coord;
		    $miRNA_star{$map_coord} = $star_coord;
		}
	    }
	}
    }
    
    ## 
    my $hp_call;
    if($star_result >= 1) {
	$hp_call = "MIRNA";
    } else {
	$hp_call = "HP";
    }
    
    # begin output
    my $output_string = '';
    my $mmmr;
    $output_string .= "$$hp_locus $$strand\n";
    $output_string .= "$displayseq\n";
    if($$strand eq "+") {
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
			    log_it($logfile, "\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
			    exit;
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
			    log_it($logfile, "\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
			    exit;
			}
			for($j = ($2 + 1); $j <= $loc_stop; ++$j) {
			    $output_string .= "<";
			}
			
			$output_string .= " ";
			# print length and abundance
			$output_string .= "l=$read_length ";
			$output_string .= "m=$antisense{$i}{$map_coord}";
			$output_string .= "\n";
			
		    }
		}
	    }
	}
    } elsif ($$strand eq "-") {
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
			    log_it($logfile, "\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
			    exit;
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
			    log_it($logfile, "\nFATAL: error in sub-routine hp_outout : failed to parse map_coord $map_coord\n");
			    exit;
			}
			for($j = ($2 - 1); $j >= $loc_start; --$j) {
			    $output_string .= "<";
			}
			
			$output_string .= " ";
			# print length and abundance
			$output_string .= "l=$read_length ";
			$output_string .= "m=$antisense{$i}{$map_coord}";
			$output_string .= "\n";
		    }
		}
	    }
	}
    }
    
    ## append fields to the input information
    my $result = $$hp_entry;
    ## hp_length, precision, duplex, and star
    $result .= "\t$hp_size_result\t$prec_result\t$duplex_result\t$star_result";
    
    ## add the hp_call, either "MIRNA" or "HP"
    $result .= "\t$hp_call";
    
    # check the output string, and suppress it if it is too big (> 400nts long or > 400 distinct small RNAs)
    # make an exception for MIRNAs, for which we always print all of the alignments.
    # this should reduce the memory footprint
    
    my @olines = split ("\n", $output_string);
	    
    unless($hp_call eq "MIRNA") {
	if((($loc_stop - $loc_start + 1) > 400) or
	   ((scalar @olines) > 403)) {  ## 3 lines are header, brackets, and sequence
	    $output_string = "$olines[0]\n$olines[1]\n$olines[2]\n xxx Alignments Not Shown - too complex xxx\n";
	}
    }
    
    ## add the alignment
    $result .= "\t$output_string";
    
    ## phew
    return $result;
}

sub dG_from_RNAeval {
    my($eval_string) = @_;
    my $dG;
    
    open(RNAEVAL, "echo \"$eval_string\" | RNAeval |");
    my @eval_line = <RNAEVAL>;
    close RNAEVAL;
    
    if($eval_line[-1] =~ /\s+\((.*\d+.*)\)/) {
	$dG = $1;
	$dG =~ s/\s//g;
    } else {
	$dG = "FAIL";
    }
    return $dG;
}

sub summarize {
    my($quant,$nohp,$dicermin,$dicermax,$logfile) = @_;  ## hash, scalar
    my %loci = ();
    my %abun = ();
    my $i;
    
    if($$nohp) {
	$abun{'ND'}{'N'} = 0;
	$loci{'ND'}{'N'} = 0;
	for($i = $$dicermin; $i <= $$dicermax; ++$i) {
	    $abun{'ND'}{$i} = 0;
	    $loci{'ND'}{$i} = 0;
	}
    } else {
	$abun{'MIRNA'}{'N'} = 0;
	$abun{'HP'}{'N'} = 0;
	$abun{'.'}{'N'} = 0;
	$loci{'MIRNA'}{'N'} = 0;
	$loci{'HP'}{'N'} = 0;
	$loci{'.'}{'N'} = 0;
	for($i = $$dicermin; $i <= $$dicermax; ++$i) {
	    $abun{'MIRNA'}{$i} = 0;
	    $abun{'HP'}{$i} = 0;
	    $abun{'.'}{$i} = 0;
	    $loci{'MIRNA'}{$i} = 0;
	    $loci{'HP'}{$i} = 0;
	    $loci{'.'}{$i} = 0;
	}
    }
    
    my @fields = ();
    my $locus;
    my $entry;
    my $rounded;
    
    while(($locus, $entry) = each %$quant) {
	@fields = split ("\t", $entry);
	
	
	# $fields[9] is the Dicer Call
	# $fields[6] is the total number of alignments
	# $fields[3] is the HP call
	++$loci{$fields[3]}{$fields[9]};
	$abun{$fields[3]}{$fields[9]} += $fields[6];
    }
    
    # report
    log_it($$logfile,"\nSummary\n\n");

    if($$nohp) {
	log_it($$logfile,"DicerCall\tLoci\tAlignments\n");
	log_it($$logfile, "N\t$loci{'ND'}{'N'}\t");
	$rounded = sprintf("%.1f",$abun{'ND'}{'N'});
	log_it($$logfile,"$rounded\n");
	for($i = $$dicermin; $i <= $$dicermax; ++$i) {
	    log_it($$logfile,"$i\t$loci{'ND'}{$i}\t");
	    $rounded = sprintf ("%.3f",$abun{'ND'}{$i});
	    log_it($$logfile,"$rounded\n");
	}
    } else {
	log_it($$logfile,"DicerCall\tNON-HP_Loci\tHP_Loci\tMIRNA_Loci\tNON-HP_Alignments\tHP_Alignments\tMIRNA_Alignments\n");
	log_it($$logfile,"N\t$loci{'.'}{'N'}\t$loci{'HP'}{'N'}\t$loci{'MIRNA'}{'N'}\t");
	$rounded = sprintf("%.1f",$abun{'.'}{'N'});
	log_it($$logfile,"$rounded\t");
	$rounded = sprintf("%.1f",$abun{'HP'}{'N'});
	log_it($$logfile,"$rounded\t");
	$rounded = sprintf("%.1f",$abun{'MIRNA'}{'N'});
	log_it($$logfile,"$rounded\n");
	for($i = $$dicermin; $i <= $$dicermax; ++$i) {
	    log_it($$logfile,"$i\t$loci{'.'}{$i}\t$loci{'HP'}{$i}\t$loci{'MIRNA'}{$i}\t");
	    $rounded = sprintf("%.1f",$abun{'.'}{$i});
	    log_it($$logfile,"$rounded\t");
	    $rounded = sprintf("%.1f",$abun{'HP'}{$i});
	    log_it($$logfile,"$rounded\t");
	    $rounded = sprintf("%.1f",$abun{'MIRNA'}{$i});
	    log_it($$logfile,"$rounded\n");
	}
    }
    log_it($$logfile,"\n");

}

sub log_it {
    my($file,$string) = @_;
    open(LOG, ">>$file");
    print LOG "$string";
    print STDERR "$string";
    close LOG;
}

sub get_folding_clusters {
    my ($clusters,$quant_hash,$minUI) = @_;  ## passed by reference .. array and hash and scalar
    my @fields = ();
    my $u_ratio;
    my @output;
    foreach my $locus (@$clusters) {
	@fields = split ("\t", $$quant_hash{$locus});
	if($fields[4] != 0) {
#	    $u_ratio = $fields[6] / $fields[4]; ## rep_total / total
	    $u_ratio = $fields[6];
	} else {
	    $u_ratio = 0;
	}
	if(($u_ratio >= $$minUI) and ($fields[7] ne "N")) { ## DicerCall cannot be "N"
	    push(@output, $locus);
	}
    }
    return @output;
}

sub requant {
    my($output,$oldclus,$finalclus,$bamfile,$dicer_min,$dicer_max,$strand_cutoff,$dicer_cutoff,$phasesize,$names,$read_group) = @_; ## passed by reference .. first one is array, hp_hash and miR_hash are hashes, others scalars

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
#    my $repnorm_norm;
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
    
    my %final_clusters = ();
    my %old_clusters = ();
    my $old;
    my $final;
    
    my $ui;
    
    my @lastfields = ();
    my $laststring;
    foreach $final (@$finalclus) {
	$final_clusters{$final} = 1;
    }
    
    foreach $old (@$oldclus) {
	$old_clusters{$old} = 1;
    }
    
    ## Delete defunct entries
    foreach $old (@$oldclus) {
	unless(exists($final_clusters{$old})) {
	    delete $$output{$old};
	}
    }
    
    foreach my $locus (@$finalclus) {
	unless(exists($$output{$locus})) {
	    
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
		log_it($logfile, "\nFATAL in sub-routine \"quant\" : could not parse coordinates from locus $locus\n");
		exit;
	    }
	    # initialize phase hash
	    for($i = $loc_start; $i <= $loc_stop; ++$i) {
		$phase_hash{$i} = 0;
	    }
	    if($$read_group) {
		open(SAM, "samtools view -F 0x4 -r $$read_group $$bamfile $locus |");
	    } else {
		open(SAM, "samtools view -F 0x4 $$bamfile $locus |");
	    }
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
		
		# determine the number of total mappings for the read from the XX:i tag.
		$nh = '';
		if($_ =~ /\tXX:i:(\d+)/) {
		    $nh = $1;
		}

		unless($nh) {
		    log_it($logfile, "\nFATAL in sub-routine \"quant\": Could not parse XX:i flag from sam line $_\n");
		    exit;
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
	    $$output{$locus} .= "$locus";
	    
	    # [1] : name from name hash
	    $$output{$locus} .= "\t$$names{$locus}";
	    
	    # [2] : STRAND .. in this sub-routine, based solely upon $$strand_cutoff.
	    
	    if($frac_watson >= $$strand_cutoff) {
		$strand = "+";
	    } elsif($frac_crick >= $$strand_cutoff) {
		$strand = "-";
	    } else {
		$strand = "\.";
	    }
	    $$output{$locus} .= "\t$strand";
	    
	    # [3] FRAC_WAT
	    $$output{$locus} .= "\t$frac_watson";
	    
	    # calculate mappings, reporting only in raw in this sub-routine
	    # [4] TOTAL
	    $$output{$locus} .= "\t$total";
	    
	    # [5] UNIQUE MAPPERS
	    $$output{$locus} .= "\t$uniques";
	    
	    # [6] REP-TOTAL
	    $ui = sprintf("%.4f",($repnorm/$total));
	    $$output{$locus} .= "\t$ui";
	    
	    # [7] DICER ... either 'N' or a number within the dicer_min to dicer_max range
	    # was calculated above
	    $$output{$locus} .= "\t$dicer_call";
	    
	    # [8] and [9]: Phase p offset and phase p-value
	    if(exists($phase_p_values{$locus})) {
		$$output{$locus} .= "\t$phase_offsets{$locus}";
		$$output{$locus} .= "\t$phase_p_values{$locus}";
	    } else {
		$$output{$locus} .= "\tND\tND";
	    }
	    
	    # [10] SHORT
	    $$output{$locus} .= "\t$internal{'short'}";
	    
	    # [11] LONG
	    $$output{$locus} .= "\t$internal{'long'}";
	    
	    # [12] through whenenver .. Dicer 
	    for($i = $$dicer_min; $i <= $$dicer_max; ++$i) {
		$$output{$locus} .= "\t$internal{$i}";
	    }
	}
	
	## ensure the correct name is inserted, for all loci
	@lastfields = split ("\t", $$output{$locus});
	$lastfields[1] = $$names{$locus};
	$laststring = join("\t", @lastfields);
	$$output{$locus} = $laststring;
	
    }
}

sub get_islands {
    my($bamfile,$mindepth,$expected_faidx,$read_group) = @_;
    # go chr by chr, using the .fai index file to get the chr names
    my @chrs = ();
    open(FAI, "$expected_faidx");
    my @fields = ();
    
    my $genome_sum = 0;
    my %fai_hash = ();
    while (<FAI>) {
	chomp;
	@fields = split ("\t", $_);
	push(@chrs, $fields[0]);
	$fai_hash{$fields[0]} = $fields[1];
	$genome_sum += $fields[1];
    }
    close FAI;
    
    my @islands = ();
    
    my $last_start;
    my $last_ok;
    
    my $island;
    
    my $mindpad;
    
    ## A progress bar based on the % of genome (in nts) analyzed, updated only at the completion of chromosomes
    log_it($logfile,"\n\tProgress in sub-routine get_islands \(dots = 5 percent\; updates only at completion of a chromosome\): ");
    my $dots_printed = 0;
    my $frac_done = 0;
    my $needed_dots = 0;
    foreach my $chr (@chrs) {
	
	$last_start = -1;
	$last_ok = -1;
	if($read_group) {
	    open(DEPTH, "samtools view -F 0x4 -r $read_group -b -u $bamfile $chr | samtools depth /dev/stdin |");
	} else {
	    open(DEPTH, "samtools view -F 0x4 -b -u $bamfile $chr | samtools depth /dev/stdin |");
	}
	while (<DEPTH>) {
	    chomp;
	    @fields = split ("\t", $_);
	    if($fields[2] >= $mindepth) {
		if($last_ok == -1) {
		    ## first start on the chr
		    $last_start = $fields[1];
		    $last_ok = $fields[1];
		} elsif ($fields[1] == $last_ok + 1) {
		    ## continuing to extend an island
		    $last_ok = $fields[1];
		} else {
		    ## close the previous island, open a new one
		    $island = "$chr" . ":" . "$last_start" . "-" . "$last_ok";
		    push(@islands,$island);
		    $last_start = $fields[1];
		    $last_ok = $fields[1];
		}
	    }
	}
	close DEPTH;
	## close the last one, if there is one
	unless($last_ok == -1) {
	    $island = "$chr" . ":" . "$last_start" . "-" . "$last_ok";
	    push(@islands,$island);
	}
	
	## update progress bar
	$frac_done += ($fai_hash{$chr} / $genome_sum);
	$needed_dots = int ($frac_done / 0.05);
	until($dots_printed == $needed_dots) {
	    log_it($logfile,"\.");
	    ++$dots_printed;
	}
	
    }
    log_it($logfile," Done\n");
    return @islands;
    
}

sub check_RNALfold {
    (open(RL, "RNALfold --version |")) || return 0;
    my $version = <RL>;
    chomp $version;
    close RL;
    if($version =~ /RNALfold \d/) {
	return $version;
    } else {
	return 0;
    }
}

sub check_RNAeval {
    (open(RE, "RNAeval --version |")) || return 0;
    my $version = <RE>;
    chomp $version;
    close RE;
    if($version =~ /RNAeval \d/) {
	return $version;
    } else {
	return 0;
    }
}

sub check_bowtie {
    (open(B, "bowtie --version |")) || return 0;
    my $version = <B>;
    chomp $version;
    close B;
    if($version =~ /bowtie version (\S+)/) {
	return $1;
    } else {
	return 0;
    }
}

sub check_bowtie_build {
    (open(B, "bowtie-build --version |")) || return 0;
    my $version = <B>;
    chomp $version;
    close B;
    if($version =~ /bowtie-build version (\S+)/) {
	return $1;
    } else {
	return 0;
    }
}   

sub trim_FA {
    my($untrimmed,$adapter) = @_;
    
    # Split comma-delimited paths...
    my @untrimmed_files = split (",",$untrimmed);
    
    # Split comma-delimited adapters...
    my @adapters = split (",",$adapter);
    
    # holder for all trimmed file names ...
    my @trimmed_files = ();
    
    my $used_adapter;
    
    # check for number of adapters...
    unless (((scalar @untrimmed_files) == (scalar @adapters)) or
	    ((scalar @adapters) == 1)) {
	log_it($logfile, "\nFATAL: The number of apdaters provided must be equal to the number of untrimmedFA files, or just one\n");
	exit;
    }
    
    foreach my $ut_file (@untrimmed_files) {
	if(@adapters) {
	    $used_adapter = shift @adapters;
	}
	log_it($logfile, "\nAdapter trimming file $ut_file with adapter $used_adapter ...");
	(open(IN, "$ut_file")) || return 0;
	my $trimmedFA = "$ut_file";
	$trimmedFA =~ s/\..*$//g;  ## strip any extension
	$trimmedFA .= "_trimmed.fasta"; ## add new extension
	(open(OUT, ">$trimmedFA")) || return 0;
	my $header;
	my $trim_len;
	my $no_insert = 0; ## includes adapter-only and no-adapter cases
	my $too_short = 0;
	my $amb = 0;
	my $ok;
	my $trim_seq;
	while (<IN>) {
	    chomp;
	    if($_ =~ /^>/) {
		$header = $_;
	    } else {
		$trim_len = 0;
		while ($_ =~ /$used_adapter/ig) {
		$trim_len = (pos $_) - (length $used_adapter);
		}
		if($trim_len == 0) {
		    ++$no_insert;
		} elsif ($trim_len < 15) {
		    ++$too_short;
		} else {
		    $trim_seq = substr($_,0,$trim_len);
		    if($trim_seq =~ /[^ATGCatcg]/) {
			++$amb;
		    } else {
			++$ok;
			print OUT "$header\n$trim_seq\n";
		    }
		}
	    }
	}
	close IN;
	close OUT;
	log_it($logfile," Done\n");
	log_it($logfile,"\tNo insert \(includes adapter-only and no-adapter cases combined\): $no_insert\n");
	log_it($logfile,"\tToo short \(less than 15nts\): $too_short\n");
	log_it($logfile,"\tAmbiguous bases after trimming: $amb\n");
	log_it($logfile,"\tOK - output: $ok\n");
	log_it($logfile,"\tResults in file $trimmedFA\n");
	push(@trimmed_files,$trimmedFA);
    }
    return @trimmed_files;
}

sub trim_FQ {
    my($untrimmed,$adapter) = @_;
    
    # Split comma-delimited paths...
    my @untrimmed_files = split (",",$untrimmed);
    
    # Split comma-delimited adapters...
    my @adapters = split (",",$adapter);
    
    # holder for all trimmed file names ...
    my @trimmed_files = ();
    
    my $used_adapter;
    
    # check for number of adapters...
    unless (((scalar @untrimmed_files) == (scalar @adapters)) or
	    ((scalar @adapters) == 1)) {
	log_it($logfile, "\nFATAL: The number of apdaters provided must be equal to the number of untrimmedFQ files, or just one\n");
	exit;
    }

    foreach my $ut_file (@untrimmed_files) {
	if(@adapters) {
	    $used_adapter = shift @adapters;
	}
	log_it($logfile, "\nAdapter trimming file $ut_file with adapter $used_adapter ...");
	(open(IN, "$ut_file")) || return 0;
	my $trimmedFQ = "$ut_file";
	$trimmedFQ =~ s/\..*$//g; ## strip extension
	$trimmedFQ .= "_trimmed.fastq"; ## add new extension
	(open(OUT, ">$trimmedFQ")) || return 0;
	my $header;
	my $trim_len;
	my $no_insert = 0; ## includes adapter-only and no-adapter cases
	my $too_short = 0;
	my $amb = 0;
	my $ok;
	my $trim_seq;
	my $plus;
	my $seq;
	my $qual;
	my $trim_qual;
	while (<IN>) {
	    chomp;
	    $header = $_;
	    $seq = <IN>;
	    chomp $seq;
	    $plus = <IN>;
	    chomp $plus;
	    $qual = <IN>;
	    chomp $qual;
	    
	    # crude validation
	    unless(($header =~ /^\@/) and
		   ($plus =~ /^\+/)) {
		log_it($logfile,"\nFATAL: FASTQ format parse error in sub-routine trim_FQ\n");
		exit;
	    }
	    $trim_len = 0;
	    while ($seq =~ /$used_adapter/ig) {
		$trim_len = (pos $seq) - (length $adapter);
	    }
	    if($trim_len == 0) {
		++$no_insert;
	    } elsif ($trim_len < 15) {
		++$too_short;
	    } else {
		$trim_seq = substr($seq,0,$trim_len);
		if($trim_seq =~ /[^ATGCatcg]/) {
		    ++$amb;
		} else {
		    ++$ok;
		    $trim_qual = substr($qual,0,$trim_len);
		    print OUT "$header\n$trim_seq\n$plus\n$trim_qual\n";
		}
	    }
	}
	
	close IN;
	close OUT;
	log_it($logfile," Done\n");
	log_it($logfile,"\tNo insert \(includes adapter-only and no-adapter cases combined\): $no_insert\n");
	log_it($logfile,"\tToo short \(less than 15nts\): $too_short\n");
	log_it($logfile,"\tAmbiguous bases after trimming: $amb\n");
	log_it($logfile,"\tOK - output: $ok\n");
	log_it($logfile,"\tResults in file $trimmedFQ\n");
	push(@trimmed_files, $trimmedFQ);
    }
    return @trimmed_files;
}

sub perform_alignment {
    my($infile,$genome,$type) = @_;
    log_it($logfile,"\nBeginning alignment of reads from $infile to genome $genome ...\n");
    log_it($logfile, "\tCheck for ebwt bowtie indices for the genome: ");
    
    # Search for the bowtie indices
    my $ebwt_check = check_ebwt($genome,$type);
    # If absent, call bowtie-build
    if($ebwt_check) {
	log_it($logfile, "PRESENT\n");
    } else {
	log_it($logfile, "ABSENT. Attempting to build with bowtie-build... ");
	my $bb_log = bowtie_build($genome,$type);
	log_it($logfile, "Done.\n\tSee $bb_log for output from bowtie-build\n");
    }
    
    # build bowtie query string
    my $bowtie = "bowtie";
    if($type eq "a") {
	$bowtie .= " -f -v 1 -a --best --strata -S $genome $infile";
    } elsif ($type eq "q") {
	$bowtie .= " -q -v 1 -a --best --strata -S $genome $infile";
    } elsif ($type eq "c") {
	$bowtie .= " -f -C";
	my $cs_index = "$genome" . "_CS";
	## check to see if 'infile' is a concatenation of a trimmed CS and trimmed QV file
	if($infile =~ /:::/) {
	    my @two_files = split (":::",$infile);
	    $bowtie .= " -Q $two_files[1] -v 1 -a --best --strata -S --col-keepends $cs_index $two_files[0]";
	} else {
	    $bowtie .= " -v 1 -a --best --strata -S --col-keepends $cs_index $infile";
	}
    }
    
    # define the outfile base name .. strip extensions from the infile
    my $outfile_base = $infile;
    $outfile_base =~ s/\..*$//g;
    
    # Open the processes
    (open(BOW, "$bowtie |")) || return 0;
    (open(OUT, "| samtools view -S -b -u - | samtools sort - $outfile_base")) || return 0;
    
    # Parsing the bowtie stream
    my $in_header = 1;
    my $last_read = "NULL";
    my @read_aligns = ();
    my @sam_fields = ();
    my $xx;
    my $index;
    my $output = 0; ## includes unmapped
    
    # Progress bar
    log_it($logfile,"Progress in bowtie \(dot = 1E5 reads processed\)\n");
    my $done = 0;
    my $one_hun_k = 0;
    my $million = 0;
    while(<BOW>) {
	if($_ =~ /^\@/) {
	    $in_header = 1;
	    # Change sort information to "coordinate"
	    $_ =~ s/SO:(\S+)/SO:coordinate/;
	    # pass it along
	    print OUT "$_";
	} else {
	    if($in_header) {
		$in_header = 0;
		# Add some other info
		print OUT "\@PG\tID:ShortStack\n";
		print OUT "\@CO\tTag XX from ShortStack indicates the total number of possible, equally likely alignments, of which just one randomly selected was reported\n";
		print OUT "\@CO\tTag XA from bowtie indicates the 'stratum' of the aligned read\n";
		print OUT "\@CO\tTag XM from bowtie indicates that the read had no reported alignments\n";
		print OUT "\@CO\tSee the SAM spec for the other, standard tags\n";
	    }
	    # data lines
	    @sam_fields = split ("\t", $_);
	    if(($sam_fields[0] ne $last_read) and ($last_read ne "NULL")) {
		$xx = scalar @read_aligns;
		##randomly select one
		$index = int(rand($xx));
		## tack on the XX tag
		$read_aligns[$index] =~ s/\n/\tXX:i:$xx\n/g;
		++$output;
		print OUT "$read_aligns[$index]";
		## reset
		@read_aligns = ();
		## Progress
		++$done;
		if($done == 100000) {
		    ++$one_hun_k;
		    log_it($logfile, ".");
		    $done = 0;
		    if($one_hun_k == 10) {
			++$million;
			log_it($logfile, " $million million reads analyzed and counting\n");
			$one_hun_k = 0;
		    }
		}
	    }
	    push(@read_aligns,$_);
	    $last_read = $sam_fields[0];
	}
    }
    # the last one
    $xx = scalar @read_aligns;
    ##randomly select one
    $index = int(rand($xx));
    ## tack on the XX tag
    $read_aligns[$index] =~ s/\n/\tXX:i:$xx\n/g;
    ++$output;
    print OUT "$read_aligns[$index]";
    
    log_it($logfile, "\n\n");
    
    # close streams
    close BOW;
    close OUT;
    
    # Report
    log_it($logfile, "Total of $output sam lines output \(includes lines for non-mapped reads\).\n");
    
    # Verify file
    my $expected_bam = "$outfile_base" . ".bam";
    unless(-r $expected_bam) {
	log_it($logfile, "\nFATAL: Expected BAM file $expected_bam was not created after bowtie run\!\n\n");
	exit;
    }
    return $expected_bam;
}

sub check_ebwt {
    my($base,$type) = @_;
    my @ebwt_suffixes = qw(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt);
    my $ebwt_count = 0;
    my $ebwt_file;
    foreach my $ebwt_suffix (@ebwt_suffixes) {
	if($type eq "c") {
	    $ebwt_file = "$base" . "_CS" . "$ebwt_suffix";
	} else {
	    $ebwt_file = "$base" . "$ebwt_suffix";
	}
	if(-r $ebwt_file) {
	    ++$ebwt_count;
	}
    }
    if($ebwt_count == 6) {
	return 1;
    } elsif ($ebwt_count > 0) {
	return 0;
    } else {
	return 0;
    }
}

sub bowtie_build {
    my($fasta,$type) = @_;
    my $log = "$fasta" . "_bowtie_build_log.txt";
    if($type eq "c" ) {
	my $out_base = "$fasta" . "_CS";
	system("bowtie-build -C $fasta $out_base > $log");
    } else {
	system("bowtie-build $fasta $fasta > $log");
    }
    return $log;
}

sub validate_bam {
    my($bam,$fai) = @_;
    log_it($logfile,"\nValidating bamfile $bam ...\n");
    
    # Is it readable?
    log_it($logfile, "\tFile readable\? ");
    unless(-r $bam) {
	log_it($logfile, "NO -- ABORT\n");
	return 0;
    }
    log_it($logfile, "YES\n");
    
    # header checks
    my %h_seqs = ();
    log_it($logfile, "\tHeader present\? ");
    (open(H, "samtools view -H $bam |")) || return 0;
    my $h;
    my $so;
    my %rg = ();
    my $rgname;
    while (<H>) {
	if($_ =~ /^\@/) {
	    unless($h) {
		$h = 1;
		log_it($logfile,"YES\n\tSorted by coordinate\? ");
	    }
	    if($_ =~ /^\@HD\t.*SO:(\S+)/) {
		if($1 eq "coordinate") {
		    $so = 1;
		} else {
		    $so = 0;
		}
	    }
	    if($_ =~ /^\@SQ\t.*SN:(\S+)/) {
		$h_seqs{$1} = 1;
	    }
	    if($_ =~ /^\@RG\t.*ID:(\S+)/) {
		$rgname = $1;
		if($_ =~ /DS:(\S+)/) {
		    $rg{$rgname} = $1;
		} else {
		    $rg{$rgname} = "None";
		}
	    }
	}
    }
    close H;
    if($so) {
	log_it($logfile, "YES\n");
    } else {
	log_it($logfile, "FAIL -- bamfile must be sorted by \'coordinate\' as indicated in the header SO:tag-- ABORT\n");
	return 0;
    }
    
    # does it match the genome?
    my %fai_names = ();
    log_it($logfile,"\tMatches genome\? ");
    (open(FAI, "$fai")) || return 0;
    my @fields = ();
    while (<FAI>) {
	@fields = split ("\t", $_);
	$fai_names{$fields[0]} = 1;
    }
    close FAI;
    my $h_seq;
    while(($h_seq) = each %h_seqs) {
	unless(exists($fai_names{$h_seq})) {
	    log_it($logfile, "FAIL -- reference name $h_seq was not found in genome -- ABORT\n");
	    return 0;
	}
    }
    log_it($logfile, "PASS\n");
    
    # Are there read groups defined, and if so, what are they?
    log_it($logfile,"\tRead Groups\? ");
    if(%rg) {
	log_it($logfile, "YES. And they are...\n");
	while(($rgname) = each %rg) {
	    log_it($logfile,"\t\t$rgname\tDescription: $rg{$rgname}\n");
	}
    } else {
	log_it($logfile, "NO\n");
	$rg{'NULL'} = 1;
    }
    
    # Does it have XX tags?  Just check the first line
    log_it($logfile, "\tHas XX tags\? ");
    (open(SAM, "samtools view $bam |")) || return 0;
    my $check = <SAM>;
    close SAM;
    if($check =~ /\tXX:i:(\d+)/) {
	log_it($logfile,"PASS\n");
    } else {
	log_it($logfile,"FAIL -- bam file needs to have XX tags per ShortStack spec -- ABORT\n");
	return 0;
    }
    
    ## Is it indexed? If not, index it.
    my $bai = "$bam" . ".bai";
    log_it($logfile, "\tIndexed file $bai readable\? ");
    if(-r $bai) {
	log_it($logfile, "YES\n");
    } else {
	log_it($logfile, "Not yet ... ");
	system "samtools index $bam";
	if(-r $bai) {
	    log_it($logfile,"Now it does\n");
	} else {
	    log_it($logfile,"FAILED to create it with samtools index -- ABORT\n");
	    return 0;
	}
    }
    
    ## Get the total number of reads
    (open(STAT, "samtools flagstat $bam |")) || return 0;
    my $n_mapped = 0;
    while (<STAT>) {
	if($_ =~ /^(\d+) \+ \d+ mapped/) {
	    $n_mapped = $1;
	}
    }
    log_it($logfile,"\tNumber of aligned reads: $n_mapped\n");

    return %rg;
}

sub check_adapter {
    my($input) = @_;
    my @entries = split (",",$input);
    foreach my $a (@entries) {
	unless($a =~ /^[ACTGactg]{8,}$/) {
	    return 0;
	}
    }
    return 1;
}

sub merge_em {
    my($bamfiles,$outdir) = @_;  ## passed by reference .. array, scalar
    log_it($logfile, "\nMerging initial alignment files...\n");
    
    # Write a header that includes all of the read groups
    # First copy the header of the first bamfile
    
    system "samtools view -H $$bamfiles[0] > ShortStack_SAM_head_temp.txt";
    
    # now add the read groups
    my $out_name = "$$outdir" . ".bam";
    my $command = "samtools merge -r -h ShortStack_SAM_head_temp.txt $out_name";
    my $rg;
    my $bamfile_part;
    foreach $bamfile_part (@$bamfiles) {
	# determine the number of reads
	(open(STAT, "samtools flagstat $bamfile_part |")) || return 0;
	my $n_mapped;
	my $n_all;
	my $n_unmapped;
	while (<STAT>) {
	    if($_ =~ /^(\d+) \+ (\d+) in total/) {
		$n_all = $1 + $2;
	    } elsif ($_ =~ /^(\d+) \+ (\d+) mapped/) {
		$n_mapped = $1 + $2;
	    }
	}
	close STAT;
	$n_unmapped = $n_all - $n_mapped;
	$rg = $bamfile_part;
	$rg =~ s/\..*$//g;
	log_it($logfile, "\tRead Group: $rg added to list .. Mapped: $n_mapped Unmapped: $n_unmapped\n");
	# add to header
	system "echo \"\@RG\tID:$rg\tDS:Mapped=$n_mapped\;Unmapped=$n_unmapped\" >> ShortStack_SAM_head_temp.txt";
	# add to command
	$command .= " $bamfile_part";
    }
    
    # Call the merger
    system "$command";
    
    # check success
    unless(-r $out_name) {
	log_it($logfile, "\nFATAL: merge of initial alignments failed\n");
	exit;
    }
    
    # clean up by deleting the component bam files and the temp header file
    foreach $bamfile_part (@$bamfiles) {
	system "rm -f $bamfile_part";
    }
    system "rm -f ShortStack_SAM_head_temp.txt";
    
    # log it
    log_it($logfile,"\n\tCompleted merge of alignments. Final file is $out_name\n");
    
    return $out_name;
}

sub trim_CS {
    my($untrimmed,$adapter,$untrimmedCSQV) = @_;
    
    # Split comma-delimited paths...
    my @untrimmed_files = split (",",$untrimmed);
    my @untrimmed_QV_files = ();
    if($untrimmedCSQV) {
	@untrimmed_QV_files = split (",",$untrimmedCSQV);
    }
	
    
    # Split comma-delimited adapters...
    my @adapters = split (",",$adapter);
    
    # holder for all trimmed file names ...
    my @trimmed_files = ();
    
    my $used_adapter;
    
    # check for number of adapters...
    unless (((scalar @untrimmed_files) == (scalar @adapters)) or
	    ((scalar @adapters) == 1)) {
	log_it($logfile, "\nFATAL: The number of apdaters provided must be equal to the number of untrimmedCS files, or just one\n");
	exit;
    }
    
    my $ut_qv_file;
    my $trimmedQV;
    foreach my $ut_file (@untrimmed_files) {
	if(@adapters) {
	    $used_adapter = shift @adapters;
	}
	my $cs_used_adapter = adapter2cs($used_adapter);
	if(@untrimmed_QV_files) {
	    $ut_qv_file = shift @untrimmed_QV_files;
	}
	log_it($logfile, "\nAdapter trimming file $ut_file with adapter $used_adapter (translated to $cs_used_adapter) ...");
	(open(IN, "$ut_file")) || return 0;
	my $trimmedCS = "$ut_file";
	$trimmedCS =~ s/\..*$//g;  ## strip any extension
	$trimmedCS .= "_trimmed.csfasta"; ## add new extension
	(open(OUT, ">$trimmedCS")) || return 0;
	
	if($ut_qv_file) {
	    (open(QIN, "$ut_qv_file")) || return 0;
	    $trimmedQV = $ut_qv_file;
	    $trimmedQV =~ s/\..*$//g;
	    $trimmedQV .= "_trimmed.qual";
	    (open(QOUT, ">$trimmedQV")) || return 0;
	}
	
	
	my $header;
	my $trim_len;
	my $no_insert = 0; ## includes adapter-only and no-adapter cases
	my $too_short = 0;
	my $amb = 0;
	my $ok;
	my $trim_seq;
	
	my $qin_line;
	
	my $cs_l_num = 0;
	my $qv_l_num = 0;
	
	my $qv_head;
	
	my @qv_vals = ();
	
	while (<IN>) {
	    chomp;
	    ++$cs_l_num;
	    if($ut_qv_file) {
		$qin_line = <QIN>;
		chomp $qin_line;
		++$qv_l_num;
	    }
	    if($_ =~ /^\#/) {
		next;
	    }
	    if($_ =~ /^>/) {
		$header = $_;
		if($ut_qv_file) {
		    $qv_head = $qin_line;
		    unless($header eq $qv_head) {
			log_it($logfile, "\nFATAL: csfasta and qual files out of sync at lines $cs_l_num (csfasta) and $qv_l_num (qual) .. headers $header and $qv_head DONT match\n");
			exit;
		    }
		}
	    } else {
		$trim_len = 0;
		while ($_ =~ /$cs_used_adapter/ig) {
		    $trim_len = (pos $_) - (length $cs_used_adapter) - 1; # subtract 1 to remove the hybrid color .. leading T is still included, so sRNA length is trim_len - 1
		}
		if($trim_len <= 1) {  ## 1 is no insert for colorspace .. leading T still there
		    ++$no_insert;
		} elsif (($trim_len - 1) < 15) {  ## accounts for leading T (or other base key)
		    ++$too_short;
		} else {
		    $trim_seq = substr($_,0,$trim_len);  ## because of above, this also chops the hybrid color
		    if($trim_seq =~ /^[ATGC][0123]+$/) {
			++$ok;
			print OUT "$header\n$trim_seq\n";
			
			if($ut_qv_file) {
			    @qv_vals = split (" ", $qin_line);
			    for(my $j = 0; $j <= ($trim_len - 2); ++$j) {
				if($j == 0) {
				    print QOUT "$qv_head\n$qv_vals[0]";
				} else {
				    print QOUT " $qv_vals[$j]";
				}
			    }
			    print QOUT "\n";
			}
			
		    } else {
			++$amb;
		    }
		}
	    }
	}
	close IN;
	close OUT;
	log_it($logfile," Done\n");
	log_it($logfile,"\tNo insert \(includes adapter-only and no-adapter cases combined\): $no_insert\n");
	log_it($logfile,"\tToo short \(less than 15nts\): $too_short\n");
	log_it($logfile,"\tAmbiguous bases after trimming: $amb\n");
	log_it($logfile,"\tOK - output: $ok\n");
	if($ut_qv_file) {
	    close QIN;
	    close QOUT;
	    log_it($logfile,"\tResults in files $trimmedCS and $trimmedQV\n");
	    my $concat = "$trimmedCS" . ":::" . "$trimmedQV";
	    push(@trimmed_files,$concat);
	} else {
	    log_it($logfile,"\tResults in file $trimmedCS\n");
	    push(@trimmed_files,$trimmedCS);
	}
    }
    return @trimmed_files;
}

sub adapter2cs {
    my($base_adapter) = @_;
    $base_adapter = uc $base_adapter;
    unless($base_adapter =~ /^[ATGC]+$/) {
        log_it($logfile, "\nFATAL: Invlaid adapter $base_adapter found in sub-routine adapter2cs\n");
        exit;
    }
    my %colors = (
        'AA' => 0,
        'CC' => 0,
        'GG' => 0,
        'TT' => 0,
        'AC' => 1,
        'CA' => 1,
        'GT' => 1,
        'TG' => 1,
        'AG' => 2,
        'CT' => 2,
        'GA' => 2,
        'TC' => 2,
        'AT' => 3,
        'CG' => 3,
        'GC' => 3,
        'TA' => 3,
        );
    my $dibase;
    my $color_adapter;
    my @letters = split('', $base_adapter);
    for(my $i = 0; $i <= ((scalar @letters) - 2); ++$i) {
        $dibase = "$letters[$i]" . "$letters[($i + 1)]";
        unless(exists($colors{$dibase})) {
            log_it($logfile, "\nFATAL: Failure to lookup the dibase $dibase in sub-routine adapter2cs\n");
            exit;
	}
        $color_adapter .= $colors{$dibase};
    }
    return $color_adapter;
}

sub check_csqv {
    my($cs,$qv) = @_;
    my @css = split (",",$cs);
    my @qvs = split (",",$qv);
    unless((scalar @css) == (scalar @qvs)) {
	log_it($logfile, "\nFATAL: Unequal number of files specified in options --untrimmedCS and --untrimmedCSQV\n");
	return 0;
    }
    foreach my $csf (@css) {
	unless(-r $csf) {
	    log_it($logfile, "\nFATAL: File $csf was not readable\n");
	    return 0;
	}
    }
    foreach my $qvf (@qvs) {
	unless(-r $qvf) {
	    log_it($logfile, "\nFATAL: File $qvf was not readable\n");
	    return 0;
	}
    }
    return 1;
}

	       

__END__

=head1 LICENSE

ShortStack.pl

Copyright (C) 2012-2013 Michael J. Axtell                                                             
                                                                                                 
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

Axtell MJ. (2013) ShortStack: Comprehensive annotation and quantification of small RNA genes.  RNA 19:740-751. doi:10.1261/rna.035279.112

=head1 VERSION

1.2.0 :: Released October 16, 2013

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 DEPENDENCIES

perl
samtools
RNAeval
RNALfold
bowtie (0.12.x OR 1.x) [only if aligning reads]
bowtie-build (0.12.x OR 1.x) [only if aligning reads]

ShortStack.pl is a perl script, so it needs perl to compile. It expects to find perl at /usr/bin/perl. If this is not where your perl is, modify line 1 of ShortStack.pl (the hashbang) accordingly. It also requires the package Getopt::Long, which I think is standard in most Perl distributions. If this package is not installed, get it from CPAN.

samtools <http://samtools.sourceforge.net/> needs to be installed in your PATH. ShortStack was developed using samtools 0.1.18. Other versions should be OK as far as I know, but let me know if not!

RNAeval and RNALfold are from the ViennaRNA package. See <http://www.tbi.univie.ac.at/~ronny/RNA/vrna2.html>. Both need to be installed in your PATH.

bowtie and bowtie-build must be version "1" .. either 0.12.x or 1.x.  These are required ONLY if you are aligning reads to the genome. Bowtie can be found at http://bowtie-bio.sourceforge.net/index.shtml. Like the other dependencies, they must be in your PATH.

=head1 INSTALL

There is no 'real' installation. After installing the dependencies (see above), you should check to make sure the ShortStack.pl is executable. It should be, but if not you can:                                                                
                                                                                                 
    chmod +x ShortStack.pl                                                         
                                                                                                 
For convenience, you can add it to your PATH .. for instance

    sudo cp ShortStack.pl /usr/bin/


=head1 USAGE
                                                                                                                             
    Shortstack.pl [options] [genome.fasta] 

There are three modes which differ in the the types of pre-analysis that are performed. Each of the modes has a different set of REQUIRED options:

Mode 1: Trim small RNA-seq reads to remove 3' adapter seqeuence, align them, and then analyze.  Required options:

    --untrimmedFA OR --untrimmedFQ OR --untrimmedCS
    
    --adapter

Mode 2: Align pre-trimmed small RNA-seq reads, and then analyze. Required options:

    --trimmedFA OR --trimmedFQ OR --trimmedCS

Mode 3: Analyze a pre-existing BAM alignment of small RNA-seq reads. Require option:

    --bamfile

Additionally, in modes 1 or 2, the option --align_only will terminate analysis after making the alignment file.

=head1 TUTORIAL

A full tutorial with sample Arabidopsis data can be found at http://axtelldata.bio.psu.edu/data/ShortStack_TestData/

=head1 OPTIONS

--help : Print a help message and then quit.

--version : Print the version number and then quit.

--outdir [string] : Name of directory to be created to receive results of the run.  Deafults to "ShortStack_[time]", where time is "UNIX time" (the number of non-leap seconds since Jan 1, 1970 UCT), if not provided

--untrimmedFA [string] : Path to untrimmed small RNA-seq data in FASTA format. Multiple datasets can be provided as a comma-delimited list.

--untrimmedFQ [string] : Path to untrimmed small RNA-seq data in FASTQ format.  Multiple datasets can be provided as a comma-delimited list.

--untrimmedCS [string] : Path to untrimmed small RNA-seq data is Colorpsace-FASTA format. Multiple datasets can be provided as a comma-delimited list. 

--untrimmedCSQV [string] : Path to untrimmed small RNA-seq quality values from SOLiD datasets, corresponding to colorspace-fasta files given in option --untrimmedCS. Multiple datasets can be provided as a comma-delimited list.

--adapter [string] : Sequence of 3' adapter to search for during adapter trimming. Must be at least 8 nts in length, and all ATGC characters. Required if either --untrimmedFA or --untrimmedFQ are specified. Multiple adapters (for when multiple input untrimmedFA/FQ files are specified) can be provided as a comma-delimited list.

--trimmedFA [string] : Path to trimmed and ready to map small RNA-seq data in FASTA format. Multiple datasets can be provided as a comma-delimited list.

--trimmedFQ [string] : Path to trimmed and ready to map small RNA-seq data in FASTQ format. Multiple datasets can be provided as a comma-delimited list.

--trimmedCS [string] : Path to trimmed and ready to map small RNA-seq data in Colorspace-FASTA format. Multiple datasets can be provided as a comma-delimited list.                                                                             
                                                                                                                                                                                                                                                
--trimmedCSQV [string] : Path to trimmed color-space quality value files, corresponding to trimmed colorspace-fasta files provided in option --trimmedCS. Multiple datasets can be provided as a comma-delimited list.              

--align_only : Exits program after completion of small RNA-seq data alignment, creating BAM file.

--bamfile [string] : Path to properly formatted and sorted BAM alignment file of small RNA-seq data.

--read_group [string] : Analyze only the indicated read-group. Read-group must be specified in the bam alignment file header. Default = [not active -- all reads analyzed]

--inv_file [string] : PATH to an einverted-produced .inv file of inverted repeats within the genome of interest.  Not required but strongly suggested for more complete annotations of hairpin-derived small RNA genes.  Default = {blank}.  Not needed for runs in "nohp" mode or runs in "count" mode (because "count" mode forces "nohp" mode as well).  A typical eniverted run uses default parameters except "-maxrepeat 10000", in order to capture long IRs.

--flag_file [string] : PATH to a simple file of genomic loci of interest.  The ShortStack-analyzed small RNA clusters will be analyzed for overlap with the loci in the flag_file .. if there is any overlap (as little as one nt), it will be reported.  Format for this file is describe below.

--mindepth [integer] : Minimum depth of mapping coverage to define an 'island'.  Default = 20.  Must be at least 2, more than 5 preferred.

--pad [integer] : Number of nucleotides upstream and downstream to extend initial islands during cluster definition.  Default = 100

--dicermin [integer] : Smallest size in the Dicer size range (or size range of interest).  Deafult = 20.  Must be between 15 and 35, and less than or equal to --dicermax

--dicermax [integer] : Largest size in the Dicer size range (or size range of interest).  Deafult = 24.  Must be between 15 and 35, and more than or equal to --dicermin

--minUI [float] : Minimum uniqueness index required to attempt RNA folding.  Must be a value between 0 and 1.  Zero forces all clusters to be folded; default: 0.1

--maxhpsep [integer] : Maximum allowed span for a base-pair during hairpin search with RNALfold; Also serves as the maximum size of genomic query to fold with RNALfold .. loci whose unpadded size is more than --maxhpsep will not be analyzed at all with RNALfold.  Default = 300.  Must be between 50 and 2000.

--minfracpaired [float] : Minimum fraction of paired nucleotides required within a valid hairpin structure.  Default = 0.67.  Allowed values are greater than 0 and less than or equal to 1.

--minntspaired [integer] : Minimum absolute number of paired nucleotides required within a valid hairpin structure.  Default = 15.  Allowed values are greater than zero and less than or equal to --maxhpsep

--maxdGperStem [float] : Maximum deltaG / stem length allowed in a valid hairpin structure.  Stem length is 0.5 * (left_stem_length + right_stem_length).  Default = -0.5

--minfrachpdepth [float] : Minimum fraction of corrected coverage within hairpin arms to keep hairpin for further analysis.  Default = 0.67.  Allowed values between 0 and 1.  See below for details.

--miRType [string] : Either "plant" or "animal".  Defaults to "plant".  This option sets --maxmiRHPPairs, --maxmiRUnpaired, and --maxLoopLength to 150, 5, and 100,000 respectively for type "plant".  For type "animal", the three are instead set to 45, 6, and 15, respectively.

--maxmiRHPPairs [integer] : Maximum number of base pairs in a valid MIRNA hairpin. default: set by --miRType "plant" to 150.  --miRType "animal" sets to 45 instead.  When provided, user settings will override miRType settings. 

--maxmiRUnpaired [integer] : Maximum number of unpaired miRNA nts in a miRNA/miRNA* duplex. default: set by --miRType "plant" to 5.  --miRType "animal" instead sets it to 6.  When provided, user settings will override miRType settings.

--maxLoopLength [integer] : maximum allowed loop length for a valid hairpin. default: set by --miRType "plant" be essentially unlimited (100,000).  --miRType "plant" sets it to 15.  When provided, user settings will override miRType settings.

--minstrandfrac [float] : Minimum fraction of mappings to one or the other strand call a polarity for non-hairpin clusters.  Also the minimum fraction of "non-dyad" mappings to the sense strand within potential hairpins/miRNAs to keep the locus annotated as a hp or miRNA.  See below for details.  Default = 0.8.  Allowed values between 0.5 and 1.

--mindicerfrac [float] : Minimum fraction of mappings within Dicer size range to annotate a locus as Dicer-derived.  Default = 0.85.  Allowed values between 0 and 1.

--phasesize [integer] : Examine phasing only for clusters dominated by the indicated size range.  Size must be within the bounds described by --dicermin and --dicermax.  Set to 'all' to examine p-values of each locus within the Dicer range, in its dominant size.  Set to 'none' to suppress all phasing analysis.  Default = 21.  Allowed values between --dicermin and --dicermax.

--count [string] : Invokes count mode, in which user-provided clusters are annotated and quantified instead of being defined de novo.  When invoked, the file provided with --count is assumed to contain a simple list of clusters.  Count mode also forces nohp mode.  Formatting details below.  Default : Not invoked.

--nohp : If "--nohp" appears on the command line, it invokes running in "no hairpin" mode.  RNA folding, hairpin annotation, and MIRNA annotation will be skipped (likely saving significant time).  Note that --count mode forces --nohp mode as well.  Default: Not invoked.


=head1 KEY FORMATTING REQUIREMENTS AND ASSUMPTIONS

=head2 Input genome.fasta file

It is critical that this be the precise genome to which the reads in the input .bam file were mapped. If it isn't, validation of the BAM alignment file will fail and the run will be aborted.

Additionally, the chromosome names in the FASTA headers must be kept SIMPLE.  Specifically, ShortStack.pl at several points parses clusters by the regex /^(\S+):(\d+)-(\d+)$/ or some variant thereof, where the first pattern is the chromosome name.  Therefore, the chromosome names must match (\S+) .. e.g. a single string of one or more non-white-space characters, with no metacharacters.  So, ">Chr1" in your reference genome is good, but ">Chr1 | XM00023 | this is a bunch of annotation blah blah blah" is bad.  This same concern applies to the input .bam file, so your chromosome names should be shortened BEFORE mapping your reads, so that they are short and they are exactly reflected in the .bam file.

If not already present, a .fai index file for the genome will be created using samtools faidx at the beginning of the run.

=head2 Small RNA-seq reads.

FASTA or FASTQ data should be devoid of comment lines and conform to FASTA or FASTQ specs. In addition, ShortStack.pl assumes that each read will occupy a single line in the file. There is no support for paired-end reads.  Colorspace-FASTA formatted data (from SOLiD) can have comment lines. Format is assumed to conform to colorspace-FASTA specifications (beginning with a nucleotide, followed by a string of colors [0,1,2,3] or ambiguity codes [.].  For SOLiD data, the quality value files can also be input (with option --untrimmedCSQV or --trimmedCSQV); these files are expected to be space delimited list of integers, in the same order as the corresponding colorspace-fasta files.

=head2 Input .bam file

It is highly recommended that you create the BAM alignment file with ShortStack.pl itself to ensure compatibility. Unfortunately, BAM files created using versions of ShortStack.pl prior to 1.0.0 (and the now-obsoleted 'Prep_bam.pl' script) will NOT be compatibile with ShortStack v 1.0.0 and higher. Sorry. There are good reasons for this, as the previously recommended alignment methods were sub-optimal.

If you do make your own BAM alignments, outside of ShortStack.pl, they must pass the following validation steps that are performed by ShortStack:

1. The header must be present

2. The sort order of the file must be 'coordinate', as indicated by the SO: tag in the header

3. All of the chromosome names found in the header MUST also be found in the genome.fasta file

4. The data lines must contain the custom XX:i:(\d+) tag, where \d+ is an integer representing the total number of valid alignment positions for that read.

5. If the option --read_group is being used, the specified read_group must be mentioned in the header in an @RG line.

The BAM file should be indexed and have the corresonding .bam.bai index file in the same path as the bamfile. However, this is not required to pass validation .. if the index is not found, it will be created during the run.

Each mapped read must have the CIGAR string set (column 6 in the SAM specification) -- ShortStack.pl determines the small RNA lengths by parsing the CIGAR string .. if any mappings (except unmapped reads, which are ignored) have "*" entered instead of a valid CIGAR string ShortStack.pl will exit and complain.

Finally, a key assumption that ShortStack.pl makes about the BAM alignment is that each original small RNA-seq read is represented by just one alignment line. For reads that had more than one valid alignment position, just one should be randomly selected, and the total number of valid alignment positions noted in the XX tag. Users of ShortStack versions prior to 1.0.0 please note that this is quite different from the earlier methods.

=head2 --count file

If running in --count mode, the user-provided file is expected to be a simple text file containing a list of coordinates in the format : [Chr]:[start]-[stop], where Chr is defined in the genome file AND in the .bam file, and start and stop are one-based, inclusive.  The same requirement for short, non-whitespaced chromosome names as discussed above holds true for input --count files.  Comment lines, that begin with '#', are ignored.  Tab-delimited files are also accepted, provided the first column has the coordinates.  The second column in tab-delimted files is assumed to be the names of the clusters, and will be used accordingly.  Any other columns in a tab-delimited input file are ignored.

Importantly, the 'Results.txt' file produced by a previous ShortStack.pl run can be used directly in subsequent runs in --count mode.  This is useful when comparing identical intervals across multiple samples.

Note that count mode also forces nohp mode.

=head2 --flag_file

Optional.  This is a list of genomic loci to scan for overlap with one or more of the small RNA loci found/analyzed by ShortStack.  Overlap of any length is reported.  The format of the file is similar to that of the --count file:  A tab-delimited text file with coordinates in the first column, and names in the second column.  Unlike for --count files, names are required to be present in the second column for --flag_file.  Coordinates must be in the format [Chr]:[start]-[stop], where Chr is defined in the genome file AND the .bam file, and start and stop are one-based, inclusive.

=head2 --inv_file 

Unless you are running in "nohp" mode, providing a inv_file will enhance the accuracy of the hairpin annotations.  RNALfold-based folding of clusters will often miss very large inverted repeats that einverted can capture.  To make the .inv file, download and install the EMBOSS package ( http://emboss.sourceforge.net/ ), then run einverted against your genome of interest.  For the purposes of ShortStack analysis, the fasta file can be ignored / deleted.  The .inv file is used for ShortStack analysis.  Note that the inv_file is not required, but ShortStack will warn you if it is missing (unless you are running in 'nohp' or 'count' mode).

=head1 OUTPUT

=head2 Results.txt

This is a simple tab-delimited text file.  The first line begins with a "#" (comment) sign, and then lists column headers.  Each subsequent line describes the key traits of a single cluster.

To import this into R, here's a tip to deal with the first line, which has the headers but begins with a "#" character.

    >results <- read.table("Results.txt", head=TRUE, sep="\t", comment.char="")

Column 1: Locus : The genome-browser-friendly coordinates of the clusters.  Coordinates are one-based, inclusive (e.g. Chr1:1-100 refers to a 100 nt interval beginning with nt 1 and ending with nt 100).

Column 2: Name : Name of cluster.  Unless the run was in --count mode and the input file of a priori clusters already had names, the names are arbitrarily designated as "Cluster_1", "Cluster_2", etc.

Column 3: FlagOverlap : Name(s) of any loci from the flag_file that overlap with the cluster are listed.  If there are two or more, they are comma-separated.  If there were none, or no flag_file was provided, than a "." is present in this column instead.

Column 4: HP : Whether this cluster appears to be hairpin-derived or not.  If not, a "." is present.  If it is a hairpin, but NOT qualified as a MIRNA, "HP" is indicated.  MIRNAs are indicated by "MIRNA".  If the run was in "--nohp" mode, than all entries in the column will be "ND" (meaning 'not determined').

Column 5: Strand : The pre-dominant strand from which the small RNA emanate.  If ".", no strand was called.  HPs and MIRNAs always have a polarity, based on the hairpin's originating strand.  Non-HP clusters have their polarity determined by the --minstrandfrac setting.

Column 6: Frac_Wat : Fraction of aligned reads to the Watson (e.g. +) strand of the cluster.  1 means all were from Watson Strand (e.g. +), 0 means all were from Crick (e.g. -) strand.

Column 7: Total : Total aligned reads within the cluster.

Column 8: Uniques : Total aligned reads derived from uniquely mapped reads .. e.g., those with XX:i:1.

Column 9: UI : Uniqueness index. This is a value >0 and <=1. To determine the UI, the ambiguity normalized counts within the locus are summed and divided by the total number of alignments. For a given read, the ambiguity-normalized value is defined as 1 / d, where d is the integer in the XX:i:[d] tag representint the total number of possible locations for that read. Thus, loci with a low UI have a large proportion of reads that are ambiguous in their origins.

Column 10: DicerCall : If "N", the cluster was not annotated as dicer-derived, per options --dicermin, --dicermax, and --mindicerfrac.  Otherwise this is a number, within the --dicermin to --dicermax size range, which indicates the most abundant small RNA size within the mappings at that cluster.

Column 11: PhaseOffset : If "ND", phasing p-value was not calculated for this cluster.  Otherwise, the offset is the one-based genomic position with which the cluster appears to be "in-phase" (based on the 5' nt of a sense-mapped small RNA).  Phasing is always in increments identicial to the Dicer size call in column 9.

Column 12: Phase_pval :  If "ND", phasing p-value was not calculated for this cluster.  Otherwise, the p-value is derived from a modified hypergeometric distribution, as described below. 

Column 13: Pairs : Number of base-pairs in the hairpin stems.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 14 : FracPaired : Fraction of the stem nucleotides that are paired.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 15 : StemLength : Defined as 0.5 * (5' arm length + 3' arm length).  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 16 : LoopLength : Number of nucleotides between the 5' and 3' stems.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 17 : dGperStem : deltaG of the stems divided by the StemLength.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 18 : FracCovHP : Fraction of the per-nucleotide coverage present in the originally found cluster that is located in the two arms of the hairpin.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 19 : HPSizeResult : "1" indicates the number of pairs in the hairpin was less than or equal to maxmiRHPPairs; 0 indicates the opposite.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 20 : PrecisionResult : The number of small RNA sequences in the stem region of the hairpin that accounted for >= 20% of the mappings.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 21 : DuplexResult : The number of possible miRNA/miRNA* duplexes in which neither partner spanned a loop and neither partner had > maxmiRUnpaired number of unpaired nucleotides in the putative duplex.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 22 : StarResult : The number of putative miRNA*'s that were actually sequenced.  Values of 1 or more here indicate a MIRNA annoation.  If this is not a HP or MIRNA locus, "NA" is entered instead.

Column 23: Short : The total mappings from reads with lengths less than --dicermin, either in raw reads (--raw mode), or mappings per million mapped.

Column 24: Long : The total mappings from reads with lengths more than --dicermax, either in raw reads (--raw mode), or mappings per million mapped.

Columns 25 - the end : The total mappings from reads with the indicated lengths.  These are the sizes within the Dicer range.

=head2 Log.txt

This is a simple log file which records the information that is also sent to STDERR during the run.

=head2 gff3 files

Two gff3-formatted files are created, one for the 'DCL' loci (those with a DicerCall that is NOT N), and the 'N' loci (those with a DicerCall of 'N'). There are NOT produced in a --count mode run.

=head2 miRNA_summary.txt

This is a tab-delimited text file that summarizes key features of the loci annotated as MIRNAs, including mature miRNA sequences, miRNA-star sequences, and the numbers of mappings for each and for the entire locus.

=head2 Hairpin and MIRNA detail files

Unless the run was done in --nohp mode, each annotated hairpin-derived and MIRNA locus will have its own simple text file to display the details of the locus.  These text files all show A) the Name and genomic coordinates of the locus, B) the sequence, in RNA form, C) the identified hairpin structure, in dot-bracket notation, and D) all mappings whose start and stop is within the interval being examined.

Reads mapped to the sense strand (sense relative to the hairpin, not necessarily relative to the genome) have "."s as placeholders and are shown in the 5'-->3' orientation.  Reads mapped to the antisense strand (antisense relative to the hairpin, not necessarily relative to the genome) have "<"s as placeholders, and are written in the 3' --> 5' orientation.  Annotated mature miRNAs have "m"s as placeholders, and annotated miRNA*'s have "*"s as placeholders.

After each read, the read length (l) and the number of mappings (m) is shown.

Note that alignments for very complex HP loci are suppressed (only the sequence and structure will be displayed in the detail file, along with a note indicating the alignment was suppressed).  This has a major impact on reducing the memory footprint of ShortStack.  In addition, since the alignments are meant for visual inspection, very complex alignments aren't really parseable by eye anyway.  Complex loci are defined as those where the hairpin is longer than 400nts AND/OR has more than 400 distinct small RNA sequences.  All MIRNA loci have the alignments presented, regardless of whether it is complex or not.

=head1 KEY METHODS

=head2 Adapter trimming (FASTA or FASTQ)

Adapter trimming by ShortStack searches for the right-most occurence of the string given in the option --adapter. If found, the adapter is chopped off. Reads where no adapter is found are discarded, as are those shorter than 15 nts after trimming. Finally, any reads that contain non-ATGC characters after trimming are also discarded. Adapter-matching is performed by a simple regex and allows for no mismatches. If the input file is in FASTQ format, the quality values are trimmed as well.

=head2 Adapter trimming (Colorspace-FASTA)

Adapter trimming of colorspace data first translates the input adapter sequence to colorspace. The right-most occurence of the color-space adapter sequence is identified. Trimming begins 1 color BEFORE the match to the color-space adapter. This removes the 'hybrid' color corresponding to the dinucleotide at the junction of the small RNA and the adapter.  Reads where no adapter is found are discarded, as are those shorter than 15 colors after trimming. Finally, any reads that contain non-0123 colors after trimming are also discarded. Adapter-matching is performed by a simple regex and allows for no mismatches. If SOLiD formatted quality values were provided (via option --untrimmedCSQV), the quality values are also trimmed at the corresponding position.

=head2 Alignments

FASTA and FASTQ alignments are performed with the bowtie settings -v 1 -a --best --strata -S. Color-space alignments use bowtie settings -C -f --col-keepends (and if quality values are available, -Q) -v 1 -a --best --strata -S.  Initial alignments are parsed by ShortStack to retain just ONE randomly selected valid alignment per read. After this selection, and the addition of the XX tag (indicating the total number of possible positions for the read), the data are then piped through samtools view and samtools sort to create a sorted BAM formatted file. The header is also modified to note ShortStack processing. Note that unmapped reads are also retained in the alignment, with the SAM Flag set as 0x4 to indicate that.

One valid question is why force bowtie to report all valid alignments, and then parse to select just one, instead of explicitly telling bowtie to report just one alignment (e.g. k 1)? The reason is that bowtie's selection of reported positions in those cases is clearly NOT random. So even though ShortStack's method incurs a performance penalty during alignment, it will guarantee a more accurate alignment.

=head2 Multiple libraries

As of version 1.1.0, ShortStack supports input of multiple FASTA or FASTQ datasets (and as of version 1.2.0, multiple colorspace-FASTA datasets as well), either trimmed or untrimmed. These are passed in via the options --untrimmedFA, --untrimmedFQ, --untrimmedCS, --trimmedFA, --trimmedFQ, or --trimmedCS as comma-delimited lists.  If the input is untrimmed, option --adapter also allows input of a corresponding comma-delimited list of adapters.  If the data are colorspace, options --untrimmedCSQV and --trimmedCSQV allow input of comma-delimited lists of quality value files corresponding to their colorspace FASTA coutnerparts given by options --untrimmedCS or --trimmedCS. The untrimmed files and the adapters are assumed to be in the same order. If all untrimmed files have the same adapter, just one sequence for option --adapter is acceptable, and will be used to guide trimming for all files.

When multiple small RNA-seq libraries are input, each one is first aligned separately, creating temporary .bam files. When this is complete, they are merged into a single alignment, which is given the name "outdir.bam" in the working directory, where "outdir" is the string given in option --outdir. During the merging, the read group information is stored, so all alignments can be de-convoluted back to their parent libraries if desired. The individual .bam alignments created intially are deleted.

=head2 Read groups

As of version 1.1.0, ShortStack incorporates the option --read_group. When specified, only the alignments from the read group specified by the option will be used for analysis. Use of this option demands that the indicated read group is specified in the header of the relevant bam alignment file.

When an analysis uses a bam alignment file that contains more than one read group (based on the bam header), and the --read_group option was NOT used in the run, the analysis will conclude with a --count mode analysis of each read group separately. This is meant to be convenient for analyses in which a de-novo small RNA gene annotation is performed using a merger of multiple libraries, followed by quantification of each locus for each small RNA-seq library separately .. this should facilitate the downstream analysis of differential expression, for instance.

=head2 de novo Cluster Discovery

Cluster discovery proceeds in two simple steps:

1. The total depth of small RNA coverage at each occupied nucleotide in the genome is examined, and initial 'islands' of coverage are defined as continuous stretches where the read depth is greater than or equal to the threshold depth specified by option --mindepth.  Note that islands could theoretically be as small as one nucleotide, since they depend on total depth of coverage, not the small RNA length per se.  Many islands may often be 20-24nts in length, corresponding to a pile of a single small RNA species.

2. The initial islands are then temporarily extended on both sides by the distance specified by option --pad.  Islands that overlap after extension are merged.  The "dangling pads" at the ends of the merged clusters are then removed.  After all extensions, resultant mergers, and end trimmings are performed, the final result is the initial clusters.  If the run is performed in --nohp mode, these are the final clusters.  If hairpins and MIRNAs are being examined, some of the clusters may be adjusted in position to fully capture the putative hairpin(s) (see below).

=head2 Hairpin and MIRNA analysis

1. Clusters are first filtered to determine eligbility for secondary structure analysis.  Only clusters whose uniqueness index* is >= [--minUI] are eligible.  In addition, clusters with a DicerCall of "N" are excluded from folding analysis.

2.  The genomic window to be subject to RNA folding is first determined.  If the locus size is > --maxhpsep, no RNA folding will take place at the locus.  Otherwise, a window with length of --maxhpsep is centered on the locus to determine the nucleotides to fold

3.  Both the top and bottom genomic strands from the window are then subjected to secondary structure prediction using RNALfold (option -L [--maxhpsep]), which returns a diverse set of often overlapping predicted structures.

4.  The structures are parsed, retaining only those that satisfy options --minfracpaired, --minntspaired, --maxLoopLength, amd --maxdGperStem.

5.  If an .inv file was provided, all inverted repeats in that file are parsed, and then filtered to also satisfy options --minfracpaired, --minntspaired, --maxLoopLength, amd --maxdGperStem.  In addition, loop lengths from einverted data are not allowed to be longer than 50% of the helix length of the putative hairpin.  Putative RNA secondary structures in dot-bracket notation are generated from the .inv alignment, not by actual RNA thermodynamic analysis.  All G-U alignments are considered paired, in addition to the standard A-U and G-C pairings.  Both strands are used, subject to passing the --minfracpaired, --minntspaired, --maxLoopLength, amd --maxdGperStem criteria.  Inverted-repeats that survive these filters are then filtered to retain only those with overlap to the original set of structure-eligible clusters (see 1 above), and the resulting set of eniverted-derived hairpins is merged with the RNALfold-derived set.

6.  Redundant hairpins are then removed.  Redundant hairpins are those whose 5' arms and 3' arms overlap.  In pairwise comparisons of redundant hairpins, the longest hairpin is retained.

7.  Hairpins that don't have overlap with the original cluster are then removed.  Because the folding window could have been extended around the cluster, there could be  putative hairpins that are not within the original cluster.  To have overlap, at least one of the hairpin's helical arms must have at least 1nt within the original cluster coordinates.

8. The pattern of small RNA expression relative to the remaining hairpins is then examined.  The per-nucleotide coverage across every base, on both strands separately, across the original locus coordinates is calculated.  If there is a single hairpin whose 5' and 3' arms contain >= [--minfrachpdepth] of the total coverage of the original locus, the hairpin is kept for futher analysis.  If more than one hairpin meets this criterion, than the one with the highest coverage fraction in the arms is retained.  Note that this step contains a correction for reads that are "dyads" .. reads that map twice to a hairpin, once in each arm, on opposite strands... this happens for perfect inverted repeat loci.  Such reads are counted towards the sense strand only for a given hairpin.

9. The pattern of small RNA expression relative to the single hairpin candidate is further scrutinized for polarity.  The fraction of all mappings in the hairpin interval must be >= [--minstrandfrac].  As in step 8, this step corrects for "dyad" reads (see above).  Hairpins that pass this step are either HP or MIRNA loci.  The coordinates of the originally determined de novo locus are discarded, and replaced with the hairpin coordinates.

10.  Each potential hairpin that remains is next analyzed to see if it qualifies as a MIRNA.  MIRNA locus annotation is designed to satisfy the criteria for de novo annotation of plant MIRNAs as described in Meyers et al. (2008) Plant Cell 20:3186-3190. PMID: 19074682.  In fact, ShortStack's criteria is a little stricter than Meyers et al., in that ShortStack has an absolute requirement for sequencing of the exact predicted miRNA* sequence for a candidate mature miRNA.  It is important to note that ShortStack's MIRNA annotation method is designed to reduce false positives at the expense of an increased rate of false negatives.  In other words, there are likely many bona fide MIRNA loci not annotated as such because they don't quite meet the strict criteria set forth below. 

- Hairpin Size:  The total number of pairs in the hairpin must be <= [--maxmiRHPPairs]

- Precision: There must be at least one candidate mature miRNA that comprises at least 20% of the total abundance of small RNAs mapped to the hairpin.  

- Duplex:  Both partners in a candidate miRNA/m, for a cluster located at Chr1:1000-2000, reads mapped to 980-1000, 1100-1123, and 2000-2021 are all counted as being within the cluster during quantification.  Note that it's possible to count the same mapping within non-overlapping clusters.

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

Note: P-values are not corrected for multiple-testing.  Consider adjustment of p-values to control for multiple testing (e.g. Bonferroni, Benjamini-Hochberg FDR, etc) if you want a defensible set of phased loci from a genome-wide analysis.

=cut


