LICENSE
    ShortStack

    Copyright (C) 2012-2018 Michael J. Axtell

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

SYNOPSIS
    Alignment of small RNA-seq data and annotation of small RNA-producing
    genes

CITATIONS
    If you use ShortStack in your work, please cite one of the following:

  VERSIONS 3.x and higher
    Johnson NR, Yeoh JM, Coruh C, Axtell MJ. (2016). G3 6:2103-2111.
    doi:10.1534/g3.116.030452

  OLDER VERSIONS
    Axtell MJ. (2013) ShortStack: Comprehensive annotation and
    quantification of small RNA genes. RNA 19:740-751.
    doi:10.1261/rna.035279.112

    Shahid S., Axtell MJ. (2013) Identification and annotation of small RNA
    genes using ShortStack. Methods doi:10.1016/j.ymeth.2013.10.004

INSTALL
  Dependencies
    All dependencies must be executable and findable in the user's PATH

    perl (version 5.x) : Generally installed in linux and mac machines by
    default. Expected to be installed at /usr/bin/perl

    samtools (version 1.x or higher) : Free from http://www.htslib.org/

    bowtie (if aligning) : Free from
    http://bowtie-bio.sourceforge.net/index.shtml .. note: requires bowtie
    ... NOT bowtie2 !

    bowtie-build (if aligning and .ebwt indices not found) : Free from
    http://bowtie-bio.sourceforge.net/index.shtml

    gzip (if aligning) : Generally installed in linux and mac machines by
    default.

    RNAfold (unless running with --nohp option to disable MIRNA search) :
    Part of the Vienna RNA package, Free from
    http://www.tbi.univie.ac.at/RNA/

    Test environment
    This release of ShortStack (3.8.4) tested on Mac OSX (10.10.5), perl
    5.18.2, samtools 1.3.1, bowtie 1.2.0, RNAfold 2.3.2. It is known that
    samtools 1.x and higher is critical (no old 0.x samtools allowed).

  Install
    There is no real installation of ShortStack. Make sure it is executable.
    For convenience, can be added to your PATH. It expects your perl
    installation to be at /usr/bin/perl.

USAGE
    Usage: ShortStack [options] {--readfile <r> | {--bamfile <b> |
    --cramfile <c>}} --genomefile <g>

    <r> : readfile must be in fasta (.fasta or .fa), colorspace-fasta
    (.csfasta), or fastq (.fastq or .fq) format, or their gzip-compressed
    versions (.fasta.gz, .fa.gz, .csfasta.gz, .fastq.gz, or .fq.gz) Can also
    be a list (seperated by spaces) of several read files.

    <b> : BAM formatted alignment file (.bam).

    <c> : CRAM formatted alignment file (.cram).

    <g> : FASTA formatted (.fa or .fasta) genome file.

TEST
    Test data and brief instructions are available at
    http://axtelldata.bio.psu.edu/data/ShortStack_TestData/

OPTIONS
    Note that we have done our best to set default settings for all options
    that are best for most users.

  General Options:
    --help : print a help message listing all options and quit

    --version : print version and quit

    --genomefile [string] : path to reference genome in .fasta or .fa
    format. Required for any run.

    --outdir [string] : name of output directory to be created for results.
    Defaults to 'ShortStack_[time]', where [time] is the current UNIX time
    according to the system. If the directory already exists, ShortStack
    exits with an error message.

  Alignment Options:
    --readfile [string] : path to readfile(s) to be aligned. valid formats:
    .fasta, .fa, .fasta.gz, .fa.gz, .fastq, .fq, .fastq.gz, .fq.gz,
    .csfasta, .csfasta.gz. Multiple files, can be specified as separate
    arguments to --readfile ... e.g. --readfile file1.fastq file2.fastq
    file3.fastq Mutually exclusive with --bamfile or --cramfile.

    --adapter [string] : sequence of 3' adapter to trim off during read-pre
    processing. Must be at least 8 bases, with only ATCG characters. If not
    specified, reads are assumed to be already trimmed.

    --bowtie_cores [integer] : Argument to be passed to bowtie's -p option,
    specifying number of processor cores to request during alignment.
    Defaults to 1. Must be an integer of 1 or more.

    --sort_mem [string] : Argument to be passed to samtools sort -m option,
    which sets the maximum memory usage during bam file sorting. If not set,
    samtools sort defaults it to 768M. Higher settings will reduce the
    overall time spent in alignment phase, at cost of more memory usage. Use
    K/M/G suffixes to specify kilobytes, megabytes, and gigabytes,
    respectively. Extremely large alignment jobs will crash (due to crash of
    samtools sort operation) if --sort_mem is not set high enough. However,
    alignment jobs will also crash if sort_mem is set too high, and all
    physical memory on your machine is exahusted.

    --mismatches [integer] : Argument to be passed to bowtie's -v option,
    specifying number of mismatches to be tolerated in a valid alignment.
    Must be either 0, 1, or 2. In cases of multiple hits, only hits with
    lowest number of mismatches kept. Default: 1.

    --cquals [string] : path(s) to color-space quality value file(s). Used
    only in conjunction with .csfasta or .csfasta.gz formatted files in
    --readfile. Compressed format for cquals is NOT allowed. Like
    --readfile, cquals can take multiple arguments for multiple files, e.g.
    --cquals file1.qual file2.qual file3.qual

    --cram : When aligning, convert final alignment to cram format instead
    of the default bam format.

    --mmap [string] : Protocol for handling multi-mapped reads. Valid
    entries are n (none), r (random), u (unique- seeded guide), or f
    (fractional-seeded guide). default: u

    --bowtie_m [string] : Setting to be passed to the -m option of bowtie.
    Over-ridden and set to 1 if option mmap is set to n. This sets the
    maximum number of multi-mappings allowed. Valid settings are integers >=
    1 OR set 'all' to disable suppression of highly multi-mapped reads.
    Default: 50

    --ranmax [string] : Reads with more than this number of possible
    alignment positions where the choice can't be guided by unequal will be
    reported as unmapped. Irrelevant if option mmap is set to n or r. Must
    be integer of 2 or greater or set to 'none' to disable. Default: 3.

    --align_only : If this switch is present, the ShortStack run will
    terminate after the alignment phase with no analysis performed.

    --show_secondaries : If this switch is present, the output alignment
    file will contain secondary alignments as well as primary alignments for
    multi-mapped reads. Secondary alignments have bit 256 set in the SAM
    FLAG field. This option can increase alignment file size, sometimes by a
    lot.

    --keep_quals : As of version 3.5, by default ShortStack alignments no
    longer store the quality values, to save space. Use of this switch will
    cause quality values to be retained. Note that this increases file size.

  Analysis Options:
    --bamfile [string] : path to input .bam alignment file of small RNAs.
    Only lines with bits 4 and 256 unset will be used. Mutually exclusive
    with --readfile or --cramfile.

    --cramfile [string] : path to input .cram alignment file of small RNAs.
    Only lines with bits 4 and 256 unset will be used. Mutually exclusive
    with --readfile or --bamfile.

    --dicermin [integer] : Minimum size of a Dicer-processed small RNA. Must
    be an integer of at least 15 and <= dicermax. Default: 20.

    --dicermax [integer] : Maximum size of a Dicer-processed small RNA. Must
    be an integer of at least 15 and >= dicermin. Deafult: 24.

    --foldsize [integer] : Size of genomic RNA segments for folding during
    MIRNA search. Any loci larger than this size will not be analyzed with
    respect for MIRNA features. Must be an integer of at least 200 and no
    larger than 1,000. Default: 300. Note that increasing this setting may
    drastically increase runtimes.

    --locifile [string] : Path to a tab-delimited plain-text file listing
    intervals to analyze. Lines starting with # are ignored. First column is
    coordinate in format Chr:start-stop, second column is names (optional),
    and any other columns are ignored. Mutually exclusive with option
    --locus.

    --locus [string] : Analyze the specified interval(s). Interval(s) is
    specified in format Chr:start-stop. Multiple intervals can be specified
    in a comma-separated list. Mutually exclusive with option --locifile.

    --nohp : Disable MIRNA search.

    --pad [integer] : Initially found clusters of small RNAs will be merged
    if the distance between them is less than or equal to the value of pad.
    Must be an integer between 0 and 50000. Default: 75.

    --mincov [string] : Clusters of small RNAs must have at least this many
    alignments. Supply an integer between 1 and 50000. Can also be a
    normalized value in reads per million (rpm) OR reads per million mapped
    (rpmm). When specifying mincov in rpm or rpmm, the mincov value must be
    a floating point number > 0 and < 500,000 followed by the string 'rpm'
    or 'rpmm'. Examples: '5' --> threshold is 5 raw reads. '3.2rpm' -->
    threshold is 3.2 reads per million mapped. '2.8rpmm' --> threshold is
    2.8 reads per million mapped. Deafult: 0.5rpm.

    --strand_cutoff [float] : Cutoff for calling the strandedness of a
    locus. Must be a floating point number between 0.5 and 1 (inclusive).
    DEFAULT: 0.8. At default of 0.8, a locus must have 80% of more of its
    reads on the top strand to be called a + strand locus, or 20% or less on
    the top strand to be a - strand locus. All others receive no strand call
    (e.g. '.'). Only stranded loci are analyzed for MIRNAs, while only
    unstranded loci are analyzed with respect to phasing. Most users
    probably want to use the default setting of 0.8.

    --total_primaries [integer] : Tell ShortStack the total number of
    primary alignments in the bam file. Specifying this value here speeds
    the analysis, since ShortStack does not need to count the reads directly
    from the bam file. Can only be specified in conjunction with --bamfile.
    This count should include all primary alignment INCLUDING unplaced ones.

SYSTEM RECOMMENDATIONS
    ShortStack was developed on Apple Mac OSX devices running 10.9 or 10.10.
    It has also been tested on Linux (CentOS and Ubuntu).

    At least 4G memory is suggested. Alignment and building bowtie indices
    tend to be the most memory-intensive portions for a given run, and
    memory usage seems to scale with genome size, but not as much with the
    number of small RNAs.

    Alignments benefit from multiple processing cores, via the
    --bowtie_cores option. All other portions are single-threaded.

    Alignment speed may also be increased using the --sort_mem option to
    increase the memory used for bam file sorting. Setting a higher
    --sort_mem will be REQUIRED for very large alignment runs to avoid
    samtools sort crashes due to too many open files.

    At least 50G of hard disk space is recommended to be available, due to
    the sometimes large size of the temporary alignment files and the final
    alignment file. Extreme settings for options --bowtie_m and --ran_max
    may cause creation of extremely huge files.

    The total time of analysis depends on several factors, including most
    prominently genome size, number of reads analyzed, whether or not bowtie
    indices need to be created, whether or not MIRNAs are being analyzed,
    and of course your equipment. Excluding building bowtie indices, we
    generally have observed run times for alignment + analysis runs to take
    between 20 minutes and 10 hours using default ShortStack settings.

ALIGNMENT METHODS
    If ShortStack is given the --readfile option, alignments of the reads
    will be performed. Specifying --readfile is mutually exclusive with both
    --bamfile or --cramfile

  Details of alignment methods and performance testing
    For full details on ShortStack's alignment methods and the results of
    performance testing, see Johnson et al. (2016) G3 6:2103-2111.
    doi:10.1534/g3.116.030452.

  Genome pre-processing
    Genome file format and naming
    All runs require a reference genome in FASTA format, specified with the
    --genomefile option. The file must end with a valid suffix .. either .fa
    or .fasta.

    Within the genome, if the name of a chromosome has whitespace
    characters, the name will be trimmed at the first whitespace character.

    Genome stitching
    If the reference genome has > 50 chromosomes/scaffolds/contigs, and the
    genome N50 length is < 1Mb, and MIRNAs are to be analyzed (e.g., --nohp
    was NOT specified), then ShortStack will 'stitch' the small chromosomes
    together to make fewer but larger chromosomes. This can drastically
    improve performance during MIRNA searching for highly fragmented genome
    assemblies. As of ShortStack 3.8, stitching has no effect on the results
    (e.g. results are reported relative to the original genome, not the
    stitched one).

    Genome indexing
    If not detected, an index of the genome will be created using samtools
    faidx.

    bowtie indices
    If not detected, bowtie-build, using all default settings, will be
    invoked to create the required six .ebwt indices of the genome. This can
    be time-consuming, and memory intensive.

  Reads pre-processing
    Reads file formats
    Small RNA reads to be aligned must be in fasta, fastq, or csfasta
    formats, or their gzip-compressed versions. File names must end with
    .fa, .fasta, .fastq, .fq, .csfasta, .fa.gz, .fasta.gz, .fastq.gz,
    .fq,gz, or .csfasta.gz.

    Multiple readfiles can be specified with option --readfile by separating
    the file names/paths with commas. Colorspace reads cannot be mixed with
    base-space reads; otherwise, mixed file formats are ok.

    If you wish to also include quality values from SOLiD data, the _QV.qual
    file(s) can be passed in through the --cquals option. Color-space
    quality values are NOT accepted in .gz compressed format.

    No paired-end support
    There is no support for paired-end reads in ShortStack. Small RNA data
    are assumed to be single-ended, and represent the 5'-->3' cDNA sequences
    of cloned RNAs.

    No condensation
    Input reads are expected to be de-condensed. That is, if a small RNA was
    sequenced 10,000 times in a run, there should be 10,000 entries, each
    with a different header name, in the input readfile. In other words,
    ShortStack is designed to take reads right off the sequencer without any
    other pre-processing (except adapter trimming .. see below).

    Unique read names required
    The small RNA reads must all have unique names within a given file. If
    this requirement is not met, alignments will be completely unreliable
    due to errors in interpreting and handling of multi-mapped reads.

    Adapter trimming
    ShortStack has a primitive 3'-adapter capability. Specify an adapter of
    at least 8nts in length with option --adapter. If nothing is given to
    --adapter, ShortStack assumes your reads are already trimmed. Trimming
    simply looks for the right-most exact match to the given apdater
    sequence, and when found, chops it off. If a read is smaller than 15nts
    after trimming, it is discarded. For more sophisticated adapter
    trimming, consider cutadapt or trimmomatic

    If quality values are present, they are trimmed as well.

  Alignment overview
    ShortStack uses bowtie to align reads. It first aligns, and processes
    the output on the fly to note how many equally good alignment positions
    were found for each read. It then uses this information in a second
    phase to 'decide' on the most likely 'correct' location for multi-mapped
    reads. The final output is a single .bam or .cram formatted alignment
    file. If multiple readfiles were input, the final bam or cram file notes
    the origin of each read with the RG tag (see sam format specification).

    mismatches
    By default, ShortStack allows up to 1 mismatch for a valid alignment.
    This helps with sequencing errors and SNPs. If a read has some
    alignments with 0 mismatches, and some with 1, only those with 0
    mismatches are kept. The option --mismatches controls this threshold,
    and can be set to 0, 1, or 2.

    *** WARNING : If the genome is large (.ebwtl bowtie indices are made,
    instead of .ebwt), there is a serious bowtie bug that has yet to be
    resolved involving the --best option.
    http://sourceforge.net/p/bowtie-bio/bugs/343/ . To get around this, when
    aligning to a 'large' reference, ShortStack forces the number of allowed
    mismatches to be 0.

    Control of multi-mappers with --bowtie_m
    In general, we find it's not worth the time or effort to deal with
    'extreme' multimapping reads. The --bowtie_m setting determines the
    threshold of 'extreme' multi-mappers. Reads that have more than
    --bowtie_m alignments are simply marked as unmapped. The default setting
    is 50. Valid settings are >=1 or set 'all' to disable any suppression of
    extreme multi-mappers (not suggested).

    Optimal placements of multi-mapped reads
    For multi-mapped reads that have between 2 and --bowtie_m number of
    equally good alignments, ShortStack has several methods to decide on the
    true origin of the read. The choice of method is specified with the
    option --mmap. The methods are:

    u: Placement guided by uniquely mapping reads. During the alignment, the
    count of uniquely mapped reads is kept in 50nt bins across the reference
    genome. The bin location is determined by the left-most coordinate of
    the uniquely mapped read. After the first phase of alignment for all
    reads (in all files) has completed, this genome-wide map of
    uniquely-mapped read counts is used to guide the decisions of the most
    likely locations of multi-mapped reads. Specifically, for a given
    multi-mapped read, the local count of uniquely mapped reads at each
    possible location is computed. The local count is that of the specific
    50nt bin the alignment lies in (again, by left-most positon) plus the
    counts of the 2 bins upstream and 2 bins downstream. All of the local
    counts are converted to fractions of the sum of all total counts. These
    fractions are then used as the probabilities of placement for the
    multi-mapped read. For instance, suppose a multi-mapped read had three
    possible positions. The read counts of uniquely mapped reads were 30,
    65, and 5. This would mean that read has a 30%, 65%, and 5% chance,
    respectively, of being assigned to each bin. The actual choice is
    probabilistic, given the computed weightings, for each read.

    f: Placement guided by all mapped reads. Like u, except that
    multi-mapped reads also contribute to the guidance densities. All reads
    contribute 1/(n of alignment positions) to each 50nt bin that the occur
    in.

    r: Placement is simply random. This is faster than u and f, but performs
    much more poorly at properly placing multi-mapped reads. Achieves high
    sensitivity, but very low precision.

    n: Multi-mapped reads are all ignored and marked as unmapped. Very fast,
    but ignores large quantities of data. Achieves high precision, but very
    low sensitivity.

    The default setting for --mmap is u

    ranmax
    When running mmap method u or f, there are some cases where no guidance
    can be given, and so the choice on where to put a multi-mapped read is
    still random. In those cases, the option ranmax will suppress any
    alignment where the choice is 'too' random. By default, --ranmax is set
    at 3, so that if a read can't be placed confidently, no placement is
    done if there are more than 3 choices.

  Alignment output format
    Final alignments are sorted bam or cram formatted alignments. bam is the
    default, while cram is created if the option --cram is set. The
    alignment file conforms to all SAM/BAM/CRAM format specifications, and
    has the following features:

    Headers contain @RG lines to describe each read-group (input readfile).

    For multi-mapped read alignments that were NOT chosen as the most likely
    alignment, bit 256 (secondary alignment) is set in the FLAG. For such
    lines, the SEQ and QUAL values are not stored, to save space. The SEQ
    and QUAL for multi-mapped alignments are kept only in the primary
    (chosen) alignment for the read.

    XX:i tags: Added by ShortStack to each line, this indicates the total
    number of valid alignments found for the read.

    XY:Z tags: Added by ShortStack to each line, this indicates how the
    reported alignment was selected: U: Uniquely mapped, P: Multi-mapped and
    placed based on probabilities calculated by mmap method u or f, R:
    Multi-mapped and randomly placed, M: Multi-mapped but marked as unmapped
    becuase the number of alignment positions exceeded --bowtie_m, O:
    Multi-mapped but marked as unmapped because no guidance possible and
    choices exceeded setting --ranmax, N: Unmapped because 0 valid
    alignments found in genome.

    XZ:f tags: Added by ShortStack to each line, this indicates the
    calculated probability of placement for the read.

BAM AND CRAM FILE REQUIREMENTS
    Existing alignments can be provided to ShortStack using the --bamfile or
    --cramfile options (for bam formatted and cram formatted alignments,
    respectively). --bamfile and --cramfile are mutually exclusive with each
    other, and with --readfile.

    Any properly formatted bam or cram file should work with ShortStack,
    subject to the requirements below. However, for best performance, it is
    recommended to use ShortStack for alignments as well.

    Requirements for input bam or cram files:

    1. Header must be present, and contain @SQ lines.

    2. File most be sorted.

    3. Read groups will not be recognized unless they are properly noted in
    the header.

    4. Paired-end, spliced, clipped, padded, or gapped alignments will be
    ignored, with a warning to the user. Reads marked as secondary
    alignments (bit 256 set in the FLAG) will not be used.

DE-NOVO CLUSTER FINDING
    Unless options --locus or --locifile are used (see below), ShortStack
    will de-novo identify clusters of small RNA accumulation genome-wide.
    Cluster definition is simple: First, all regions containing at least one
    primary alignment are found where the maximum distance between the ends
    of the alignments is <= option --pad (default: 75). Second, if the
    number of alignments in the cluster is >= option --mincov (default:
    0.5rpm), the cluster is kept. The mincov threshold can also be specified
    in terms of reads per million by using a value such as 3.2rpm (which
    specifies the threshold to be 3.2 reads per million). Using a rpm
    threshold allows the sensitivity of cluster discovery to be normalized
    between libraries of different sizes. Alternatively, reads per million
    mapped (rpmm) can be specified: A --mincov of 1.2rpmm indicates 1.2
    reads per million mapped is the threshold. rpm is a fraction of total
    library size, while rpmm is a fraction of only the aligned & placed
    fraction of the library.

UNPLACED SMALL RNAS
    As of version 3.6 and higher, de-novo identification of small RNA
    clusters also will include reporting of unplaced small RNAs ... small
    RNAs that were not placed on the reference genome. Only small RNAs with
    an abundance higher than the limit set by option --mincov are reported.
    These small RNAs typically inlcude RNAs that could not be aligned
    anywhere on the reference, as well as multi-mapped RNAs where ShortStack
    did not choose a alignment position for (see alignment methods).

USER-SPECIFIED CLUSTERS
    Users can supply specific loci to analyze in two ways. For just one or a
    few loci, the option --locus can be used. The argument should be a
    coordinate in the format Chr:Start-Stop. Multiple loci can be specified
    in option --locus by using commas to delimit them.

    For larger lists of user-defined loci, and external file can be used
    instead, specified with option --locifile. The file is a plain-text ,
    tab-delimited format. The first column should list the coordinates in
    Chr:Start-Stop format. An optional second column can contain names of
    the loci. Any other columns will be ignored. Also, lines starting with #
    will be ignored.

    Options --locus and --locifile are mutually exclusive. Also, if either
    is used, no de- novo cluster finding occurs.

ANALYSIS METHODS
    Regardless of whether the small RNA clusters were de-novo discovered or
    user-defined, the analysis methods of each cluster are the same. The
    major methods are described below:

  Read-group-specific counts
    Quantification occurs separately for each read-group in the alignment.
    The results are in the 'Counts.txt' file, which has the observed number
    of reads, the mean number of reads for the ten re-samplings, and the
    standard deviations. When there are multiple read-groups, this is
    helpful to gather the raw data for differential expression analysis.

    There is always a read-group called 'main', which is all read-groups
    combined.

  Strandedness of loci
    Loci where >= 80% of the primary alignments are on the top genomic
    strand are noted with a strand of +. Loci where <=20% of the primary
    alignment are on the top genomic strand are noted with a strand of -.
    All other loci are given a strand of .

  Major RNA
    For each locus analyzed, the single most abundant RNA is noted and the
    sequence reported. In cases where there is a tie, the reported major RNA
    is chosen arbitrarily from among the tied RNAs.

  Complexity
    Complexity is a metric that varies from >0 to 1. It is calculated as (n
    distinct alignments) / (abundance of alignments), thus lower values
    indicate loci dominated by just a few dominant RNAs, while higher values
    indicate loci with more diverse sets of small RNAs.

  DicerCall
    The 'DicerCall' reflects the predominant RNA size observed in the locus.
    However, if < 80% of the reads in a locus are NOT within the bounds
    described by the options --dicermin and --dicermax, then the DicerCall
    is 'N' instead. DicerCalls of N usually reflect loci where the small
    RNAs are NOT related to an RNAi process ... most often, breakdown
    products of abundant RNAs.

    Note that the predominant RNA size need not be a majority .. for
    instance a locus with 40% 21 mers, 39% 22 mers, and 21% 23 mers would
    have a DicerCall of 21.

  MIRNAs
    ShortStack's MIRNA analysis is meant to eliminate false positives. It
    therefore probably allows some degree of false negatives (e.g., loci
    that really are MIRNAs but are not annotated as such). MIRNA analysis in
    ShortStack version 3 is a step-wise process. If a locus fails a certain
    step, it is removed from consideration and given a specific code
    indicate what step it failed. The codes are below in step-wise order.
    The Major RNA is always hypothesized to be the mature miRNA in the
    locus.

    Note that MIRNA analysis is limited to loci that are <= the length
    specified by option --foldsize ... the default setting is 300 nts.
    Increasing this size may allow you to find more MIRNAs, but will also
    slow down the runtime of the process.

    MIRNA analysis codes:

    N0: not analyzed due to run in --nohp mode.

    N1: no reads at all aligned in locus

    N2: DicerCall was invalid (< 80% of reads in the Dicer size range
    defined by --dicermin and --dicermax).

    N3: Major RNA abundance was less than 2 reads.

    N4: Major RNA length is not in the Dicer size range defined by
    --dicermin and --dicermax.

    N5: Locus size is > than maximum allowed for RNA folding per option
    --foldsize (default is 300 nts).

    N6: Locus is not stranded (>20% and <80% of reads aligned to top strand)

    N7: RNA folding attempt failed at locus (if occurs, possible bug?)

    N8: Strand of possible mature miRNA is opposite to that of the locus

    N9: Retrieval of possible mature miRNA position failed (if occurs,
    possible bug?)

    N10: General failure to compute miRNA-star position (if occurs, possible
    bug?)

    N11: Possible mature miRNA had > 5 unpaired bases in predicted precursor
    secondary structure.

    N12: Possible mature miRNA was not contained in a single predicted
    hairpin

    N13: Possible miRNA/miRNA* duplex had >2 bulges and/or >3 bulged nts

    N14: Imprecise processing: Reads for possible miRNA, miRNA-star, and
    their 3p variants added up to less than 50% of the total reads at the
    locus.

    N15: Maybe. Passed all tests EXCEPT that the miRNA-star was not
    sequenced. INSUFFICIENT evidence to support a de novo annotation of a
    new miRNA family.

    Y: Yes. Passed all tests INCLUDING sequencing of the exact miRNA-star.
    Can support a de novo annotation of a new miRNA family.

    For loci where MIRNA analysis returns a Y (yes) result, a plain-text
    summary of the locus and its secondary structure is found in the MIRNAs
    directory.

    Users should be aware that sometimes ShortStack will annotate known
    miRNA-stars as miRNAs, if the abundance of the miRNA-star in the
    analyzed library is higher.

    MIRNA analysis can be disabled with the --nohp option. This may save
    significant analysis time.

    As of ShortStack version 3.x, MIRNA analysis is geared toward plant
    MIRNAs. It probably is just fine for animal MIRNAs too, but has not been
    tested on them.

  Phasing
    Phasing here refers to a signature of periodicity of small RNA abundance
    that reflects dsRNA processing from a defined terminus. Phased siRNA
    clusters often are triggered by an upstream small RNA-mediated clevage
    event which causes RNA-dependent RNA polymerase activity and subsequent
    siRNA production from the terminus defined by the cleavage event.

    For valid loci, ShortStack 3.7 and above uses a modified version of the
    formula described by Guo et al. (2015) (doi:
    10.1093/bioinformatics/btu628), S = PR * PN * ln(1 + PArpm), where S is
    the phase score, PR is the phase ratio (see Axtell 2010 doi:
    10.1007/978-1-60327-005-2_5), PN is the number of distinct sequences
    that are phased, and PArpm is the abundance of the phased reads in units
    of reads per million.

    ShortStack calculates the phase score in a 21 nt phase size for loci
    with a DicerCall of 21, or in a 24 nt phase size for loci with a
    DicerCall of 24, and returns the score. Higher phasing scores indicate
    more phasing signature. Phase scores range from very near 0 (worst) up.

    The modification of the Guo et al. formula, first implemented in
    ShortStack version 3.7, makes the PhaseScore numbers comparable between
    different libraries. A score of ~30 or more indicates a well-phased
    locus.

    Not all loci are subject to phasing analysis. Loci with no reads at all
    aligned, a DicerCall of anything except 21 or 24, a Locus Size of < 3 *
    DicerCall, and stranded loci (>= 80% of reads on top strand OR <= 20% of
    reads on top strand) are not analyzed. These are assigned a PhaseScore
    of -1.

OUTPUT FILES
    All output files are in the directory created by ShortStack, whose name
    is specified by option --outdir The exceptions are the .fai file (genome
    index file) created if it is not present and the six ebwt bowtie index
    files that are created if not present ... these are all put in the same
    location as the input genome file.

  Log file
    A log of the run messages is created and written to Log.txt. It is the
    same as the messages printed to STDERR during the run.

  ErrorLogs
    For debugging. Most users won't need to look at this. It stores the
    verbose outputs of bowtie-build, bowtie, samtools sort, and samtools
    merge commands that are not kep in the main Log. Sometimes these are
    helpful in diagnosing bugs, particularly samtools sort and merge bugs
    due to memory issues.

  Stitched genome file
    If the input genome was 'stitched' (see above), the stitched genome file
    will be put in the ShortStack outdir, along with its fai index
    temporarily during the run. Both files will be deleted at the end of the
    run so you won't see them unless your run was killed before completion
    for some reason.

  Results file
    The file Results.txt is a plain-text tab-delimited file that contains
    the core results of the analysis. The columns are labeled in the first
    row, and are:

    1. Locus: Coordinates of the locus in format Chr:Start-Stop

    2. Name: Name of the locus

    3. Length: Length of the locus (nts)

    4. Reads: Total number of primary alignments in the locus

    5. RPM: Total number of primary alignments normalized to reads per
    million. Note the the normalization factor includes all primary
    alignments .. both mapped and unmapped.

    6. UniqueReads: Number of uniquely aligned primary alignments in locus.

    7. FracTop: Fraction of primary alignments aligned to the top genomic
    strand

    8. Strand: Strand call for the locus

    9. MajorRNA: Most abundant RNA at locus. In cases of tie, MajorRNA is
    arbitrarily chosen from the tied entries.

    10. MajorRNAReads: Number of primary alignments for the MajorRNA.

    11. Complexity: A number >0 and <= 1 that reflects the complexity of
    small RNA production from the locus. Defined by
    (n_distinct_read_sequences) / (abundance of all reads). Lower numbers
    indicate loci that are more dominated by a single highly abundant RNA.

    12. DicerCall: If >= 80% of the primary alignments were reads >=
    dicermin and <= dicermax, DicerCall is a number that indicates the
    predominant size of the RNA population from the locus. If the 80%
    threshold was not met, then DicerCall is N instead. Can also be NA if
    the locus had no aligned reads.

    13. MIRNA: Results of MIRNA analysis. Codes starting with N indicate not
    a MIRNA, Y means yes. See above for full description of codes.

    14. PhaseScore: Phasing score for a phase size of 21 or 24nts according
    to a modified version of equation 3 of Guo et al (2015) doi:
    10.1093/bioinformatics/btu628. If the locus had a DicerCall of 21, phase
    score is for a 21 nt phasing register. If the locus had a DicerCall of
    24, the phase score is for a 24 nt phasing register. See above for full
    description of phasing analysis.

    15. Short: Number of primary alignments that were shorter than
    --dicermin

    16. Long: Number of primary alignments that were longer than --dicermax

    17-end: Number of primary alignments of the indicated RNA size.

  Unplaced file
    The Unplaced.txt file is plain text, tab-delimited. It shows each
    unplaced small RNA whose abundance was >= the limit set by --minocv.
    This file is only created in a de-novo run. RNAs are sorted first by
    length (ascending) and then by ASCII (ascending). Columns show the
    sequence, its length, the total number of reads, the reads per million
    (RPM) and the number of equally valid genome alignments.

    In some cases the number of genome alignments may not be able to be
    found. An entry of '?' indicates that the number of hits is unknown.
    This will occur if the BAM file used was not created by ShortStack. An
    entry of '??' indicates that conflicting information about the number of
    hits is stored in the bam file.

  Counts file
    The Counts.txt file is plain text, tab-delimited. For each locus, it
    shows the total raw read counts. Each read-group is broken out
    seperately, and the sum of all read groups is also shown (termed
    'main'). Data from unplaced small RNAs, if present, are also included in
    Counts.txt

  MIRNAs directory
    This directory contains plain-text descriptions of each locus that was
    judged 'M' or 'Y' in MIRNA analysis. The files show the sequence of the
    locus, the predicted RNA secondary structure in dot-bracket notation,
    and the locations of the miRNA and miRNA-star. If the miRNA-star was not
    sequenced, its sequence is shown as 'X's instead of the real sequence.

    Below this top line, all other small RNAs aligned to the locus are
    shown. Those aligned to the opposite strand have '<' as delimiters
    instead of '.'.

    Lower-case nts in the displayed small RNA sequences indicate positions
    where the small RNA sequence differs from the reference sequence. Note
    that the reference sequence, not the small RNA sequences, are used to
    compute predicted secondary structures.

    l: length of RNA, a: number of alignments.

  GFF3 files
    If the run was a de-novo analysis, three gff3 files are created to
    indicate the positions of the discovered loci.

    ShortStack_N.gff3 has the loci with DicerCalls of 'N' (e.g., those that
    are unlikely RNAi-related).

    ShortStack_D.gff3 has the loci with DicerCalls that were not 'N' (e.g.,
    those that ARE likely RNAi-related).

    ShortStack_All.gff3 has ALL loci (it is the merger of the other two gff3
    files).

  bam or cram alignment file
    If aligment was performed, the final bam or cram formatted alignment
    will also be in the ShortStack outdir. The ShortStack-specific tags of
    these files are described above (section Alignment output format).

FAQ
  bowtie2 is newer than bowtie. Why do you still require bowtie but disallow bowtie2?
    Answer: Three reasons. 1) unlike bowtie2, bowtie has support for
    colorspace data, and 2) According to the manuals for both programs,
    bowtie2 is optimized for longer (>50 nts) reads, while bowtie is
    optimized for shorter reads. 3) Time. Despite the above comments, we
    will explore this transition in future versions of ShortStack.

  Why does ShortStack say that are known MIRNA loci are NOT MIRNAs?
    Answer: MIRNA annotation by ShortStack is, by design, meant to strongly
    reduce, perhaps eliminate, false positives. Any locus given a MIRNA
    result of 'Y' by ShortStack has sufficiently strong evidence to support
    its annotation of a miRNA. However, reduction of false positives comes
    at a price .. there will be some false negatives .. true MIRNAs that are
    not reported as such by ShortStack. Users should consider a 'No' MIRNA
    result by ShortStack to mean that the evidence in that particular small
    RNA-seq run did not offer 100% proof that the locus was a MIRNA.

  I ran the same analysis, with the same reads, and the same settings a second time, and received slightly different results. Is this a bug?
    Answer: No. This is caused by the treatment of multi-mapped reads.
    Because the decisions on which of the possible alignment positions are
    probabilistic, some small number of the reads will be differ in their
    selected primary positions when alignments are repeated. This is normal,
    and typically the differences are minor.

  I get different numbers of MIRNAs with ShortStack 3 relative to earlier versions. Why?
    Answer: The MIRNA detection methods have changed significantly. You may
    find ShortStack 3 to be more strict (find fewer MIRNAs) relative to
    earlier versions. This is because false positives are really minimized
    with ShortStack 3, potentially at the expense of some false negatives.

  What happened to the flagfile option from earlier versions of ShortStack?
    Answer: It is gone. This was a rather crude way to assess overlaps
    between ShortStack-discovered clusters and loci of a user's interest.
    Use bedtools instead (using the gff3 output from ShortStack).

  Are the read counts reported by ShortStack normalized in any way?
    Answer: The column 'Reads' in the Results is just the raw reads. The
    column 'RPM' in the Results is the reads per million, calculated on the
    basis of all aligned + unaligned reads in the library.

  ShortStack seems to be slow. Why? And how to make it go faster?
    Answer: To make alignments go faster, use the --bowtie_cores and
    --sort_mem options to make full use of your system. Their default
    settings (1 core, and sort_mem of 768M) are quite low to ensure success
    on low-powered machines, but if you have more cores and memory
    available, raising these will speed alignments along quite a bit.
    Another way to make alignments go faster is to specifiy r or n for
    option --mmap. But there is a trade-off there .. r causes multi-mapped
    reads to be just placed randomly instead of more intelligently, while n
    causes all multi-mappers to be marked as unmapped. So, if you use --mmap
    of r or n, you will get a much faster alignment, but a much less
    sensitive and precise one. There are also some tricks to increase the
    speed of analysis, but they all also involve some down-side. You can set
    option --nohp, which means that MIRNAs will not be tested for. This will
    increase the speed of analysis but of course you won't be able to
    annotate MIRNAs. You can also adjust the option --mincov to have a
    higher threshold. This will cause fewer loci to be discovered (only
    those with higher expression levels), so analysis time will be reduced.
    But of course the trade-off there is that you will not discover loci
    with lower expression levels. The current default of mincov 0.5rpm
    should be a good balance of sensitivity and speed for most applications.
