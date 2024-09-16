# ShortStack
Alignment of small RNA-seq data and annotation of small RNA-producing genes

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/shortstack/README.html) ![Conda](https://img.shields.io/conda/dn/bioconda/shortstack)

# Author
Michael J. Axtell, Penn State University, mja18@psu.edu

# IMPORTANT CHANGE - Condensed Reads
As of version 4.1.0 ShortStack creates and uses "condensed reads" for alignments. This speeds things up and keeps file sizes smaller. But it requires some upgrades and some understanding. Ensure that `strucVis` version is >= 0.9 and `ShortTracks` version is >= 1.2. This is enforced by bioconda but users who manually install will need to upgrade. See [Outputs](#outputs) and [Overview of Methods](#overview-of-methods) for more about condensed reads.

# Table of Contents
- [Citations](#citations)
- [Installation](#installation)
- [Usage](#usage)
- [Resources](#resources)
- [Testing and Examples](#testing-and-examples)
- [Vignettes](#vignettes)
- [Outputs](#outputs)
- [Visualizing Results](#visualizing-results)
- [Overview of Methods](#overview-of-methods)
- [How to go FAST](#how-to-go-fast)
- [ShortStack Version 4 Major Changes](#shortstack-version-4-major-changes)
- [Issues](#issues)
- [FAQ](#faq)

# Citations
If you use `ShortStack` in support of your work, please cite one or more of the following:

- Johnson NR, Yeoh JM, Coruh C, Axtell MJ. (2016). G3 6:2103-2111.
    doi:10.1534/g3.116.030452
- Shahid S., Axtell MJ. (2013) Identification and annotation of small RNA genes using ShortStack. Methods doi:10.1016/j.ymeth.2013.10.004
- Axtell MJ. (2013) ShortStack: Comprehensive annotation and quantification of small RNA genes. RNA 19:740-751. doi:10.1261/rna.035279.112

# Installation

You can either use the `conda` package manager to install from the bioconda channel, or manually set up an environment. Use of `conda`/bioconda is *highly* recommended!

## Install using conda (recommended)

First, install `conda`, and then set it up to use the bioconda channel following the instructions at <https://bioconda.github.io>

Then, follow instructions below based on your system to install the dependencies to a new environment and activate the environment.

### Linux or Intel-Mac
```
conda create --name ShortStack4 shortstack 
conda activate ShortStack4
```

### Silicon-Mac
Some dependencies have not been compiled for the newer Silicon-based Macs on bioconda, so we need to force conda to install the osx-64 (Intel) versions instead. Silicon Macs can run Intel code using built-in Rosetta translation.
```
conda create --name ShortStack4
conda activate ShortStack4
conda config --env --set subdir osx-64
conda install shortstack
```

## Manual installation

Create an environment that contains the following packages / tools compiled and installed. *Make note of required versions! (last updated for ShortStack release 4.1.0)*

- `python` >= 3.12.3 <https://www.python.org>
- `samtools` >= 1.20 <https://www.htslib.org>
- `bowtie` >= 1.3.1 <https://bowtie-bio.sourceforge.net/index.shtml>
- viennarna 2.* <https://www.tbi.univie.ac.at/RNA/documentation.html>
- `tqdm` <https://tqdm.github.io>
- `numpy` <https://numpy.org>
- biopython <https://biopython.org>
- `strucVis` >= 0.9 <https://github.com/MikeAxtell/strucVis>
- `ShortTracks` >= 1.2 <https://github.com/MikeAxtell/ShortTracks>
- `bedtools` >= 2.31.1 <https://bedtools.readthedocs.io/en/latest/>
- `cutadapt` >= 4.9 <https://cutadapt.readthedocs.io/en/stable/>

Then, download the `ShortStack` script from this github repo. Make it executable `chmod +x ShortStack` and then copy it into your environment's PATH.

# Usage
```
ShortStack [-h] [--version] --genomefile GENOMEFILE [--known_miRNAs KNOWN_MIRNAS] (--readfile [READFILE ...] | --bamfile [BAMFILE ...]) [--outdir OUTDIR] [--adapter ADAPTER | --autotrim] [--autotrim_key AUTOTRIM_KEY] [--threads THREADS] [--mmap {u,f,r}] [--align_only] [--dicermin DICERMIN] [--dicermax DICERMAX] [--locifile LOCIFILE | --locus LOCUS] [--nohp] [--dn_mirna] [--strand_cutoff STRAND_CUTOFF] [--mincov MINCOV] [--pad PAD] [--make_bigwigs]
```

## Required
- `--genomefile GENOMEFILE` : Path to the reference genome in FASTA format. Must be indexable by both `samtools faidx` and `bowtie-build`, or already indexed.
- `(--readfile [READFILE ...] | --bamfile [BAMFILE ...])` : *Either* `--readfile` or `--bamfile` is required.
    - `--readfile [READFILE ...]` : Path(s) to one or more files of reads in `fastq` or `fasta` format. May be `gzip` compressed. Multiple files are separated by spaces. Inputting reads triggers alignments to be performed.
    - `--bamfile [BAMFILE ...]` : Path(s) to one or more files of aligned sRNA-seq data in BAM format. Multiple files are separated by spaces. BAM files must match the reference genome given in `--genomefile`.

## Recommended
- `--known_miRNAs KNOWN_MIRNAS` : Path to FASTA-formatted file of known mature miRNAs. FASTA must be formatted such that a single RNA sequence is on one line only. ATCGUatcgu characters are acceptable. These RNAs are typically the sequences of known microRNAs; for instance, a FASTA file of mature miRNAs pulled from <https://www.mirbase.org>. These known miRNA sequences are aligned to the genome and used to nucleate searches for loci that meet all expression-based and secondary structure-based requirements for *MIRNA* locus identification. See also option `--dn_mirna`.
- `--outdir OUTDIR` : Specify the name of the directory that will be created for the results.
    - default: `ShortStack_[time]`, where `[time]` is the Unix time stamp according to the system when the run began.
- `--autotrim` : This is strongly recommended **when supplying untrimmed reads via `--readfile`**. The `autotrim` method automatically infers the 3' adapter sequence of the untrimmed reads, and the uses that to coordinate read trimming. However, do **not** use `--autotrim` if your input reads have already been trimmed!
    - Note: `autotrim` currently assumes your library strategy generated reads where nucleotide 1 of the read is the first biological / sRNA-derived nucleotide, and the 3' adapter starts immediately after the last sRNA nucleotide. It further assumes there are no random nucleotides (Ns) in the 3' adapter sequence. If your data do not meet these assumptions you cannot use `--autotrim`. Instead, remove your adapters by other appropriate methods and input the trimmed reads using `--readfile` *without* option `--autotrim`.
    - Note: mutually exclusive with `--adapter`.
- `--threads THREADS` : Set the number of threads to use. More threads = faster completion.
    - default: 1

## Other options
- `-h` : Print a help message and then quit.
- `--version` : Print the version and then quit.
- `--adapter ADAPTER` : Manually specify a 3' adapter sequence to use during read trimming. Mutually exclusive with `--autotrim`. The `--adapter` option will apply the same adapter sequence to trim **all** given readfiles.
    - Note: Use of `--adapter` is discouraged. In nearly all cases, `--autotrim` is a better bet for read trimming.
- `--autotrim_key AUTOTRIM_KEY` : A DNA sequence to use as a known suffix during the `--autotrim` procedure. ShortStack's autotrim discovers the 3' adapter by scanning for reads that begin with the sequence given by `AUTOTRIM_KEY`. This should be the sequence of a small RNA that is known to be highly abundant in all of the libraries. The default sequence is for miR166, a microRNA that is present in nearly all plants at high levels. For non-plant experiments, or if the default is not working well, consider providing an alternative to the default.
    - default: `TCGGACCAGGCTTCATTCCCC` (miR166)
- `--mmap {u,f,r}` : Sets the mode by which multi-mapped reads are handled. These modes are described in [Johnson et al. (2016)](https://doi.org/10.1534/g3.116.030452). The default `u` mode has the best performance.
    - `u` : (Default) Only uniquely-aligned reads are used as weights for placement of multi-mapped reads.
    - `f` : Fractional weighting scheme for placement of multi-mapped reads.
    - `r` : Multi-mapped read placement is random.
- `--align_only` : This switch will cause ShortStack to terminate after the alignment phase; no analysis occurs.
- `--dicermin DICERMIN` : An integer setting the minimum size (in nucleotides) of a valid small RNA. Together with `--dicermax`, this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. >= 80% of the reads in a given cluster must be in the range indicated by `--dicermin` and `--dicermax`.
    - default: 21
- `--dicermax DICERMAX` : An integer setting the minimum size (in nucleotides) of a valid small RNA. Together with `--dicermin`, this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. >= 80% of the reads in a given cluster must be in the range indicated by `--dicermin` and `--dicermax`.
    - default: 24
- `--locifile LOCIFILE` : Path to a file of pre-determined loci to analyze. This will prevent *de novo* discovery of small RNA loci. The file may be in gff3, bed, or simple tab-delimited format (Chr:Start-Stop[tab]Name). Mutually exclusive with `--locus`.
- `--locus LOCUS` : A single locus to analyze, given as a string in the format Chr:Start-Stop (using one-based, inclusive numbering). This will prevent *de novo* discovery of small RNA loci. Mutually exclusive with `--locifile`.
- `--nohp` : Switch that prevents search for microRNAs. This saves computational time, but *MIRNA* loci will not be differentiated from other types of small RNA clusters.
- `--dn_mirna` : Switch that activates a *de novo* search for *MIRNA* loci. By default ShortStack will confine *MIRNA* analysis to loci where one or more queries from the `--known_miRNAs` file are aligned to the genome. Activating *de novo* searching with `--dn_mirna` does a more comprehensive genome-wide scan for *MIRNA* loci. Loci discovered with `--dn_mirna` that do not overlap already known microRNAs should be treated with caution.
- `--strand_cutoff STRAND_CUTOFF` : Floating point number that sets the cutoff for standedness. Must be > 0.5 and < 1.
    - default: 0.8. Loci with >80% reads on the top genomic strand are '+' stranded, loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'
- `--mincov MINCOV` : Minimum alignment depth, in units of reads per million, required to nucleate a small RNA cluster during *de novo* cluster search. Must be an floating point number > 0. 
    - default: 1
- `--pad PAD` : Initial peaks (continuous regions with depth exceeding argument `--mincov`) are merged if they are this distance or less from each other. Must be an integer >= 1. 
    - default: 200
- `--make_bigwigs` : Enable creation of bigwig files from sRNA-seq alignments. Applies only when performing alignments via input to `--readfile`.

# Resources
## Memory
Peak memory (RSS) primarily scales with genome size, and especially with the size of the largest chromosomes in the genome assembly. Generally between 4-10GB memory per thread is more than enough, but very large genomes should be monitored for memory overruns and requested memory adjusted accordingly. Because of the multi-threading, it is important to scale available memory with the number of threads being used.

## Disk
During the alignment phase ShortStack will potentially write many large, but temporary, files to disk. These can easily exceed 100GB disk space especially when there are a lot of multi-mapping reads. Ensure that the location of `--outdir` has ample free disk space (plan on 200GB minimum to be safe).

## CPU and running time
All compute-intensive parts of ShortStack are now multi-threaded. Providing more threads via `--threads` generally will decrease run-times in near-linear fashion. Be sure to scale memory with thread use though .. 4-10GB RAM per thread, depending on genome size, seems usually sufficient.

Read alignment is often the most time-consuming portion of the analysis. *MIRNA* identification (triggered by `--known_miRNAs` and/or `--dn_mirna`) is also time-consuming. Larger genomes generally run slower compared to smaller genomes. Highly fragmented genome assemblies (*e.g.* very high numbers of chromosomes/scaffolds) can be particularly slow because of the index lookup costs associated with thousands of entries. Consider obtaining and using better genome assemblies, or removing very short scaffolds from highly fragmented genome assemblies.

# Testing and Examples
## Gather Test Data
### Reference Genome
For testing we will use the *Arabidopsis thaliana* TAIR10 genome assembly, including plastid and mitochondrial genomes. The Axtell Lab website is serving a copy of this that is easy to obtain:

```
curl https://plantsmallrnagenes.science.psu.edu/ath-b10/Arabidopsis_thalianaTAIR10.fa -ko Arabidopsis_thalianaTAIR10.fa
```
### Raw sRNA-seq Reads
We can obtain some example *A. thaliana* sRNA-seq data from SRA using the SRA toolkit. If you don't already have the sra toolkit installed you can obtain it at https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

After installing the sra toolkit you should configure it .. see <https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration>

*tip* I like to configure sra tools to prefetch and download to the working directory. That way I can remember to get rid of the pre-fetched directories when I have dumped the FASTQ.

Then you can use the standard `prefetch` and `fasterq-dump` method to retrieve the FASTQ files.

```
prefetch SRR3222443 SRR3222444 SRR3222445
fasterq-dump SRR3222443 SRR3222444 SRR3222445
```
You will now have 3 `.fastq` files of raw (untrimmed) sRNA-seq reads. These data are derived from Col-0 *Arabidopsis thaliana* immature inflorescence tissues (see Wang et al. 2017 <https://doi.org/10.1111/tpj.13463>)

#### Alternative using SRA run browser
Sometimes sra-tools is a pain to use. Alternatively, you can use the web interface at <https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&display=metadata> to grab the example data, searching by SRR accession number, and downloading FASTQ.

### Known miRNAs
To get a list of known miRNAs, we will use all [miRBase](https://www.mirbase.org) annotated mature miRNAs from miRBase. First, download the `mature.fa` file from miRBase at <https://mirbase.org/download/>. Then filter it to get only the `ath` ones (*e.g.* the ones from *A. thaliana*).

```
grep -A 1 '>ath' mature.fa | grep -v '\-\-' > ath_known_miRNAs.fasta
```

## Example Run
This example is a full run. It takes 3 raw (untrimmed) readfiles, identifies the adapters, trims the reads, indexes the genome, aligns the reads, discovers small RNA loci, and annotates high-confidence *MIRNA* loci. The example uses 6 threads; this can be adjusted up or down depending on your system's configuration; more threads decrease execution time but the response is non-linear (diminishing returns with very high thread numbers). The examples uses the test data described above.

```
ShortStack --genomefile Arabidopsis_thalianaTAIR10.fa --readfile SRR3222443.fastq SRR3222444.fastq SRR3222445.fastq --autotrim --threads 6 --outdir ExampleShortStackRun --known_miRNAs ath_known_miRNAs.fasta
```

On my compute cluster using 6 threads this completes in about 9 minutes. All results are in the directory specified by `--outdir`, "ExampleShortStackRun". The outputs are described in the section below called "Outputs".

# Vignettes
More examples of common tasks can be found on the wiki at <https://github.com/MikeAxtell/ShortStack/wiki>

# Outputs
## Results.txt
A tab-delimited text file giving key information for all small RNA clusters. Columns:
1. *Locus*: Coordinates of the locus in Chrom:Start-Stop format, one-based, inclusive.
2. *Name*: Name for the Locus. *De novo* loci are named like `Cluster_1`, etc. Loci provided by the user with option `--locifile` preserve the user-given names.
3. *Chrom*: Name of the chromosome
4. *Start*: One-based start position (left-most) of the locus, inclusive.
5. *End*: One-based end position (right-most) of the locus, inclusive.
6. *Length*: Length of the locus (base-pairs)
7. *Reads*: Number of aligned sRNA-seq reads that overlap this locus.
8. *DistinctSequences*: Number of distinct sRNA sequences that overlap this locus. A single DistinctSequence can have one or more reads.
9. *FracTop*: Fraction of Reads aligned to the top genomic strand.
10. *Strand*: Inferred strandednes of the locus, inferred from FracTop and the `--strand_cutoff` setting.
11. *MajorRNA*: Sequence of the single most abundant RNA sequence at the locus.
12. *MajorRNAReads*: Number of reads for the MajorRNA aligned to this locus.
13. *Short*: Number of reads with a size less than `--dicermin` aligned to the locus. By default these are reads < 21 nucleotides in length.
14. *Long*: Number of reads with a size greater than `--dicermax` aligned to the locus. By default these are reads > 24 nucleotides in length.
15. *21*: Number of 21 nucleotide reads aligned to the locus.
16. *22*: Number of 22 nucleotide reads aligned to the locus.
17. *23*: Number of 23 nucleotide reads aligned to the locus.
18. *24*: Number of 24 nucleotide reads aligned to the locus.
19. *DicerCall*: If >= 80% of all aligned reads are within the boundaries of `--dicermin` and `--dicermax`, than the DicerCall gives the size of most abundant small RNA size. If < 80% of the aligned reads are in the `--dicermin` and `--dicermax` boundaries, DicerCall is set to 'N'. Loci with a DicerCall of 'N' are unlikely to be small RNAs related to the Dicer-Like/Argonaute system of gene regulation.
20. *MIRNA*: Did the locus pass all criteria to be called a *MIRNA* locus? If so, 'Y'. If not, 'N'.
21. *Known_miRNAs*: Semicolon delimited list of user-provided known RNAs that aligned to the locus. If none, 'NA'.

## Counts.txt
A tab-delimited text file giving the raw alignment counts for each locus in each separate sample. Only produced if there was more than one sRNA-seq file used to create alignments. This file is useful for downstream analyses, especially differential expression analysis.

## autotrim_details.tsv
If ShortStack trimmed the reads, details about the trim process are stored in this file. Input_reads and Output_reads tally the total number of reads in those two categories. Reads that were NOT output were either too short after trimming or no adapter was found ("Too_short_reads" and "No_adapter_reads", respectively).

## alignment_details.tsv
A tab-delimited text file that gives details about small RNA-seq alignments as a function of mapping type, read length, and sample. This can be useful for plotting purposes and subsequent quality control of small RNA-seq data. It is only produced if alignments are performed.
- *mapping_type*
    - U: Uniquely mapped (not a multimapper).
    - P: Multimapper placed using the method set by option `--mmap`.
    - R: Multimapper placed at random.
    - H: Very highly multi-mapped read (>=50 hits).
    - N: Unmapped reads.

## Results.gff3
Small RNA loci in the gff3 format. Suitable for use on genome browsers. For loci that are annotated as *MIRNAs* there will be an additional entry for the mature microRNA position. The 'score' column in the gff3 format stores the number of sRNA-seq aligned reads at that locus.

## known_miRNAs.gff3
When the user provides known RNA sequences via `--known_miRNAs`, they are aligned to the reference genome. Every (perfect) alignment to the reference is stored and reported in the known_miRNAs.gff3 file. The score column shows the number of alignments that start and end at the exact coordinates and strand.

**Important** : known_miRNAs are aligned and shown in the `known_miRNAs.gff3` file regardless of whether any empirical small RNA-seq data are found. Thus, expect to find entries with a score of 0; these are cases where no instances of the given knownRNA were aligned to that location in the genome.

## strucVis/
The directory `strucVis/` contains visualizations of each locus that was annotated as a *MIRNA* locus. These are made by the script [strucVis](https://github.com/MikeAxtell/strucVis). For each locus there is a postscript file and a plain-text file. Both show the coverage of aligned small RNA-seq data as a function of position, aligned with the predicted RNA secondary structure of the inferred *MIRNA* hairpin precursor. These files are meant for manual inspection of *MIRNA* loci.

## mir.fasta
This is a FASTA formatted file containing hairpin, mature miRNA, and miRNA* sequences derived from ShortStack's identification of *MIRNA* loci. These are genomic sequences, and the genomic coordinates are noted in the FASTA header. ShortStack's determination of mature miRNA vs. miRNA* strands is based on abundance of alignments at that particular locus. These designations may not always be accurate for an entire *MIRNA* family .. sometimes one paralog can attract most of the true mature miRNA alignments, leaving the other paralogs with mostly true miRNA* alignments. Take care when performing annotations.

## Trimmed and condensed readfiles (fasta formatted)
If raw reads were trimmed by ShortStack, they will also have been "condensed" and then written to FASTA formatted files. Read condensation writes a single unique sequence just once, regardless of how many reads had that sequence in the original input. The read depth is noted in the condensed sequence FASTA header. The names will have a lower-cased 't' appended to the front to signify "trimmed", and "_condensed" to indicate that they were condensed such that each unique sequence is only written once. If the input files were .gz compressed, then the trimmed files will be too.

Example of condensed format:
```
>tSRR3222443_Cd1482943_509066
TCGGACCAGGCTTCATTCCCC
```

- `tSRR3222443` : Prefix specific for the input read file. The `t` indicates trimmed.
- `Cd1482943` : `Cd` means "condensed". The integer after is a unique integer within that file.
- `_509066` : The trailing integer that represents the total number of reads that have this sequence. In this case, 509,066 reads with this sequence.

## .bam file(s)
If the ShortStack run was aligning reads, one or more `.bam` files will be found. The bam format stores large scale alignment data. Corresponding bam index files (`.bam.csi`) will also be found. 
- *IMPORTANT* : As of ShortStack release 4.1.0, bam files will be written using condensed reads. This saves file space and has other advantages during processing. HOWEVER, the bam files containing condensed reads CANNOT be used naively to calculate read depth. Each alignment of a condensed sequence has one line in the BAM file. However, that condensed sequence may have had more than one read placed at that position. See the format details below:

### bam/sam format for ShortStack >= 4.1.0 condensed reads
ShortStack and bowtie both add several optional SAM tags to alignment lines. 

| Tag | Type | Source | Meaning |
| --- | ---- | ----- | ------- | 
| XA:i:\<int> | integer | bowtie | Aligned sequence belongs to stratum \<int> |
| MD:Z:\<S> | string | bowtie | String representation of the mismatched reference bases in the alignment. See SAM format specification for details |
| NM:i:\<int> | integer | bowtie | Edit distance to the reference |
| XM:i:\<int> | integer | bowtie | For a read with no reported alignments, <N> is 0 if the read had no alignments. If -m was specified and the read’s alignments were suppressed because the -m ceiling was exceeded, <N> equals the -m ceiling + 1, to indicate that there were at least that many valid alignments (but all were suppressed). In -M mode, if the alignment was randomly selected because the -M ceiling was exceeded, <N> equals the -M ceiling + 1, to indicate that there were at least that many valid alignments (of which one was reported at random). ShortStack uses bowtie `-k` of 12 and does not set `-m`.
| XX:i:\<int> | integer | ShortStack | Total number of possible alignment positions in this genome |
| XY:Z:\<S> | string | ShortStack | A single letter that notes the type of placement that ShortStack made: `U`: Sequence was uniquely aligned to this location in the genome. `P`: Sequence was multimapped to the genome and at least some of it's reads were placed here using ShortStack's multi-mapping read algorithm. `R`: Sequence was multimapped to the genome and it's placements on the genome were randomly picked. `H`: Sequence was very highly multimapped to the genome (>= 20 hits). Read allocation was random and only 20 of the possible locations were even considered. `N`: Sequence was not mapped at all to the reference genome.
| XZ:f:\<F>| floating point number | ShortStack | Fraction of the reads for this sequence that were allocated to this alignment position. For sequences with XY:Z type `P`, this reflects the multi-mapping algorithm's judgement on read allocation. For sequences with XY:Z type `R` or `H`, this is an equal split with all other alignment positions. For sequences with XY:Z type `U` or `N`, this will be 1, because 100% of the reads are allocated to this position.
| YS:Z:\<S> | string | ShortStack | String describing the sequence length category. Possible categories are `<21` , `21`, `22`, `23`, `24`, and `>24`.
| XW:i:\<int> | integer | ShortStack | Count of reads of this sequence that were allocated to this alignment position |
| RG:\<Z> | string | bowtie | Read-group name. This will be set when more than one file of reads was input to alignments. |

### Understanding multi-mapped condensed sequences using bam tags

Here's a condensed sequence that was multimapped and positions assigned using ShortStack's multi-mapping algorithm: `tSRR3222443_Cd1482943_509066`. As described above, this sequence was observed for 509,066 of the original reads. There are 7 total lines in the bam file that list this sequence. The key information is below:

| Chromosome | Position | XX:i | XY:Z | XZ:f | XW:i |
| --- | --- | --- | --- | --- | --- | 
| 2 | 19176241 | 7 | P | 0.995 | 506593 |
| 3 | 22922301 | 7 | P | 0 | 73 |
| 5 | 2838737 | 7 | P | 0 | 73 |
| 5 | 2840708 | 7 | P | 0 | 72 |
| 5 | 16775524 | 7 | P | 0 | 73 |
| 5 | 17516378 | 7 | P | 0.004 | 2109 |
| 5 | 25504882 | 7 | P | 0 | 73 |

We see that this sequence had seven possbile alignment positions (XX:i = 7) and that ShortStack used it's multi-mapper algorithm to allocate the reads (XY:Z = P). ShortStack thought that most of these reads, 99.5% of them (XZ:f = 0.995), belonged at the chromosome 2 location. This is 506,593 of the reads (XW:i = 506,593). 0.04% of the reads (2109) were placed at chromosome 5:17516378, and small numbers of reads (72 to 73) allocated to the others. Note that the XZ:f numbers are rounded to three decimal places, so a value of '0' is actually a quite small fraction in this example.


## .bw (bigwig) files
If the ShortStack run was aligning reads AND option `--make_bigwigs` was set, multiple `.bw` files in the [bigwig format](https://genome.ucsc.edu/goldenpath/help/bigWig.html) will be produced. These bigwig files are produced by the program [ShortTracks](https://github.com/MikeAxtell/ShortTracks) and are useful for visualizing data on genome browsers (in particular [JBrowse2](https://jbrowse.org/jb2/)). The scales of the bigwig files are normalized to reads per million; thus, multiple tracks can be compared.

Types of bigwig files produced by ShortStack:
- readlength/stranded: A set of 8 `.bw` files with suffixes in the format `_x_y.bw`
    - x : A number ('21', '22', '23-24') or the string 'other', indicating the sizes of sRNAs tracked in that file.
    - y : Either 'p' or 'm' for the plus or minus genomic strand.
- readgroup: Only produced when there are multiple samples in the alignment. A set of n `.bw` files, with one per read-group. Because values are normalized to reads-per-million, these tracks are directly comparable to each other.

The README for [ShortTracks](https://github.com/MikeAxtell/ShortTracks) has details about how to load these data onto [JBrowse2](https://jbrowse.org/jb2/) using "multi-wiggle" tracks for nice visualization of sRNA-seq data.

# Visualizing Results

Loci annotated as *MIRNA* can be visualized from the `srucVis/` files. These show the predicted RNA secondary structures with the small RNA-seq read depth coverage.

## Genome Browsers
The output of ShortStack is designed to work with genome browsers. Specifically, the files `Results.gff3`, `known_miRNAs.gff3`, the `.bam` files, and the optional `.bw` files can be directly visualized on either major genome browser (IGV, JBrowse).

[JBrowse2](https://jbrowse.org/jb2/) has the ability to create "multi-wiggle" tracks. These tracks show multiple quantitative data tracks at once, bound to a common quantitative axis. The `.bw` bigwig files created by ShortStack & ShortTracks are normalized to reads-per-million, allowing direct comparisons in a multi-wiggle track. This allows visualization of size, coverage, and strandedness of the data. See the README for [ShortTracks](https://github.com/MikeAxtell/ShortTracks) for details. I recommend using the Desktop version of [JBrowse2](https://jbrowse.org/jb2/).

# Overview of Methods
## Indexing
Two types of genome indices are required. The first is an `.fai` index, created by `samtools faidx`. The second is a `bowtie` index of the genome file, which comprises several files ending with `.ebwt` or `.ebwtl`. If either index is missing, ShortStack will attempt to create it. The bowtie index is created using the `bowtie-build` executable from `bowtie` using default settings and the number of threads given by the `--threads` argument. Any index files created by ShortStack are written to the same location as the `--genomefile`.  Therefore, the `--genomefile` should be located in a writable directory if the required indices are not already built at the start of the run.

## Read trimming
Read trimming using `--autotrim` assumes that the input FASTQ data comes from forward strand sequence reads where the first letter of the read corresponds to the first nucleotide of the sense-strand small RNA. Read trimming seeks to first identify the 3' adapter sequence used in the FASTQ file, and then to use the tool `cutadapt` to remove the adapters. Using the recommended `--autotrim` setting, the sequence of the 3' adapter is automatically detected. This is accomplished by looking for reads that begin with an `--autotrim_key`, which by default is the miR166 sequence 'TCGGACCAGGCTTCATTCCCC'. The adapter is inferred by looking at the sequences that follow the `--autotrim_key`. The `cutadapt` settings used for read trimming are `cutadapt -a [adapter] -o [tInput.fastq] -j [--threads] --discard-untrimmed -m 12 --report minimal [Input.fastq]`. These settings discard untrimmed reads (where no adapter was found), and only retained reads that, after trimming, are at least 12 nucleotides long.

## Read condensation
*New as of ShortStack 4.1.0* :: Small RNA-seq data often contain multiple reads with the same sequence. These are often microRNAs or siRNAs that were highly abundant in the sample. "Read condensation" refers to the process of condensing the input for alignment such that each unique sequence is only represented once in the condensed FASTA file. Suppose the sequence TCGGACCAGGCTTCATTCCCC was the sequence for 509,066 reads. Instead of repetitively writing, and aligning, the same sequence over and over we can instead write a single FASTA entry. Example of condensed format:
```
>tSRR3222443_Cd1482943_509066
TCGGACCAGGCTTCATTCCCC
```

- `tSRR3222443` : Prefix specific for the input read file. The `t` indicates trimmed.
- `Cd1482943` : `Cd` means "condensed". The integer after is a unique integer within that file.
- `_509066` : The trailing integer that represents the total number of reads that have this sequence. In this case, 509,066 reads with this sequence.

Using read condensation has several advantages:
- Speeds up alignments
- Reduces file sizes of temporary files and of the final BAM files.

However, there is a pitfall:
- BAM files are *no longer* one alignment = one read. Read depth must be inferred from the XW:i tags in the bam file. Use extreme caution with analysis of these bam files made with 'condensed' reads. 

## Alignments
Alignment of trimmed fastq data uses `bowtie`. There are usually two phases to alignment. In the first phase reads are mapped to the genome using bowtie settings `bowtie -p [--threads] -v 1 -k 50 -S --best --strata -x [--genomefile] [trimmedFASTQ]`. This allows zero or one mismatch hits, keeping only the zeroes if both zero and one-hit cases exist. All hits up to a maximum of 50 are stored for multi-mapping reads. In the second phase a single location for each of the multimapping reads is decided upon. This decision is made by one of the three possible methods set by `--mmap` : u, f, or r. Setting u uses a local weighting scheme set by only the uniquely mapping alignments in an area. Setting f uses a local weighting scheme that includes both the unique and multi-mapped possibilities in an area. Setting r just randomly takes one of the possible location as the reported one. These methods are fully described in [Johnson et al., 2016](https://doi.org/10.1534/g3.116.030452). Each individual FASTQ file is processed and converted to a sorted BAM file. If more than one FASTQ file was input, a single sorted merged file, `merged_alignments.bam` is created.

## Cluster discovery
Genomic intervals where the depth of small RNA coverage, in reads-per-million, is greater than or equal to `--mincov` (1 by default) are identified. Each of these intervals are then extended in both directions by the length given by setting `--pad`. Regions that overlap after extension are merged. Note that if one or more *MIRNA* loci are later found to overlap the intial cluster, the initial cluster is removed from the output (and only the refined, trimmed *MIRNA* region(s) are reported).

## MIRNA annotation
*MIRNA* annotation has two entry points for initial searches: Locations of aligned user-provided sequences from `--known_miRNAs` and, if option `--dn_mirna` is True, any 21 or 22 nt read whose abundance exceeds the depth of `--mincov`. From these initial starting points, ShortStack first examines the local region to find miR/miR-star-like patterns of read accumulation (essentially, "two-peaks" of read coverage on the same genomic strand that might correspond to the miR/miR-star pair). If such a pattern is found, the RNA secondary structure in the local area is predicted. The sRNA-seq alignments in conjunction with the predicted RNA secondary structure are analyzed with respect to the criteria in [Axtell and Meyers, 2018](https://doi.org/10.1105/tpc.17.00851). If the criteria are met, the locus is annotated as a *MIRNA*.

# How to go FAST
ShortStack is already pretty fast. But if you want more speed, try:
- If performing alignments, set `mmap` to `r`. Setting `mmap r` means that multi-mapped reads will just be placed using a random-guess instead of the weighted method described by [Johnson et al. (2016)](https://doi.org/10.1534/g3.116.030452). The slower, weighted methods are slower but more accurate. 
- Use more `--threads` as your system allows. Up to 10 or 20. Be sure to allocate enough total memory (about 8-10GB RAM per thread, maybe higher for large genomes).
- Use a better reference genome: Highly fragmented genome assemblies are much slower than well-assembled genomes.

# ShortStack version 4 Major Changes
ShortStack version 4 is a major update. The major changes are:
## Major Changes

- Completeley re-written in `python3`.
- Streamlined installation using a `conda` recipe hosted on bioconda.
- All compute-intensive processes are now multi-threaded, so execution times are faster when the user specifies higher values of `--threads`.
- Much more reliance on other tools (`bedtools`, `cutadapt` for instance) .. less re-inventing of wheels.
- Output of hairpin structure visualizations using [strucVis](https://github.com/MikeAxtell/strucVis).
- Output of genome-browser-ready quantitative coverage tracks of aligned small RNAs using [ShortTracks](https://github.com/MikeAxtell/ShortTracks).
- *MIRNA* locus identification has been thoroughly changed to increase sensitivity while maintaining specificity.
- *MIRNA* locus identification can now be guided by user-provided 'known RNAs'. In contrast, truly *de novo* annotation of *MIRNA* loci, in the absence of matching the sequence of a 'known RNA' is disabled by default. This change in philosophy acknowledges that, in most well-studied organisms, most high-confidence microRNA families are already known.
- Change the license to MIT from GPL3.
- (As of version 4.1.0) Implement use of condensed reads for global speed up and reduced file sizes.

## Option changes:
- Drop support for cram format (options `--cram`, `--cramfile`  eliminated)
- Drop support for colorspace (option `--cquals` eliminated)
- Replace option `--bowtie_cores` with `--threads`
- Eliminate option `--bowtie_m`. Now -k 50 is always used.
- Eliminate option `--ranmax`. Now mmappers will always be placed (except mode u)
- Eliminate SAM tags XY:Z:O and XY:Z:M .. no more suppression of mmap reads
- Add SAM tag XY:Z:H .. highly repetitive read (50 or more hits, not all known).
- Add SAM tag YS:Z .. small RNA size information
- Eliminate option `--keep_quals`. Quality values will always be stored in the bam file if input was fastq.
- Modify option `--locus` so that it only accepts a single locus query.
- Eliminate option `--total_primaries` .. instead use a fast hack to rapidly calculate this.
- Option `--locifile` now understands .bed and .gff3 formats, as well as the original simple tab-delimited format.
- Added options `--autotrim` and `--autotrim_key`. This allows automatic detection of 3' adapters by tallying the most common sequence that occurs after a known, highly abundant small RNA (given by `autotrim_key`).
- Add option `--known_miRNAs`. Provide a FASTA file of known mature small RNA sequences to search for and to nucleate searches for qualifying *MIRNA* loci.
- Add option `--dn_mirna`. The `--dn_mirna` activates a *de novo* search for *MIRNA* loci independent of those that align to the 'known RNAs' provided by the user. By default, `--dn_mirna` is not active.
- Remove option `--show_secondaries` (as of 4.1.0)
- Remove option `--no_bigwigs` (as of 4.1.0)
- Add option `--make_bigwigs` (as of 4.1.0)

# Issues
Please post issues, comments, bug reports, questions, etc. to the project github page at <https://github.com/MikeAxtell/ShortStack>.

# FAQ

- **I ran an analysis and found no loci annotated as *MIRNA* loci!** 
    - By default, ShortStack will not do a *de novo* search for loci that qualify as *MIRNA* loci. To search for *MIRNA* loci the user has to explicitly request it, using either or both of the options `--known_miRNAs` and `--dn_mirna`. `known_miRNAs` provides a list of known mature miRNA sequences. Places where these sequences align to the reference genome are examined to see if the small RNA alignment pattern and predicted RNA secondary structure qualifies as a *MIRNA* locus. The switch `--dn_mirna` turns on a *de novo MIRNA* search. The *de novo MIRNA* search is turned off by default to reduce false annotations.  The idea is that most mature miRNAs are known in most species by now.
- **What happened to the phasing scores?**
    - I decided to omit phasing scores as of ShortStack version 4.0. This is because I gradually have lost confidence the accuracy of genome-wide scans to provide acceptable sensitivity *and* specificity for scoring phasing. For a detailed analysis of the challenges of calling phasing of siRNA clusters in genome-wide analyses, see [Polydore et al. (2018)](https://doi.org/10.1002/pld3.101). I am considering bringing phasing scores back, but just for 21-22 nt siRNA loci, in a future release.
- **Installation fails with conda**
    - Many problems I've seen with installation using conda stem from incorrectly configured conda installations. It is critical to have conda configured exactly as specified by the [bionconda project](https://bioconda.github.io). Your `.condarc` file (usually found in the user's home directory) must have the channels configured with `conda-forge` as the highest priority channel, `bioconda` as the second-highest priority, and `defaults` as the lowest priority. Additionally, the 'strict' channel priority option must be enabled. See instructions at [bionconda project](https://bioconda.github.io). A correctly configured conda should have a `.condarc` file that looks like:
```
channels:
    - conda-forge
    - bioconda
    - defaults
channel_priority: strict
repodata_fns:
    - repodata.json
```
- **An older version of ShortStack is installed by conda**
    - This can occur because of quirks in the way that conda searches for package updates .. sometimes it relies on the local cache of available packages instead of scanning the remotes. To overcome this, explicitly request the desired version. For instance, to get ShortStack version 4.0.1:
```
conda create --name ShortStack4
conda activate ShortStack4
conda install shortstack=4.0.1
```
- **Strange errors**
    - Check that all dependencies are up to date. In particular, look at `samtools`. ShortStack requires a modern `samtools` version. Use of old versions of `samtools` can cause cryptic errors. Currently, `samtools` >= 1.20 is required. This is noted here in the README. It is also enforced in the bioconda recipe. You can check the versions of dependencies that ShortStack requires at the [bioconda recipe here](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/shortstack/meta.yaml).

