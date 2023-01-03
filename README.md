# ShortStack

Alignment of small RNA-seq data and annotation of small RNA-producing genes

# Author
Michael J. Axtell, Penn State University, mja18@psu.edu

# Alpha testing
*Important* : The ShortStack 4 branch is in alpha testing mode! Not ready for production use!

# Citations
If you use `ShortStack` in support of your work, please cite one or more of the following:

- Johnson NR, Yeoh JM, Coruh C, Axtell MJ. (2016). G3 6:2103-2111.
    doi:10.1534/g3.116.030452
- Axtell MJ. (2013) ShortStack: Comprehensive annotation and
    quantification of small RNA genes. RNA 19:740-751.
    doi:10.1261/rna.035279.112
- Shahid S., Axtell MJ. (2013) Identification and annotation of small RNA
    genes using ShortStack. Methods doi:10.1016/j.ymeth.2013.10.004

# Installation

**NOTE** This is for alpha testing! After first release there will instead be a full `conda` recipe available from bioconda which will simplify installation!

## Create a suitable environment

First, install `conda`, and then set it up to use the bioconda channel following the instructions at <https://bioconda.github.io>

Then, follow instructions below based on your system to install the dependencies to a new environment and activate the environment.

### Linux or Intel-Mac
```
conda create --name ShortStack4 strucvis shorttracks bedtools biopython bowtie cutadapt numpy tqdm 
conda activate ShortStack4
```

### Silicon-Mac
Some dependencies have not been compilied for the newer Silicon-based Macs on bioconda, so we need to force conda to install the osx-64 (Intel) versions instead. Silicon Macs can run Intel code using built-in Rosetta translation.
```
conda create --name ShortStack4
conda activate ShortStack4
conda config --env --set subdir osx-64
conda install strucvis shorttracks bedtools biopython bowtie cutadapt numpy tqdm
```

## Manually install script
Once the environment is prepared, manually install the `ShortStack` script. Download the file (or archive) from the ShortStack4 branch on gitub at <https://github.com/MikeAxtell/ShortStack/tree/ShortStack4>. *Be sure you are using the ShortStack4 branch, not the master branch.*

Place the script `ShortStack` into the correct location in your conda environment. You can find this by issuing a `conda info --envs` command to find the correct path. See example from my laptop below.

```
conda info --envs
```
The above command shows me the directory for my conda environment:
```
ShortStack4           *  /Users/mja18/miniconda3/envs/ShortStack4
```
Knowing this, we can install the `ShortStack` script within the `/bin/` directory of that location. 
```
chmod +x ShortStack
mv ShortStack /Users/mja18/miniconda3/envs/ShortStack4/bin/
```

# Usage
```
ShortStack [-h] [--version] --genomefile GENOMEFILE [--knownRNAs KNOWNRNAS] (--readfile [READFILE ...] | --bamfile [BAMFILE ...]) [--outdir OUTDIR] [--adapter ADAPTER | --autotrim] [--autotrim_key AUTOTRIM_KEY] [--threads THREADS] [--mmap {u,f,r}] [--align_only] [--dicermin DICERMIN] [--dicermax DICERMAX] [--locifile LOCIFILE | --locus LOCUS] [--nohp] [--dn_mirna] [--strand_cutoff STRAND_CUTOFF] [--mincov MINCOV] [--pad PAD]
```

## Required
- `--genomefile GENOMEFILE` : Path to the reference genome in FASTA format. Must be indexable by `samtools faidx`, or already indexed.
- `(--readfile [READFILE ...] | --bamfile [BAMFILE ...])` : *Either* `--readfile` or `--bamfile` is required.
    - `--readfile [READFILE ...]` : Path(s) to one or more files of reads in `fastq` or `fasta` format. May be `gzip` compressed. Multiple files are separated by spaces. Inputting reads triggers alignments to be performed.
    - `--bamfile [BAMFILE ...]` : Path(s) to one or more files of aligned sRNA-seq data in BAM format. Multiple files are separated by spaces. BAM files must match the reference genome given in `--genomefile`.

## Recommended
- `--knownRNAs KNOWNRNAS` : Path to FASTA-formatted file of known small RNAs. FASTA must be formatted such that a single RNA sequence is on one line only. ATCGUatcgu characters are acceptable. These RNAs are typically the sequences of known microRNAs; for instance, a FASTA file of mature miRNAs pulled from <https://www.mirbase.org>. Providing these data increases the accuracy of *MIRNA* locus identification.
- `--outdir OUTDIR` : Specify the name of the directory that will be created for the results.
    - default: `ShortStack_[time]`, where `[time]` is the Unix time stamp according to the system when the run began.
- `--autotrim` : This is strongly recommended **when supplying untrimmed reads via `--readfile`**. The `autotrim` method automatically infers the 3' adapter sequence of the untrimmed reads, and the uses that to coordinate read trimming. However, do **not** use `--autotrim` if your input reads have already been trimmed!
    - Note: mutually exlcusive with `--adapter`.
- `--threads THREADS` : Set the number of threads to use. More threads = faster completion.
    - default: 1

## Other options
- `-h` : Print a help message and then quit.
- `--version` : Print the version and then quit.
- `--adapter ADAPTER` : Manually specify a 3' adapter sequence to use during read trimming. Mutually exclusive with `--autotrim`. The `--adapter` option will apply the same adapter sequence to trim **all** given readfiles.
    - Note: Use of `--adapter` is discouraged. In nearly all cases, `--autotrim` is a better bet for read trimming.
- `--autotrim_key AUTOTRIM_KEY` : A DNA sequence to use as a known suffix during the `--autotrim` procedure. ShortStack's autotrim discovers the 3' adapter by scanning for reads that begin with the sequence given by `AUTOTRIM_KEY`. This should be the sequence of a small RNA that is known to e highly abundant in all of the libraries. The default sequence is for miR166, a microRNA that is present in nearly all plants at high levels. For non-plant experiments, or if the default is not working well, consider providing an alternative to the default.
    - default: `TCGGACCAGGCTTCATTCCCC` (miR166)
- `--mmap {u,f,r}` : Sets the mode by which multi-mapped reads are handled. These modes are described in [Johnson et al. (2016)](https://doi.org/10.1534/g3.116.030452). The default `f` mode has the best performance.
    - `u` : Only uniquely-aligned reads are used as weights for placement of multi-mapped reads.
    - `f` : (Default) Fractional weighting scheme for placement of multi-mapped reads.
    - `r` : Multi-mapped read placement is random.
- `--align_only` : This switch will cause ShortStack to terminate after the alignment phase; no analysis occurs.
- `--show_secondaries` : If this switch is set, ShortStack will retain secondary alignments for multimapped reads. This will increase bam file size, possibly by a lot.
- `--dicermin DICERMIN` : An integer setting the minimum size (in nucleotides) of a valid small RNA. Together with `--dicermax`, this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. >= 80% of the reads in a given cluster must be in the range indicated by `--dicermin` and `--dicermax`.
    - default: 21
- `--dicermax DICERMAX` : An integer setting the minimum size (in nucleotides) of a valid small RNA. Together with `--dicermin`, this option sets the bounds to discriminate Dicer-derived small RNA loci from other loci. >= 80% of the reads in a given cluster must be in the range indicated by `--dicermin` and `--dicermax`.
    - default: 24
- `--locifile LOCIFILE` : Path to a file of pre-determined loci to analyze. This will prevent *de novo* discovery of small RNA loci. The file may be in gff3, bed, or simple tab-delimited format (Chr:Start-Stop[tab]Name). Mutually exclusive with `--locus`.
- `--locus LOCUS` : A single locus to analyze, given as a string in the format Chr:Start-Stop (using one-based, inclusive numbering). This will prevent *de novo* discovery of small RNA loci. Mutually exclusive with `--locifile`.
- `--nohp` : Switch that prevents search for microRNAs. This saves computational time, but *MIRNA* loci will not be differentiated from other types of small RNA clusters.
- `--dn_mirna` : Switch that activates a *de novo* search for *MIRNA* loci. By default ShortStack will confine *MIRNA* analysis to loci where one or more queries from the `--knownRNAs` file are aligned to the genome. Activating *de novo* searching with `--dn_mirna` does a more comprehensive genome-wide scan for *MIRNA* loci. Loci discovered with `--dn_mirna` that do not overlap already known microRNAs should be treated with caution.
- `--strand_cutoff STRAND_CUTOFF` : Floating point number that sets the cutoff for standedness. Must be > 0.5 and < 1.
    - default: 0.8. Loci with >80% reads on the top genomic strand are '+' stranded, loci with <20% reads on the top genomic strand are '-' stranded, and all others are unstranded '.'
- `--mincov MINCOV` : Minimum alignment depth required to nucleate a small RNA cluster during *de novo* cluster search. Must be an integer >=1. 
    - default: 5
- `--pad PAD` : Initial peaks (continuous regions with depth exceeding argument mincov are merged if they are this distance or less from each other. Must be an integer >= 1. 
    - default: 75

# System Recommendations

# Overview of Methods

# Outputs

# Visualizing Results

# ShortStack version 4 Major Changes

# Issues

# FAQ
