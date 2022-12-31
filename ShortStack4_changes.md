## Changes

- Re-written in `python3`. Requires python >= 3.7
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
- Eliminate option `--total_primaries` .. instead use a hack to rapidly calculate this.
- Option `--locifile` now understands .bed and .gff3 formats, as well as the original simple tab-delimited format.
- Added options `--autotrim` and `--autotrim_key`. This allows automatic detection of 3' adapters by tallying the most common sequence that occurs after a known, highly abundant small RNA (given by `autotrim_key`).
- Most (all) compute-intensive processes are now multi-threaded.

