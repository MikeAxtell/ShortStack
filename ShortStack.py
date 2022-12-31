#!/usr/bin/env python

version = '4.0'

# Main script / control is at the BOTTOM of this file

# Imports from Python Standard Library
import argparse
import time
import os
import re
import pprint
import shutil
import sys
import csv
import multiprocessing as mp
import gzip
import subprocess
import random

from collections import Counter

# Imports of other modules that are NOT in the Standard Library. Must check.
### TO DO .. WRITE A CHECK HERE FOR REQUIRED NON STANDARD MODULES!

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from tqdm import tqdm

def make_outdir(odir):
    """Create the output directory for a ShortStack run.

    Will exit if unsuccessful.
    Requires modules os.
    """
    os.makedirs(odir)

def loc_check(args, fai):
    """Entry point to input and check values for --locifile or --locus.

    Inputs:
    args - the arguments namespace for the run.
    fai - a dictionary from the genome fai index file mapping chr names
     to their lengths.
    Output:
    If neither --locus nor --locifile is set, returns None.
    If --locus is set, returns a dictionary if check is successful.
    If --locifile is set, sends filehandle to a parser function. Returns
    dictionary if successful.
    Requires modules sys, os.
    """
    if (args.locus is None) and (args.locifile is None):
        return None

    # for --locus cases:
    if args.locus is not None:
        failmsg = loc_check_core(args.locus, fai)
        if failmsg is not None:
            sys.exit('Invalid --locus : ' + failmsg)
        else:
            loci = {args.locus: 'User_locus'}
            return loci

    # for --locifile cases need to send it to the appropriate parser
    elif args.locifile is not None:
        # locifile can be .bed, .gff3/.gff, or else assumed to be tab-delimited.
        base, lf_ext = os.path.splitext(args.locifile)
        loci = locifile_parse(args.locifile, fai, lf_ext)
        return loci

def locifile_parse(locifile, fai, ext):
    """Parse and validate a locifile.

    Input:
    locifile - filepath to the loci file (args.locifile)
    fai - dictionary of chromosome names mapped to lengths
    ext - file extension, used to determine gff/gff3 , bed, or other.
    Output:
    loci - dictionary of valid loci with coordinates mapped to names.

    Requires modules csv, re, sys.
    """
    # comment lines, track lines, and browser lines will be silently ignored.
    badline = re.compile('^#|^track\s|^browser\s')
    # patterns used for name parsing from gff files.
    gffID = re.compile('ID=([^;]+)')
    gffName = re.compile('Name=([^;]+)')
    i = 0
    loci = {}
    multiples = {}
    filename = locifile
    with open(locifile) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            # silently ignore empty lines, or badlines
            if not row:
                continue
            if badline.match(row[0]) is not None:
                continue

            # i is a counter for valid-ish rows in the file.
            i += 1

            # determine locus and name. This depends on the file format.
            if(ext == '.bed'):
                # bed minimally must have 3 columns set.
                # bed is 0-based, half-open. We must convert to our standard
                # 1-based, fully-closed coordinate system.
                if len(row) < 3:
                    print('WARNING: row {0}'.format(row) +
                    ' in bed file {0}'.format(filename) +
                    ' is invalid bed format because it has less than 3 columns.')
                    continue
                elif len(row) < 4:
                    # If there are only three columns, need to make our own name.
                    name = 'User_locus_' + str(i)
                    new_start = 1 + int(row[1])
                    locus = str(row[0]) + ':' + str(new_start) + '-' + str(row[2])
                else:
                    # Bed had 4 or more columns, so name is 4th column.
                    name = str(row[3])
                    new_start = 1 + int(row[1])
                    locus = str(row[0]) + ':' + str(new_start) + '-' + str(row[2])
            elif ((ext == '.gff') or (ext == '.gff3')):
                # gff must have 9 columns, no more, no less.
                if len(row) != 9:
                    print('WARNING: row{0}'.format(row) +
                    ' in gff file {0}'.format(filename) +
                    ' is invalid gff because it does not have 9 columns.')
                    continue
                else:
                    # try to find a name from gff3 column 9
                    if gffID.search(str(row[8])) is not None:
                        name = gffID.search(str(row[8]))[1]
                    elif gffName.search(str(row[8])) is not None:
                        name = gffName.search(str(row[8]))[1]
                    else:
                        # defaults to automatic name if ID= or Name= not found.
                        name = 'User_locus_' + str(i)
                    locus = str(row[0]) + ':' + str(row[3]) + '-' + str(row[4])
            else:
                # if not bed or gff, assume a standard ShortStack-y locifile.
                locus = str(row[0])
                if(len(row) > 1):
                    name = str(row[1])
                else:
                    i_str = str(i)
                    name = 'User_locus_' + i_str

            # Check validity of the locus.
            failmsg = loc_check_core(locus, fai)
            if failmsg is not None:
                print('WARNING: locus {0}'.format(locus) +
                ' from locifile {0}'.format(filename) +
                ' is invalid: ' + failmsg)
            elif locus in loci:
                # duplicate entries are only warned at the first duplicate encountered.
                if not locus in multiples:
                    print('WARNING: locus {0}'.format(locus) +
                    ' from locifile {0}'.format(filename) +
                    ' is present more than once. It will only be analyzed once' +
                    ' with first name found, which is {0}'.format(loci[locus]))
                    multiples[locus] = 1
            else:
                loci[locus] = name
    if(loci):
        return loci
    else:
        msg = 'File {0} had no valid loci! Aborting run.'.format(filename)
        sys.exit(msg)

def loc_check_core(locus, fai):
    """Check validity of a locus.

    Checks that a locus is in valid format, that it matches a
    chromosome in the fai index, and that it fits on that chromosome.
    Input:
    locus - string of the locus to be checked
    fai - dictionary of chr names mapped to lengths
    Output:
    msg - An error message if the check fails. If passed, msg is None
    Requires module re.
    """
    loc_pattern = re.compile('(^\S+):(\d+)-(\d+)$')
    if not loc_pattern.match(locus):
        msg = 'Invalid format'
        return msg
    else:
        chr = loc_pattern.match(locus).group(1)
        start = int(loc_pattern.match(locus).group(2))
        stop = int(loc_pattern.match(locus).group(3))

    # chr must be defined in the genome fai
    if chr not in fai:
        msg = 'Chromosome {0} is not defined in genome'.format(chr)
        return msg
    # start must be <= l_stop
    if start > stop:
        msg = 'Start of {0} is greater than stop of {1}'.format(start,stop)
        return msg
    # start must be >= 1
    if start < 1:
        msg = 'Start of {0} is not >= 1'.format(start)
        return msg
    # stop must fit onto the relevant chromosome.
    if stop > fai[chr]:
        msg = ('Stop of {0} is longer than chromosome {1}'.format(stop, chr) +
        ' length of {0}'.format(fai[chr]))
        return msg

def genome_prep(args):
    """Retrieve index information for a genome FASTA.fai file.

    Input:
    args - args namespace object containing parsed arguments.
    Output:
    fai_dict - A dictionary of chr names mapped to lengths.
    Requires modules os, sys, csv, subprocess
    Requires python >= 3.7 to ensure insertion order of the dictionary is
    stable.
    """
    # check for .fai index file. Attempt to create it if it is not there.
    fai = args.genomefile + '.fai'
    if not os.path.isfile(fai):
        # attempt to create fai file
        #pysam.faidx(args.genomefile)
        subprocess.run(f'samtools faidx {args.genomefile}', shell=True)
    if not os.path.isfile(fai):
        msg = ('Could neither find nor create' +
        ' required genomefile index {0}'.format(fai))
        sys.exit(msg)

    # Store the chr names and lengths in a dictionary.
    # Take advantage of stable dictionary key order in python >= 3.7
    fai_dict = {}
    with open(fai) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            chr = str(row[0])
            len = int(row[1])
            fai_dict[chr] = len
    return fai_dict

def check_executables(args):
    """Determine and check for existence of required executables.

    Input:
    args - The args namespace object derived from parsing command line.

    If one or more required executables are not found, exits program.
    Uses modules shutil and sys.
    """
    rqd = []

    # RNAfold and strucVis are needed unless --nohp or --align_only is set
    if not ((args.nohp is True) or (args.align_only is True)):
        rqd.append('RNAfold')
        rqd.append('strucVis')

    # bowtie and bowtie-build are needed if --readfile is is not an empty list.
    # OR if --knownRNAs is present
    if (args.readfile is not None) or (args.knownRNAs is not None):
        rqd.append('bowtie')
        rqd.append('bowtie-build')
    
    # ShortTracks is required if user is aligning reads. So, if --readfile is not empty
    if args.readfile is not None:
        rqd.append('ShortTracks')

    # cutadapt is needed if --autotrim or --adapter is set
    if (args.autotrim is True) or (args.adapter is not None):
        rqd.append('cutadapt')

    # samtools is always required
    rqd.append('samtools')

    for excbl in rqd:
        expath = shutil.which(excbl)
        if expath is not None:
            print('Required executable {0} : {1}'.format(excbl, expath))
        else:
            msg = 'Required executable {0} : Not found!'.format(excbl)
            sys.exit(msg)

def valid_genomefile(x):
    """Tests existence and extension of --genomefile

    x : string input from option --genomefile
    Extension must be .fa or .fasta.
    Uses modules os and sys.
    """
    if not os.path.isfile(x):
        msg = 'Invalid genomefile {0} : is not a file'.format(x)
        sys.exit(msg)
    base,ext = os.path.splitext(x)
    if not (ext == '.fa' or ext == '.fasta'):
        msg = ('Invalid genomefile {0}'.format(x) +
               ' : does not end in .fa or .fasta'.format(x))
        sys.exit(msg)
    return x

def valid_knownRNAs(x):
    """Tests existence and extension of --knownRNAs

    x : string input from option --knownRNAs
    Extension must be .fa or .fasta.
    Uses modules os and sys.
    """
    if not os.path.isfile(x):
        msg = 'Invalid knownRNAs {0} : is not a file'.format(x)
        sys.exit(msg)
    base,ext = os.path.splitext(x)
    if not (ext == '.fa' or ext == '.fasta'):
        msg = ('Invalid knownRNAs {0}'.format(x) +
               ' : does not end in .fa or .fasta'.format(x))
        sys.exit(msg)
    return x

def valid_mincov(x):
    mc = int(x)
    if mc < 1:
        msg = 'Invalid mincov {0} : Must be integer of 1 or more'.format(mc)
        sys.exit(msg)
    return mc

def valid_pad(x):
    pad = int(x)
    if pad < 1:
        msg = 'Invalid mincov {0} : Must be integer of 1 or more'.format(pad)
        sys.exit(msg)
    return pad

def valid_adapter(x):
    """Verifies option --adapter.

    x : string input from option --adapter.
    Adapter must be all ATGC and >= 8 characters long.
    Requires modules re and sys.
    Used during argument parsing.
    """
    dna = re.compile('^[ATCG]+$')
    if not dna.match(x):
        msg = 'Invalid adapter {0} : Non ATGC characters'.format(x)
        sys.exit(msg)
    if len(x) < 8:
        msg = 'Invalid adapter {0} : Must be at least 8 letters'.format(x)
        sys.exit(msg)
    return x

def valid_autotrim_key(x):
    """Verifies option --autotrim_key.

    x : string input from option --autotrim_key.
    Adapter must be 20-30 bases long, and all ATUGCatugc.
    Returns upper-cased, U to T version.
    Requires modules re and sys.
    Used during argument parsing.
    """
    y = x.upper().replace('U', 'T')
    dna = re.compile('^[ATCG]+$')
    if not dna.match(y):
        msg = 'Invalid autotrim_key {0} : Non ATGC characters'.format(y)
        sys.exit(msg)
    if (len(y) < 20) or (len(y) > 30):
        msg = 'Invalid autotrim_key {0} : Must be 20-30 letters long'.format(y)
        sys.exit(msg)
    return y

def valid_threads(x):
    """Verifies option --thread.

    x: string input from option --valid_threads.
    It must be an integer of 1 or higher.
    Used during argument parsing.
    Uses module sys.
    """
    thr = int(x)
    if thr < 1:
        msg = 'Invalid threads {0} : Must be integer of 1 or more'.format(thr)
        sys.exit(msg)
    return thr

def valid_placemax(x):
    """Verifies option --placemax.

    x: string input from option --valid_placemax.
    Must be an integer of 0 or more.
    Used during argument parsing.
    Uses module sys.
    """
    pmx = int(x)
    if not pmx >= 0:
        msg = 'Invalid placemax {0} : Must be integer of 0 or more'.format(x)
        sys.exit(msg)
    return pmx


def valid_dicer(x):
    """Verifies options --dicermin or --dicermax.


    x : string input from options --dicermin or --dicermax.
    Must be an integer of >= 15. Used during argument parsing.
    Uses module sys.
    """
    d = int(x)
    if not d >= 15:
        msg = 'Invalid dicermin or dicermax entry {0} :Must be >= 15.'.format(x)
        sys.exit(msg)
    return d

def valid_locus(x):
    """Verifies option --locus.

    x : string input from option --locus.
    Checks for correct format, and that stop is
    > then start. Used during argument parsing.
    Requires modules re and sys.
    """
    loc_pattern = re.compile('(^\S+):(\d+)-(\d+)$')
    if not loc_pattern.match(x):
        msg='Invalid locus {0} - Must be in format Chr:start-stop'.format(x)
        sys.exit(msg)
    l_start = loc_pattern.match(x).group(2)
    l_stop = loc_pattern.match(x).group(3)
    if not l_stop >= l_start:
        msg = 'Invalid locus {0} - stop must be >= start'.format(x)
        sys.exit(msg)
    return x

def valid_strand_cutoff(x):
    """Verifies option --strand_cutoff.

    x : string input from option --strand_cutoff.
    Must be a float with value
    > 0.5 and < 1. Used during argument parsing.
    Uses module sys.
    """
    sco = float(x)
    if not (sco > 0.5 and sco < 1):
        msg = 'Invalid strand_cutoff {0}: must be float > 0.5 and < 1'.format(x)
        sys.exit(msg)
    return sco

def check_file_list_ext(file_list, exts_list):
    """Check a file list for isfile and proper extensions.

    file_list : list of paths for files.
    exts_list : list of valid extensions.

    Given a list of paths and a list of valid suffixes,
    check to see that a) each file exists and b) each has
    a valid suffix. Requires modules os, sys, and re.
    """

    # compile a regex of the valid exts
    pattern = '$|'.join(exts_list) + '$'
    pattern = pattern.replace('.', '\.')
    ok_exts = re.compile(pattern)

    # check existence of each file, and then its extension
    for fpath in file_list:
        if not os.path.isfile(fpath):
            msg = 'file {0} not found'.format(fpath)
            sys.exit(msg)
        if not ok_exts.search(fpath):
            msg = ('file {0} has an invalid extension. '.format(fpath) +
                   'Valid extensions are {0}.'.format(exts_list))
            sys.exit(msg)

def ShortStack_argparse(ss_version):
    """Command line parsing for ShortStack.

    ss_version : string giving current ShortStack version name.
    returns arguments as a namepsace, if it is successful.
    Requires modules argparse and sys.
    """
    parser = argparse.ArgumentParser(prog='ShortStack')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {0}'.format(ss_version))
    parser.add_argument('--genomefile',
                        required=True,
                        type=valid_genomefile,
                        help='FASTA file of the reference genome (required)')
    parser.add_argument('--knownRNAs',
                        type=valid_knownRNAs,
                        help='FASTA file of known/suspected mature microRNAs')
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument('--readfile',
                              nargs='*',
                              help='One or more files of reads' +
                              ' (fq, fa, gzip OK)')
    parser_group.add_argument('--bamfile',
                              nargs='*',
                              help='One or more BAM alignment files')
    parser.add_argument('--outdir',
                        default = 'ShortStack_' + str(int(time.time())),
                        help='Output directory name. Defaults to' +
                        ' ShortStack_time')
    parser_group3 = parser.add_mutually_exclusive_group()
    parser_group3.add_argument('--adapter',
                        type=valid_adapter,
                        help='3-primer adapter sequence to trim off of' +
                        'reads. If given applies to all input fastq files.' +
                        ' Mutually exclusive with --autotrim.')
    parser_group3.add_argument('--autotrim',
                        action='store_true',
                        help='If this switch is set, automatically discover' +
                        ' the 3-prime adapter from each input readfile, and ' +
                        'trim it. This uses the sequence from --autotrim_key' +
                        ' to discover the adapter sequence. Mutually' +
                        ' exlcusive with --adapter.')
    parser.add_argument('--autotrim_key',
                        type = valid_autotrim_key,
                        default = 'TCGGACCAGGCTTCATTCCCC',
                        help = 'Sequence of an abundant, known small RNA' +
                        ' to be used to discover the 3-prime adapter' +
                        ' sequence. Has no effect unless --autotrim is' +
                        ' specified. Defaults to TCGGACCAGGCTTCATTCCCC' +
                        ' (miR166). Can be upper or lower-case, T or U' +
                        ' and must be 20-30 bases long.' )
    parser.add_argument('--threads',
                        type=valid_threads,
                        default=1,
                        help='Number of threads to use (integer) - default: 1')
    parser.add_argument('--mmap',
                        choices=['u', 'f', 'r'],
                        default='u',
                        help='Protocol for multi-mapped reads: u, f, or r' +
                        ' - default: u')
    parser.add_argument('--align_only',
                        action='store_true',
                        help='If this switch is set, ShortStack quits after' +
                        ' performing alignments without any analyses' +
                        ' performed.')
    parser.add_argument('--show_secondaries',
                        action='store_true',
                        help='If this switch is set, ShortStack will retain' +
                        ' secondary alignments for multimapped reads.' +
                        ' This will' +
                        ' increase bam file size, possibly by a lot.')
    parser.add_argument('--dicermin',
                        default=21,
                        type=valid_dicer,
                        help='Minimum size of a valid Dicer-processed ' +
                        'small RNA.' +
                        ' Must be integer >= 15 and <= dicermax. Default: 20.')
    parser.add_argument('--dicermax',
                        default=24,
                        type=valid_dicer,
                        help='Maximum size of a valid Dicer-processed' +
                        ' small RNA.' +
                        ' Must be integer >= 15 and <= dicermax. Default: 24.')
    parser_group2 = parser.add_mutually_exclusive_group()
    parser_group2.add_argument('--locifile',
                               help='File listing intervals to analyze.' +
                               ' Can be simple tab-delimited, .bed, or .gff3.' +
                               ' Tab-delimited format is column 1 with' +
                               ' coordinates Chr:start-stop, column 2 with' +
                               ' names. Input file assumed to be simple tab-' +
                               'delimited unless file name ends in .bed or' +
                               ' .gff3.' +
                               ' Mutually exclusive with --locus.')
    parser_group2.add_argument('--locus',
                               type=valid_locus,
                               help='Analyze the specified interval,' +
                               ' given in' +
                               ' format Chr:start-stop. Mutually ' +
                               ' exclusive with' +
                               ' --locifile.')
    parser.add_argument('--nohp',
                        action='store_true',
                        help='If this switch is set, RNA folding will' +
                        ' not take' +
                        ' place, thus MIRNA loci cannot be ' +
                        ' annotated. This' +
                        ' does however save CPU time.')
    parser.add_argument('--dn_mirna',
                        action='store_true',
                        help='If this switch is set, a de novo' +
                        ' search for new MIRNA loci will be performed.' +
                        ' By default, de novo MIRNA finding is not performed ' +
                        ' and MIRNA searches are limited to loci matching' +
                        ' RNAs from --knownRNAs that align to the genome')
    parser.add_argument('--strand_cutoff',
                        type=valid_strand_cutoff,
                        default = 0.8,
                        help='Cutoff for calling the strandedness of a ' +
                        'small RNA' +
                        ' locus. Must be a floating point > 0.5 and < 1.' +
                        ' Default: 0.8.')
    parser.add_argument('--mincov',
        type=valid_mincov,
        default=5,
        help='Minimum alignment depth required to nucleate a ' +
        'small RNA cluster during de novo cluster search.' +
        ' Must be an integer >=1. Default: 5')
    parser.add_argument('--pad',
        type=valid_pad,
        default=200,
        help='Initial peaks (continuous regions with depth '+
        'exceeding argument mincov are merged if they are '+
        'this distance or less from each other. Must be an '+
        'integer >= 1. Default: 75')

    ### Options in perl ShortStack that are not yet implemented here

    # --sort_mem  ... ignored for now. Is it needed?
    # --mismatches ... re-visit after writing alignment section
    # --foldsize ... re-visit when analysis is constructed

    # Parse the arguments
    args = parser.parse_args()

    # check bamfile(s) or readfiles(s)
    if args.readfile is not None:
        readfile_exts = ['.fq', '.fastq', 
                         '.fq.gz', '.fastq.gz']
        check_file_list_ext(args.readfile, readfile_exts)

    if args.bamfile is not None:
        bamfile_exts = ['.bam']
        check_file_list_ext(args.bamfile, bamfile_exts)

    # check dicermin relative to dicermax
    if not args.dicermax >= args.dicermin:
        msg = ('dicermax of {0}'.format(args.dicermax) +
               ' is not >= dicermin of {0}'.format(args.dicermin))
        sys.exit(msg)
    return args

def start_message(args, version):
    """Report settings to user at start of ShortStack run.

    args : The arguments namespace object output by ShortStack_argparse.
    version : string giving current ShortStack version name.
    Requires modules time and pprint.
    """
    print('\nShortStack version {0}\n'.format(version))
    print('Beginning run')
    print('Options:')
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(vars(args))

def align(args, fai, trimmed_files):
    """Entry point for ShortStack's alignment procedures.

    Input:
    args - namespace object from startup.
    fai - Dictionary of chr names and lengths from the user's genome.
    trimmed_files - list of trimmed files. May be empty if user provided
        bamfiles instead

    Output:
    bam - PATH to the final bamfile to be used for analysis.

    Requires:
    multiprocessing as mp
    from collections import Counter
    subprocess
    """
    if args.bamfile is not None:
        # No alignment. But we need to check the user's provided file(s)
        bam = check_user_bamfile(args, fai)
        # Could still need bowtie indices if knownRNAs is provided
        if args.knownRNAs is not None:
            bowtie_indices(args)
        return bam
    else:
        # Alignments imminent. Need to check for bowtie indices.
        bowtie_indices(args)
    
    print('')
    print('Beginning alignment phase')  
    final_bams = []
    # Initiate a tsv output file with details of mappings
    map_rpt_file = init_map_rpt(args)

    for t_file in trimmed_files:
        print('')
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print(f'Aligning {t_file}')
        print(f'First pass alignment with bowtie using {str(args.threads)} threads')
        (rs1_sam_files, al_dens) = al_pass_1(t_file, args)

        # pass 2 - place multimappers using multiprocessing.
        print(f'Second pass - placing multimappers using {str(args.threads)} threads' +
        f' to process {len(rs1_sam_files)} chunks')
        #print(' Chunk numbers in process / finished:', end='')
        if __name__ == '__main__':
            # Use sam chunks of <= 1M reads for multiprocessing
            # make list of inputs for multiprocessing.starmap
            rs2_iter = list(zip(rs1_sam_files, [al_dens] * len(rs1_sam_files),
                [args] * len(rs1_sam_files), 
                list(range(1, len(rs1_sam_files) + 1))))
            with mp.Pool(args.threads) as pool:
                rs2_results = pool.starmap(al_pass_2, 
                tqdm(rs2_iter, total=len(rs1_sam_files), desc='Chunks sent', 
                unit=' Chunk', leave=None))
            rs2_files, all_counters = zip(*rs2_results)
            # merge counters
            lib_counts = Counter()
            for counter in all_counters:
                lib_counts = lib_counts + counter
            # concatenate rs2 sam files
            cat_args = ['cat']
            cat_args.extend(list(rs2_files))
            #rs2_sam = rs2_files[0].replace('_1_rs2.sam', '_rs2.sam')
            rs2_sam = '_rs2.sam'.join(rs2_files[0].rsplit('_1_rs2.sam', 1))
            rs2_sam_fh = open(rs2_sam, "w")
            subprocess.run(cat_args, stdout=rs2_sam_fh)
            rs2_sam_fh.close()
            # clean up files
            rm_args = ['rm', '-f']
            rm_args.extend(list(rs2_files))
            subprocess.run(rm_args)
        else:
            sys.exit('FATAL: function align not called as __main__') 
        
        print('')
        print ('Converting to sorted bam format')
        # Convert to bam
        # rs2_bam = rs2_sam.replace('.sam', '.bam') 
        rs2_bam = '.bam'.join(rs2_sam.rsplit('.sam', 1))
        bam_conv_args = (['samtools', 'view', '-b', '--threads',
           str(args.threads), '--output', rs2_bam, rs2_sam])
        rs2_bam_sp = subprocess.run(bam_conv_args)
        if rs2_bam_sp.returncode != 0:
            sys.exit('FATAL: samtools view encountered an error!')
        rm_args = ['rm', '-f', rs2_sam]
        subprocess.run(rm_args)

        # sort
        final_bam = '.bam'.join(rs2_bam.rsplit('_rs2.bam', 1)) 
        bsort_args = (['samtools', 'sort', '--threads', str(args.threads),
            '-o', final_bam, rs2_bam])
        bsort_sp = subprocess.run(bsort_args)
        if bsort_sp.returncode != 0:
            sys.exit('FATAL: samtools sort encountered an error!')
        rm_args = ['rm', '-f', rs2_bam]
        subprocess.run(rm_args)
        
        # index
        idx_args = (['samtools', 'index', '-@', str(args.threads), 
            final_bam])
        idx_sp = subprocess.run(idx_args)
        if idx_sp.returncode != 0:
            sys.exit('FATAL: samtools index encountered an error!')
        
        # Report, both a simple report to stdout and a detailed one to a file.
        mapping_rpt(lib_counts, map_rpt_file, final_bam)

        # add file name to list
        final_bams.append(final_bam)

    # If there is more than one entry in final_bams, need to merge and index
    if len(final_bams) > 1:
        # do stuff
        print('')
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print('Merging and indexing alignments') 
        mb = args.outdir + '/merged_alignments.bam'
        merge_args = (['samtools', 'merge', '--threads', str(args.threads), 
            '-o', mb, '-r'])
        merge_args.extend(final_bams)
        subprocess.run(merge_args)
        idx_args = ['samtools', 'index', '-@', str(args.threads), mb]
        subprocess.run(idx_args)

        # ShortTracks x 2 : mode readgroup and mode readlength
        print('')
        print('Creating browser tracks by readgroup using ShortTracks')
        callShortTracks(mb, 'readgroup', False)
        print('')
        print('Creating browser tracks by readlength and strand using ShortTracks')
        callShortTracks(mb, 'readlength', True)

        return mb
    else:
        # ShortTracks x 1 : mode readlength
        print('')
        print('Creating browser tracks by readlength and strand using ShortTracks')
        callShortTracks(final_bams[0], 'readlength', True)
        return final_bams[0]

def callShortTracks(bamfile, mode, stranded):
    cmd = f'ShortTracks --mode {mode} --bamfile {bamfile}'
    if stranded is True:
        cmd = cmd + ' --stranded'
    subprocess.run(cmd, shell=True)

def mapping_rpt(lib_counts, map_rpt_file, bamfile):
    """Produce a mapping report from a lib_counts counter object

    A simple report is sent to stdout
    Details in tsv form are appended to the map_rpt_file
    
    Requires: Counter from collections
    """

    prefixes = ['U', 'P', 'R', 'H', 'N']
    suffixes = ['<21', str(21), str(22), str(23), str(24), '>24']
    ks = []
    read_total = sum(lib_counts.values()) 
    with open(map_rpt_file, "a") as fh:
        for p in prefixes:
            for s in suffixes:
                k = p + s
                o = (bamfile + '\t' + p + '\t' + s + '\t' + 
                    str(lib_counts[k]) + '\n')
                fh.write(o)
    u_sum = 0
    for s in suffixes:
        k = 'U' + s
        u_sum += lib_counts[k]
    u_perc = round(100 * (u_sum / read_total), 1)
    print(f'Uniquely mapped (U) reads: {u_sum}/{read_total} ({u_perc}%)')

    p_sum = 0
    for s in suffixes:
        k = 'P' + s
        p_sum += lib_counts[k]
    p_perc = round(100 * (p_sum / read_total), 1)
    print('Multi-mapped reads placed (P) with guidance:' +
        f' {p_sum}/{read_total} ({p_perc}%)')
    
    r_sum = 0
    for s in suffixes:
        k = 'R' + s
        r_sum += lib_counts[k]
    r_perc = round(100 * (r_sum / read_total), 1)
    print('Multi-mapped reads randomly (R) placed:' +
        f' {r_sum}/{read_total} ({r_perc}%)')
    
    h_sum = 0
    for s in suffixes:
        k = 'H' + s
        h_sum += lib_counts[k]
    h_perc = round(100 * (h_sum / read_total), 1)
    print('Very highly (H) multi-mapped reads (>=50 hits):' +
        f' {h_sum}/{read_total} ({h_perc}%)')
 
    n_sum = 0
    for s in suffixes:
        k = 'N' + s
        n_sum += lib_counts[k]
    n_perc = round(100 * (n_sum / read_total), 1)
    print('Not mapped (N) reads (no hits):' +
        f' {n_sum}/{read_total} ({n_perc}%)')
 

def init_map_rpt(args):
    """Create and add headers to a file reporting the results of mapping

    """
    map_rpt_file = args.outdir + '/alignment_details.tsv'
    topline = 'readfile\tmapping_type\tread_length\tcount'
    fh = open(map_rpt_file, "w")
    print(topline, file=fh)
    fh.close()
    return map_rpt_file

def al_pass_2(rs1_file, al_dens, args, n):
    """Placement of multimappers from a SAM file.
    
    Requires Counter from collections
    """
    #rs2_sam = rs1_file.replace('_rs1.sam', '_rs2.sam')

    #print(f' {str(n)}', end='')

    rs2_sam = '_rs2.sam'.join(rs1_file.rsplit('_rs1.sam', 1))
    rs2_sam_fh = open(rs2_sam, "w")
    
    # Get an empty counter for this sam file
    lib_counters = Counter()

    # Count actual lines in the rs1 file, for progress tracker
    #wc = subprocess.run(["wc", "-l", rs1_file],
    #capture_output = True, text = True)
    #wc_pat = re.compile('(^\d+)')
    #wc_num = int(wc_pat.match(wc.stdout).group(1))

    # Read in contents of file to memory
    with open(rs1_file) as fh:
        samlines = fh.readlines()
    
    ## line by line, accumulating alignment blocks
    last_read = None
    alignments = []
    #with open(rs1_file) as rs1_fh:
    for samline in samlines:
        if samline.startswith('@'):
            # headers, just print
            print(samline, end='', file=rs2_sam_fh)
        else:
            # not a header.
            samfields = samline.rstrip().split("\t")
            readname = samfields[0]
            # Have we passed to a new readname?
            if readname != last_read and last_read is not None:
                # Process that last block
                lib_counters = al_block_process2(lib_counters,
                    al_dens, alignments, args, rs2_sam_fh)
                    
                # reset
                alignments = []

            # set last_read and add
            last_read = readname
            alignments.append(samline)

    # Process last block of alignment
    lib_counters = al_block_process2(lib_counters,
        al_dens, alignments, args, rs2_sam_fh)
    rs2_sam_fh.close()

    # clean up the rs1 file
    rm_rs1_args = ['rm', '-f', rs1_file]
    # Uncomment the below after testing
    subprocess.run(rm_rs1_args)

    return(rs2_sam, lib_counters)
 
def al_pass_1(t_file, args):
    """First pass aligment using bowtie.

    Input:
    - t_file : file path to trimmed fastq(gz) file of reads to align.
    - args : The ShortStack arguments object.

    Output:
    - rs1_sam_files : a list of read-sorted SAM files.
    - al_dens : alignment density dictionary
    
    Requires:
    os
    sys
    subprocess
    re
    """
    # file paths and handles
    (head, tail) = os.path.split(t_file)
    (root, ext) = os.path.splitext(tail)

    # guard against extensions like .fastq.gz or .fq.gz ...
    #  we don't want to name files with .fastq or .fq if they
    #  aren't that format!

    fqlike = re.compile('\.fq$|\.fastq$')
    if fqlike.search(root):
        (root2, ext2) = os.path.splitext(root)
        final_root = root2
    else:
        final_root = root

    rs1_counter = 1
    rs1_sam = args.outdir + '/' + final_root + '_' + str(rs1_counter) + '_rs1.sam'
    rs1_sam_fh = open(rs1_sam, "w")

    # bowtie arguments
    bt_args = (['bowtie', '-p', str(args.threads), '-v', '1',
        '-k', '50', '-S', '--best', '--strata', '-x',
        args.genomefile, t_file])
    
    # bowtie subprocess, piping stdout 
    bt1 = subprocess.Popen(bt_args, stdout=subprocess.PIPE, 
        text=True, bufsize=1)
        
    # holder for alignment(s) of a single read
    alignments = []

    # Dictionary storing read densities in 50nt bins
    al_dens = {}
    
    # Other holder and counters
    last_read = None
    read_count = 0
    rs1_sam_files = [rs1_sam]

    # parse alignments as they come, and print to SAM
    for samline in tqdm(bt1.stdout, unit=' bowtie output lines',
    leave=None):
        # Processing alignment lines
        if not samline.startswith('@'):
            # split the line
            samfields = samline.rstrip().split("\t")
            readname = samfields[0]
            # Have we passed to a new read? 
            if readname != last_read:
                # Process the last block of alignments
                al_dens = al_block_process(al_dens, alignments, args, 
                    rs1_sam_fh)

                # reset
                alignments = []

                #check read count. Move to new SAM file if it is 1 million
                read_count += 1
                if read_count >= 1000000:
                    rs1_sam_fh.close()
                    rs1_counter += 1
                    rs1_sam = (args.outdir + '/' + root + '_' + 
                        str(rs1_counter) + '_rs1.sam')
                    rs1_sam_fh = open(rs1_sam, "w")
                    read_count = 0
                    rs1_sam_files.append(rs1_sam)

            #set last read name encountered and add to alignments
            last_read = readname
            alignments.append(samline)

        # header lines fall here, just print them
        else:
            print(samline, end='', file=rs1_sam_fh)

    # Process the last block of alignments
    al_dens = al_block_process(al_dens, alignments, args, 
        rs1_sam_fh)

    # ensure subprocess completed
    bt1.wait()

    # clean up and check returncode
    rs1_sam_fh.close()
    if bt1.returncode != 0:
        sys.exit('FATAL: bowtie encountered an error!')

    # return data
    return(rs1_sam_files, al_dens)

def check_user_bamfile(args, fai):
    """Merging (if needed) and validation of user-provided bamfile(s).

    Input:
    args - argument namespace from startup.
    fai - dictionary of chr names and lengths from startup.

    Output:
    bam - validated single bamfile path.

    Requires subprocess, os, sys.
    """
    # Did the user provide multiple bamfiles? If so, merge them, keeping
    #  read groups inferred from file names.
    if len(args.bamfile) > 1:
        # in order to use pysam.merge() , have to write file of bamfile paths.
        # also for -b option of samtools merge
        bflist = args.outdir + '/user_bamfiles.txt'
        f = open(bflist, 'w')
        for bf in args.bamfile:
            f.write(os.path.abspath(bf) + "\n")
        f.close()

        # Merge the files, then make sure file exists.
        bam = args.outdir + '/merged.bam'
        #pysam.merge("-r", "-b", bflist, bam)
        subprocess.run(f'samtools merge -b {bflist} {bam}')
        if os.path.exists(bam) is False:
            sys.exit("Merging of the bamfiles provided by user failed!")
    else:
        # here the user provided a single bamfile.
        bam = args.bamfile[0]

    # index the file. samtools index won't complain if an index already exists.
    #pysam.index(bam)
    subprocess.run(f'samtools index {bam}', shell=True)

    return bam

######## end of old align.py functions

def al_block_process2(lib_counters, al_dens, alignments, args, rs_sam2_fh):
    """Decide on multimappers and output final SAM lines

    Input:
    - lib_counters : Dictionary of counters for XY and YS tags
    - al_dens : read density for guidance .. can be empty for method r, n
    - alignments : list of samlines all from same read
    - args : the ShortStack argument object
    - rs_sam2_fh : file handle for SAM output

    Output:
    - lib_counters : updated
    
    The XX, XY, XZ, and YS tags are added to the samlines as well

    Requires module random

    """
    hits = len(alignments)
    xx_i = 'XX:i:' + str(hits)
    # One hits are easy no matter what
    if hits == 1:
        # Could be mapped or unmapped
        samfields = alignments[0].rstrip().split('\t')
        ys_z_v = get_ys_tag_v(samfields[9])
        ys_z = 'YS:Z:' + ys_z_v
        xz_f = 'XZ:f:' + '1.000'
        if int(samfields[1]) & 4:
            # Unmapped
            xy_z_v = 'N'
            xy_z = 'XY:Z:' + xy_z_v
            xx_i = 'XX:i:' + str(0)
        else:
            # Mapped
            xy_z_v = 'U'
            xy_z = 'XY:Z:' + xy_z_v
        # increment counters
        # testing
        k = xy_z_v + ys_z_v
        #print(k)
        #print(lib_counters)
        lib_counters[k] += 1
        # modify the alignment and output to SAM
        al_out = (alignments[0].rstrip() + '\t' + xx_i + '\t' +
            xy_z + '\t' + xz_f + '\t' + ys_z + '\n')
        # print to file
        print(al_out, end='', file=rs_sam2_fh)
    else:
        # More than one alignment. Time to decide!
        ## Calc probability of each read
        probs = []
        if args.mmap == 'r':
            # easy
            for x in range(0, len(alignments)):
                probs.append(round(1/hits, 3))
            xy_z_v = 'R'
            xy_z = 'XY:Z:' + xy_z_v
        else:
            (probs, xy_z_v) = get_mmap_probs(alignments, al_dens)
            xy_z = 'XY:Z:' + xy_z_v

        ## Make a choice, using random.choices() with weighted choice
        chosen = random.choices(alignments, weights=probs)[0]

        ## Output .. are secondaries being output?
        for idx, al in enumerate(alignments):
            if (al == chosen) or (args.show_secondaries is True):
                # output this one
                samfields = al.rstrip().split('\t')
                ## set remaining tags
                xz_f = 'XZ:f:' + str(probs[idx])
                ys_z_v = get_ys_tag_v(samfields[9])
                ys_z = 'YS:Z:' + ys_z_v
                # increment counters, unless secondary
                if al == chosen:
                    k = xy_z_v + ys_z_v
                    lib_counters[k] += 1

                # add flag for secondaries
                if al != chosen:
                    old_flag = int(samfields[1])
                    new_flag = str(old_flag + 256)
                    samfields[1] = new_flag
                # Build the output
                al_out = ('\t'.join(samfields) + '\t' + xx_i + '\t' +
                    xy_z + '\t' + xz_f + '\t' + ys_z +'\n')
                # print
                print(al_out, end = '', file=rs_sam2_fh)

    return lib_counters    


 
def get_mmap_probs(alignments, al_dens):
    """Compute multi-mapped read probabilities at each of the possible sites

    Uses the alignment densities (al_dens) calculated earlier.

    Returns a list of probabilities  and an xy_z value or either R or P.

    """
    counts = []
    total = 0
    for al in alignments:
        samfields = al.rstrip().split('\t')
        center_bin = int(samfields[3]) // 50
        count = 0
        for x in range(center_bin - 2, center_bin + 2):
            dkey = samfields[2] + ':' + str(x)
            if dkey in al_dens:
                count += al_dens[dkey]
        counts.append(count)
        total += count
    if total > 0:
        probs = [round(x / total, 3) for x in counts] 
        if len(alignments) < 50:
            xy_z_v = 'P'
        else:
            xy_z_v = 'H'
    else:
        # No guidance possible. It is a random choice
        probs = []
        for x in range(0, len(alignments)):
            probs.append(round(1/len(alignments), 3))
        if len(alignments) < 50:
            xy_z_v = 'R'
        else:
            xy_z_v = 'H'
    return (probs, xy_z_v)

def get_ys_tag_v(seq):
    """Given SEQ, return the entry for a YS tag.

    YS:Z tags indicate small RNA length.
    Can be '<21', '21', '22', '23', '24', or '>24'.
    All strings.

    Returns the string
    """
    seqlen = len(seq)
    if seqlen < 21:
        ys_tag_v = '<21'
    elif seqlen == 21:
        ys_tag_v = str(21)
    elif seqlen == 22:
        ys_tag_v = str(22)
    elif seqlen == 23:
        ys_tag_v = str(23)
    elif seqlen == 24:
        ys_tag_v = str(24)
    else:
        ys_tag_v = '>24'
    return ys_tag_v

def al_block_process(al_dens, alignments, args, rs_sam_fh):
    """Analyze and print a block of alignments for one read

    Input: 
    - al_dens : A dictionary of read densities in 50nt bins
    - alignments : A list of samlines, including the newline character
    - args : The ShortStack argument object
    - rs_sam_fh : An open file handle for printing SAM lines to

    Output:
    - al_dens : Updated dictionary

    """
    hits = len(alignments)
    xx_alignments = []
    for al in alignments:
        samfields = al.rstrip().split('\t')
        # Store in density dictionary if 1 hit or if method f so long as
        #  read was truly mapped
        if ((hits == 1) or (args.mmap == 'f')) and not int(samfields[1]) & 4:
            readbin = int(samfields[3]) // 50
            dkey = samfields[2] + ':' + str(readbin)
            to_add = 1 / hits
            al_dens[dkey] = al_dens.get(dkey, 0) + to_add

        print(al, end='', file=rs_sam_fh)

    return al_dens


def find_adapter(fqfile, prefix, n):
    """Given path to fastq/fastq.gz and a prefix, find best 3' adapter.

    Inputs:
    - fqfile : path to fastq or gzip-compressed fastq file
    - prefix : sequence to use a 'prefix' .. a common microRNA
    - n : Integer for use by tqdm for positioning progress tracker

    Output:
    - adapter : Up to 20nt long best adapter sequence.
    
    Depends on FastqGeneralIterator from Bio.SeqIO.QualityIO
    """
    prefixcount = 0
    seqcount = 0
    aseqs = {}
    if(re.search("gz$", fqfile)):
        fqhandle = gzip.open(fqfile, "rt")
    else:
        fqhandle = open(fqfile, "rt")
    for title, seq, qual in tqdm(FastqGeneralIterator(fqhandle),
        position=n, leave=None, unit_scale=0.25, desc=fqfile,
        unit='reads', mininterval=1):
        seqcount += 1
        if(seq.startswith(prefix)) :
            prefixcount += 1
            aseq = seq[len(prefix):]
            if(len(aseq) > 20) :
                aseq = aseq[0:20]
            if aseq in aseqs:
                aseqs[aseq] += 1
            else:
                aseqs[aseq] = 1
    fqhandle.close()
    #print("\t\tPrefix", prefix, "found", prefixcount, "times out of", seqcount,
    #    " records.")
    if prefixcount > 0:
        bestaseq = sorted(aseqs, key=aseqs.get, reverse=True)[0]
        bestperc = round((aseqs[bestaseq] / prefixcount) * 100, 1)
        #print("\t\tMost frequent adapter (first 20nts):", bestaseq, "found",
        #    aseqs[bestaseq], "times (", bestperc, '%)')
        ###print (f'For {fqfile} adapter {bestaseq} will be used for trimming.')
        return (fqfile, bestaseq) 
        
    else:
        #print(f'Autotrim failed to find an adapter for {fqfile}!')
        #print(f' Consider a different autotrim_key!') 
        # We don't sys.exit() here because this function may be running
        #  as one of several threads. Main script should abort instead.
        return (fqfile, None)

def trimmer(fqfile, adapter, args):
    """Trim fastq files using cutadapt with standardized settings

    Inputs:
    - fqfile : path to the untrimmed fastq(gz) file
    - adapter : 3' adapter sequence to be used by cutadapt
    - args : ShortStack arguments

    Output: filepath to trimmed fastq(gz) file

    Requires : 
    - cutadapt (external program on PATH)
    - os
    - subprocess
    - sys

    """
    # 1. determine trimmed file path (in the SS outdir)
    (head, tail) = os.path.split(fqfile)
    tfile = args.outdir + '/t' + tail
    
    # 2. build the cutadapt call
    ca_args = (['cutadapt', '-a', adapter, '-o', tfile, '-j', str(args.threads),
    '--discard-untrimmed', '-m', '12', '--report', 'minimal', fqfile
    ])
    
    # 3. execute the cutadapt job
    print('')
    print(f'Trimming readfile {fqfile} ...')
    #ca = subprocess.run(args=ca_args, capture_output=True, text=True)
    ca = subprocess.run(args=ca_args)
    if ca.stdout is not None:
        print(ca.stdout)
    if ca.stderr is not None:
        print(ca.stderr)
    if ca.returncode == 0:
        print(f'  Completed. Trimmed file is {tfile}.')
        return tfile
    else:
        sys.exit('  FATAL: adapter trimming with cutadapt failed!')    

def bowtie_indices(args):
    """Check for and if necessary build bowtie indices

    Inputs:
    - args : The ShortStack argument object

    Output:
    - nothing, but dies if a failure occurs

    Requires:
    - bowtie-build (external)
    - os
    - subprocess
    - sys

    """
    # Expect .1.ebwt through .4.ebwt and .rev.1.ebwt and .rev.2.ebwt
    #  Could also all be ebwtl, with an 'l' on the end. But always 6 files.

    set1 = ([args.genomefile + '.1.ebwt', args.genomefile + '.2.ebwt',
        args.genomefile + '.3.ebwt', args.genomefile + '.4.ebwt',
        args.genomefile + '.rev.1.ebwt', args.genomefile + '.rev.2.ebwt'
    ])
    
    set2 = ([args.genomefile + '.1.ebwtl', args.genomefile + '.2.ebwtl',
        args.genomefile + '.3.ebwtl', args.genomefile + '.4.ebwtl',
        args.genomefile + '.rev.1.ebwtl', args.genomefile + '.rev.2.ebwtl'
    ])
    set1ok = True
    set2ok = True
    for file1 in set1:
        if not os.path.exists(file1):
            set1ok = False
    for file2 in set2:
        if not os.path.exists(file2):
            set2ok = False
    if (not set1ok) and (not set2ok):
        # indices need to be built
        # report to user
        print('')
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print('Required bowtie indices not found. Building them ...')
         
        bb_args = (['bowtie-build', '--threads', str(args.threads), '-q',
            args.genomefile, args.genomefile  
        ])
        bb = subprocess.run(args=bb_args, capture_output=True, text=True)
        if bb.stdout is not None:
            print(bb.stdout)
        if bb.stderr is not None:
            print(bb.stderr)
        if(bb.returncode == 0):
            print('    Completed')
        else:
            sys.exit('    FATAL: bowtie-build did not complete successfully')

def get_trimmed(args):
    """Read trimming including finding adapter sequences

    Inputs:
    - args : The ShortStack argument object

    Output:
    - trimmed_files : List of trimmed read files

    Requires:
    - time
    - multiprocessing as mp
    - sys

    Notes:
    Must be run as __main__ because of multiprocessing

    """
    if __name__ == '__main__':
        found_adapters = []
        if (args.readfile is not None) and (args.autotrim is True):
            # auto-detect the 3' adapters in each readfile, using as many
            # threads as the user lets us.
            print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
            print('')
            print(f'Initiating autotrim with ' +
                f'an autotrim_key of {args.autotrim_key}')
            
            # tqdm set up
            tqdm.set_lock(mp.RLock())

            with mp.Pool(processes=args.threads,
                initializer=tqdm.set_lock, 
                initargs=(tqdm.get_lock(),)) as pool:
                # need to iterate over readfile but also send args.autotrim_key
                # Numbers are for tqdm
                autotrim_iter = list(zip(args.readfile,
                    [args.autotrim_key] * len(args.readfile),
                    list(range(1, len(args.readfile) + 1))))
                for found_adapter in pool.starmap(find_adapter, autotrim_iter):
                    if(all(found_adapter)):
                        found_adapters.append(found_adapter)
                    else:
                        msg = ('FATAL: Failed to infer an adapter from one or' +
                            ' more fastq files!. Try a different autotrim_key' +
                            ' or trim your reads outside of ShortStack!')
                        sys.exit(msg)
        print('')
        print('Adapters:')
        pprint.pprint(found_adapters)
        
        # Adapter trimming
        trimmed_files = []

        if (args.readfile is not None) and (found_adapters):
            print('')
            print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
            print('Initiating read trimming with cutadapt')
            # method using found_adapters .. from autotrim 
            # found_adapters is a list of tuples; each tuple is fqfile, adapter
            for file_adapter in found_adapters:
                trimmed_file = trimmer(file_adapter[0], file_adapter[1], args)
                trimmed_files.append(trimmed_file)
        elif (args.readfile is not None) and (args.adapter is not None):
            print('')
            print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
            print('Initiating read trimming with cutadapt')
            # method using a single user-provided adapter
            for read_file in args.readfile:
                trimmed_file = trimmer(read_file, args.adapter, args)
                trimmed_files.append(trimmed_file)
        elif args.readfile is not None:
            # no trimming at all .. trimmed_files are the user's readfile
            trimmed_files = args.readfile
        return trimmed_files
    else:
        sys.exit('Fatal: function get_trimmed must be called as __main__')

def denovo_chrom(chrom, bamfile, args):
    #print(f'chrom {chrom} bamfile {bamfile}')
    #print(f'Starting to work on chrom {chrom}')
    st_depth_args = (['samtools', 'depth', '-a', '-r', chrom,
        '-g', '256', bamfile])
    st_depth = subprocess.Popen(st_depth_args, stdout=subprocess.PIPE,
        text=True, bufsize=1)
    reader = csv.reader(st_depth.stdout, delimiter="\t")
    depth = []
    for row in reader:
        depth.append(int(row[2]))
    #print(f'Ingested raw data for chrom {chrom}')
    depth_np = np.array(depth)
    #print('Created np array')
    
    # Get the array indices that correspond to depth positions
    #  above the threshold. These correspond to 0-based chromosome
    #  coordinates.
    high_pos = np.where(depth_np >= args.mincov)[0]

    # Find the consecutive groups
    # https://stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-in-a-numpy-array
    consec = np.split(high_pos, np.where(np.diff(high_pos) != 1)[0]+1)
    
    # Write as bed3 file
    peaksbed = args.outdir + '/Ref_' + chrom + '_peaks.bed'
    with open(peaksbed, "w") as fh:
        for x in consec:
            bstart = x[0] # bed start zero-indexed
            bend = x[-1] + 1 # bed end does NOT include interval (weird bed spec)
            o = str(chrom) + '\t' + str(bstart) + '\t' + str(bend) + '\n'
            fh.write(o)

    # return file name
    return peaksbed
    
def denovo_clusters(bamfile, args, fai, loci):
    """Find clusters based on depth threshold and merging
    Inputs:
    - bamfile : path to indexed bamfile
    - args : the shortstack args object
    - fai : dictionary of chromosome names and lengths
    - loci : dictionary of user-provided loci (might be None)

    Requires:
    - mutliprocessing as mp
    - re

    """

    if loci is not None:
        # Write user-provided loci in bed4 format
        user_bed = args.outdir + '/user.bed'
        ufh = open(user_bed, "w")
        for loc in loci:
            loc_pattern = re.compile('(^\S+):(\d+)-(\d+)$')
            chr = loc_pattern.match(loc).group(1)
            start = int(loc_pattern.match(loc).group(2))
            stop = int(loc_pattern.match(loc).group(3))
            bstart = start - 1
            ufh.write(str(chr) + '\t' + str(bstart) + '\t' + str(stop) + '\t' + str(loci[loc]) + '\n')
        ufh.close()
        return user_bed
    elif __name__ == '__main__':
        # Report to user
        print('')
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print('Defining small RNA clusters de novo')

        # loop through each chromosome using multithreading
        chroms = []
        for chrom in fai:
            chroms.append(chrom)
        #chroms = fai.keys()
        dn_iter = list(zip(chroms, [bamfile] * len(chroms),
            [args] * len(chroms)))
        #print(f'chroms: {chroms} dn_iter: {dn_iter}')
        with mp.Pool(args.threads) as pool:
            Ref_peak_beds = pool.starmap(denovo_chrom, dn_iter)
    else:
        sys.exit('FATAL: function align not called as __main__')
    
    # cat together the peak files from each chrom
    all_peaks_bed = args.outdir + '/All_peaks.bed'
    all_peaks_bed_fh = open(all_peaks_bed, "w")
    Ref_peak_beds.insert(0, 'cat')
    #print(Ref_peak_beds)
    subprocess.run(Ref_peak_beds, stdout=all_peaks_bed_fh)
    all_peaks_bed_fh.close()

    # Merge adjacent intervals with bedtools merge, according to argument pad
    cluster_bed = args.outdir + '/deNovoClusters.bed'
    cluster_bed_fh = open(cluster_bed, "w")
    bt_merge_args = ['bedtools', 'merge', '-d', str(args.pad), '-i', all_peaks_bed]
    subprocess.run(bt_merge_args, stdout=cluster_bed_fh)
    cluster_bed_fh.close()

    # clean up
    rm_args = ['rm', '-rf', all_peaks_bed]
    rm_args.extend(list(Ref_peak_beds[1:]))  ## removes the 'cat' that was tacked on above
    subprocess.run(rm_args)

    # return file name
    return cluster_bed

def quant(args, merged_bam, cluster_bed):
    """Calculate ShortStack quantification metrics from bed file

    Input:
    - args : Main argument object 
    - merged_bam : Filepath to the master .bam file
    - cluster_bed : Filepath to the bedfile of loci to quantify

    Output:
    - quant_data : Formatted like quant_data['Chr1:1234-2345'][~~item~~]. 
    Items are start, end, chrom,
    Length, Reads, UniqueReads, FracTop, Strand, MajorRNA, MajorRNAReads,
    Short, Long, ~~inclusive sizes~~, DicerCall
    - possible_mirs : Set of 21 and 22nt reads that originated from
    within a cluster and were >= mincov. A list, with entries like
    'AUGCC,Chr1,1234,+,23' (Seq,Chrom,Pos,Strand,NumberOfReads)

    Requires:
    - samtools (in PATH)
    - subprocess
    - csv
    """
    possible_mirs = []
    quant_data = {}
    lines_done = 0
    lines_chunk = 0

    # tqdm
    # get 'n' from the bed file name
    prefix = args.outdir + '/bedchunk'
    suffix = '.bed'
    no_prefix = re.sub(prefix, '', cluster_bed)
    n = int(re.sub(suffix, '', no_prefix))

    # for read group parsing
    rg_pattern = re.compile('RG:Z:(.*)')
    rgs = get_read_groups(merged_bam)

    # Get lines of the bedfile
    with open(cluster_bed) as fh:
        bedlines = fh.readlines()
    
    for bedline in tqdm(bedlines, position=n, leave=None):
    

    #time_file = args.outdir + '/times.csv'
    #time_fh = open(time_file, "w")
    #with open(cluster_bed) as fh:
        #for bedline in fh:
            #start_time = time.perf_counter_ns()
            #if lines_done >= 16300:
            #    print(f' Working on {bedline}')
            # construct samtools coordinates
        bedfields = bedline.rstrip().split("\t")
        sstart = int(bedfields[1]) + 1
        stcoords = bedfields[0] + ':' + str(sstart) + '-' + bedfields[2]

            # Counters
            # Sequences
            # strand_count
            # Lengths

        seqs = Counter()
        mir_count = Counter()
        strand_count = Counter()
        seq_lens = Counter()
        rg_counts = Counter()

        # samtools view command
        svcmd = ['samtools', 'view', '-F', '256', merged_bam, stcoords]
        sv = subprocess.Popen(svcmd, stdout=subprocess.PIPE,
            text=True, bufsize=1)
        reader = csv.reader(sv.stdout, delimiter="\t")
        for row in reader:
            flag = int(row[1])
            pos = row[3]  # keep as string
            if flag & 16:
                samseq = revcomp(row[9]).replace("T", "U")
                strand = '-'
            else:
                samseq = row[9].replace("T", "U")
                strand = '+'
            strand_count[strand] += 1
            seqs[samseq] += 1
            if len(samseq) < args.dicermin:
                seq_lens['Short'] += 1
            elif len(samseq) > args.dicermax:
                seq_lens['Long'] += 1
            else:
                seq_lens[str(len(samseq))] += 1
            
            # Try to find a read group
            if rgs:
                for field in reversed(row):
                    if rg_pattern.match(field):
                        this_rg = rg_pattern.match(field).group(1)
                        rg_counts[this_rg] += 1
                        break

            # tracking 21 and 22 nt reads for possible mir search
            if len(samseq) in range(21, 23, 1):
                    # process for bed format
                chrom = str(row[2])
                bstarti = int(pos) - 1
                bstart = str(bstarti)
                bend = str(bstarti + len(samseq))
                name = '.'
                score = str(0)
                k = (chrom, bstart, bend, name, score, strand)
                    #k = str(samseq) + ',' + row[2] + ',' + pos + ',' + strand
                mir_count[k] += 1
            #view_done_time = time.perf_counter_ns()
            #view_elapsed_time = view_done_time - start_time
            # filter the possible mirs based on mincov
        for m in mir_count:
            if mir_count[m] >= args.mincov:
                    # add the read count into the 'score' location
                newm = (m[0], m[1], m[2], m[3], str(mir_count[m]), m[5])
                newm_string = '\t'.join(newm)
                    #mx = m + ',' + str(mir_count[m])
                possible_mirs.append(newm_string)
            
            # Create a dictionary with unique Seqs => counts
            #seqdict = {}
            #for uniqseq in set(seqs):
            #    seqdict[uniqseq] = seqs.count(uniqseq)
            
            # Locus, Length, Reads, UniqueReads, FracTop, Strand, MajorRNA, MajorRNAReads, DicerCall, Short, Long, [..sizes..]
        quant_data[stcoords] = {}
        quant_data[stcoords]['start'] = sstart
        quant_data[stcoords]['end'] = int(bedfields[2])
        quant_data[stcoords]['chrom'] = str(bedfields[0])
        quant_data[stcoords]['Length'] = int(bedfields[2]) - int(bedfields[1])
        quant_data[stcoords]['Reads'] = strand_count['+'] + strand_count['-']

        if len(seqs) == 0:
                # no reads in the interval
            quant_data[stcoords]['UniqueReads'] = 0
            quant_data[stcoords]['FracTop'] = 'NA'
            quant_data[stcoords]['Strand'] = '.'
            quant_data[stcoords]['MajorRNA'] = 'NA'
            quant_data[stcoords]['MajorRNAReads'] = 0
            quant_data[stcoords]['Short'] = 0
            quant_data[stcoords]['Long'] = 0
            for rlen in range(args.dicermin, args.dicermax + 1, 1):
                quant_data[stcoords][str(rlen)] = 0
            quant_data[stcoords]['DicerCall'] = 'NA'
            if rgs:
                for rg in rgs:
                    quant_data[stcoords][rg] = 0
        else:
            quant_data[stcoords]['UniqueReads'] = len(seqs)
            quant_data[stcoords]['FracTop'] = strand_count['+'] / quant_data[stcoords]['Reads']
            if quant_data[stcoords]['FracTop'] >= args.strand_cutoff:
                quant_data[stcoords]['Strand'] = '+'
            elif quant_data[stcoords]['FracTop'] <= 1 - args.strand_cutoff:
                quant_data[stcoords]['Strand'] = '-'
            else:
                quant_data[stcoords]['Strand'] = '.'
                #value_key_pairs = ((value, key) for (key,value) in seqdict.items())
                #sorted_value_key_pairs = sorted(value_key_pairs, reverse=True)
            quant_data[stcoords]['MajorRNA'] = seqs.most_common(1)[0][0]
            quant_data[stcoords]['MajorRNAReads'] = seqs.most_common(1)[0][1]
                #len_counts = Counter()
                #for useq in seqdict:
                #    if len(useq) < args.dicermin:
                #        len_counts['Short'] += seqdict[useq]
                #    elif len(useq) > args.dicermin:
                #        len_counts['Long'] += seqdict[useq]
                #    else:
                #        len_counts[str(len(useq))] += seqdict[useq]
            quant_data[stcoords]['Short'] = seq_lens['Short']
            quant_data[stcoords]['Long'] = seq_lens['Long']
            predom_len = str(args.dicermin)
            predom_reads = 0
            for rlen in range(args.dicermin, args.dicermax + 1, 1):
                quant_data[stcoords][str(rlen)] = seq_lens[str(rlen)]
                if seq_lens[str(rlen)] > predom_reads:
                    predom_len = str(rlen)
                    predom_reads = seq_lens[str(rlen)]
            frac_ok = (quant_data[stcoords]['Reads'] - (quant_data[stcoords]['Short'] + 
                quant_data[stcoords]['Long'])) / quant_data[stcoords]['Reads']

                #print(frac_ok)
                
            if frac_ok >= 0.8:
                quant_data[stcoords]['DicerCall'] = predom_len
            else:
                quant_data[stcoords]['DicerCall'] = 'N'
            if rgs:
                for rg in rgs:
                    quant_data[stcoords][rg] = rg_counts[rg]
        lines_done += 1
        lines_chunk += 1
        if lines_chunk == 100:
            #print(f'File {cluster_bed} lines done {lines_done}')
            lines_chunk = 0
            #loc_done_time = time.perf_counter_ns()
            #loc_elapsed_time = loc_done_time - view_done_time
            #total_reads = len(seqs)
            #time_fh.write(str(total_reads) + ',' + str(view_elapsed_time) + ',' + str(loc_elapsed_time) + '\n')
            #print(f'{total_reads},{view_elapsed_time},{loc_elapsed_time}')
            #if str(bedfields[0]) == '2':
            #    time_fh.close()
            #    return None

    return (quant_data, possible_mirs)

def revcomp(seq):
    """Reverse complements a string

    It is very dumb. Only works on capital ATGC. All bet are off
    if input is otherwise
    https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
    """

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def quant_controller (args, merged_bam, all_bed):
    """Controls quant process using multithreading

    Input:
    - args : The main ShortStack argument object
    - merged_bam : The master bam file of sRNA-seq alignments
    - all_bed : Filepath to bed file containing all bed intervals

    Output:
    - qdata : Formatted like quant_data['Chr1:1234-2345'][~~item~~]. 
    Items are start, end, chrom,
    Length, Reads, UniqueReads, FracTop, Strand, MajorRNA, MajorRNAReads,
    Short, Long, ~~inclusive sizes~~, DicerCall
    - final_bed_path : Path to bedfile of possible de novo mature miRNA locations

    Requires:
    - Must be run as __main__
    - time
    - multiprocessing as mp
    - bedtools, rm, sort (all in PATH)
    - subprocess

    """
    # Split the bed file into chunks that parallelize analysis
    if __name__ == '__main__':
        # Report to user
        print('')
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print(f'Analyzing cluster properties using {args.threads} threads')

        bed_lines = 0
        with open(all_bed) as fh:
            for bl in fh:
                bed_lines += 1
        chunk_size = int(bed_lines / args.threads)
        chunk_files = []
        bed_base = args.outdir + '/bedchunk'
        chunk_n = 1
        chunk_file = bed_base + str(chunk_n) + '.bed'
        cfh = open(chunk_file, "w")
        chunk_files.append(chunk_file)
        blines = 0
        with open(all_bed) as bfh:
            for bl in bfh:
                cfh.write(bl)
                blines += 1
                if blines > chunk_size:
                    cfh.close()
                    chunk_n += 1
                    chunk_file = bed_base + str(chunk_n) + '.bed'
                    cfh = open(chunk_file, "w")
                    chunk_files.append(chunk_file)
                    blines = 0
        cfh.close()

        #tqdm setup for mp
        tqdm.set_lock(mp.RLock())

        q_iter = list(zip([args] * len(chunk_files), 
            [merged_bam] * len(chunk_files), chunk_files))
        #print(q_iter)
        #sys.exit()
        with mp.Pool(processes=args.threads, 
            initializer=tqdm.set_lock, initargs=(tqdm.get_lock(),)) as pool:
            q_results = pool.starmap(quant, q_iter)
        q_data_t, p_mirs_l = zip(*q_results)

        #p_mirs = p_mirs_l[0] # ???

        # Boundary condition: no possible mirs found
        if len(p_mirs_l) < 1:

            # clean up intermediate files
            rm_args = ['rm', '-f']
            for chunk_file in chunk_files:
                rm_args.append(chunk_file)
            subprocess.run(args=rm_args)
            final_bedpath = None
            return(qdata, final_bedpath)

        # Write the possible mirs to a temporary bed file
        pmir_path = args.outdir + '/possible_mirs.bed'
        with open(pmir_path, "w") as pmir_bedfh:
            for pmir_entry in p_mirs_l:
                for line in pmir_entry:
                    pmir_bedfh.write(line)
                    pmir_bedfh.write('\n')
        
        # process the qdata tuple into a single dictionary
        q_data = {}
        for q_data_entry in q_data_t:
            q_data = merge_two_dicts(q_data, q_data_entry)
        
        # bedtools cluster -s -i possible_mirs.bed | sort -k7,7n -k5,5nr
        # cluster and sort the possible mirs.
        # The idea is to assume that overlapping ones (on same genomic strand)
        # are part of same query. Keep only the most abundant one.
        bt_args = ['bedtools', 'cluster', '-s', '-i', pmir_path]
        clus_bed_path = args.outdir + '/possible_mirs_clustered.bed'
        clus_bed_fh = open(clus_bed_path, "w")
        subprocess.run(args=bt_args, stdout=clus_bed_fh, text=True)
        clus_bed_fh.close()

        sort_args = ['sort', '-k7,7n', '-k5,5nr', clus_bed_path]
        so_clus_bed_path = args.outdir + '/possible_mirs_clustured_sorted.bed'
        so_clus_bed_fh = open(so_clus_bed_path, "w")
        subprocess.run(args=sort_args, stdout=so_clus_bed_fh, text=True)
        so_clus_bed_fh.close()

        clus_n = 'foo'
        final_bed_path = args.outdir + '/possible_mirs_clean.bed'
        final_bed_fh = open(final_bed_path, "w")
        with open(so_clus_bed_path) as so_clus_bed_fh:
            for bedline in so_clus_bed_fh:
                f = bedline.rstrip().split("\t")
                if f[6] != clus_n:
                    # keep it, dropping the cluster number
                    f2 = f[0:6]
                    outline = '\t'.join(f2) + '\n'
                    final_bed_fh.write(outline)
                clus_n = f[6]

        final_bed_fh.close()

        # clean up intermediate files
        rm_args = ['rm', '-f', pmir_path, clus_bed_path, so_clus_bed_path]
        for chunk_file in chunk_files:
            rm_args.append(chunk_file)
        subprocess.run(args=rm_args)

        print()
        print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
        print(f' Completed')

        return(q_data, final_bed_path)
    else:
        sys.exit('Fatal: function quant_controller needs to be __main__ !')

def mirna(args, merged_bam, fai, pmir_bedfile):
    """Entry point into microRNA locus annotation

    Input:
    - args :  ShortStack argument object
    - merged_bam : Master file of sRNA-seq alignments
    - fai : Dictionary of ChromosomeNames => Lengths
    - pmir_bedfile : path to bed formatted file of possible mirs from de novo

    Output: ??

    Requires:
    - time
    - subprocess
    - bowtie (in PATH)
    - samtools (in PATH)
    - bedtools (in PATH)
    """
    if args.nohp is True:
        return None
    
    # if no knownRNAs AND no de novo places to look, return nothing
    if (args.knownRNAs is None) and (pmir_bedfile is None):
        return None
    
    # similarly, if no knownRNAs and de novo is turned off, return nothing
    if (args.knownRNAs is None) and (args.dn_mirna is False):
        return None
    
    # report to user
    print('')
    print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
    print('Searching for valid microRNA loci')

    if args.knownRNAs is not None:
        # Align user's sequences to the genome
        print('Aligning knownRNAs sequences to genome')
        # we must sanitize user input .. capitalize, U to T
        sane_knownRNAs = sanitize_knownRNAs(args)

        # align using bowtie
        # "score" set to zero, reflecting unknown alignments at this point
        mir_samfile = args.outdir + '/user_mir.sam'
        un_file = args.outdir + '/knownRNAs_unaligned.fasta'
        bt_args = ['bowtie', '-f', '-a', '-S', '-v', '0', '-p',
            str(args.threads), '--mapq', '0', '--un', un_file, '-x',  
            args.genomefile, sane_knownRNAs, mir_samfile]
        subprocess.run(args = bt_args)

        # convert to bam, sort, and convert to bed
        mir_unsortedbam = args.outdir + '/user_mir_unsorted.bam'
        sv_args = ['samtools', 'view', '--threads', str(args.threads),
            '-b', mir_samfile, '-o', mir_unsortedbam]
        subprocess.run(args = sv_args)
        mir_bam = args.outdir + '/user_mir.bam'
        sort_args = ['samtools', 'sort', '--threads', str(args.threads),
            '-o', mir_bam, mir_unsortedbam]
        subprocess.run(args = sort_args)
        
        bam2bed_args = ['bedtools', 'bamtobed', '-i', mir_bam]
        user_mir_bed_f = args.outdir + '/user_mir.bed'
        user_mir_bed_fh = open(user_mir_bed_f, "w")
        subprocess.run(args=bam2bed_args, stdout=user_mir_bed_fh)
        user_mir_bed_fh.close()

        # cleanup
        rm_args = ['rm', '-f', sane_knownRNAs, mir_samfile, 
            mir_unsortedbam, mir_bam]
        subprocess.run(args = rm_args)
    else:
        user_mir_bed_f = None

    # Let's not combine the user's file and the de novo file.
    # Instead, process them sequentially.

    # A file for the combined, passed mir loci in bed format
    passed_mir_bed = args.outdir + '/passed_mir_redundant.bed'
    passed_mir_bed_fh = open(passed_mir_bed, "w")

    # counter to tell if nothing is passed
    n_passed = 0
    n_returned = 0

    if user_mir_bed_f is not None:
        print('')
        print('Screening of possible microRNAs from user provided knownRNAs')
        if __name__ == '__main__':
            with open(user_mir_bed_f) as f:
                bedlines = f.read().splitlines()
            mir_iter = list(zip(bedlines, [args] * len(bedlines),
            [merged_bam] * len(bedlines), [fai] * len(bedlines)))
            with mp.Pool(args.threads) as pool:
                user_mloci1 = pool.starmap(mir_analysis,
                tqdm(mir_iter, total=len(bedlines),
                desc = 'Candidates examined', leave=None),
                chunksize=1)
            # user_mloci1 is returned from above
            # Unpack it
            user_q_bedlines, user_locus_bedlines = zip(*user_mloci1)

            # Write the revised (with scores added) user_q_bedlines
            u_q_file = args.outdir + '/knownRNAs.bed'
            u_q_fh = open(u_q_file, "w")
            for line in user_q_bedlines:
                u_q_fh.write(line)
                u_q_fh.write('\n')
            u_q_fh.close()

            # Add the passed mirna loci to the final bed file
            for line in user_locus_bedlines:
                n_returned += 1
                if line is not None:
                    n_passed += 1
                    passed_mir_bed_fh.write(line)
                    passed_mir_bed_fh.write('\n')
            
    # next, the de novo analyses, if they exist, and were requested
    if args.dn_mirna is True:
        pmir_filt_bedfile = args.outdir + '/pmir_filt.bed'
        if (user_mir_bed_f is not None) and (pmir_bedfile is not None):
            # remove any denovos that are exactly in the user file keep the rest
            user_locations = []
            with open(user_mir_bed_f) as user_mir_bed_fh:
                for user_bed_line in user_mir_bed_fh:
                    bf = user_bed_line.rstrip().split("\t")
                    loc = bf[0] + ':' + bf[1] + ':' + bf[2] + ':' + bf[5]
                    user_locations.append(loc)

            # Then the de novos that are NEW are written
            with open(pmir_filt_bedfile, "w") as pmir_filt_bedfh:
                with open(pmir_bedfile) as pmir_bedfh:
                    for pmir_bedline in pmir_bedfh:
                        bf = pmir_bedline.rstrip().split("\t")
                        loc = bf[0] + ':' + bf[1] + ':' + bf[2] + ':' + bf[5]
                        if loc not in user_locations:
                            pmir_filt_bedfh.write(pmir_bedline)
            pmir_bedfh.close()
            # the original de novo file is no longer needed
            # test
            #print(f'-1 attempting to remove {pmir_bedfile}')
            #rm_args = ['rm', '-f', pmir_bedfile]
            #subprocess.run(args=rm_args)
        elif pmir_bedfile is not None:
            # it becomes the "filtered" file
            # the original de novo file is no longer needed so mv instead of cp
            # test
            #print(f'0 attempting to mv {pmir_bedfile} to {pmir_filt_bedfile}')
            mv_args = ['mv', pmir_bedfile, pmir_filt_bedfile]
            subprocess.run(args=mv_args)
        
        # Now, the de novo can be run
        if pmir_bedfile is not None:
            print('')
            print('Screening of possible de novo microRNAs')
            if __name__ == '__main__':
                with open(pmir_filt_bedfile) as f:
                    bedlines = f.read().splitlines()
                mir_iter = list(zip(bedlines, [args] * len(bedlines),
                [merged_bam] * len(bedlines), [fai] * len(bedlines)))
                with mp.Pool(args.threads) as pool:
                    denovo_mloci1 = pool.starmap(mir_analysis,
                    tqdm(mir_iter, total=len(bedlines),
                    desc = 'Candidates examined', leave=None),
                    chunksize=1)
                # denovo_mloci1 is returned from above
                # Unpack it
                dn_q_bedlines, dn_locus_bedlines = zip(*denovo_mloci1)

                # Don't care about dn_q_bedlines
                # Write the dn_locus_bedlines to final file
                n_returned = 0
                for line in dn_locus_bedlines:
                    n_returned += 1
                    if line is not None:
                        n_passed += 1
                        passed_mir_bed_fh.write(line)
                        passed_mir_bed_fh.write('\n')

                # Testing
                #debug = args.outdir + '/debug.bed'
                #debugf = open(debug, "w")
                #for line in dn_q_bedlines:
                #    debugf.write(line)
                #    debugf.write('\n')
                #debugf.close()

                # Delete the pmir_filt_bedfile and the pmir_bedfile
                # test
                # to delete {pmir_bedfile}')
                subprocess.run(f'rm -f {pmir_filt_bedfile} {pmir_bedfile}', shell=True)

    passed_mir_bed_fh.close()
    
    # Delete the user_mir_bed_f
    if user_mir_bed_f is not None:
        subprocess.run(f'rm -f {user_mir_bed_f}', shell=True)
    
    # Stop now if none were passed
    if n_passed == 0:
        print('')
        print('No microRNA loci were found!')
        return None

    # Remove overlapping passed microRNA loci
    nr_passed_mir_bed = collapse_mir_loci(args, passed_mir_bed)

    # Call strucVis
    # To do: strucVis on conda!
    call_strucVis(nr_passed_mir_bed, args, merged_bam)

    # call quant on the mirna loci
    mir_qdata, junk = quant_controller(args, merged_bam, nr_passed_mir_bed)

    # Report n found to user
    n_found = 0
    with open(nr_passed_mir_bed) as bed:
        for line in bed:
            n_found += 1
    print('')
    print(f'Found a total of {n_found} MIRNA loci')
    print('')
    print('<<< WARNING >>>')
    print('Do not rely on these results alone to annotate new MIRNA loci!')
    print('The false positive rate for de novo MIRNA identidication is low, but NOT ZERO')
    print('Insepct each locus, especially the strucVis output, and see')
    print('https://doi.org/10.1105/tpc.17.00851 , https://doi.org/10.1093/nar/gky1141')

    # Look for file of denovo, if found, filter it
    # .. this removes any de-novo loci that overlap with 
    # .. any of the passed microRNAs.
    dnfile = args.outdir + '/deNovoClusters.bed'
    if os.path.exists(dnfile) is True:
        filt_dnfile = args.outdir + '/filtDeNovoClusters.bed'
        # Use bedtools subtract
        cmd = f'bedtools subtract -A -a {dnfile} -b {nr_passed_mir_bed}'
        cmd = cmd + f' > {filt_dnfile}'
        subprocess.run(cmd, shell=True)
        # remove the first file
        subprocess.run(f'rm {dnfile}', shell=True)

    # return the mir_qdata
    return mir_qdata

def mir_analysis(bedline, args, merged_bam, fai):
    """Finding valid microRNA loci

    Inputs:
    - combined_mir_bedfile : Filepath to bedfile of possible microRNA locations
    - args : ShortStack's arguments object
    - merged_bam : Master sRNA-seq alignment file
    - fai : Dictionary of Chrom => Lengths

    Outputs: ?

    Requires: 
    """
    # Along the way, we will re-write the user_mir.bed file to put in
    # read counts in the score column


    bedfields = bedline.split("\t")

    # Is the exact query actuallu aligned at this location?
    # We only need to ask this if it is as user query,
    #  in which case the score column from bed format will be 0.
    if int(bedfields[4]) == 0:
        q_exact_count = count_exact_location(bedfields, merged_bam)
    else: 
        q_exact_count = int(bedfields[4])

    # bedfield score field is replaced with the exact query count
    bedfields[4] = str(q_exact_count)
    new_bedline = '\t'.join(bedfields)
    
    if q_exact_count == 0:
        #print('TEST failed check_is_usermir_there')
        #new_bedline = new_bedline + '\tnoReads'
        return (new_bedline, None)

    # Get the locus to fold, if any
    # Screens rapidly for strand cutoff and two peaks
    mir_locus = get_mir_locus(bedfields, fai, args, merged_bam)
    if mir_locus == None:
        #print('TEST failed get_mir_locus')
        #new_bedline = new_bedline + '\tfailed_mir_locus'
        return (new_bedline, None)
    
    # Get the dotbracket
    cmd = 'samtools faidx'
    if bedfields[5] == '-':
        cmd = cmd + ' -i'
    cmd = cmd + f' {args.genomefile} {mir_locus}'
    cmd = cmd + f' | RNAfold --noPS'
    foldjob = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    foldlines = foldjob.stdout.split('\n')
    dotbracket = foldlines[2].rstrip().split(' ')[0]

    # Analyze the fold
    new_locus = analyze_fold(mir_locus, dotbracket, bedfields, merged_bam, args)

    # Will be None or not
    if new_locus is not None:
        # construct and return a bedline of the mir_locus
        # It needs to somehow store the miRNA location too?
        # Can use the thickStart and thickEnd fields (6 and 7 in zero-based)
        chr, start, stop = parse_locus(new_locus)
        mirbed = []
        mirbed.append(chr) # Chrom
        mirbed.append(str(int(start) - 1)) # Start in bed indexing system
        mirbed.append(str(stop)) # End
        mirbed.append(bedfields[3]) # Name, which has little/no meaning right now
        mirbed.append(str(bedfields[4])) # Score, which here is the n of reads for the mature only
        mirbed.append(bedfields[5]) # Strand
        mirbed.append(str(bedfields[1])) # Start of the mature miRNA, in bed indexing system
        mirbed.append(str(bedfields[2])) # Stop of the mature miRNA, in bed indexing system
        mirbedline = '\t'.join(mirbed)
        #new_bedline = new_bedline + '\tpassed'
        return (new_bedline, mirbedline)
    else:
        #new_bedline = new_bedline + '\tfailed_analyze_fold'
        return (new_bedline, None)

def precision_test(peaks, merged_bam, args, fai):
    # peaks{'querySAM': [], 'upSAM': [], 'downSAM': []}
    # both up and down could be empty
    if (peaks['upSAM'] is None) and (peaks['downSAM'] is None):
        return None, None

    # determine range of SAM alignments
    hiSAM = 1
    loSAM = 1000000000000000 # assume no chromosome is > 1 quadrillion bps!
    precise_aligns = len(peaks['querySAM'])
    for qsam in peaks['querySAM']:
        start = int(qsam[3])
        end = int(qsam[3]) + len(qsam[9]) - 1
        if start < loSAM:
            loSAM = start
        if end > hiSAM:
            hiSAM = end
    if peaks['upSAM'] is not None:
        precise_aligns += len(peaks['upSAM'])
        for qsam in peaks['upSAM']:
            start = int(qsam[3])
            end = int(qsam[3]) + len(qsam[9]) - 1
            if start < loSAM:
                loSAM = start
            if end > hiSAM:
                hiSAM = end
    if peaks['downSAM'] is not None:
        precise_aligns += len(peaks['downSAM'])
        for qsam in peaks['downSAM']:
            start = int(qsam[3])
            end = int(qsam[3]) + len(qsam[9]) - 1
            if start < loSAM:
                loSAM = start
            if end > hiSAM:
                hiSAM = end
    chrom = peaks['querySAM'][0][2] # first alignment is fine, all chrom same
    flag = peaks['querySAM'][0][1] # likewise all aligns will have same strand
    if int(flag) & 16:
        f = '-f' # samtools view require
    else:
        f = '-F' # samtools view exclude
    
    # Determine the locus boundaries... just add 24nts on each side
    #  subject to chromosome edges
    loc_start = max(loSAM - 24, 1)
    loc_end = min(hiSAM + 24, fai[chrom])

    # Count alignments, on the same strand, in the interval
    st_loc = chrom + ':' + str(loc_start) + '-' + str(loc_end)
    st_args = ['samtools', 'view', '-c', f, str(16), merged_bam, st_loc]
    st = subprocess.run(args=st_args, capture_output=True, text=True)
    total_aligns = int(st.stdout.rstrip())

    if total_aligns > 0:
        precision = precise_aligns / total_aligns
    else:
        precision = 0
    
    #print(f'total_aligns: {total_aligns}')
    #print(f'precise_aligns: {precise_aligns}')
    #print(f'precision: {precision}')

    return precision, st_loc


def get_mir_locus(bedfields, fai, args, merged_bam):
    """Given bedfields, compute upstream and downstream intervals for miRNA folding
    Scans for strand and peaks, based on area-under-curve counts from samtools depth

    Returns either None, or a mir_locus in format Chr:Start-Stop (1-based inclusive)
    """
    # testing
    #print(f'Working on bedfields {bedfields}')
    query_start = int(bedfields[1]) + 1 # in one-based
    query_stop = int(bedfields[2])
    query_strand = bedfields[5]
    chrom = bedfields[0]

    # -200 to +200, relative to query_peak_start
    # samtools view gracefully ignores queries that go over the chrom end.
    # But it won't gracefully ignore negative numbers
    region_start = max(query_start - 200, 1)
    region_stop = min(query_stop + 200, fai[chrom])
    reg = str(chrom) + ':' + str(region_start) + '-' + str(region_stop)

    # upstream window definition 
    up_reg_start = region_start
    up_reg_stop = min(query_stop + 20, fai[chrom])
    up_reg = str(chrom) + ':' + str(up_reg_start) + '-' + str(up_reg_stop)

    # downstream window definition
    down_reg_start = max(query_start - 20, 1)
    down_reg_stop = region_stop
    down_reg = str(chrom) + ':' + str(down_reg_start) + '-' + str(down_reg_stop)

    # Only two calls against the bam file
    # Retrieve depths across whole region irrespective of strand
    depths_all = subprocess.run(f'samtools depth -a -r {reg} {merged_bam}', shell=True, 
    text=True, capture_output=True)

    # Retrieve depths across whole region for only plus strand
    depths_plus = subprocess.run(f'samtools depth -a -G 16 -r {reg} {merged_bam}', shell=True, 
    text=True, capture_output=True)

    # Parse into dictionaries
    da_lines = depths_all.stdout.split('\n')
    dp_lines = depths_plus.stdout.split('\n')
    da = {} # position => depth
    dp = {} # position => depth
    for da_line in da_lines:
        if da_line:
            fields = da_line.split('\t')
            da[int(fields[1])] = int(fields[2])
    for dp_line in dp_lines:
        if dp_line:
            fields = dp_line.split('\t')
            dp[int(fields[1])] = int(fields[2])
    #print('da')
    #print(da)
    #print('')
    #print('dp')
    #print(dp)
    
    # Infer the minus genomic strand depths
    dm = {}
    for pos in range(region_start, region_stop + 1):
        m = da[pos] - dp[pos]
        dm[pos] = m
    #print('dm')
    #print(dm)
    
    # Up or down is determined by largest auc on the query's genomic strand
    up_auc = 0
    for pos in range(up_reg_start, up_reg_stop + 1):
        if query_strand == '+':
            up_auc += dp[pos]
        else:
            up_auc += dm[pos]
    
    down_auc = 0
    for pos in range(down_reg_start, down_reg_stop + 1):
        if query_strand == '+':
            down_auc += dp[pos]
        else:
            down_auc += dm[pos]
    
    if up_auc >= down_auc:
        # Take the up
        side_taken = 'up'
        mir_locus = up_reg
        mir_start = up_reg_start
        mir_stop = up_reg_stop
        stranded_auc = up_auc
    else:
        # Take the down
        side_taken = 'down'
        mir_locus = down_reg
        mir_start = down_reg_start
        mir_stop = down_reg_stop
        stranded_auc = down_auc

    # Screen based on fraction on correct strand.
    # Need to compute auc in the interval for both strands
    unstranded_auc = 0
    for pos in range(mir_start, mir_stop + 1):
        unstranded_auc += da[pos]
    
    frac_correct = stranded_auc / unstranded_auc
    if frac_correct < float(args.strand_cutoff):
        #print('Failed frac correct strand')
        return None
    
    # Screen based on two peaks at a "precision" of 0.65
    # Note this is a fuzzy precision estiamate, not the final one
    # The idea is to make sure there is at least one plausible set of peaks
    # that may end up being a miR/miR*
    ## Determine the auc in the query window
    auc_q = 0
    for pos in range(query_start, query_stop + 1):
        if query_strand == '+':
            auc_q += dp[pos]
        else:
            auc_q += dm[pos]
    
    if side_taken == 'up':
        # query is on the right, we scan on it's left
        for i in range(mir_start, query_start - (query_stop - query_start)):
            this_auc = 0
            for j in range(i, i + (query_stop - query_start)):
                if query_strand == '+':
                    this_auc += dp[j]
                else:
                    this_auc += dm[j]
            if this_auc == 0:
                continue
            precision = (auc_q + this_auc) / stranded_auc
            if precision >= 0.65:
                return mir_locus
    else:
        # query is on the left, we scan on its right
        for i in range(query_stop + 10, mir_stop - (query_stop - query_start)):
            this_auc = 0
            for j in range(i, i + (query_stop - query_start)):
                if query_strand == '+':
                    this_auc += dp[j]
                else:
                    this_auc += dm[j]
            if this_auc == 0:
                continue
            precision = (auc_q + this_auc) / stranded_auc
            if precision >= 0.65:
                return mir_locus
    # If you get here, failed the two peaks check
    #print('FAILED 2 peaks')
    return None

def count_exact_location(bedfields, merged_bam):
    """Checks a bed position to count exact reads

    Inputs:
    - bedfields : list of data from a bed file
    - merged_bam : The master sRNA-seq alignment file

    Output:
    - readcount : integer

    Requires:
    - samtools (in PATH)
    - subprocess
    """

    mir_name = bedfields[3]
    mir_chrom = bedfields[0]
    mir_pos = int(bedfields[1]) + 1
    mir_end = int(bedfields[2])
    mir_strand = bedfields[5]

    st_expr = f'\"pos == {mir_pos} && endpos == {mir_end}\"'
    if mir_strand == '+':
        flag_filt = '-F 16'
    else:
        flag_filt = '-f 16'
    
    mb_coords = (mir_chrom + ':' + str(mir_pos) + '-' +
        str(mir_pos))
    mb = subprocess.run(f'samtools view -c {flag_filt} -e {st_expr} {merged_bam} {mb_coords}',
    shell=True, capture_output=True, text=True)

    read_count = int(mb.stdout.rstrip())
    return read_count


def sanitize_knownRNAs(args):
    # sane means upper case, with Us converted to Ts
    # This is required for bowtie mapping
    sane_filepath = args.outdir + '/user_knownRNAs.fasta'
    sane_fh = open(sane_filepath, "w")
    for seq_record in SeqIO.parse(args.knownRNAs, "fasta"):
        out = ('>' + str(seq_record.id) + '\n' + 
            str(seq_record.seq.upper().back_transcribe()) + '\n')
        sane_fh.write(out)
    sane_fh.close()
    return(sane_filepath)

def fold_mir_locus(mir_locus, args, bedfields, precision):
    # RNAfold a given locus, from the strand specified in bedfields
    # return the dot-bracket line
    
    # Determine strand from bedfields
    strand = bedfields[5] # will be +, -, or possibly .

    # Retrieve genomic sequence as a simple string
    #chr, start, stop = parse_locus(mir_locus)
    #genome = Fasta(args.genomefile)
    #zstart = start - 1 # pyfadix ranges are bed-like
    #if strand == '-':
    #    foldSeq = str(genome[chr][zstart:stop].reverse.complement)
    #else:
    #    foldSeq = str(genome[chr][zstart:stop])

    # RNAfold, capturing the dot-bracket line
    if strand == '-':
        rf = subprocess.run(f'samtools faidx -i {args.genomefile} {mir_locus} | RNAfold', shell=True,
        capture_output=True, text=True)
    else:
        rf = subprocess.run(f'samtools faidx {args.genomefile} {mir_locus} | RNAfold', shell=True,
        capture_output=True, text=True)
    db = rf.stdout.split("\n")[1]
    # remove the trailing deltaG part
    db = db.split(' ')[0]
    return db

def parse_locus(locus):
    # parse a locus in the samtools / ShortStack format
    # like chr:start-stop
    loc_pattern = re.compile('(^\S+):(\d+)-(\d+)$')
    chr = loc_pattern.match(locus).group(1)
    start = int(loc_pattern.match(locus).group(2))
    stop = int(loc_pattern.match(locus).group(3))
    return chr, start, stop

def get_q_rels(mir_locus, bedfields, dotbracket):
    # Given a locus like Chr1:1230-2233 and
    #  a list from a bed line,
    #  convert genomic coordinates of the interval 
    #  from the bed entry to 1-based relative coordinates
    #  of the sub-sequence

    # First check for validity of query
    if mir_locus is None:
        return None, None
    if dotbracket is None:
        return None, None

    #print(f'Sending mir_locus {mir_locus}')
    foldchr, foldstart, foldstop = parse_locus(mir_locus)
    qstart = int(bedfields[1]) + 1 # plus one because of bed weirdness
    qstop = int(bedfields[2])
    strand = bedfields[5]

    if strand == '-':
        # its harder
        q_rel_start = foldstop - qstop + 1
        q_rel_stop = foldstop - qstart + 1
    else:
        q_rel_start = qstart - foldstart + 1
        q_rel_stop = qstop - foldstart + 1
    
    # testing
    #print(f'foldchr: {foldchr} foldstart: {foldstart} foldstop: {foldstop}')
    #print(f'qstart: {qstart} qstop: {qstop} strand: {strand}')
    #print(f'q_rel_start: {q_rel_start} q_rel_stop: {q_rel_stop}')
    #print("\n")
    return q_rel_start, q_rel_stop

def test_structure(dotbracket, q_rel_start, q_rel_stop):
    # Takes a dotbracket structure, a start and and stop,
    # and evaulates the given interval for 
    # a) all pairs in same direction
    # b) <= 5 unpaired nts
    # Pass returns True
    # Fail returns False
    
    # adjust for 0-based half-open index slicing numbers
    a = q_rel_start - 1

    # get structure
    struc = dotbracket[a:q_rel_stop]

    n_unp = struc.count('.')
    n_left = struc.count('(')
    n_right = struc.count(')')

    if (n_unp <= 5) and ((n_left == 0) or (n_right == 0)):
        return True
    else:
        return False

def get_s_rels(q_rel_start, q_rel_stop, dotbracket):
    # Compute position of miR* given miR location on a dotbracket

    # Get a lookup table of paired positions
    pairs = get_pairs(dotbracket)

    # determine the 3' end (stop) of miR* from q_rel_start
    s_rel_stop = None
    for attempts, qpos in enumerate(range(q_rel_start, q_rel_stop)):
        # testing
        #print ('finding s_rel_stop')
        #print(f'attempts: {attempts}')
        #print(f'qpos: {qpos}')
        #print('')
        if pairs[qpos] is not None:
            s_rel_stop = pairs[qpos] + attempts + 2
            break
    
    # determine the 5' end (start) or miR* from q_rel_stop - 2
    s_rel_start = None
    for attempts, qpos in enumerate(range((q_rel_stop - 2), q_rel_start, -1)):
        # testing
        #print ('finding s_rel_start')
        #print(f'attempts: {attempts}')
        #print(f'qpos: {qpos}')
        #print('')
        if pairs[qpos] is not None:
            s_rel_start = pairs[qpos] - attempts
            break
    
    # testing
    #print(f'dotbracket: {dotbracket}')
    #print(f's_rel_start: {s_rel_start}')
    #print(f's_rel_stop: {s_rel_stop}')
    #print('')

    # Sanity check. Length of computed miR* needs to be close to query
    qlen = q_rel_stop - q_rel_start + 1
    slen = s_rel_stop - s_rel_start + 1
    diff = slen - qlen
    if diff not in range(-2, 2):
        return None, None
    
    return s_rel_start, s_rel_stop
    
def get_pairs(dotbracket):
    # Compute positional lookups for a dotbracket structure
    # Note conversion from enumerate's zero-based system to one-based
    pairs = {}
    for i, j in enumerate(dotbracket):
        pairs[(i + 1)] = None
    
    stack = []
    for i, j in enumerate(dotbracket):
        if j == '(':
            stack.append((i + 1))
        elif j == ')':
            bp1 = stack.pop()
            bp2 = i + 1
            pairs[bp1] = bp2
            pairs[bp2] = bp1
    # testing
    #print(dotbracket)
    #print(pairs)
    #sys.exit()
    return pairs

#s_gen_start, s_gen_stop = get_star_genomic(mir_locus, bedfields, s_rel_start, s_rel_stop)
def get_star_genomic(mir_locus, bedfields, s_rel_start, s_rel_stop):
    # given relative coordinates, return genomic coordinates for a miR* 
    strand = bedfields[5]
    chr, start, stop = parse_locus(mir_locus)
    if strand == '+' :
        s_g_start = start + s_rel_start - 1
        s_g_stop = start + s_rel_stop - 1
    else:
        # - strand stuff
        s_g_stop = stop - s_rel_start + 1
        s_g_start = stop - s_rel_stop + 1
    return s_g_start, s_g_stop

def compute_final_precision(merged_bam, s_gen_start, s_gen_stop, mir_locus, bedfields):
    q_win_start = int(bedfields[1])  # zero to one base conversion, add 1nt slop
    q_win_end = int(bedfields[2]) + 1 # one nt slop on the other side
    q_exp = f'\"pos >= {str(q_win_start)} && endpos <= {str(q_win_end)}\"'

    s_win_start = s_gen_start - 1
    s_win_stop = s_gen_stop + 1
    s_exp = f'\"pos >= {str(s_win_start)} && endpos <= {str(s_win_stop)}\"'

    if bedfields[5] == '+':
        flag_filt = '-F 16'
    else:
        flag_filt = '-f 16'

    # count reads in query window
    q_count = subprocess.run(f'samtools view -c -e {q_exp} {flag_filt} {merged_bam} {mir_locus}', 
    shell=True, capture_output=True, text=True)

    # count reads in star window
    s_count = subprocess.run(f'samtools view -c -e {s_exp} {flag_filt} {merged_bam} {mir_locus}', 
    shell=True, capture_output=True, text=True)

    # count total reads
    t_count = subprocess.run(f'samtools view -c {merged_bam} {mir_locus}',
    shell=True, capture_output=True, text=True)

    # Do the math
    t_n = int(t_count.stdout.rstrip())
    if t_n == 0:
        return 0
    
    # has to be at least one star read, regardless
    s_n = int(s_count.stdout.rstrip())
    if s_n == 0:
        return 0

    q_n = int(q_count.stdout.rstrip())
    precision = (q_n + s_n) / t_n
    return precision

def analyze_fold(mir_locus, dotbracket, bedfields, merged_bam, args):
    """Assess a folded microRNA hairpin candidate  
    """
    # Get relative position of exact query on the dot brackets
    q_rel_start, q_rel_stop = get_q_rels(mir_locus, bedfields, dotbracket)
    if q_rel_start is None:
        return None
    
    # Test the structure of the exact query
    # True if OK, False if failed
    q_struc = test_structure(dotbracket, q_rel_start, q_rel_stop)
    if q_struc is False:
        return None
    
    # Compute relative position of the predicted microRNA*
    s_rel_start, s_rel_stop = get_s_rels(q_rel_start, q_rel_stop, dotbracket)
    if s_rel_start is None:
        return None

    # Compute position of miR* in genomic coordinates given relative coordinates
    s_gen_start, s_gen_stop = get_star_genomic(mir_locus, bedfields, s_rel_start, s_rel_stop)
    if s_gen_start is None:
        return None
    
    # Now that the miR* position is known, the locus needs to be trimmed
    #  in size.
    new_start = min(int(bedfields[1]) - 21, s_gen_start - 20)
    new_end = max(int(bedfields[2]) + 20, s_gen_stop + 20)
    new_locus = bedfields[0] + ':' + str(new_start) + '-' + str(new_end)
    
    # It is possible that the new locus folds differently. So, re-analyze!
    # Get the dotbracket
    cmd = 'samtools faidx'
    if bedfields[5] == '-':
        cmd = cmd + ' -i'
    cmd = cmd + f' {args.genomefile} {new_locus}'
    cmd = cmd + f' | RNAfold --noPS'
    foldjob = subprocess.run(cmd, shell=True, text=True, capture_output=True)
    foldlines = foldjob.stdout.split('\n')
    new_dotbracket = foldlines[2].rstrip().split(' ')[0]

    # Get relative position of exact query on the dot brackets
    q_rel_start, q_rel_stop = get_q_rels(new_locus, bedfields, new_dotbracket)
    if q_rel_start is None:
        return None
    
    # Test the structure of the exact query
    # True if OK, False if failed
    q_struc = test_structure(new_dotbracket, q_rel_start, q_rel_stop)
    if q_struc is False:
        return None
    
    # Compute relative position of the predicted microRNA*
    s_rel_start, s_rel_stop = get_s_rels(q_rel_start, q_rel_stop, new_dotbracket)
    if s_rel_start is None:
        return None

    # Compute position of miR* in genomic coordinates given relative coordinates
    s_gen_start, s_gen_stop = get_star_genomic(new_locus, bedfields, s_rel_start, s_rel_stop)
    if s_gen_start is None:
        return None

    # Compute the actual precision (based on miR* and the final locus coords)
    final_precision = compute_final_precision(merged_bam, s_gen_start, 
    s_gen_stop, new_locus, bedfields)
    if final_precision < 0.75:
        return None
    else:
        # Testing        
        return new_locus

def collapse_mir_loci(args, passed_mir_bed):

    # Testing
    #print(f'In collapse_mir_loci received passed_mir_bed as {passed_mir_bed}')
    #sys.exit()

    # Sort the bed file and call bedtools cluster
    cmd = f'sort -k1,1 -k2,2n {passed_mir_bed} | '
    cmd = cmd + f'bedtools cluster'
    clust = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    # Open an output file for the non-redundant bed file
    outfile = args.outdir + '/passed_mir_final.bed' 
    out = open(outfile, "w")
    
    # Scan results
    lines = clust.stdout.split('\n')
    current_cluster = 0
    best_score = 0
    for line in lines:
        if line:
            fields = line.split('\t')
            this_clus_n = int(fields.pop())
            if this_clus_n != current_cluster:
                if current_cluster > 0:
                    # Write results
                    out.write(best_line)
                    out.write('\n')
                    # Reset
                    best_score = 0
            this_score = int(fields[4])
            if this_score > best_score:
                best_line = '\t'.join(fields)
                best_score = this_score

            # increment current_cluster every time
            current_cluster = this_clus_n
    
    # last one
    out.write(best_line)
    out.write('\n')

    # clean up
    subprocess.run(f'rm -f {passed_mir_bed}', shell=True)

    out.close()
    return outfile

def call_strucVis(bed, args, bam):
    # Create directories for the outputs
    psdir = args.outdir + '/strucVis'
    subprocess.run(f'mkdir {psdir}', shell=True)
    
    # get locus information into memory
    # make an interable for use in multiprocessing
    print('')
    print ('Creating visualizations of microRNA loci with strucVis')
    bedlines = []

    with open(bed) as f:
        for line in f:
            bedlines.append(line.rstrip())
    
    for line in (tqdm(bedlines, desc='Loci processed', leave=None)):
        bedfields = line.rstrip().split('\t')

        # Build strucVis command
        cmd = f'strucVis -b {bam} -g {args.genomefile}'
        onestart = int(bedfields[1]) + 1 # zero-based to one-based adjustment
        coords = str(bedfields[0]) + ':' + str(onestart) + '-' + str(bedfields[2])
        cmd = cmd + f' -c {coords}'
        if bedfields[5] == '-':
            strand = 'minus'
        else:
            strand = 'plus'
        cmd = cmd + f' -s {strand}'
        psfilebase = psdir + '/'
        locname = ''
        if bedfields[3] != '.':
            psfilebase = psfilebase + bedfields[3] + '__'
            locname = bedfields[3] + '__'
        san_coords = coords.replace(':', '_')
        psfile = psfilebase + san_coords + '.ps'
        txtfile = psfilebase + san_coords + '.txt'
        locname = locname + san_coords
        cmd = cmd + f' -p {psfile} -n {locname}'
        cmd = cmd + f' > {txtfile}'

        # Execute the strucVis command
        subprocess.run(cmd, shell=True)
    
def write_files(args, qdata, mir_qdata):
    print('')
    print('Writing final files')
    rgs = get_read_groups(merged_bam)
    # Counts.txt file, only if there are 2 or more read groups
    if rgs:
        if len(rgs) > 1:
            # Coords Name MIRNA RG1 ... RGX
            countsf = args.outdir + '/Counts.txt'
            countsfh = open(countsf, "w")
            header = 'Coords\tName\tMIRNA\t' + '\t'.join(rgs) + '\n'
            countsfh.write(header)

    # Results.txt file
    resultsf = args.outdir + '/Results.txt'
    resultsfh = open(resultsf, "w")
    res_header = 'Locus\tName\tChrom\tStart\tEnd\tLength\tReads\tUniqueReads\t'
    res_header = res_header + 'FracTop\tStrand\tMajorRNA\tMajorRNAReads\t'
    res_header = res_header + 'Short\tLong\t'
    for rlen in range(args.dicermin, args.dicermax + 1, 1):
        res_header = res_header + str(rlen) + '\t'
    res_header = res_header + 'DicerCall\tMIRNA\n'
    resultsfh.write(res_header)

    # Gather bedlines to examine
    # bed files to draw from:
    # deNovoClusters.bed if it exists; 
    # .  if not, filtDeNovoClusters.bed if it exists
    # .    if not, user.bed (user loci)
    # passed_mir_final.bed if it exists
    nm_bedlines = [] # not MIRNA
    m_bedlines = [] # MIRNA
    
    dn = args.outdir + '/deNovoClusters.bed'
    fdn = args.outdir + '/filtDeNovoClusters.bed'
    u = args.outdir + '/user.bed'
    m = args.outdir + '/passed_mir_final.bed'
    
    if os.path.exists(dn):
        with open(dn) as dnh:
            bedlines = dnh.readlines()
        nm_bedlines = nm_bedlines + bedlines
    elif os.path.exists(fdn):
        with open(fdn) as fdnh:
            bedlines = fdnh.readlines()
        nm_bedlines = nm_bedlines + bedlines
    elif os.path.exists(u):
        with open(u) as uh:
            bedlines = uh.readlines()
        nm_bedlines = nm_bedlines + bedlines
    
    if os.path.exists(m):
        with open(m) as mh:
            bedlines = mh.readlines()
        m_bedlines = m_bedlines + bedlines
    
    # Now report out
    if nm_bedlines:
        for line in nm_bedlines:
            f = line.rstrip().split('\t')
            onestart = int(f[1]) + 1
            coords = f[0] + ':' + str(onestart) + '-' + f[2]
            if len(f) > 3:
                name = f[3]
            else:
                name = 'NA'
            
            # Results.txt
            rline = f'{coords}\t{name}\t'
            rline = rline + qdata[coords]['chrom'] + '\t'
            rline = rline + str(qdata[coords]['start']) + '\t'
            rline = rline + str(qdata[coords]['end']) + '\t'
            rline = rline + str(qdata[coords]['Length']) + '\t'
            rline = rline + str(qdata[coords]['Reads']) + '\t'
            rline = rline + str(qdata[coords]['UniqueReads']) + '\t'
            rline = rline + str(qdata[coords]['FracTop']) + '\t'
            rline = rline + qdata[coords]['Strand'] + '\t'
            rline = rline + qdata[coords]['MajorRNA'] + '\t'
            rline = rline + str(qdata[coords]['MajorRNAReads']) + '\t'
            rline = rline + str(qdata[coords]['Short']) + '\t'
            rline = rline + str(qdata[coords]['Long']) + '\t'
            for rlen in range(args.dicermin, args.dicermax + 1, 1):
                rline = rline + str(qdata[coords][str(rlen)]) + '\t'
            rline = rline + qdata[coords]['DicerCall'] + '\t'
            rline = rline + 'N\n'
            resultsfh.write(rline)

            # Counts.txt
            if rgs:
                if len(rgs) > 1:
                    cline = f'{coords}\t{name}\tN'
                    for rg in rgs:
                        cline = cline + '\t' + str(qdata[coords][rg])
                    cline = cline + '\n'
                    countsfh.write(cline)
    if m_bedlines:
        for line in m_bedlines:
            f = line.rstrip().split('\t')
            onestart = int(f[1]) + 1
            coords = f[0] + ':' + str(onestart) + '-' + f[2]
            if len(f) > 3:
                name = f[3]
            else:
                name = 'NA'
            
            # Results.txt
            rline = f'{coords}\t{name}\t'
            rline = rline + mir_qdata[coords]['chrom'] + '\t'
            rline = rline + str(mir_qdata[coords]['start']) + '\t'
            rline = rline + str(mir_qdata[coords]['end']) + '\t'
            rline = rline + str(mir_qdata[coords]['Length']) + '\t'
            rline = rline + str(mir_qdata[coords]['Reads']) + '\t'
            rline = rline + str(mir_qdata[coords]['UniqueReads']) + '\t'
            rline = rline + str(mir_qdata[coords]['FracTop']) + '\t'
            rline = rline + mir_qdata[coords]['Strand'] + '\t'
            rline = rline + mir_qdata[coords]['MajorRNA'] + '\t'
            rline = rline + str(mir_qdata[coords]['MajorRNAReads']) + '\t'
            rline = rline + str(mir_qdata[coords]['Short']) + '\t'
            rline = rline + str(mir_qdata[coords]['Long']) + '\t'
            for rlen in range(args.dicermin, args.dicermax + 1, 1):
                rline = rline + str(mir_qdata[coords][str(rlen)]) + '\t'
            rline = rline + mir_qdata[coords]['DicerCall'] + '\t'
            rline = rline + 'Y\n'
            resultsfh.write(rline)

            # Counts.txt
            if rgs:
                if len(rgs) > 1:
                    cline = f'{coords}\t{name}\tY'
                    for rg in rgs:
                        cline = cline + '\t' + str(mir_qdata[coords][rg])
                    cline = cline + '\n'
                    countsfh.write(cline)

    resultsfh.close()
    if rgs:
        if len(rgs) > 1:
            countsfh.close()
    

def get_read_groups(bam):
    rgs = []
    rg = subprocess.run(f'samtools view -H {bam} | grep \'@RG\'', 
    shell=True, text=True, capture_output=True)
    rgpattern = re.compile('ID:(\S+)$')
    rglines = rg.stdout.split('\n')
    for rgline in rglines:
        if rgline:
            if rgpattern.search(rgline):
                this_rg = rgpattern.search(rgline).group(1)
                rgs.append(this_rg)
    return rgs

def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z

## Main control / script
if __name__ == '__main__':
    # Process command line, setting and checking options.
    args = ShortStack_argparse(version)
    
    # Report settings to user.
    start_message(args, version)

    # Check for required executables.
    check_executables(args)

    # Genome preparation. Gives dictionary fai, which has chr names and lengths
    fai = genome_prep(args)

    # locus / locifile check. If neither present, the loc object is None.
    # otherwise, loci is a dictionary of loci mapped to names.
    loci = loc_check(args, fai)

    # Create the output directory (or die tryin)
    make_outdir(args.outdir)

    # Adapter discovery and adapter trimming
    trimmed_files = get_trimmed(args)
    
    # Alignment
    merged_bam = align(args, fai, trimmed_files)
    if args.align_only is True:
        sys.exit()
    
    # Analysis
    cluster_bed = denovo_clusters(merged_bam, args, fai, loci) 
    qdata, pmir_bedfile = quant_controller(args, merged_bam, cluster_bed)
    #pmir_bedfile = None
    
    # MIRNA search
    mir_qdata = mirna(args, merged_bam, fai, pmir_bedfile)

    # Clean up a possible intermediate file that for some cryptic reason
    #  is not able to be deleted earlier in the script
    todel = args.outdir + '/possible_mirs_clean.bed'
    if os.path.exists(todel):
        subprocess.run(f'rm {todel}', shell=True)

    # Write Results.txt and Counts.txt
    write_files(args, qdata, mir_qdata)

    # You're done!
    print('')
    print(time.strftime("%a %d %b %Y %H:%M:%S %z %Z", time.localtime()))
    print('Run Completed!')







    
    


