#!/usr/bin/perl -w
use Getopt::Long;

$version = "0.1.1-DEV";

# Define a usage statement
$usage = "$0 $version
Input a sam, sam\.gz, or \/bam file sorted by read names,
Tally mappings per read and note with NH:i: tags,
Output a \.bam file sorted by chromosomal location, and a corresponding \.bai index file\.

NOTES:
A\. If your input already has NH:i: tags, the orignal tags will be over-written\!
B\. Ensure that your input is sorted by read name\!
C\. Ensure that your input has CIGAR strings entered \(column 5 of the SAM format\) if you want to analyze the data with ShortStack\.pl

USAGE:
1\. For \.bam files : samtools view [input\.bam\] | $0 --prefix [prefix] --genome [genome\.fasta]\n
2\. For \.sam\.gz files : gzip -d -c [input\.sam\.gz] | $0 --prefix [prefix] --genome [genome\.fasta]\n
3\. For uncompressed \.sam files : $0 --prefix [prefix] --genome [genome\.fasta] < [input\.sam]


OPTIONS \(both required\)
-- prefix [string] : base name for the final \.bam and \.bam\.bai files that will be created
-- genome [string] : Path to the fasta formatted genome file used for the mapping

Type perldoc $0 for more details or see the README

";


# initial option definitions
$prefix = '';
$genome = '';

# get user options from command line
GetOptions ('prefix=s' => \$prefix,
	    'genome=s' => \$genome);

# Verify existence of prefix
unless($prefix) {
    die "FATAL: Please provide a prefix for your file name with option --prefix and try again\!\n\n$usage\n";
}

# Verify genome file is readable
unless(-r $genome) {
    die "FATAL: Please provide the genome fasta file used for mapping with option --genome and try again\!\n\n$usage\n";
}

# Verify that samtools is present and executable
# check for required installation of samtools and get version, or quit and complain
($samtools_version,$full_samtools_test_text) = get_samtools_version();
if($samtools_version eq "not found") {
    die "samtools not found\n$full_samtools_test_text\n\n$usage\n";
}


# Initialize temp file
$tempfile = "$prefix" . "\.sam\.gz";
if(-e $tempfile) {
    die "FATAL: File $tempfile already exists\.  Please choose another prefix and try again\n\n$usage\n";
}
open(TEMP, "| gzip > $tempfile");  ## compress on the fly

# check on other temp files that will need to be created
$tempfile2 = "$prefix" . "_temp" . "\.bam";
if(-e $tempfile2) {
    die "FATAL: File $tempfile2 already exists\.  Please pick another prefix and try again\n\n$usage\n";
}

$final_file = "$prefix" . "\.bam";
$final_file_index = "$final_file" . "\.bai";
if((-e $final_file) or
   (-e $final_file_index)) {
    die "FATAL: $final_file or $final_file_index already exists\.  Please choose another prefix and try again\n\n$usage\n";
}

# Add NH:i: tags and output to the temp .sam.gz file, still will be sorted by read name at this point

$x = 0;
$n = 0;
$j = 0;

%freqs = ();
%sizes = ();  # if CIGAR strings available, size distribution of mapped reads will be reported

$last_read = "not_a_read";
@storage = (); 
while (<STDIN>) {
    chomp;
    if($_ =~ /^\@/) {
	print TEMP "$_\n";
	next;  ## ignore header lines
    }
    if($_ =~ /^\#/) {
	print TEMP "$_\n";
	next;  ## also ignore comment lines if present (not expected to be)
    }
    
    @fields = split ("\t", $_);
    ++$x;
    if ($x == 100000) {
	print STDERR "\.";
	$x = 0;
	++$n;
	if($n == 10) {
	    $n = 0;
	    ++$j;
	    print STDERR " $j Million Mappings and Counting\n";
	}
    }
    
    # check for proper CIGAR string, and if not present, make a note and warn the user at the end
    if($fields[5] eq "\*") {
	$cigar_warning = 1;
    }
    
    if($last_read ne $fields[0]) {  ## you have passed into the next entry, or the first one, so evaluate if there is any data
	# get read size of CURRENT mapping, if its mapped
	unless($fields[1] & 4) {
	    $read_size = parse_cigar($fields[5]);
	    ++$sizes{$read_size};
	}
	if((scalar @storage) > 0) {
	    ++$freqs{$mappings};
	    ## process the previous read, output
	    ## need to add NH:i tag unless modified output suppressed
	    $nh_i = "NH:i:" . "$mappings";
	    foreach $stored_line (@storage) {
		chomp $stored_line;
		$stored_line =~ s/\tNH:i:\d+//g;  ## remove any existing NH:i: flags
		$stored_line .= "\t$nh_i";  ## add the manually found one
		print TEMP "$stored_line\n";  ## output
	    }
	}
	# reset 
	@storage = ();
	$mappings = 0;
    }
    ## store the current line
    push(@storage,$_);
    
    ## check the flag to see if read is unmapped
    ## examines the flag-sum to see if either '4' or '8' is present (see SAM specification)
    $mapped = check_mapped($fields[1]);
    ## above $mapped is zero if not mapped, 1 if mapped
    
    if($mapped) {
	++$mappings;
	## for mapped reads, get the flags, if present, for the current mapping
    }
    $last_read = $fields[0];
}

## don't forget to analyze the last read in the file
if((scalar @storage) > 0) {
    ++$freqs{$mappings};
    ## process the previous read, output
    ## need to add NH:i tag
    $nh_i = "NH:i:" . "$mappings";
    foreach $stored_line (@storage) {
	chomp $stored_line;
	$stored_line =~ s/\tNH:i:\d+//g;  ## remove any existing NH:i: flags
	$stored_line .= "\t$nh_i";  ## add the manually found one
	print TEMP "$stored_line\n";  ## output
    }
}
close TEMP;

print STDERR "done with phase 1\n";
print STDERR "\nSmall RNA Length Distribution of Mapped Reads\n";
print STDERR "Length\tFrequency\n";
# first show the size distribution of the mapped reads
@sorted_sizes = sort {$a <=> $b} (keys %sizes);
for($i = $sorted_sizes[0]; $i <= $sorted_sizes[-1]; ++$i) {
    print STDERR "$i\t";
    if(exists($sizes{$i})) {
	print STDERR "$sizes{$i}\n";
    } else {
	print STDERR "0\n";
    }
}
print STDERR "\n";

@f = sort {$a <=> $b} (keys %freqs);
$total_mapped = 0;
$total_mappings = 0;
print STDERR "Mappings-per-read\tFrequency\n";
foreach $fre (@f) {
    print STDERR "$fre\t$freqs{$fre}\n";
    unless($fre == 0) {  ## unmapped reads don't count towards total mappings!
	$total_mapped += $freqs{$fre};
	$total_mappings += ($freqs{$fre} * $fre);
    }
}
print STDERR "TOTAL READS MAPPED: $total_mapped\n";
print STDERR "TOTAL MAPPINGS: $total_mappings\n";

if($cigar_warning) {
    print STDERR "\n\*\*\* WARNING \*\*\*\n";
    print STDERR "CIGAR strings are missing from at least some of your mapped reads\.  This alignment therefore cannot be analyzed by ShortStack\.pl\.  Use an aligner that outputs proper CIGAR strings in the SAM/BAM format\.\n";
}

# use samtools to create a location-sorted .bam file from the read-sorted .sam.gz temp file
print STDERR "\nCreating initial \.bam file\n";
system "samtools view -bT $genome $tempfile > $tempfile2";
print STDERR "\nSorting bamfile by chromosomal location\n";
system "samtools sort $tempfile2 $prefix";
print STDERR "\nIndexing the sorted file\n";
system "samtools index $final_file";

# clean up
print STDERR "Deleting temp files\n";
system "rm -f $tempfile";
system "rm -f $tempfile2";
print STDERR "Completed\n";
print STDERR "Results in files $final_file and its index $final_file_index\n";
exit;


sub check_mapped { ## given the flag field from a sam data line
    my($digit) = @_;
    my $mapped;
    
    if(($digit & 4) or
       ($digit & 8)) {
	$mapped = 0;
    } else {
	$mapped = 1;
    }
    return $mapped;
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

__END__

=head1 LICENSE

Prep_bam.pl

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

Input a sam, sam.gz, or bam file sorted by read names,                                                                            
Tally mappings per read and note with NH:i: tags,                                                                                    
Output a .bam file sorted by chromosomal location, and a corresponding .bai index file.

=head1 INSTALL

install samtools from <http://samtools.sourceforge.net/> and ensure that samtools is in your PATH

ensure the script is executable                                                                  
                                                                                                 
    chmod +x Prep_bam.pl                                                         
                                                                                                 
ensure the script is in your PATH (examples):                                                    
                                                                                                 
    sudo cp Prep_bam.pl /usr/bin/                                                
                                                                                                 
OR just for one session assuming script is in your working directory:                            
                                                                                                 
    PATH=$PATH:.                                                                                 
                                                                                                 
ensure 'perl' is located in /usr/bin/ .. if not, edit line 1 of script accordingly                 

=head1 USAGE
                                                                                                                             
1. For .bam files : samtools view [input.bam] | $0 --prefix [prefix] --genome [genome.fasta]                                  

2. For .sam.gz files : gzip -d -c [input.sam.gz] | $0 --prefix [prefix] --genome [genome.fasta]                          

3. For uncompressed .sam files : $0 --prefix [prefix] --genome [genome.fasta] < [input.sam]      

=head1 OPTIONS

-- prefix [string] : base name for the final .bam and .bam.bai files that will be created    
                                      
-- genome [string] : Path to the fasta formatted genome file used for the mapping           

=head1 NOTES

=head2 Input format and assumptions

It is critical that the input sam, sam.gz, or bam file be sorted by read name. Additionally, it is assumed that read names are unique within the file.  If either of these assumptions are not true, the results of this script are not reliable.

=head2 Exisiting NH:i: tags

If there are already NH:i: tags in the input data, they will be discarded and over-written by the script.

=head2 CIGAR strings

The presence of proper CIGAR strings on all mappings is checked, and a warning is issued if any are missing.  The CIGAR string is used by ShortStack.pl to determine the length of the mapped small RNA, so if any are missing, the alignment cannot be analyzed by ShortStack.pl.  If that happens, consider using another aligner, or change your alignment protocol, so that proper CIGAR strings are present in column 6 of the SAM format (see SAM specification for details).

=head1 VERSIONS

0.1.0 : Initial release. April 29, 2012

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=cut


