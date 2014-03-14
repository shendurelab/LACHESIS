#!/usr/bin/perl -w
use strict;

# make_bed_around_restriction_site.pl: Make a BED file representing the regions around all occurrences of a restriction site.
#
# For syntax, run with no arguments.
#
# The output BED file is designed for use with bedtools intersect, as follows:
# bedtools intersect -abam [SRR.bam] -b [$BED_out] > [SRR.REduced.bam]
# samtools view -h [SRR.REduced.bam] > [SRR.REduced.sam]
# This restricts a SAM/BAM file to only include reads close to a restriction site, which is a good way to filter Hi-C data, according to Fig. 1b of this paper:
# http://www.nature.com/ng/journal/v43/n11/full/ng.947.html
# Also see PreprocessSAM.pl, which uses the output file.
#
# Josh Burton
# April 2013




if ( scalar @ARGV != 3 ) {
    
    # Report syntax.
    print "$0: Make a BED file representing the regions around all occurrences of a restriction site.\n\n";
    print "Input:\n";
    print "1. A fasta file representing a genome (reference or draft assembly.)\n";
    print "2. A restriction site sequence (e.g., HindIII = AAGCTT, NcoI = CCATGG, Dpn1 = GATC).\n";
    print "3. A range, representing how much space around the sequence to include.  Recommend 500 basedon Yaffe & Tanay, Nat. Genetics 2011.\n";
    print "\n";
    exit;
}



# Get command-line arguments.
my ( $FASTA_in, $motif_seq, $range ) = @ARGV;

my $verbose = 0;



# Derive an output filename.
my $BED_out = "$FASTA_in.near_$motif_seq.$range.bed";


# Determine how many letters needed to be added to each line in order to find instances of the sequence that bridge lines in the fasta.
my $N_prev_chars = length($motif_seq) - 1;


my $contig_name = '';
my $offset = 0;
my $prev_chars;
my @motif_positions;


# Open the input fasta file and read through it line-by-line.
print localtime() . ": Opening file $FASTA_in (this may take a while)\n";
open IN, '<', $FASTA_in or die "Can't find file `$FASTA_in'";
open BED, '>', $BED_out or die;

while (<IN>) {
    my $line = $_;
    chomp $line;
    
    # If this is a header line, we're done with this contig/chromosome (unless we just started), and start a new contig/chromosome.
    if ( $line =~ /^\>(\S+)/ ) {
	
	# The hash %motif_positions contains all positions on the (now complete) old contig at which this motif appears.
	# Convert this list of positions to a set of BED lines, as necessary.
	my ( $prev_start, $prev_end ) = (-1,-1);
	foreach my $pos ( @motif_positions ) {
	    if ( $prev_end == -1 ) {
		$prev_start = $pos;
		$prev_end   = $pos;
	    }
	    if ( $prev_end + 2*$range < $pos ) {
		$prev_start =         $range if $prev_start < $range;
		$prev_end = $offset - $range if $prev_end > $offset - $range; # prevent overflow past the end of the contig/chromosome
		print BED "$contig_name\t", $prev_start - $range, "\t", $prev_end + $range, "\n";
		$prev_start = $pos;
	    }
	    #print "pos = $pos\n";
	    $prev_end = $pos;
	}
	
	# Print the final BED line for this contig/chromosome.
	if (@motif_positions) {
	    $prev_start =         $range if $prev_start < $range;
	    $prev_end = $offset - $range if $prev_end > $offset - $range; # prevent overflow past the end of the contig/chromosome
	    print BED "$contig_name\t", $prev_start - $range, "\t", $prev_end + $range, "\n";
	}
	
	# Get the new contig's name.
	$contig_name = $1;
	print localtime() . ": $contig_name\n" if $verbose;
	
	# Reset other contig-related variables.
	$offset = 0;
	$prev_chars = '';
	@motif_positions = ();
    }
    
    # Otherwise, read through this contig/chromosome.
    else {
	if ( $offset != 0 ) { die unless $prev_chars; }
	
	my $verbose = 0;
	
	# Look for instances of this motif in this line of the fasta (including the overlap characters from the previous line, tacked on at the beginning.)
	my $motif_loc = -1;
	while (1) {
	    $motif_loc = index "$prev_chars$line", $motif_seq, $motif_loc + 1;
	    last if ( $motif_loc == -1 ); # no more instances found
	    
	    # Found a motif location!  Add this index to the list of contig positions at which the motif has been seen.
	    # Adjust the index so it properly describes the zero-indexed motif position in this contig.
	    push @motif_positions, $motif_loc + $offset - length $prev_chars;
	    
	    print "$contig_name\t$offset\t$prev_chars\t->\t$motif_loc\n" if $verbose;
	}
	
	
	# Save the last few characters of this line, so that they can be appended onto the next line in a search for the sequence.
	my $line_len = length $line;
	$prev_chars = substr( $line, $line_len - $N_prev_chars );
	$offset += $line_len;
    }
    
    
}

close IN;
close BED;


print localtime() . ": Done!\n";

