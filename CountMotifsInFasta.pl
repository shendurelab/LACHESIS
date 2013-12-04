#!/usr/bin/perl -w
use strict;


# CountMotifsInFasta.pl
#
# Input: (1) a fasta file; (2) a motif (e.g., a restriction site).
#
# Find the number of instances of the motif in each contig.
#
# Report a list of contig names with their motif counts.  This file should be entered into the Lachesis INI file under the key "RE_COUNTS_FILE".
#
#
# Josh Burton
# April Fools' Day, 2013






sub add_contig_data($$);



my $verbose = 0;


# Get input arguments.
die "Syntax: $0 <assembly-fasta> <restriction site motif>\ne.g.: $0 mouse/assembly/assembly.fasta AAGCTT\n" unless ( scalar @ARGV == 2 );
my ($fasta,$motif) = @ARGV;

# Check that input arguments make sense.
die "Can't find assembly fasta file $fasta: $!" unless -e $fasta;

# AAGCTT is HindIII restriction site (most Hi-C datasets including human, mouse); GATC is DpnII restriction site (fly)
die "Nonsensical motif sequence '$motif': Must contain only A,C,G,T" if $motif =~ /[^ACGT]/;
my $motif_rc = reverse $motif;
$motif_rc =~ tr/ACGT/TGCA/;
warn "WARNING: Motif '$motif' is not its own reverse complement.  This script doesn't check for the sequence in rc.\n" unless $motif eq $motif_rc;








print localtime(). ": CountMotifsInFasta!\n";


# These arrays will contain the output data
my @contig_names;
my @contig_lens;
my @contig_counts;


my $count = 0;
my $contig_name = '';
my $contig_seq = '';
my $N_lines = 0;

# Read the FASTA file line-by-line.
print localtime() . ": Reading file $fasta...\n";
open IN, '<', $fasta or die;

while (<IN>) {
    chomp;
    my $line = $_;
    
    print localtime() . ": Each dot . represents 1M lines\n" if ( $verbose && !$N_lines );
    $N_lines++;
    
    # If this is the beginning of a new contig: count motifs in the previous contig *if there is one), then start over.
    if ( $line =~ /^\>([\w\|\.]+)/ ) {
	my $new_name = $1;
	
	if ( $N_lines != 1 ) {
	    add_contig_data($contig_name,\$contig_seq);
	}
	$contig_name = $new_name;
	$contig_seq = '';
    }
    
    # Otherwise, just add to the growing sequence of this contig.
    else {
	$contig_seq .= $line;
    }
    
    if ( $verbose && $N_lines % 100000 == 0 ) { print '.'; $_++; }
}



close IN;


# Report the final contig!
add_contig_data($contig_name,\$contig_seq);


print "\n" if $verbose;




# Open output file for writing.
my $outfile = "$fasta.counts_$motif.txt";
print localtime(). ": Done reading.  $N_lines lines read.  Now writing to $outfile\n";




# Print results to the output file.
open OUT, '>', $outfile or die;

my $N_contigs = scalar @contig_names;
foreach (0..$N_contigs-1) {
    my $ratio = ( $contig_counts[$_] ? $contig_lens[$_] / $contig_counts[$_] : "inf" );
    if ( $verbose ) { print "$contig_names[$_]\t$contig_lens[$_]\t$contig_counts[$_]\t$ratio\n"; }
    print OUT "$contig_names[$_]\t$contig_counts[$_]\n";
}

close OUT;

print localtime() . ": Done!\n";





sub add_contig_data($$)
{
    my ($name,$rseq) = @_;
    
    # TEMP: apparently at some point I wanted to only include the 10Mb on either end of the contig.  I'm not sure why this was ever a thing, but it should
    # have had no effect because most draft assemblies don't have contigs/scaffolds bigger than 20 Mb.
    #my $max_contig_length = 10000000;
    #if ( length $$rseq > 2 * $max_contig_length ) {
	#$$rseq = substr( $$rseq, 0, $max_contig_length ) . substr( $$rseq, -$max_contig_length );
    #}
    
    my $len = length $$rseq;
    
    #print localtime() . ": add_contig_data\t$name\tseq of length $len\n" if $verbose;
    
    my $count = @{[$$rseq =~ /$motif/g]}; # one-liner taken from: http://www.webmasterworld.com/forum13/3853.htm
    
    # Add to global data structures.
    push @contig_names,  $name;
    push @contig_lens,   $len;
    push @contig_counts, $count;
}
