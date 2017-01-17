#!/usr/bin/perl -w
use strict;

#############################################################################
#                                                                           #
# This software and its documentation are copyright (c) 2014-2015 by Joshua #
# N. Burton and the University of Washington.  All rights are reserved.     #
#                                                                           #
# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  #
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                #
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  #
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      #
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT #
# OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR  #
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                #
#                                                                           #
#############################################################################



# CountMotifsInFasta.pl
#
# Find the number of instances of the motif in each contig.
#
# Input: (1) a fasta file; (2) a motif (e.g., a restriction site).
#
# Report a list of contig names with their motif counts.  This file should be entered into the Lachesis INI file under the key "RE_COUNTS_FILE".
#
# Protip: If you want the combined number of instances of two motifs in each contig, do the following:
# CountMotifsInFasta.pl assembly.fasta motif1
# CountMotifsInFasta.pl assembly.fasta motif2
# paste assembly.fasta.counts_{motif1,motif2}.txt | awk 'OFS="\t" {print $1,$2+$4}' > assembly.fasta.counts_motif1_motif2.txt
#
#
# Josh Burton
# April Fools' Day, 2013






sub add_contig_data($$);



my $verbose = 0;


# Get input arguments.
die "\nCountMotifsInFasta: Count the number of occurrences of a motif in each contig of a fasta file.\nSyntax: $0 <assembly-fasta> <restriction site motif>\ne.g.: $0 mouse/assembly/assembly.fasta AAGCTT\nTo simultaneously count the occurrences of multiple motifs at once, separate them with an underscore, e.g.: 'AAGCTT_CCATGG'\n\nCreates a file at <assembly-fasta>.counts_<motif>.txt\n\n" unless ( scalar @ARGV == 2 );
my ($fasta,$motif) = @ARGV;

# Check that input arguments make sense.
die "Can't find assembly fasta file $fasta: $!" unless -e $fasta;

# Conver the single motif token into multiple motifs.
my @motifs = split /_/, $motif;


# Convert motifs into regexes, to account for RC'ing sequences and generalized IUPAC codes (e.g., N = any of ACGT, R = A or G, etc.)
my @motif_regexes;
foreach my $motif (@motifs) {
    die "Nonsensical motif sequence '$motif': Must contain only A,C,G,T and other IUPAC codes" if $motif =~ /[^ACGTRYSWKMBDHVN]/;
    
# Make sure the motifs are all self-RCs.
    my $motif_rc = reverse $motif;
    $motif_rc =~ tr/ACGTRYSWKMBDHV/TGCAYRWSMKVHDB/;
    warn "WARNING: Motif '$motif' is not its own reverse complement.  This script doesn't check for the sequence in rc.\n" unless $motif eq $motif_rc;
    
    # Make a reverse-complement copy of each motif.
    my $regex = $motif;
    my $regex_rc = reverse $regex;
    $regex_rc =~ tr/ACGTRYSWKMBDHV/TGCAYRWSMKVHDB/;
    
    # Unroll the IUPAC codes from single letters into Perl-parseable regular expressions.
    foreach ($regex,$regex_rc) {
	s/R/\[AG\]/g;
	s/Y/\[CT\]/g;
	s/S/\[CG\]/g;
	s/W/\[AT\]/g;
	s/K/\[GT\]/g;
	s/M/\[AC\]/g;
	s/B/\[CGT\]/g;
	s/D/\[AGT\]/g;
	s/H/\[ACT\]/g;
	s/V/\[ACG\]/g;
	s/N/\[ACGT\]/g;
    }
    
    print "$motif\t->\t$regex\t$regex_rc\n";
    push @motif_regexes, $regex;
}








print localtime(). ": CountMotifsInFasta!\n";


# These arrays will contain the output data
my @contig_names;
my @contig_lens;
my @contig_counts;


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
    
    # If this is the beginning of a new contig: count motifs in the previous contig (if there is one), then start over.
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
    
    my $len = length $$rseq;
    
    #print localtime() . ": add_contig_data\t$name\tseq of length $len\n" if $verbose;
    
    my $count = 0;
    map { $count += @{[$$rseq =~ /$_/g]} } @motif_regexes; # one-liner taken from: http://www.webmasterworld.com/forum13/3853.htm
    
    # Add to global data structures.
    push @contig_names,  $name;
    push @contig_lens,   $len;
    push @contig_counts, $count;
}
