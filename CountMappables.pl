#!/net/gs/vol3/software/modules-sw/perl/5.14.2/Linux/RHEL6/x86_64/bin/perl -w
use strict;
use Memory::Usage;
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This software and its documentation are copyright (c) 2014-2015 by Joshua //
// N. Burton and the University of Washington.  All rights are reserved.     //
//                                                                           //
// THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  //
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      //
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT //
// OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR  //
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



# CountMappables.pl
#
# Estimate the 'mappability' of a set of contigs.  For each contig, mappability is defined as the fraction of reads aligning to that contig which have MQ>0.
#
# Input: a set of one or more BAM files.
# Output: a two-column file.  The first column gives a contig name (in the same order as the BAM header(s)), and the second gives the mappability fraction.
#
# Josh Burton
# November 2013




my @BAMs_in = glob "human/to_prefosmid/SRR??????.bam";
map { die unless -e } @BAMs_in;
#@BAMs_in = (@BAMs_in[1..3]);
#@BAMs_in = ( 'human/to_prefosmid/SRR400263.bam' );

# First get a list of the contig names.  Do this by parsing one of the BAM files' headers.
my @contig_names;
open HEADER, "samtools view -H $BAMs_in[0] |" or die;

foreach (<HEADER>) {
    next unless /^\@SQ\tSN\:(.+)\tLN\:/; # grab contig names from the header
    push @contig_names, $1;
}

close HEADER;


print localtime() . ": ", scalar @contig_names, " contigs found.\n";

print localtime() . ": Reading ", scalar @BAMs_in, " BAMs...\n";


# Determine how many read alignments to each contig have mapping quality = 0.
my ( %MQ0_aligns, %MQ1_aligns ); # "MQ1" means mapping quality > 0

# Read each BAM file line-by-line.
foreach my $BAM (@BAMs_in) {
    
    print localtime() . ": Reading file: $BAM\n";
    $|++;
    #open IN, "samtools view $BAM |" or die "Can't open file $BAM: $!"; # samtools view includes reads but excludes header
    #open IN, "samtools view $BAM | head -n100000 |" or die; # samtools view includes reads but excludes header
    open IN, '<', "$BAM.sam" or die "Can't open file $BAM: $!"; # use SAMs
    
    
    while (<IN>) {
	next if /^\@/; # skip header lines
	
	# We only care about the third token (contig name) and fifth token (mapping quality).
	my ($t1,$t2,$t3,$t4,$t5,@more) = split;
	next if $t3 eq '*'; # read doesn't map
	
	# Fill the data structure.
	if ( $t5 > 0 ) { $MQ1_aligns{$t3}++; }
	else           { $MQ0_aligns{$t3}++; }
	
	#print "$t3\t$t5\n";
    }
    
    close IN;
    
    print "Hash sizes:\t", scalar keys %MQ1_aligns, "\t", scalar keys %MQ0_aligns, "\n";
}



# Calculate and print results for each contig.
print localtime() . ": Results!\n";

foreach my $contig (@contig_names) {
    
    my $MQ0 = $MQ0_aligns{$contig} || 0;
    my $MQ1 = $MQ1_aligns{$contig} || 0;
    my $frac_MQ1 = $MQ1 / ( ($MQ0 + $MQ1) || 1 );
    
    
    print "DATA\t$contig\t$MQ0\t$MQ1\t$frac_MQ1\n";
    
}


print localtime() . ": Done!\n";
