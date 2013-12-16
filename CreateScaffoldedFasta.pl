#!/net/gs/vol3/software/modules-sw/perl/5.14.2/Linux/RHEL6/x86_64/bin/perl -w
use strict;

# CreateScaffoldedFasta.pl
#
# Post-Lachesis step: Take the Lachesis output orderings/orientations and use them to create a new fasta file.
#
# Syntax:
# CreateScaffoldedFasta.pl <assembly-fasta-file> <OUTPUT_DIR>
#
# Inputs:
# -- Assembly fasta file.
# -- Lachesis OUTPUT_DIR; must have a subdirectory main_results/ containig a file clusters.by_name.txt and a file group*.ordering for each group
#
# Outputs:
# -- A file "Lachesis_assembly.fasta" in OUTPUT_DIR
#
# 
#
# Josh Burton
# November 2013




sub LoadFasta( $ );
sub WriteFasta( $$$$ );
sub RC( $ );


# Defaults in the absence of parameters.
my $assembly_fasta = "../fly/assembly/assembly.fasta";
my $output_dir = "out/Kc167";
my $unknown_gap_size = 10000; # default size to apply to gaps between contigs/scaffolds when the gap size is not known from Lachesis

# Read parameters from command-line arguments.
if ( scalar @ARGV == 2 ) {
    ( $assembly_fasta, $output_dir ) = @ARGV;
}



print localtime(). ": CreateScaffoldedFasta.pl with input fasta = $assembly_fasta, OUTPUT_DIR = $output_dir\n";

my $fasta_out = "$output_dir/Lachesis_assembly.fasta";



# Find the clusters.by_name.txt file.
my $clusters_file = "$output_dir/main_results/clusters.by_name.txt";
die "ERROR: Couldn't find clusters file ('clusters.by_name.txt') in $output_dir/main_results/.  Is this the right directory?  Was Lachesis successfully run?\n"
    unless -e $clusters_file;

# Find the group*.ordering files.  The number of these files indicates the number of groups.  They should appear with increasing ID values, starting at 0.
my $N_orderings = 0;
while ( 1 ) {
    last unless -e "$output_dir/main_results/group$N_orderings.ordering";
    $N_orderings++;
};

# There should be at least one ordering file!
unless ( $N_orderings ) {
    print "ERROR: Couldn't find any ordering files ('group*.ordering') in $output_dir/main_results/.  Is this the right directory?  Was Lachesis successfully run?\n";
    exit;
}



print localtime() . ": Found $N_orderings ordering files ('group*.ordering' in $output_dir/main_results/).\n";




# Read in the assembly fasta file.
print localtime() . ": Reading in sequences from assembly file $assembly_fasta\n";
my ( $fasta_names, $fasta_seqs ) = &LoadFasta( $assembly_fasta );
my $N_contigs = scalar keys %$fasta_seqs;
die unless $N_contigs == scalar @$fasta_names;
print localtime() . ": Found $N_contigs contigs/scaffolds in assembly.\n";



# Now, write to the output fasta file.  There will be three phases of writing::
# PHASE 1: For each ordering, create a single scaffold and write it to file.
# PHASE 2: For each contig that is clustered into a group but not ordered within that group, write it to file with a name indicating its group.
# PHASE 3: For each contig that is not clustered, write it to file.
open OUT, '>', $fasta_out or die;


# PHASE 1: For each ordering, create a single scaffold and write it to file.

my %contigs_used; # hash indicating which contigs have already been written.

# Read through each group*.ordering file in turn.  From each file, cobble together a scaffold sequence.
foreach my $i ( 0..$N_orderings-1 ) {
    my $ordering_file = "$output_dir/main_results/group$i.ordering";
    
    print localtime(). ": Creating a scaffold from file $ordering_file...\t";
    
    open IN, '< ', $ordering_file or die "Can't find file $ordering_file: $!";
    
    my $scaffold = '';
    my $N_contigs = 0;
    
    # Each noncommented line in an ordering file has 5 tokens: contig ID, contig name, contig orientation (1=rc), orientation score, gap size ('.' if unknown.)
    while (<IN>) {
	next if /^\#/; # skip commented lines
	my ( $ID, $name, $rc, $q, $gap_size ) = split;
	
	# Create the sequence of N's representing the gap.
	$gap_size = $unknown_gap_size if $gap_size eq '.';
	my $gap = 'N' x $gap_size;
	
	die "ERROR: Ordering file $ordering_file includes contig named '$name', not found in fasta file $assembly_fasta\n" unless exists $fasta_seqs->{$name};
	
	$contigs_used{$name}++;
	
	# Reverse-complement the sequence if necessary.
	$fasta_seqs->{$name} = &RC( $fasta_seqs->{$name} ) if $rc;
	
	# Append the sequence, and its subsequent gap, to the scaffold in progress.
	$scaffold .= $fasta_seqs->{$name};
	$scaffold .= $gap;
	$N_contigs++;
    }
    
    
    # Remove the last gap.  It didn't mean anything anyway.
    $scaffold =~ s/N+$//;
    my $len = length $scaffold;
    print "length = $len\n";
    
    next unless $N_contigs; # we may get an empty scaffold if we have a group in which no contigs were ordered (e.g., a singleton)
    
    # Set the scaffold name.
    my $scaffold_name = "Lachesis_group${i}__${N_contigs}_contigs__length_$len";
    
    # Write the scaffold to file.
    &WriteFasta( \*OUT, $scaffold_name, $scaffold, 80 );
    
    close IN;
}




# PHASE 2: For each contig that is clustered into a group but not ordered within that group, write it to file with a name indicating its group.

# Open the clusters file and look through it for contigs that have been clustered but did not appear in the orderings.
print localtime() . ": Reading clusters file $clusters_file and finding contigs that have been clustered but not ordered...\n";
open IN, '<', $clusters_file or die "Can't find file $clusters_file: $!";
my $ID = 0;
my $N_unordered_contigs = 0;

while (<IN>) {
    
    # Each (non-commented) line simply contains a tab-delimited list of contig names in this cluster.  The cluster ID is implied by the line order.
    next if /^\#/; # skip commented lines
    my @names = split;
    
    foreach my $name (@names) {
	next if $contigs_used{$name};
	$contigs_used{$name}++;
	
	die "ERROR: Clusters file $clusters_file includes contig named '$name', not found in fasta file $assembly_fasta\n" unless exists $fasta_seqs->{$name};
	
	&WriteFasta( \*OUT, "${name}__unordered_in_group$ID", $fasta_seqs->{$name}, 80 );
	$N_unordered_contigs++;
    }
    
    $ID++;
}

close IN;

print localtime(). ": Wrote $N_unordered_contigs contigs that were clustered but not ordered by Lachesis.\n";



# PHASE 3: For each contig that is not clustered, write it to file.
print localtime() . ": Placing contigs that were not clustered\n";

my $N_unclustered_contigs = 0;

foreach my $name ( @$fasta_names ) {
    next if $contigs_used{$name};
    die "ERROR: Clusters file $clusters_file includes contig named '$name', not found in fasta file $assembly_fasta\n" unless exists $fasta_seqs->{$name};
    
    &WriteFasta( \*OUT, "${name}__unclustered", $fasta_seqs->{$name}, 80 );
    $N_unclustered_contigs++;
} 

print localtime(). ": Wrote $N_unclustered_contigs contigs that were not clustered by Lachesis.\n";


close OUT;


print localtime(). ": Done!\n";









# SUBROUTINES







# LoadFasta: Convert a fasta file to contigs.
# Outputs:
# 1. An array of contig names.
# 2. A hash of contig name to contig sequence.
sub LoadFasta( $ ) {
    
    #print localtime() . ": LoadFasta: $_[0]\n";
    
    open IN, '<', $_[0] or die;
    
    my $contig_name;
    my @contig_names;
    my %contig_seqs;
    while (<IN>) {
	chomp;
	if ( /^\>(.+)/ ) {
	    $contig_name = $1;
	    push @contig_names, $contig_name;
	}
	else {
	    $contig_seqs{$contig_name} .= $_;
	}
    }
    
    close IN;
    
    die "ERROR: LoadFasta: Couldn't parse file $_[0] properly.  Are you sure this is a FASTA file?" unless scalar @contig_names >= 1 && scalar keys %contig_seqs >= 1;
    
    return ( \@contig_names, \%contig_seqs );
}







# WriteFasta: Write a contig/scaffold to a fasta file.
# Arguments: output filehandle; contig/scaffold name; contig/scaffold sequence; number of sequence characters per line (70 or 80 is standard)
sub WriteFasta( $$$$ ) {
    
    my ( $fh, $name, $seq, $N_chars_per_line ) = @_;
    
    # Write the name.
    print $fh ">$name\n";
    
    my $seq_len = length $seq;
    my $i = 0;
    
    # Write the sequence, breaking it up into separate lines.
    while ( $i < $seq_len ) {
	print $fh substr( $seq, $i, $N_chars_per_line ), "\n";
	$i += $N_chars_per_line;
    }
    
}



# RC: Reverse-complement a sequence.
sub RC($ ) {
    my $seq = reverse $_[0];
    $seq =~ tr/ACGTacgtRYrySWswKMkmBDHVbdhv/TGCAtgcaYRyrWSwsMKmkVHDBvhdb/; # consider all possible IUPAC codes (except Ns, those can stay the same)
    return $seq;
}
