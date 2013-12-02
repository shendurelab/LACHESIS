#!/usr/bin/perl -w
use strict;

# commentify_INI.pl
#
# A simple script to add comments (and blank lines) to an INI file to make it readable.  The comments are hard-wired to accurately describe the contents of 
# an INI file if it has all of the right fields in the right order.  The counterpart of this script is decommentify_INIs.pl, which does the reverse.
#
# This is a convenience script that allows me to easily maintain consistent commenting across a set of INI files.  I will make all of my modifications to the
# comments in this file.  Changes in the number of key-value pairs should also be reflected here (in the variable $expected_N_lines.)
#
# Josh Burton
# August 2013


# Check for input arguments.
if ( @ARGV != 1 ) {
    print STDERR "\nSyntax: $0 <INI-file>   (writes to stdout)\n\n";
    exit;
}

my $expected_N_lines = 26;


# Open the filename provided, and get the lines out of it.
die "Can't find input file `$ARGV[0]`: $!\n" unless -e $ARGV[0];
open IN, '<', $ARGV[0] or die;
my @lines = <IN>;
map {chomp} @lines;
close IN;

# Check that each line describes a key-value pair.
map { die "ERROR: The following line in $ARGV[0] doesn't seem to describe a key-value pair:\n$_" unless /^\w+\s+\=\s/ } @lines;

die "ERROR: Wrong number of lines in $ARGV[0].  It should be $expected_N_lines, but saw ", scalar @lines, " instead.\n" unless scalar @lines == $expected_N_lines;


# Now, print the entire commented file.

print <<FILE_END;
###############################################################################################################################################################
#
#
#   Lachesis.ini - initialization file for Lachesis runs
#
#   This file controls how the Lachesis program is run.  Changing this file will change the inputs to Lachesis, which algorithms it runs, the heuristic
#   parameters it uses, and so forth.
#
#   When modifying Lachesis.ini, be careful about formatting!  Each line must be of the form "key = value", where key and value are text strings, with no
#   internal spaces.  There must be some space around the equals sign, and there should be no leading whitespace.
#   In some cases there can be multiple values, which should be separated by tabs or spaces, without commas.
#   Commented lines (such as this one) are ignored.  Do not append comments at the end of an otherwise non-commented line.  To add or remove comments en masse
#   from this file, use the scripts commentify_INI.pl and decommentify_INI.pl in INIs/.
#   You may add or remove empty lines and commented lines, but DO NOT remove any parameters, or change the order in which the parameters appear in the file.
#   Lachesis.ini is parsed in the module RunParams.cc, in the function ParseIniFile().  If this file is formatted incorrectly, this function throws an error.
#
#   Note that values for boolean keys must be either '0' or '1'.  Any other value (including e.g., "true") will fail.
#
#
#   Josh Burton
#   June 2013
#
#
###############################################################################################################################################################



# Species.  If assembling a species for which there is already a reference, use "human" or "mouse" or "fly"; these names are hard-wired into the code.
$lines[0]




#################################################
#
#   PARAMETERS FOR THE INPUT ASSEMBLY
#

# Directory containing the input SAM/BAM files that describe Hi-C read alignments to the draft assembly.
# This directory path can be absolute, or relative to the directory in which Lachesis is run.  This should be listed before SAM_FILES.
$lines[1]

# The input SAM/BAM files that describe Hi-C read alignments to the draft assembly.
# These files should exist in SAM_DIR, above, and they should already be pre-processed (e.g., to remove PCR duplicates) and sorted by read name (not by order).
# If any of these files fails to exist, Lachesis won't run.
$lines[2]

# A file containing the number of restriction enzyme (RE) sites per contig/scaffold in the draft assembly.
# This file should contain a row for each contig in the draft assembly, and two columns: one with the contig name, and one with an integer.
# For example, you can create the first column alone by running `grep "^\\>" [fasta-filename] | cut -c2-`.
# The second column will take more work.  Try running the script CountMotifsInFasta.pl.
$lines[3]






#################################################
#
#   PARAMETERS FOR THE REFERENCE GENOME
#
#   If you are assembling a genome for which a reference already exists, you can use these parameters to inform Lachesis of the reference genome sequence.
#   You must align your draft assembly to the reference genome.  Lachesis will then evaluate its performance by comparing its clustering, ordering, and
#   orienting results to the results implied by the alignments to reference.
#   If USE_REFERENCE = 0, none of these parameters are examined or used.
#

# Use a reference genome?  Options: 0 (false), 1 (true).  If 0, then SIM_BIN_SIZE, BLAST_FILE_HEAD, and REF_ASSEMBLY_NAMES_FILE are ignored.
$lines[4]

# If the draft assembly is just the reference genome chopped up into simulated bins, indicate the bin size here.  Otherwise set to 0.
$lines[5]

# File head for BLAST alignments.  You must align the draft assembly to the reference genome using BLAST (UNIX command: `blastn -outfmt 7 ...`).
# The output should go into a set of one or more files called <BLAST_FILE_HEAD>.*.blast.out, where * = 1,2,3,...
# Lachesis will create a file at <OUTPUT_DIR>/cached_data/TrueMapping.assembly.txt.  Once this file exists, you no longer need the BLAST files.
# Alternatively, if SIM_BIN_SIZE > 0 (above), BLAST_FILE_HEAD is ignored because no alignments to reference are needed.
$lines[6]

# 'Names files' for fastas.  These files simply list the names of the contigs in a fasta file.  They can be generated by running `grep "^\\>" [fasta-filename]`.
$lines[7]
$lines[8]





#################################################
#
#   WHICH ALGORITHMIC STEPS TO PERFORM?
#
#   For each of these keys, options are 0 (false), 1 (true).
#   The orientation step is considered part of the ordering step.
#   If USE_REFERENCE = 1 (above), there will be reference-assisted evaluation in each of clustering, ordering, and reporting.
#

$lines[9]
$lines[10]
$lines[11]

# At the beginning of clustering, the Hi-C links are loaded from the SAM files, and then the cache file <OUTPUT_DIR>/cached_data/all.GLM is created.
# If this cache file already exists, and if OVERWRITE_GLM=0, the links are loaded from cache, which saves time.  Set to 1 if the content of SAM_FILES changes.
$lines[12]

# At the beginning of ordering, the links are loaded from the SAM files, and then the cache files <OUTPUT_DIR>/cached_data/group*.CLM are created.
# If these cache files already exist, and if OVERWRITE_CLMS=0, the links are loaded from cache, which saves time.
# Set to 1 if you change anything about the clustering, so that the change will propagate to the ordering.  Otherwise Lachesis will throw an error.
$lines[13]





#################################################
#
#   HEURISTIC PARAMETERS FOR CLUSTERING AND ORDERING
#
#   For these keys, the accepted values are integers and decimal numbers.
#

# Number of clusters.  Set this to the number of chromosomes in the input assembly.
$lines[14]
# Only use contigs as informative for clustering if they have at least this many restriction enzyme (RE) sites.
$lines[15]
# Only use contigs as informative for clustering if they have LESS than this much times the average density of Hi-C links.
# Contigs with too many Hi-C links tend to be in heterochromatin or other repeat-rich regions.
$lines[16]
# Boolean (0/1).  If 1, then non-informative contigs (the ones that fail the CLUSTER_MIN_RE_SITES, CLUSTER_MAX_LINK_DENSITY tests) may be added to groups after
# clustering is over, if they fit cleanly into one group.
$lines[17]
# Boolean (0/1).  Draw a 2-D heatmap of the entire Hi-C link dataset before clustering.
$lines[18]
# Boolean (0/1).  Draw a 2-D dotplot of the clustering result, compared to truth.  This is time-consuming and eats up file I/O.  Ignored if USE_REFERENCE=0.
# The dotplots go to ~/public_html/dotplot.SKY.*.jpg
$lines[19]


# Minimum number of RE sites in contigs allowed into the initial trunk.
$lines[20]
# Minimum number of RE sites in shreds considered for reinsertion.
$lines[21]



#################################################
#
#   OUTPUT PARAMETERS
#
#   These parameters do not change the results of Lachesis, but they change where and how the results are reported.
#

# Directory which will contain the output.  This directory will have subdirectories called main_results/ and cached_data/.
$lines[22]

# Output file for the ReportChart, the main summarized piece of text output.  This file will be created in OUTPUT_DIR (if DO_REPORTING = 1).
$lines[23]

# IDs of groups to exclude from the reporting numbers (e.g. groups determined to be small, crappy, and/or chimeric.)  If not excluding any groups, set to "-1".
$lines[24]
# Quality filter.  Contigs whose orientation quality scores (differential log-likelihoods) are this or greater are considered high-quality.
# The quality scores depend on Hi-C read coverage, so you'll want to try out some values in order to achieve an informative differentiation in the ReportChart.
$lines[25]
FILE_END
