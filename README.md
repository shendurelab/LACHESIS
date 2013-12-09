             |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
             |                                                  |
             |                                                  |
             |                     LACHESIS                     |
             |                                                  |
             |  A software tool to measure the thread of life.  |
             |                                                  |
             |                                                  |
             |     Created by Josh Burton (jnburton@uw.edu)     |
             |       in the Department of Genome Sciences       |
             |   at the University of Washington, Seattle, WA   |
             |                                                  |
             |    Paper:  http://dx.doi.org/10.1038/nbt.2727    |
             |                                                  |
             |                                                  |
             |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
             |                                                  |
             |                                                  |
             |       USER'S MANUAL and QUICK START GUIDE        |
             |                                                  |
             |                                                  |
             |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|



###############################################################################


                          |~~~~~~~~~~~~~~~~~~~~~~~~|
                          |                        |
                          |    TABLE OF CONTENTS   |
                          |                        |
                          |~~~~~~~~~~~~~~~~~~~~~~~~|


INSTALLING LACHESIS
1. System requirements
2. Downloading the LACHESIS package
3. Compiling the LACHESIS package
4. Testing out LACHESIS on a sample dataset

RUNNING LACHESIS
1. Input requirements
2. Creating an INI file
3. Aligning the Hi-C reads to the draft assembly
4. Filtering the Hi-C reads
5. Aligning the draft assembly to the reference assembly, if there is one
6. Counting the restriction sites in the draft assembly contigs/scaffolds
7. Extracting the assembly contig names
8. Running LACHESIS
9. Interpreting the LACHESIS results

TROUBLESHOOTING

COPYRIGHT AND DISCLAIMER

ACKNOWLEDGMENTS



###############################################################################



INSTALLING LACHESIS



1. System requirements

To setup and run LACHESIS, you will need a computer running in a UNIX environment with at least 16GB of memory, with the following software installed:
-- gcc, the C++ compiler (http://gcc.gnu.org/)
-- The boost C++ libraries (http://www.boost.org/)
-- The SAMtools toolkit for handling SAM/BAM files (http://samtools.sourceforge.net/)
You may also need the following software:
-- The short-read aligner BWA (http://bio-bwa.sourceforge.net/) or another such aligner
-- The BLAST aligner in command-line form (http://www.ncbi.nlm.nih.gov/books/NBK1763/)
-- The bedtools toolkit for handling genomic intervals (https://code.google.com/p/bedtools/)



2. Downloading the LACHESIS package

Download the LACHESIS package from http://shendurelab.github.io/LACHESIS/ into a UNIX filesystem.  If you download the tarball (LACHESIS.tar.gz), unpacki it using the following UNIX commands:

tar xzvf LACHESIS.tar.gz
cd LACHESIS/



3. Compiling the LACHESIS package

To compile LACHESIS, you must first download and install two other libraries: boost (available at www.boost.org) and SAMtools (available at samtools.sourceforge.net).  Once these are installed, set the shell environment variables $LACHESIS_BOOST_DIR and $LACHESIS_SAMTOOLS_DIR to point to the directories containing these packages.  The command for setting an environment variable will depend on what Unix shell you are using.  Then, to compile LACHESIS, simply type `make` in the main LACHESIS directory.



4. Testing out LACHESIS on a sample dataset

This LACHESIS distribution includes a small sample dataset - specifically, a reduced version of the hESC dataset - that allows you to try out running LACHESIS.  To run it, go into the LACHESIS directory and type `Lachesis INIs/test_case.ini`.  The file test_case.ini gives the parameters that control how LACHESIS is run, in this case including the input files from the sample dataset.  On my computer, LACHESIS takes 3 minutes to run on this dataset, most of which is in file I/O.

A directory called out_test/ will be created and will contain the results from this run.  A summary of the results is in the file REPORT.txt; the main output files are in the subdirectory main_results/; other intermediate results are in the subdirectory cached_data.  The results from this test case won't be very good because the dataset of Hi-C links is so small, but they should give you an idea of what to expect from LACHESIS.



###############################################################################



RUNNING LACHESIS



1. Input requirements

So you want to use LACHESIS to scaffold your de novo assembly.  You must start with two (or three) input files:

-- Your Hi-C reads, in an alignable format
-- Your draft (de novo) assembly, in fasta format
-- (Optional) A reference assembly, in fasta format

If you are going to use BWA to align your Hi-C reads, then the reads must be in a set of fastq files, and there should be two files for each library: one for the first read in each pair, and one for the second.  If you are using another aligner, you may have your Hi-C reads in another format as required by that aligner.

Obviously you can't supply a reference assembly if no such assembly exists yet.  But if there is one, the reference assembly can be input to LACHESIS in order to provide some very useful reference-based evaluation of the result.  Even a partial sequence assembly can be helpful.



2. Creating an INI file

Before running LACHESIS, you will need to create an initialization file, or INI file (*.ini).  LACHESIS will parse this file and determine everything it needs to know in order to run, including the locations of all other input files.  To create an INI file, make a copy of any of the INI files in the subdirectory INIs/ and edit that copy.  The INI file contains several parameters; you will want to read through its documentation and set all parameters as you see fit.

Five of the INI parameters describe files that LACHESIS uses: SAM_FILES (in SAM_DIR), RE_COUNTS_FILE, BLAST_FILE_HEAD, DRAFT_ASSEMBLY_NAMES_FILE, and REF_ASSEMBLY_NAMES_FILE.  These all refer to files that you will need to create before running LACHESIS, by processing your input datasets.  See the included file "Data Prep Flowchart.png" for a visual guide to this process.



3. Aligning the Hi-C reads to the draft assembly

LACHESIS scaffolds de novo assemblies using the locations of the Hi-C reads on the assembly contigs/scaffolds.  Hence you must align the Hi-C reads to the assembly contigs/scaffolds.  You must ultimately produce files containing paired-read alignments in the SAM/BAM format (either SAMs or BAMs are acceptable; see http://samtools.sourceforge.net/).  For each library of Hi-C reads (represented by a pair of fastq files) you should create a single SAM/BAM file describing the read pairs.  DO NOT create a pair of SAM/BAM files, one for each read in the pair!

Aligning is a computationally intensive process and will take far more CPU time than LACHESIS itself will.  You can use any aligner that produces SAM or BAM files, but we recommend BWA (http://bio-bwa.sourceforge.net/), which we used during the development of LACHESIS.  If you use BWA, you will need to have your Hi-C reads in fastq format as described above, and you will need to perform the following steps:

-- Use `bwa index` to index your draft assembly.  If your genome is large (>10 Mb), you will need to use the flag `-a bwtsw`.  This will create a file "<assembly>.fasta.bwt".
-- Use `bwa aln` to align the Hi-C reads to the indexed assembly.  Each individual fastq file is aligned separately.
-- Use `bwa sampe` to determine the optimal placement of each read pair.  This is where the two fastq files from each library are combined.  The output of bwa sampe will be a SAM file (*.sam).  To save disk space, you can convert this file to a BAM file (*.bam) using samtools (`samtools view -bS <sam-file>`).



4. Filtering the Hi-C reads

Many of the Hi-C reads can be determined to be noise, rather than signal.  It is not strictly necessary to remove these reads from your BAM files before running LACHESIS.  But it will reduce the files' size, reduce the I/O time required by LACHESIS to read them in, and may reduce LACHESIS' error rate.

We have included the scripts PreprocessSAMs.pl and PreprocessSAMs.sh, which perform the filtering methods we used.  PreprocessSAMs.pl removes all reads that do not map within 500bp of a restriction enzyme site, following the suggestion of Yaffe & Tanay (Nature Genetics, 2011).  To do this, it uses the script make_BED_around_RE_site.pl (also included) and also the bedtools library.  PreprocessSAMs.pl also removes unpaired reads, which LACHESIS cannot use.  Make sure to set the variable $RE_site in PreprocessSAMs.pl.  You can also use the script PreprocessSAMs.sh as a batch script to run PreprocessSAMs.pl on a set of SAM/BAM files (set the variables SAMs and ASSEMBLY).  The final set of SAM/BAM files is specified in the INI file as the parameter SAM_FILES, which names the files, and SAM_DIR, which names the directory they are in.



5. Aligning the draft assembly to the reference assembly, if there is one

You can skip this step if you aren't using a reference assembly.  LACHESIS needs to know where the contigs/scaffolds from the draft assembly belong on the reference chromosomes.  It can process the text output from the blastn command-line program, which you can read about and download at http://www.ncbi.nlm.nih.gov/books/NBK1763/.  You will need to perform the following steps:

-- Create a BLAST database for the reference assembly: `makeblastdb -in <ref-fasta> -dbtype=nucl -out <ref-fasta>.blastdb`
-- Run blastn.  Make sure you use the argument 'outfmt -7' so that blastn produces output in the format that LACHESIS expects.
-- If you have a large (mammalian-scale) draft assembly, blastn will take a long time to align it.  You may want to split it up into separate fasta files and run blastn on each file separately.  See the included scripts blast.sh and blast.qsub.sh to see how we did this on the human assembly.  The blastn output files from each chunk should be named <BLAST_FILE_HEAD>.NNN.blast.out, where NNN = 1, 2, 3, ...
-- The blastn output file(s) are specified in the INI file by the parameter BLAST_FILE_HEAD.  Note that LACHESIS parses the BLAST output the first time it runs, then creates a cached file so that it doesn't have to take the time to parse the BLAST output again.  This means that if you change the set of BLAST files, you'll need to delete the file at <OUTPUT_DIR>/cached_data/TrueMapping.assembly.txt in order to get LACHESIS to load your new set of BLAST files.



6. Counting the restriction sites in the draft assembly contigs/scaffolds

To normalize the Hi-C data properly, LACHESIS needs to know how many instances of the RE site are in each contig/scaffold in the draft assembly.  You can create a file for this using the script CountMotifsInFasta.pl, which we have included.  In CountMotifsInFasta.pl, set the variable $motif to the restriction enzyme sequence (e.g., "AAGCTT") and set the variable $fasta to your <draft-assembly>.fasta file.  CountMotifsInFasta.pl will produce a file at <draft-assembly>.counts_<motif>.txt.  This filename is specified in the INI file as the parameter RE_COUNTS_FILE.



7. Extracting the assembly contig names

LACHESIS likes to know the names of the contigs/scaffolds in the draft assembly, and also the reference assembly, if there is one.  The simplest way to do this is with the UNIX command `grep`, as follows:

grep "^\>" <assembly-fasta>  >  <assembly-names-fasta>

The output filename <assembly-names-fasta> is specified in the INI file as DRAFT_ASSEMBLY_NAMES_FILE.  If you are using a reference assembly, you should run the same command on it, and its output filename is specified in the INI file as REF_ASSEMBLY_NAMES_FILE.



8. Running LACHESIS

If you've done steps 2-7 above, and you've created an INI file with the appropriate input file names, you should be good to go.  Try running LACHESIS with the command: `Lachesis <your-INI-file>`.  If there are any problems with the INI file, including missing input files, LACHESIS will immediately abort and will explain what went wrong.  If not, then LACHESIS will begin loading in files, and all you need to do is wait and see what it does.

There are many parameters in the INI file beyond the inputs to LACHESIS.  You may need to tweak some of these parameters in order to get an optimal assembly.  Be sure to read the documentation in the INI file carefully so that you know what you're doing.



9. Interpreting the LACHESIS results

LACHESIS will create a set of output files in the directory OUTPUT_DIR that you have specified in the INI file.  This directory contains a file "REPORT.txt" that will give a topline summary of LACHESIS' performance.  It also contains two subdirectories, main_results/ and cached_data/.  The main_results/ directory will contain the following files:
I. clusters.txt and clusters.by_name.txt):  These files indicate the clustering results.  Each of LACHESIS' chromosome group is shown as a line, and the input contigs/scaffolds in that group are listed on the line, either by ID (clusters.txt) or by contig name (clusters.by_name.txt).
II. group*.ordering: These files indicate the ordering and orienting results.  There is one file for each group.  In each file is a list of input contigs/scaffolds in order, with their orientations, orientation quality scores, and gap sizes.

To create the final assembly fasta, run the included script CreateScaffoldedFasta.pl.  After you have run this, OUTPUT_DIR will contain the file "Lachesis_assembly.fasta".  This fasta file will contain your output assembly, in three sections of successive contigs:
I. One large scaffold for each of the ordered and oriented chromosome groups predicted by LACHESIS (i.e., each of the group*.ordering files.)
II. All of the input contigs/scaffolds that have been clustered into a chromosome group by LACHESIS, but were not ordered within that group.  They will be given names indicating what chromosome group they belong in.
III. All of the input contigs/scaffolds that were not clustered at all by LACHESIS.



###############################################################################



TROUBLESHOOTING



LACHESIS is a good piece of software, but it isn't perfect.  You may run it and get a result you weren't expecting.  You may also run it and get no result at all because it crashes.

If LACHESIS crashes, the first thing you should do is look carefully at its output.  It might give a verbose explanation of what went wrong and give you a good idea for how to fix it.  You may also receive an "assertion error", which looks like this: "Assertion ... failed."  That means that at some stage of the algorithm, LACHESIS encountered something specific that it wasn't expecting.  An assertion error will come with a reference to the file (*.cc or *.h) and the line number where the error occurred.  Try looking at that line in the file, which should contain the function assert().  There should be some comments around that line that explain what might be causing the assertion error.

In general, we've made a strong effort to make LACHESIS a well-designed and well-commented piece of code.  If you're familiar with C++, you should be able to poke around in the source code and get an idea of what's going on.  We recommend starting with the top-level module, Lachesis.cc, and working from there.

If all else fails, and you still need help running LACHESIS, please e-mail Josh Burton at jnburton@uw.edu.



###############################################################################



COPYRIGHT AND DISCLAIMER

The LACHESIS software package and all software and documentation contained with it are copyright © 2012-2013 by Josh Burton and the University of Washington.  All rights are reserved.  This software is supplied 'as is' without any warranty or guarantee of support.  The University of Washington is not responsible for its use, misuse, or functionality.



###############################################################################



ACKNOWLEDGEMENTS

Thanks to Jay Shendure for leadership, management, and many ideas.

Thanks to Aaron McKenna for making LACHESIS available over GitHub.

Thanks to Aaron McKenna, Qi Zhou, and Christopher Beitel for patiently helping me test LACHESIS for bugs and compatibility.



When using Lachesis, please cite Burton et al., Nature Biotechnology 31:1119-1125 (2013).  Link here: dx.doi.org/10.1038/nbt.2727
