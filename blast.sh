#!/bin/bash -f

# Align a set of contigs (a draft assembly in progress) to the human reference genome.
# BLAST (not cross_match, not BWA) is the best aligner for this.
# The output file will be used by AnalyzeHiC for validating the scaffolding of the draft assembly with Hi-C data.
#
#
# Josh Burton
# February 2013

# 
# Use the bash shell to interpret this job script 
#$ -S /bin/bash
# 
# Run in parallel on as many as 16 processors.
#################$ -pe serial 1-16
# 
# Only submit this job to nodes that have at least 8GB of RAM free.
##$ -l mem_requested=8G
#########$ -l mem_requested=8G, h_vmem=8G

# Make a cool header and print it to the SGE stdout and stderr files.
# Having this information in your output file will help track down any errors which might occur during the job's run.
START_TIME=`/bin/date`
echo "********************************************************************"
echo "*"
echo "*    COMMAND = $0 $@"
echo "*    \$PWD = $PWD"
echo "*    \$USER = $USER"
echo "*    \$HOSTNAME = $HOSTNAME"
echo "*    START TIME = $START_TIME"
echo "*"
echo "********************************************************************"
echo

# SGE jobs do not run in your login environment, so you'll need to 
# load your environment, or atleast the modules your script needs to run
source $HOME/.bashrc

# Directories for use in the command
ASSEMBLY_DIR=$HOME/vol10/src/HiC/human/assembly
REFDB=$HOME/vol10/hg19/all/Homo_sapiens_assembly19.fasta.blastdb
#ASSEMBLY_DIR=$HOME/vol10/src/HiC/Chlamydomonas/assembly
#REFDB=$HOME/vol10/src/HiC/Chlamydomonas/ref/Creinhardtii_236.fa.blastdb
CHUNK=$1

# Script or command(s) to run via SGE
tm blastn -query $ASSEMBLY_DIR/assembly.$CHUNK.fasta -db $REFDB -perc_identity 99 -evalue 100 -word_size 50 -out $ASSEMBLY_DIR/assembly.$CHUNK.blast.out -outfmt 7 #-num_threads 16 # TEMP
#blastn -query $ASSEMBLY_DIR/assembly.$CHUNK.fasta -db $REFDB -out $ASSEMBLY_DIR/assembly.$CHUNK.blast.manual.out -outfmt 7 -perc_identity 99 -evalue 100 -word_size 100


# RUNTIMES (human)
# For 190 chunks of 100 each
# Chunk 100: 9m39s
# Chunk 120: 2m52s, 4.5GB
# Chunk 150: 51s, 4GB
# Chunk 188: 1m15s, 4GB (1m59s if -word_size 45)
# Chunk 189: ~30s
