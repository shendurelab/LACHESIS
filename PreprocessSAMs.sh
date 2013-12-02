#!/bin/bash -f


# PreprocessSAMs.sh
#
# This is a wrapper script for PreprocessSAMs.pl, which pre-processes SAM/BAM files so they can be used with Lachesis.
# Modify the variables "DIR", "SAMs", and "ASSEMBLY" and below, and you're good to go.
# You can submit this to a cluster via qsub.
# See PreprocessSAMs.pl for more documentation on what the pre-processing entails.
#
# Josh Burton
# July 2013



### STUFF FOR QSUBBING

# Use the bash shell to interpret this job script 
#$ -S /bin/bash
# Reserve space for this job.  This helps it compete against, for example, a barrage of tiny jobs that require 1 CPU each.
#$ -R y
# Only submit this job to nodes that have at least this much RAM free.
# Both of these parameters are needed to get the JVM (which Picard uses) to run.
##$ -l h_vmem=8G, mem_requested=8G

source $HOME/.bashrc # make bwa, samtools, picard available

### STUFF FOR QSUBBING - end

DIR=$HOME/vol10/src/Lachesis # This directory should contain PreprocessSAMs.pl.
# SAMs: A set of SAM/BAM files (both are allowable, as long as the extensions accurately describe the files).  Path is relative to $DIR.
SAMs="human/to_postfosmid/SRR400260.bam human/to_postfosmid/SRR400261.bam human/to_postfosmid/SRR400262.bam human/to_postfosmid/SRR400263.bam human/to_postfosmid/SRR442155.bam human/to_postfosmid/SRR442156.bam human/to_postfosmid/SRR442157.bam"
ASSEMBLY=human/postfosmid/assembly.fasta # the fasta file representing the draft assembly to which the SAMs are aligned.  Pathsis relative to $DIR.


echo "SAMs = $SAMs"

for SAM in $SAMs ; do
    $DIR/PreprocessSAMs.pl $DIR/$SAM $DIR/$ASSEMBLY &
done
