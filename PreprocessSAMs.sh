#!/usr/bin/env bash
##///////////////////////////////////////////////////////////////////////////////
##//                                                                           //
##// This software and its documentation are copyright (c) 2014-2015 by Joshua //
##// N. Burton and the University of Washington.  All rights are reserved.     //
##//                                                                           //
##// THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  //
##// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF                //
##// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  //
##// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY      //
##// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT //
##// OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR  //
##// THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                //
##//                                                                           //
##///////////////////////////////////////////////////////////////////////////////

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

## One should never need to do this in a script.
##source $HOME/.bashrc # make bwa, samtools, picard available

### STUFF FOR QSUBBING - end

##DIR=$HOME/vol10/src/Lachesis # This directory should contain PreprocessSAMs.pl.
DIR=$(dirname $(which Lachesis))
# SAMs: A set of SAM/BAM files (both are allowable, as long as the extensions accurately describe the files).  Path is relative to $DIR.
SAMs=${0:-"human/to_postfosmid/SRR400260.bam human/to_postfosmid/SRR400261.bam human/to_postfosmid/SRR400262.bam human/to_postfosmid/SRR400263.bam human/to_postfosmid/SRR442155.bam human/to_postfosmid/SRR442156.bam human/to_postfosmid/SRR442157.bam"}
ASSEMBLY=${1:-"human/postfosmid/assembly.fasta"} # the fasta file representing the draft assembly to which the SAMs are aligned.  Pathsis relative to $DIR.

echo "SAMs = $SAMs"

for SAM in $SAMs ; do
    $DIR/PreprocessSAMs.pl $DIR/$SAM $DIR/$ASSEMBLY &
done
