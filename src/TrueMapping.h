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


/**************************************************************************************************************************************************************
 *
 * TrueMapping.h
 *
 * A TrueMapping is a list of all the contigs in a draft assembly, with an indication of where they map onto a reference assembly.  If you are assembling a
 * genome for which a reference assembly already exists, you'll want to set up a TrueMapping object so you can evaluate Lachesis' results.  On the other hand,
 * if you're assembling a genome truly de novo, you'll never need to instantiate a TrueMapping.
 *
 * The motivation for a TrueMapping is as follows.  The Lachesis algorithm, as implemented in the GenomeLinkMatrix and ChromLinkMatrix classes, is designed
 * to take a set of contigs and scaffold them together into one or more scaffolds using Hi-C links that have been mapped onto them.  The best way to evaluate
 * the success of Lachesis is to scaffold together a draft assembly for an organism that also has a completed reference assembly, and then compare the
 * Lachesis-created scaffold with the "true" mapping of the draft assembly's contigs.  The true mapping is stored in the TrueMapping class, which is simply a
 * data structure designed for this purpose.
 *
 * The input to the TrueMapping class is a text file of BLAST output representing the alignment of the draft assembly to the true reference.  To create this
 * file, run blastn with -outfmt 7.  It may be a computational hassle to get blastn to run with a reasonable amount of time and memory usage.  To speed it up,
 * restrict the space of allowable alignments (here's a set of options that was used successfully on human: -perc_identity 99 -evalue 100 -word_size 50).
 *
 * Alternatively, a TrueMapping class can be created for a non-de novo assembly that consists of the reference genome split up into "bins" of equal length.
 *
 * Note that the contigs in the draft assembly may technically be scaffolds themselves, if they are given as contigs with long strings of interior N's to
 * represent gaps.  This is how ALLPATHS reports its draft scaffolds.  They are still treated as contigs by Lachesis.
 *
 * Also note that the TrueMapping currently assumes only one "true" location of each draft contig on the reference.  Alignments in BLAST output are sorted by
 * bit score, and the TrueMapping just takes the first-listed (highest bit score) alignment and treats it as gospel.  This may cause misleading results,
 * especially for small contigs that may lie in repetitive regions.
 *
 *
 * Josh Burton
 * February 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _TRUE_MAPPING__H
#define _TRUE_MAPPING__H

#include <assert.h>
#include <string>
#include <vector>
#include <iostream>
#include <map>
using namespace std;




struct BlastAlignmentVec; // defined in the .cc





class TrueMapping
{
 public:

  /* CONSTRUCTORS */

  // Default constructor.  This is only invoked when another object that contains a TrueMapping object (e.g., Reporter) is instantiated.
  TrueMapping() {}
  // Load in a set of BLAST alignments for a de novo GLM.  Also load in query and target names.  The dummy_SAM_file is for getting contig lengths.
  // Files should exist: <BLAST_file_head>.blast.out, <BLAST_file_head>.*.blast.out (for * = 1,2,3,...)
  TrueMapping( const string & species, const vector<string> & query_names, const vector<string> & target_names, const string & BLAST_file_head, const string & out_dir, const string & dummy_SAM_file );
  // Create a TrueMapping for a de novo GLM made by chopping up a reference genome into bins of size BIN_SIZE.  The species must be human (for now).
  TrueMapping( const string & species, const int BIN_SIZE, const vector<string> & query_names, const vector<string> & target_names );
  // Create a TrueMapping for a non-de novo GLM.
  TrueMapping( const string & species, const int & bin_size, const vector<string> & chrom_names, const map<string,int> & chrom_lengths );

  /* QUERY FUNCTIONS */
  string species() const { return _species; }
  int NQueries() const { return _query_names.size(); }
  int NTargets() const { return _target_names.size(); }
  string TargetName ( const int target_ID ) const { return _target_names.at(target_ID); }
  int TargetID( const string & target_name ) const;
  // Functions that input a query ID start with the letter "Q".
  bool   QMaps       ( const int query_ID ) const { return _target.at(query_ID) != -1; } // return true iff this query maps at all
  int    QTargetID   ( const int query_ID ) const { assert( query_ID < NQueries() ); return _target.at(query_ID); } // return -1 if !QMaps()
  string QTargetName ( const int query_ID ) const { if ( QTargetID(query_ID) == -1 ) return ""; return TargetName( QTargetID(query_ID) ); }
  int    QTargetStart( const int query_ID ) const { return min( _start.at(query_ID), _stop. at(query_ID) ); } // Start() <= Stop() always
  int    QTargetStop ( const int query_ID ) const { return max( _stop .at(query_ID), _start.at(query_ID) ); }
  bool   QRC         ( const int query_ID ) const { return _start.at(query_ID) > _stop.at(query_ID); } // false if !QMaps()
  // Quality scores.  The extra "Q" stands for "Quality".
  double QQAlignability( const int query_ID ) const { return _qual_alignability.at(query_ID); }
  double QQSpecificity ( const int query_ID ) const { return _qual_specificity .at(query_ID); }
  bool   QHighQuality  ( const int query_ID ) const { return QQAlignability(query_ID) >= 0.5 && QQSpecificity(query_ID) >= 0.9; } // HEUR: quality threshold
  // Print: Print human-readable output of all alignments.  If genome_order == true, sort output lines in reference order; otherwise, sort by query ID.
  void Print( const bool genome_order ) const;
  // PrintSeqLengthOnTargets: Print a chart report about the total length of query sequence aligned to each target sequence.
  void PrintSeqLengthOnTargets( const string & dummy_SAM_file, ostream & out = cout ) const;

  /* Query mapping functions */

  // Mapping of contig ID to chromosome ID in the fasta order
  vector<int> QueriesToChromIDs() const;
  // Mapping of contig ID to genome ordering (e.g., the first contig on the first chromosome in the fasta maps to 0, etc.)
  // If you call this, consider calling ReorderQueries() afterward.
  vector<int> QueriesToGenomeOrder() const;
  // Find the first contig on a chromosome.  Useful for non-de novo assemblies.
  int FirstContigOnChrom( const int chrom_ID ) const;

  /* Modification functions */

  // Remove a target sequence from consideration by marking all query sequences that map to it as unmapped.  (Don't remove its name from the lists though.)
  // This function is useful if a reference fasta contains unannotated or unplaced scaffolds (e.g., the "het" and "U" reference contigs in Drosophila.)
  void RemoveTarget( const int target_ID );
  void RemoveTarget( const string & target_name ) { RemoveTarget( TargetID( target_name ) ); }
  // Merge two target sequences - that is, mark them as both the same sequence in actuality.  Assign a new name to the combined sequence.
  // This function is useful if a reference fasta contains two arms of a chromosome as separate target sequences (e.g., 2L/2R and 3L/3R in Drosophila.)
  void MergeTargets( const int target_ID_1, const int target_ID_2, const string & merged_name );
  void MergeTargets( const string & target_name_1, const string & target_name_2, const string & merged_name )
  { MergeTargets( TargetID( target_name_1 ), TargetID( target_name_2 ), merged_name ); }



  // Reorder the query sequences.  This accepts output from QueriesToHumanGenome().
  void ReorderQueries( const vector<int> & new_ordering );


 private:

  // Helper function for the TrueMapping constructor.
  // Read alignment info from the BLAST files and write them in a simple format to TrueMapping_file.  If TrueMapping_file already exists, just read from it
  // directly, to save runtime.  Either way, load the alignment data into this TrueMapping object.
  void ReadBlastAlignsFromFileSet( const string & species, const string & dummy_SAM_file, const vector<string> & BLAST_files, const string & TrueMapping_file );

  /* The following three functions assume that species() == "human" and that the chromosomes are named in accordance with the standard in HumanGenome.h. */

  // Return a mapping of contig ID to chromosome ID in the human genome.
  vector<int> QueriesToHumanChromIDs() const;
  // Return a mapping of contig ID to overall ordering in the human genome (e.g., the first contig on chr1 maps to 0, and so on.)
  // If you call this, consider calling ReorderQueries() afterward.
  vector<int> QueriesToHumanGenome() const;
  // Make a mapping of *target* ID to chromosome ID in the human genome.
  vector<int> TargetToHumanChromIDs() const;



  /* DATA */
  string _species;


  // Lists of contig names.
  vector<string> _query_names; // query names (draft assembly)
  vector<string> _target_names; // target names (reference assembly)

  // The alignments themselves.  Query and target contigs here are indicated by an integer, which is their index in the above name lists.
  vector<int> _target; // query contig -> target contig
  vector<int> _start, _stop; // query contig -> alignment start/stop on target (iff start > stop, then the alignment is RC)

  // Alignment quality scores, as calculated in TabulateAlignsToTarget.
  vector<double> _qual_alignability, _qual_specificity;

};




#endif
