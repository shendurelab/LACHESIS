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
 * SAMStepper
 *
 * The SAMStepper class is a helper class that allows the user to load in a set of SAM/BAM files and then conveniently step through all of the alignments.
 *
 * Usage example:
 *
 * SAMStepper stepper( SAM_files_vector );
 * stepper.FilterAligned(); // call this and other Filter() functions before calling next_read() or next_pair()
 *
 * for ( bam1_t * align = stepper.next_read(); align != NULL; align = stepper.next_read() ) {
 * -- or --
 * for ( pair< bam1_t *, bam1_t *> aligns = stepper.next_pair(); aligns.first != NULL; aligns = stepper.next_pair() ) {
 *
 *   // do something with align
 *
 * }
 *
 *
 *
 *
 * Notes:
 * You can input any combination of SAM or BAM files, as long as they have the suffix "sam" or "bam" (case-insensitive).
 * The SAMStepper will figure out which is which.
 * Note that the alignments returned by next_read() and next_pair are only pointers; the alignment objects themselves are stored in the SAMStepper object.
 * Hence a subsequent call to next_read() or next_pair() will change the pointer, as well as any bam1_core_t &'s (but not bam1_core_t's) pointing inside the
 * alignment.
 *
 *
 * This class was originally a helper class in the AlgorithmOctopus module of SimCancer (July 2011).
 *
 * Josh Burton
 * December 2012
 *
 *************************************************************************************************************************************************************/



#ifndef __SAM_STEPPER_H
#define __SAM_STEPPER_H

#include <vector>
#include <string>
#include <assert.h>
using namespace std;

// To compile using this, you must include -I<samtools dir>
#include <sam.h>


class SAMStepper
{

 public:
  // Main constructor: Load in a set of files.
  SAMStepper( const vector<string> & SAM_files );
  // Constructor: Load in only one file.  This is a wrapper to the above constructor.
  SAMStepper( const string & SAM_file );

  // Destructor.
  ~SAMStepper();

  /* FILTERING FUNCTIONS
   * These filter the set of reads that are reported via next_read() and next_pair().  These functions must be called before calling next().
   * NOTE: These functions sometimes act wonky because the read's FLAG doesn't always match its alignment - e.g., I've seen a lot of reads that are aligned
   * yet have the 0x4 FLAG set - wtf, bwa. */

  // FilterAligned: Only accept reads that are aligned (i.e., FLAG & 0x4 == 0.)
  void FilterAligned();
  // FilterAlignedPairs: Only accept reads in pairs that both aligned (i.e., FLAG & 0x4 == 0 and FLAG & 0x8 == 0).
  void FilterAlignedPairs();
  // FilterRegion: Only accept reads that align to a specific chromosome and/or range.
  void FilterRegion( const int chrom, const int start = -1, const int end = -1 );


  // next_read(): Main function to get an alignment.  Return NULL if there are no more alignments to get.
  bam1_t * next_read();
  // next_pair(): Get a pair of alignments.  Assume paired reads appear in consecutive order in the file, and skip over unpaired reads.
  // Return pair<NULL,NULL> if there are no more paired alignments to get.
  pair< bam1_t *, bam1_t * > next_pair();

  /* QUERY FUNCTIONS */

  int N_targets() const { assert( _SAM != NULL ); return _SAM->header->n_targets; }

  // as_SAM_line: A wrapper to bam_format1, which formats a bam1_t object as a string in the format of a line in a SAM file (with no newline).
  char * as_SAM_line( const bam1_t * align ) const { return bam_format1( _SAM->header, align ); }

  // N_aligns_read: The number of alignments read so far - i.e., the number of calls to next_read() (or 2x calls to next_pair()) that did not return NULL.
  // Alignments that are filtered out, including unpaired alignments skipped by next_pair(), are not counted.
  int64_t N_aligns_read() const { return _N_aligns_read; }
  int64_t N_pairs_read()  const { return _N_pairs_read; }

 private:


  /* PRIVATE FUNCTIONS */

  // Init(): Set up the SAMStepper object.  This is immediately called by any constructor.
  void Init();

  // Open the next SAM/BAM file in the set, via open_SAM().
  samfile_t * open_next();


  /* PRIVATE DATA */

  const vector<string> _SAM_files; // list of SAM filenames

  int _ID; // index (in list) of the currently open SAM/BAM file
  samfile_t * _SAM; // the currently open SAM/BAM file
  bam1_t * _align, * _align2; // the alignment object(s) to return - the actual objects are stored here!
  bool _init; // has Init() been called yet?
  int64_t _N_aligns_read; // number of alignments returned so far via next_read() *or* next_pair().  This does NOT include alignments filtered out.
  int64_t _N_pairs_read; // number of pairs of alignments returned so far via next_pair().

  // Optional filtering variables.  These can be set by the Filter() functions, and they control which alignments are returned by next().
  bool _filter_aligned;
  bool _filter_aligned_pair;
  int _filter_chrom, _filter_start, _filter_end;

  static const bool _verbose = false;
};




// Wrapper to samopen() which figures out the open mode (SAM vs. BAM.)
// This function is separate from the SAMStepper class and is designed to be usable outside it.
// To avoid memory leaks, be sure to eventually call samclose() on all pointers returned from open_SAM().
samfile_t * open_SAM( const string & SAM_file );

// The following 5 functions all use open_SAM.

// NTargetsInSAM(): Return the number of target sequences in this SAM file.
int NTargetsInSAM( const string & SAM_file );
// TargetNames(): Return a vector of the names of the target sequences in this SAM file.
vector<string> TargetNames    ( const string & SAM_file );
// TargetLengths(): Return a vector of the lengths of the target sequences in this SAM file.
vector<int>  TargetLengths    ( const string & SAM_file );
// TargetNHits(): Return a vector of the number of query sequences aligning to each target in this SAM file (note: this counts ALL reads, without filtering.)
// This function is time-consuming.
vector<int>    TargetNHits    ( const string & SAM_file );
// TargetCoverages(): Return a vector of the coverage of target sequences by query sequence hits.  Calls TargetLengths() and TargetNHits().
// Note that the "coverage" here is defined as the NUMBER of query sequences aligning to each target, per target bp.  Read lengths are NOT taken into account.
vector<double> TargetCoverages( const string & SAM_file );


#endif
