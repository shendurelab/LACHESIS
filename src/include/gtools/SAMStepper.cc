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


// For documentation, see SAMStepper.h
#include "SAMStepper.h"

#include <string.h>
#include <assert.h>
#include <vector>
#include <string>
#include <iostream>
#include "../TimeMem.h"

// Boost libraries
#include <boost/algorithm/string.hpp> // to_upper
#include <boost/filesystem.hpp>

// To compile using this, you must include -I<samtools dir>
#include <sam.h>


// Constructor for multiple files.
SAMStepper::SAMStepper( const vector<string> & SAM_files )
  : _SAM_files( SAM_files ),
    _init( false )
{
  Init();
}


// Constructor for one file.
SAMStepper::SAMStepper( const string & SAM_file )
  : _SAM_files( vector<string>( 1, SAM_file ) ),
    _init( false )
{
  Init();
}




// Destructor: cleanup.
SAMStepper::~SAMStepper()
{
  if ( _N_aligns_read > 0 ) free(_align ->data);
  if ( _N_pairs_read  > 0 ) free(_align2->data);
  delete _align;
  delete _align2;
  samclose(_SAM);
}



// Init(): Set up the SAMStepper object.  This is immediately called by any constructor.
void
SAMStepper::Init()
{
  assert( !_init ); // this will fail if Init() has already been called.

  // Check that all SAM files exist.
  assert( !_SAM_files.empty() );
  for ( size_t i = 0; i < _SAM_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( _SAM_files[i] ) );

  // Initialize local variables.
  _ID = 0;
  _align  = new bam1_t();
  _align2 = new bam1_t();
  _init = true;
  _N_aligns_read = 0;
  _N_pairs_read  = 0;

  // All filters are off by default.
  _filter_aligned = false;
  _filter_aligned_pair = false;
  _filter_chrom = -1;
  _filter_start = -1;
  _filter_end = -1;

  // Open the first SAM/BAM file.  (Don't yet read any alignments from it!)
  assert( !_SAM_files.empty() );
  _SAM = open_next();
}






void
SAMStepper::FilterAligned()
{
  assert( _N_aligns_read == 0 ); // can't turn on a filter after beginning to go through alignments
  _filter_aligned = true;
}


void
SAMStepper::FilterAlignedPairs()
{
  assert( _N_aligns_read == 0 ); // can't turn on a filter after beginning to go through alignments
  _filter_aligned = true;
  _filter_aligned_pair = true;
}


void
SAMStepper::FilterRegion( const int chrom, const int start, const int end )
{
  assert( _N_aligns_read == 0 && _N_pairs_read == 0 ); // can't turn on a filter after beginning to go through alignments

  _filter_aligned = true; // a read can't align to a specific region if it doesn't align at all!
  _filter_chrom = chrom;
  _filter_start = start;
  _filter_end = end;
}




// Open the next SAM/BAM file in the set, via open_SAM().  Do not actually read any alignments from it.
samfile_t *
SAMStepper::open_next()
{
  if ( _verbose ) cout << Time() << ": SAMStepper is opening file: " << _SAM_files[_ID] << endl;
  return open_SAM( _SAM_files[_ID] );
}




// next_read(): Main function to get an alignment.  Fills the variable _align and also returns it.  Return NULL if there are no more alignments to get.
bam1_t *
SAMStepper::next_read()
{
  // Get the next alignment from the currently open SAM/BAM file.
  while ( samread( _SAM, _align ) != -1 ) {
    const bam1_core_t & c = _align->core;

    // If alignment and/or range filters have been set, skip over any read that doesn't meet them.
    bool aligned = !( c.flag & 0x4 );
    if ( _filter_aligned && !aligned ) continue;

    if ( _filter_aligned_pair && ( c.flag & 0x8 ) ) continue;

    if ( _filter_chrom != -1 && _filter_chrom != c.tid ) continue;
    if ( ( _filter_start != -1 && _filter_end != -1 ) && ( c.pos < _filter_start || c.pos > _filter_end ) ) continue;

    // Sanity check.  Aligned reads should satisfy all of these.
    if ( aligned ) {
      assert( c.tid != -1 ); // note that c.pos == -1 can occur because for some reason c.pos is 1 less than the number in the POS column of the samfile
      assert( c.n_cigar > 0 );
    }

    _N_aligns_read++;
    return _align;
  }


  // If control reaches here, we've exhausted this SAM/BAM file.
  // If there are no more files, we're done: return NULL.
  if ( _ID + 1 == (int)_SAM_files.size() ) return NULL;

  // If there's another SAM/BAM file to open, open it.
  _ID++;
  samclose(_SAM);
  _SAM = open_next();

  // Now recurse on next_read() so that we can actually find and return an alignment.
  return next_read();
}




// next_pair(): Get a pair of alignments.  Assume paired reads appear in consecutive order in the file, and skip over unpaired reads.
// Return pair<NULL,NULL> if there are no more paired alignments to get.
pair< bam1_t *, bam1_t * >
SAMStepper::next_pair()
{
  pair< bam1_t *, bam1_t * > null = make_pair( (bam1_t *) NULL, (bam1_t *) NULL );

  // Call next_read() twice to get two new reads.  Each call to next_read() fills the local variable _align.
  if ( next_read() == NULL ) return null;
  bam_copy1( _align2, _align );
  if ( next_read() == NULL ) return null;

  // If necessary, call next_read() repeatedly until two consecutive reads have the same fragment name (or one is a substring of the other).;
  // This can happen if one read passes the filters but the others does not, including times when the two reads have inconsistent 0x4 and 0x8 FLAGs
  // (which oddly happens sometimes with MT-aligned reads in this one SAM file.)
  while ( 1 ) {

    // Find the first non-matching character between the names.
    // If it's punctuation (i.e., not alphanumeric), or if we reach the end of the strings, the names are considered to match.
    //cout << "NAMES: " << bam1_qname(_align) << "\t" << bam1_qname(_align2) << endl;
    const char * name1 = bam1_qname(_align);
    const char * name2 = bam1_qname(_align2);
    int len1 = strlen(name1), len2 = strlen(name2);
    //cout << "LENGTHS = " << len1 << '\t' << len2 << '\t' << name1[len2] << '\t' << boolalpha << isalpha(name1[len2]) << endl;
    if ( strncmp( name1, name2, min( len1, len2 ) ) == 0 ) {
      if ( len1 == len2 ) break;
      else if ( len1 > len2 ) { if ( !isalnum( name1[len2] ) ) break; }
      else if ( len2 > len1 ) { if ( !isalnum( name2[len1] ) ) break; }
    }

    // No match found.  Move on.
    _N_aligns_read--; // offset the increase due to the extra next_read() call
    bam_copy1( _align2, _align );
    if ( next_read() == NULL ) return null;
  }


  // Return _align2 before _align because it came first in the file.
  _N_pairs_read++;
  return make_pair( _align2, _align );
}





// open_SAM(): A wrapper to samopen() which figures out the open mode (SAM vs. BAM.)
samfile_t *
open_SAM( const string & file )
{
  const string suffix = boost::to_upper_copy( file.substr( file.size() - 3 ) );
  assert( suffix == "BAM" || suffix == "SAM" );
  const string mode = ( suffix == "BAM" ? "rb" : "r" );
  return samopen( file.c_str(), mode.c_str(), 0 );
}



// NTargetsInSAM(): Return the number of target sequences in this SAM file.
int
NTargetsInSAM( const string & SAM_file )
{
  if ( !boost::filesystem::is_regular_file( SAM_file ) ) {
    cout << "ERROR: SAMStepper::NTargetsInSAM: Can't find file `" << SAM_file << "'" << endl;
    assert(0);
  }

  samfile_t * sam = open_SAM( SAM_file );
  int N_targets = sam->header->n_targets;
  samclose(sam);

  return N_targets;
}




// TargetNames(): Return a vector of the names of the target sequences in this SAM file.
vector<string>
TargetNames( const string & SAM_file )
{
  if ( !boost::filesystem::is_regular_file( SAM_file ) ) {
    cout << "ERROR: SAMStepper::TargetNames: Can't find file `" << SAM_file << "'" << endl;
    assert(0);
  }

  vector<string> target_names;

  samfile_t * sam = open_SAM( SAM_file );
  for ( int i = 0; i < sam->header->n_targets; i++ )
    target_names.push_back( sam->header->target_name[i] );
  samclose(sam);

  return target_names;
}


// TargetLengths(): Return a vector of the lengths of the target sequences in this SAM file.
vector<int>
TargetLengths( const string & SAM_file )
{
  if ( !boost::filesystem::is_regular_file( SAM_file ) ) {
    cout << "ERROR: SAMStepper::TargetLengths: Can't find file `" << SAM_file << "'" << endl;
    assert(0);
  }

  vector<int> target_lengths;

  samfile_t * sam = open_SAM( SAM_file );
  for ( int i = 0; i < sam->header->n_targets; i++ )
    target_lengths.push_back( sam->header->target_len[i] );
  samclose(sam);

  return target_lengths;
}



// TargetNHits(): Return a vector of the number of query sequences aligning to each target in this SAM file (note: this counts ALL reads, without filtering.)
// This function is time-consuming.
vector<int>
TargetNHits( const string & SAM_file )
{
  if ( !boost::filesystem::is_regular_file( SAM_file ) ) {
    cout << "ERROR: SAMStepper::TargetNHits: Can't find file `" << SAM_file << "'" << endl;
    assert(0);
  }

  // Set up a SAMStepper object to read in the reads.
  SAMStepper stepper( SAM_file );

  vector<int> target_N_hits( stepper.N_targets(), 0 );

  // For each read that aligned, look at what target sequence it aligned to, and tally the number of hits for that target.
  // ALT METHOD: Call FilterAligned() on the SAMStepper, so no need to check the tids.  This is faster but buggy because the 0x4 FLAG doesn't always match.
  for ( bam1_t * align = stepper.next_read(); align != NULL; align = stepper.next_read() ) {

    const int & tid = align->core.tid;
    if ( tid != -1 ) target_N_hits[tid]++; // if tid == -1, the read didn't align
  }

  return target_N_hits;
}



// TargetCoverages(): Return a vector of the coverage of target sequences by query sequence hits.  Calls TargetLengths() and TargetNHits().
// Note that the "coverage" here is defined as the NUMBER of query sequences aligning to each target, per target bp.  Read lengths are NOT taken into account.
vector<double>
TargetCoverages( const string & SAM_file )
{
  if ( !boost::filesystem::is_regular_file( SAM_file ) ) {
    cout << "ERROR: SAMStepper::TargetCoverages: Can't find file `" << SAM_file << "'" << endl;
    assert(0);
  }

  vector<int> target_N_hits = TargetNHits( SAM_file );
  vector<int> target_lens = TargetLengths( SAM_file );

  assert( target_N_hits.size() == target_lens.size() );

  // Calcualate each contig's coverage.
  vector<double> target_covs( target_N_hits.size(), 0 );
  for ( size_t i = 0; i < target_N_hits.size(); i++ )
    target_covs[i] = double( target_N_hits[i] ) / target_lens[i];

  return target_covs;
}
