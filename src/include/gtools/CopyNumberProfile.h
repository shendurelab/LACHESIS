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


/******************************************************************************
 *
 * CopyNumberProfile
 *
 * A simple data structure representing a set of copy-number calls across
 * the human genome.  A CopyNumberProfile is loaded from a BEDgraph file.
 *
 *
 * The main query function is CN(), which can be called on either a single
 * position or a chrom_interval, and returns a called copy number.  If you have
 * a sorted set of positions or chrom_intervals, using the CopyNumberProfile
 * isn't as quite as fast as the stepping-through-both-sorted-sets-at-once
 * method (as demonstrated in HapSplitter::CallHapResolvedCopyNumbers) because
 * each individual call requires CopyNumberProfile to look through the whole
 * set of copy-number calls.
 *
 *
 * The CopyNumberProfile can also store HRCNs (Haplotype-Resolved Copy Numbers)
 * in the same fashion as CNs.  To load HRCNs into a CopyNumberProfile object,
 * call AddHRCN().  To query the HRCN of an interval, call HRCN().
 * Note that the HRCNs are stored in a distinct set of chrom_intervals from the
 * CNs, which is necessarily more restrictive.  Hence a position or interval
 * may have a defined CN() but not HRCN().
 *
 *
 *
 * Josh Burton
 * August 2012
 *
 *****************************************************************************/


#ifndef _COPY_NUMBER_PROFILE__H
#define _COPY_NUMBER_PROFILE__H



#include "ChromInterval.h"

#include <string>
#include <map>
#include <vector>
using namespace std;


class CopyNumberProfile
{
 public:

  /* CONSTRUCTORS */

  CopyNumberProfile();
  // Load data in from a BEDgraph file.
  // TODO: Make this constructor load HRCNs as well, if the BED file has 2 extra fields.
  CopyNumberProfile( const string & BEDgraph_file, const int gap_size = 0 );

  /* FUNCTIONS TO ADD CN/HRCN DATA */

  // You should call SetGapSize before calling AddCN, AddHRCN, or AddHRCNData.
  void SetGapSize( const int s ) { _gap_size = s; }

  // AddCN, AddHRCN: Add a CN or HRCN interval.
  // Merge with existing intervals if possible (depending on gap_size.)
  // Note that this merging will only work properly if these functions are
  // called in ascending order of chrom_intervals.
  // A more general solution is easy to write, but would take longer to run.
  void AddCN  ( const chrom_interval & interval, const int CN );
  void AddHRCN( const chrom_interval & interval, const pair<int,int> & HRCN );
  void AddHRCN( const chrom_interval & interval, const int major_CN, const int minor_CN );

  // Load data from a set of HRCN calls.  This also creates CN calls.
  void AddHRCNData( const map< chrom_interval, pair<int,int> > & HRCN_calls );

  /* MAIN QUERY FUNCTIONS */

  // CN: Copy number.  A return value of -1 indicates no called data: e.g.,
  // position/interval is inside a centromere or between windows, or interval
  // is not wholly inside a single copy-number call.
  int CN( const string & chrom, const int pos ) const;
  int CN( const chrom_interval & query ) const;
  // CN_region: Return a contiguous region of common copy number.  If there
  // isn't one, return an empty chrom_interval().
  chrom_interval CN_region( const string & chrom, const int pos ) const;
  chrom_interval CN_region( const chrom_interval & query ) const;
  // CN_regions: Return all CN regions that overlap this interval (though not their CNs.)
  vector<chrom_interval> CN_regions( const chrom_interval & query ) const;

  // HRCN: Haplotype-resolved copy number.
  // A return value of <-1,-1> indicates no called data.
  pair<int,int> HRCN( const string & chrom, const int pos ) const;
  pair<int,int> HRCN( const chrom_interval & query ) const;


  // Output functions.
  void PrintCNStats  ( ostream & out = cout ) const;
  void PrintHRCNStats( ostream & out = cout ) const;
  void PrintHRCNs    ( ostream & out = cout ) const;
  void WriteHRCNs( const string & BEDfile, const bool append = false ) const;


 private:

  // Helper functions for CN() and HRCN().
  pair<chrom_interval,int> CN_and_interval( const string & chrom, const int pos ) const;
  pair<chrom_interval, pair<int,int> > HRCN_and_interval( const string & chrom, const int pos ) const;


  // Copy-number calls, and haplotype-resolved copy-number (HRCN) calls.
  // These vectors are indexed by chromosome.  The chrom_intervals on which CNs
  // are called are guaranteed not to overlap with one another.
  vector< map< chrom_interval, int > > _CNs;
  vector< map< chrom_interval, pair<int,int> > > _HRCNs;

  // Acceptable size of gaps between CNs and HRCNs.
  // Gaps of size <= _gap_size will be merged by AddCN and AddHRCN.
  int _gap_size;

};





#endif
