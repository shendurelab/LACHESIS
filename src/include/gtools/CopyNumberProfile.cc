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


// For documentation, see CopyNumberProfile.cc.
#include "CopyNumberProfile.h"

#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>


// Modules in ~/include/gtools.
#include "HumanGenome.h"
#include "ChromInterval.h" // chrom_interval
#include "FileParsers.h" // ParseBEDgraph

// Boost libraries.
#include <boost/algorithm/string.hpp> // split
#include <boost/lexical_cast.hpp>




const vector<string>  chrID_to_chrom = HumanGenome_chroms();
const map<string,int> chrom_to_chrID = HumanGenome_chrom_IDs();



CopyNumberProfile::CopyNumberProfile()
  : _gap_size( 0 )
{
  _CNs.resize  ( HumanGenome_n_chroms );
  _HRCNs.resize( HumanGenome_n_chroms );
}



// Load data in from a BEDgraph file.
// Files written with CopyNumberProfile::Write() can be loaded here (CNs only!)
CopyNumberProfile::CopyNumberProfile( const string & BEDgraph_file, const int gap_size )
  : _gap_size( gap_size )
{
  _CNs.resize  ( HumanGenome_n_chroms );
  _HRCNs.resize( HumanGenome_n_chroms );

  // First, parse the BEDgraph and get a set of copy-number calls in intervals.
  vector< pair<chrom_interval,double> > CN_calls = ParseBEDgraph( BEDgraph_file );

  // Now, put these calls into the local data structure.
  // The AddCN function merges the intervals if necessary, if they're in order.
  for ( size_t i = 0; i < CN_calls.size(); i++ )
    AddCN( CN_calls[i].first, lround( CN_calls[i].second ) );
}




/************************************
 *                                  *
 *     FUNCTIONS TO ADD IN DATA     *
 *                                  *
 ************************************/



// Add a CN interval.  Merge it with existing intervals if necessary.
// This code is copied from AddHRCN, below.
void
CopyNumberProfile::AddCN( const chrom_interval & interval, const int CN )
{
  int chrID = interval.chrID;
  assert( chrID != -1 );

  // Consider the already-existing HRCN intervals on this chromosome.
  // Determine whether this HRCN interval can be merged into the existing final
  // interval (e.g. it's close to the existing interval and has the same HRCN.)
  // This logic assumes AddHRCN() will be called in increasing order,
  // and hence we only need to consider the final chrom_interval in the map.
  bool merged = false;
  if ( !_CNs[chrID].empty() ) {
    map<chrom_interval,int>::iterator it = _CNs[chrID].end();
    it--;

    assert( interval > it->first ); // this assert will fail if chrom_intervals aren't added in increasing order; see comment above

    if ( interval.distance( it->first ) <= _gap_size && CN == it->second ) {

      // Merge the new region into the existing final region.
      merged = true;
      chrom_interval merger = it->first;
      _CNs[chrID].erase( it );
      merger = merger.union_with_gap( interval );
      _CNs[chrID].insert( make_pair( merger, CN ) );

    }
  }


  // If this interval wasn't merged in, add it to the _CNs data structure.
  if ( !merged ) _CNs[chrID].insert( make_pair( interval, CN ) );
}



// Wrapper to the other version of AddHRCN.
void
CopyNumberProfile::AddHRCN( const chrom_interval & interval, const pair<int,int> & HRCN )
{
  AddHRCN( interval, HRCN.first, HRCN.second );
}



// Add an HRCN interval.  Merge it with existing intervals if necessary.
void
CopyNumberProfile::AddHRCN( const chrom_interval & interval, const int major_CN, const int minor_CN )
{
  int chrID = interval.chrID;
  assert( chrID != -1 );
  assert( major_CN >= minor_CN );
  pair<int,int> HRCN = make_pair( major_CN, minor_CN );

  // Verify that this region's HRCN is consistent with its already-recorded CN.
  assert( CN( interval ) == major_CN + minor_CN );


  // Consider the already-existing HRCN intervals on this chromosome.
  // Determine whether this HRCN interval can be merged into the existing final
  // interval (e.g. it's close to the existing interval and has the same HRCN.)
  // This logic assumes AddHRCN() will be called in increasing order,
  // and hence we only need to consider the final chrom_interval in the map.
  bool merged = false;
  if ( !_HRCNs[chrID].empty() ) {
    map< chrom_interval, pair<int,int> >::iterator it = _HRCNs[chrID].end();
    it--;

    assert( interval > it->first ); // this assert will fail if chrom_intervals aren't added in increasing order; see comment above

    if ( interval.distance( it->first ) <= _gap_size && HRCN == it->second ) {

      // Merge the new region into the existing final region.
      merged = true;
      chrom_interval merger = it->first;
      _HRCNs[chrID].erase( it );
      merger = merger.union_with_gap( interval );
      _HRCNs[chrID].insert( make_pair( merger, HRCN ) );

    }
  }


  // If this interval wasn't merged in, add it to the _HRCNs data structure.
  if ( !merged ) _HRCNs[chrID].insert( make_pair( interval, HRCN ) );
}





// Load data from a set of HRCN calls.  This also creates CN calls.
// In order to merge CNs/HRCNs properly, you should call SetGapSize() first.
void
CopyNumberProfile::AddHRCNData( const map< chrom_interval, pair<int,int> > & HRCN_calls )
{
  // Loop over the input chrom_intervals.  Note that they are in sorted order
  // in the map<> and thus are suitable for calling AddHRCN().
  map< chrom_interval, pair<int,int> >::const_iterator it;
  for ( it = HRCN_calls.begin(); it != HRCN_calls.end(); ++it ) {
    AddCN  ( it->first, it->second.first + it->second.second );
    AddHRCN( it->first, it->second );
  }
}




/************************************
 *                                  *
 *       MAIN QUERY FUNCTIONS       *
 *                                  *
 ************************************/




int
CopyNumberProfile::CN( const string & chrom, const int pos ) const
{
  return CN_and_interval( chrom, pos ).second;
}



int
CopyNumberProfile::CN( const chrom_interval & interval ) const
{
  if ( interval.empty() ) return -1;
  // Only return a non-null answer if there's a single chrom_interval that
  // contains this entire chrom_interval.
  pair<chrom_interval,int> answer = CN_and_interval( interval.chrom(), interval.start );
  if ( answer.first.contains( interval ) ) return answer.second;
  else return -1;
}


chrom_interval
CopyNumberProfile::CN_region( const string & chrom, const int pos ) const
{
  return CN_and_interval( chrom, pos ).first;
}


chrom_interval
CopyNumberProfile::CN_region( const chrom_interval & query ) const
{
  if ( query.empty() ) return chrom_interval();
  // Only return a non-null answer if there's a single chrom_interval that
  // contains this entire chrom_interval.
  pair<chrom_interval,int> answer = CN_and_interval( query.chrom(), query.start );
  if ( answer.first.contains( query ) ) return answer.first;
  else return chrom_interval();
}




// CN_regions: Return all CN regions that overlap this interval (though not their CNs.)
vector<chrom_interval>
CopyNumberProfile::CN_regions( const chrom_interval & query ) const
{
  vector<chrom_interval> regions;
  if ( query.empty() ) return regions; // empty in -> empty out

  // Get the set of copy-number calls (i.e., chrom_intervals) on this chrom.
  const map<chrom_interval,int> & intervals = _CNs[ query.chrID ];

  // Loop over all the copy-number calls and search for intervals that overlap the query.
  map<chrom_interval,int>::const_iterator iter;
  for ( iter = intervals.begin(); iter != intervals.end(); iter++ ) {
    chrom_interval overlap = iter->first.intersection( query );
    if ( !overlap.empty() ) regions.push_back( overlap );
  }

  return regions;
}



pair<int,int>
CopyNumberProfile::HRCN( const string & chrom, const int pos ) const
{
  return HRCN_and_interval( chrom, pos ).second;
}



pair<int,int>
CopyNumberProfile::HRCN( const chrom_interval & query ) const
{
  if ( query.empty() ) return make_pair( -1, -1 );
  // Only return a non-null answer if there's a single chrom_interval that
  // contains this entire chrom_interval.
  pair<chrom_interval, pair<int,int> > answer = HRCN_and_interval( query.chrom(), query.start );
  if ( answer.first.contains( query ) ) return answer.second;
  else return make_pair( -1, -1 );
}





/* PRIVATE HELPER FUNCTIONS */


// Input a chromosomal position.
// Return the copy number at that position, as well as the interval over which
// that copy number holds.  If none found, return chrom_interval() and -1.
pair<chrom_interval,int>
CopyNumberProfile::CN_and_interval( const string & chrom, const int pos ) const
{
  // Get the set of copy-number calls (i.e., chrom_intervals) on this chrom.
  const map<chrom_interval,int> & intervals = _CNs[ chrom_to_chrID.at(chrom) ];

  // Loop over all the copy-number calls and search for an interval that
  // contains this position.
  map<chrom_interval,int>::const_iterator iter;
  for ( iter = intervals.begin(); iter != intervals.end(); iter++ )
    if ( iter->first.contains(pos) )
      return *iter;

  // If no interval contains this position, return the dummy answer.
  return make_pair( chrom_interval(), -1 );
}



// Input a chromosomal position.
// Return the HRCN at that position, as well as the interval over which
// that HRCN holds.  If none found, return chrom_interval() and -1.
pair<chrom_interval, pair<int,int> >
CopyNumberProfile::HRCN_and_interval( const string & chrom, const int pos ) const
{
  // Get the set of HRCN (i.e., chrom_intervals) on this chrom.
  const map<chrom_interval, pair<int,int> > & intervals = _HRCNs[ chrom_to_chrID.at(chrom) ];

  // Loop over all the copy-number calls and search for an interval that
  // contains this position.
  map<chrom_interval, pair<int,int> >::const_iterator iter;
  for ( iter = intervals.begin(); iter != intervals.end(); iter++ )
    if ( iter->first.contains(pos) )
      return *iter;

  // If no interval contains this position, return the dummy answer.
  return make_pair( chrom_interval(), make_pair(-1,-1) );
}



/************************************
 *                                  *
 *         OUTPUT FUNCTIONS         *
 *                                  *
 ************************************/



void
CopyNumberProfile::PrintCNStats( ostream & out ) const
{
  int n_calls = 0;
  int64_t total_len = 0;

  out << "Number of copy-number calls per chromosome:" << endl;
  for ( size_t i = 0; i < HumanGenome_n_chroms; i++ )
    if ( !_CNs[i].empty() ) {
      out << chrID_to_chrom[i] << "\t" << _CNs[i].size() << endl;
      n_calls += _CNs[i].size();
      for ( map<chrom_interval,int>::const_iterator it = _CNs[i].begin(); it != _CNs[i].end(); ++it )
	total_len += it->first.len();
    }
  out << "Total number of copy-number calls: " << n_calls << endl;
  out << "Total length of copy-number calls: " << total_len << endl;
}


void
CopyNumberProfile::PrintHRCNStats( ostream & out ) const
{
  int n_calls = 0;
  int64_t total_len = 0;

  out << "Number of HRCN calls per chromosome:" << endl;
  for ( size_t i = 0; i < HumanGenome_n_chroms; i++ )
    if ( !_HRCNs[i].empty() ) {
      n_calls += _HRCNs[i].size();
      out << chrID_to_chrom[i] << "\t" << _HRCNs[i].size() << endl;
      for ( map<chrom_interval,pair<int,int> >::const_iterator it = _HRCNs[i].begin(); it != _HRCNs[i].end(); ++it )
	total_len += it->first.len();
    }
  out << "TOTAL\t" << n_calls << endl;
  out << endl;
}



void
CopyNumberProfile::PrintHRCNs( ostream & out ) const
{
  int n_calls = 0;
  out << "CopyNumberProfile::PrintHRCNs" << endl;

  // Print the number of HRCN calls per chromosome.
  for ( size_t i = 0; i < HumanGenome_n_chroms; i++ )
    if ( !_HRCNs[i].empty() ) {
      n_calls += _HRCNs[i].size();
      map< chrom_interval, pair<int,int> >::const_iterator it;
      for ( it = _HRCNs[i].begin(); it != _HRCNs[i].end(); it++ )
	out << "\t" << it->first << "\t(" << it->second.first << ":" << it->second.second << ")" << endl;
    }

  out << "TOTAL CALLS:\t" << n_calls << endl;
  out << endl;
}




// Write the HRCN calls in this CopyNumberProfile to a BED file.
// BED files written with Write() can be read in the constructor
// CopyNumberProfile(string) (CNs only!)
void
CopyNumberProfile::WriteHRCNs( const string & BEDfile, const bool append ) const
{
  ofstream out( BEDfile.c_str(), append ? ios::app : ios::out );
  map<chrom_interval, pair<int,int> >::const_iterator iter;

  // Output the chromosomes in HumanGenome order.
  for ( size_t i = 0; i < HumanGenome_n_chroms; i++ )
    for ( iter = _HRCNs[i].begin(); iter != _HRCNs[i].end(); iter++ ) {
      pair<int,int> HRCN = iter->second;
      int CN = HRCN.first + HRCN.second;
      out << iter->first.str_BED() << '\t' << CN << '\t' << HRCN.first << '\t' << HRCN.second << endl;
    }

  out.close();
}
