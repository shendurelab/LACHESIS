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
 * ChromInterval
 *
 * A simple data structure representing an interval on a chromosome.
 * Each chrom_interval contains 4 ints: ID, chrID, start, stop
 *
 *
 * One-letter notations:
 * x: a chrom_interval [a,b)
 * a,b: start and stop (the interval contains a, not b)
 * p: a chromosomal position; may or may not be in [a,b)
 *
 *
 * Josh Burton
 * January 2012
 *
 *****************************************************************************/


#ifndef _CHROM_INTERVAL__H
#define _CHROM_INTERVAL__H

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

// Local modules in the gtools library.
#include "VCF_variant_info.h"


class chrom_interval {

  /* MAIN ELEMENTS */
 public:
  int ID; // this can be used as an index in a vector of chrom_intervals; it has no effect on the object's behavior, not even on operator==
  int chrID; // chromosome ID; converted to name with _chrom_order and _chrom_name.
  int start, stop; // the interval is [start,stop)


  /* STATIC ELEMENTS */
 private:
  // Maps to convert chromosome name (e.g., "chrX") to order (e.g., 23).
  // These maps allow us to store chrom_intervals' chromosomes internally as
  // ints, which is more space-efficient than strings and simplifies sorting.
  static const vector<string> _chrom_name;
  static const map<string,int> _chrom_order;
  static int chrom_order( const string & chrom_name );




  /* CONSTRUCTORS */
 public:
  chrom_interval() : ID(-1), chrID(-1), start(INT_MAX), stop(-INT_MAX) {}
  chrom_interval( const string & s ); // parse form "chr:start-stop"
  chrom_interval( const string & chrom_name, const int & a, const int & b );
 private:
  chrom_interval( const int & chrom_ID, const int & a, const int & b );
 public:

  /* MODIFICATION FUNCTIONS */
  void add( const int & pos ); // add a position to the interval
  void add( const VCF_variant_info & var );
  void clear() { *this = chrom_interval(); }


  /* QUERY FUNCTIONS */

  // Basic queries for data.
  string chrom() const { return _chrom_name.at(chrID); }
  int len()  const { if ( empty() ) return 0; return stop - start; }
  int size() const { if ( empty() ) return 0; return stop - start; }
  bool empty() const { return chrID == -1 && start == INT_MAX && stop == -INT_MAX; }
  bool proper() const { return stop >= start; } // should be true if !empty

  // Test whether a chrom_interval "contains" something.
  bool contains( const int p ) const { return p >= start && p < stop; }
  bool contains( const VCF_variant_info & var ) const { return contains( var.pos ) && ( _chrom_name.at(chrID) == var.chrom ); }
  bool contains( const chrom_interval & x ) const;


  // Two chrom_intervals "overlap" if they have at least one point in common.
  // If the end of one coincides with the beginning of the other, they "abut"
  // but do not overlap.
  bool overlaps( const chrom_interval & x ) const;
  bool abuts( const chrom_interval & x ) const;

  // The "distance" between two chrom_intervals is the amount of empty space
  // between them.  It is 0 if the intervals overlap or abut, and -1 if they
  // are on different chromosomes.
  int distance( const chrom_interval & x ) const;

  // The "overlap size" of two chrom_intervals is the size of the intersection
  // of the intervals.  If is 0 if they do not overlap.
  int overlap_size( const chrom_interval & x ) const { return intersection(x).size(); }

  // Return the intersection/union (in the set theory sense) of two intervals.
  // If there is no overlap, intersection() returns an empty chrom_interval.
  // union_no_gap() and union_with_gap() return the union of two intervals,
  // except that if there is no intersection or abutment, union_no_gap()
  // returns an empty chrom_interval while union_with_gap() returns both
  // intervals plus the gap (unless they're on different chromosomes.)
  chrom_interval intersection  ( const chrom_interval & x ) const;
  chrom_interval union_no_gap  ( const chrom_interval & x ) const;
  chrom_interval union_with_gap( const chrom_interval & x ) const;


  /* OPERATORS */

  // Comparison operator for chrom_intervals.  Used for sorting.
  friend bool operator< ( const chrom_interval &x1, const chrom_interval &x2 );
  friend bool operator> ( const chrom_interval &x1, const chrom_interval &x2 );
  friend bool operator==( const chrom_interval &x1, const chrom_interval &x2 );
  friend bool operator!=( const chrom_interval &x1, const chrom_interval &x2 )
  { return !(x1==x2); }


  // Comparison operator between ints and chrom_intervals.
  // Let p be a position and [a,b) be an interval.  If p >= b, then p > [a,b];
  // if p < a, then p < [a,b); otherwise, p !< [a,b) and p !> [a,b).
  friend bool operator<( const int & p, const chrom_interval & x ) { return p <  x.start; }
  friend bool operator>( const int & p, const chrom_interval & x ) { return p >= x.stop; }
  friend bool operator<( const chrom_interval & x, const int & p ) { return p >= x.stop; }
  friend bool operator>( const chrom_interval & x, const int & p ) { return p <  x.start; }

  /* OUTPUT FUNCTIONS */

  // Output format: "chrom:start-stop".
  string str() const;
  friend ostream & operator<<( ostream & out, const chrom_interval x );
  // Output format: "chrom-start-stop" (necessary e.g., for Windows filenames.)
  string str_dashed() const;
  // Output format: "chrom\tstart\tstop" (useful for BED files)
  string str_BED() const;
};


// Convert a line in a BED or BEDgraph file to a chrom_interval.
chrom_interval
chrom_interval_from_BED_line( const string & line );





#endif
