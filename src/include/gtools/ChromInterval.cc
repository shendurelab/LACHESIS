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


// For documentation, see ChromInterval.h.
#include "ChromInterval.h"

#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <assert.h>

// Local modules in the gtools library.
#include "HumanGenome.h"
#include "VCF_variant_info.h"

// Boost libraries.
#include <boost/algorithm/string.hpp> // split
#include <boost/lexical_cast.hpp>


// Initialize the static data structures _chrom_name and _chrom_order.
// These allow us to convert between chromosome name strings and an internal
// representation of chromosomes as integers.
// The conversion defines the canonical ordering of chroms in the human genome.
const vector<string>  chrom_interval::_chrom_name  = HumanGenome_chroms();
const map<string,int> chrom_interval::_chrom_order = HumanGenome_chrom_IDs();





// Convert a chromosome name to an integer chromosome ID.
// This function includes automatic error-checking (for unrecognized inputs)
// and is the preferred method for internal conversions.
int
chrom_interval::chrom_order( const string & chrom_name )
{
  map<string,int>::const_iterator it = _chrom_order.find( chrom_name );
  // If this assert fails, the chrom name is unrecognized.
  assert ( it != _chrom_order.end() );
  return it->second;
}




/* CONSTRUCTORS */



// Parse a string of the form "chr:start-stop".
// If the input string is not of this form, this function will fail an assert
// or a lexical_cast.
chrom_interval::chrom_interval( const string & s )
  : ID(-1)
{
  // Parse the chromosome name.
  size_t p0 = s.find(':');
  assert( p0 != string::npos ); // input string must contain ':'
  chrID = chrom_order( s.substr( 0, p0 ) );

  // Parse the start and stop integers.
  size_t p1 = s.find( '-', p0+1 );
  assert( p1 != string::npos ); // string must contain '-' after ':'
  start = boost::lexical_cast<int>( s.substr( p0+1, p1-p0-1 ) );
  stop  = boost::lexical_cast<int>( s.substr( p1+1 ) );

  assert( start <= stop );
}


// Constructor from chrom name and start/stop.
chrom_interval::chrom_interval( const string & chrom_name, const int & a, const int & b )
  : ID(-1), chrID( chrom_order(chrom_name) ), start(a), stop(b)
{
  assert( start <= stop );
}


// Constructor from chrom internal numbering and start/stop.
// This constructor is private.
chrom_interval::chrom_interval( const int & chrom_ID, const int & a, const int & b )
  : ID(-1), chrID(chrom_ID), start(a), stop(b)
{
  assert( start <= stop );
}







bool
chrom_interval::contains( const chrom_interval & x ) const
{
  if ( chrID != x.chrID ) return false;
  if ( start >  x.start ) return false;
  if ( stop  <  x.stop  ) return false;
  return true;
}




// Two chrom_intervals "overlap" if they have at least one point in common.
bool
chrom_interval::overlaps( const chrom_interval & x ) const
{
  if ( chrID != x.chrID ) return false;
  if ( start >= x.stop  ) return false;
  if ( stop  <= x.start ) return false;
  return true;
}


// If the end of one chrom_interval coincides with the beginning of the other,
// they "abut" but do not overlap.
bool
chrom_interval::abuts( const chrom_interval & x ) const
{
  if ( chrID != x.chrID ) return false;
  if ( start == x.stop  ) return true;
  if ( stop  == x.start ) return true;
  return false;
}


// The "distance" between two chrom_intervals is the amount of empty space
// between them.  It is 0 if the intervals overlap or abut, and -1 if they
// are on different chromosomes.
int
chrom_interval::distance( const chrom_interval & x ) const
{
  if ( chrID != x.chrID ) return -1;
  if ( start >= x.stop  ) return start - x.stop;
  if ( stop  <= x.start ) return x.start - stop;
  return 0;
}




// Return the intersection (shared region) between these two chrom_intervals.
// If there is no intersection, return an empty chrom_interval().
chrom_interval
chrom_interval::intersection( const chrom_interval & x ) const
{
  if ( !overlaps(x) ) return chrom_interval();
  return chrom_interval( chrID, max( start, x.start ), min( stop, x.stop ) );
}


// union_no_gap: Return the set-theoretical union of two chrom_intervals, if it
// exists as a single chrom_interval.  That is, if the two chrom_intervals
// intersect or abut, return a single chrom_interval encompassing both; if not,
// return an empty chrom_interval().
chrom_interval
chrom_interval::union_no_gap( const chrom_interval & x ) const
{
  if ( !overlaps(x) && !abuts(x) ) return chrom_interval();
  return chrom_interval( chrID, min( start, x.start ), max( stop, x.stop ) );
}


// union_with_gap: Return the set-theoretical union of two chrom_intervals,
// adding in the gap between the intervals if necessary.  If they're on
// different chromosomes, return an empty chrom_interval().
chrom_interval
chrom_interval::union_with_gap( const chrom_interval & x ) const
{
  if ( chrID != x.chrID ) return chrom_interval();
  return chrom_interval( chrID, min( start, x.start ), max( stop, x.stop ) );
}




// Add this position (presumably on the same chromosome) to the interval.
void
chrom_interval::add( const int & pos )
{
  assert( !empty() ); // if empty, we have no way of knowing what chromosome
  start = min( start, pos );
  stop  = max( stop , pos );
}

void
chrom_interval::add( const VCF_variant_info & var )
{
  if ( empty() ) chrID = chrom_order( var.chrom );
  else assert( _chrom_name.at(chrID) == var.chrom );
  add( var.pos );
}






string
chrom_interval::str() const
{
  if ( empty() ) return "<empty chrom_interval>";
  return _chrom_name.at(chrID) + ":" + boost::lexical_cast<string>(start) + "-" + boost::lexical_cast<string>(stop);
}




ostream &
operator<<( ostream & out, const chrom_interval x )
{
  out << x.str();
  return out;
}



string
chrom_interval::str_dashed() const
{
  if ( empty() ) return "";
  return _chrom_name.at(chrID) + "-" + boost::lexical_cast<string>(start) + "-" + boost::lexical_cast<string>(stop);
}


string
chrom_interval::str_BED() const
{
  if ( empty() ) return "";
  return _chrom_name.at(chrID) + "\t" + boost::lexical_cast<string>(start) + "\t" + boost::lexical_cast<string>(stop);
}




bool
operator<( const chrom_interval & x1, const chrom_interval & x2 )
{
  // Compare chrom_intervals by chromosome (internal numbering!) and then by
  // interval location.
  if ( x1.chrID != x2.chrID ) return x1.chrID < x2.chrID;
  if ( x1.start != x2.start ) return x1.start < x2.start;
  return x1.stop < x2.stop;
}


bool
operator>( const chrom_interval & x1, const chrom_interval & x2 )
{ return ( x2 < x1 ); }


// Note that the equality operator does NOT consider the variable ID.
bool
operator==( const chrom_interval & x1, const chrom_interval & x2 )
{
  if ( x1.start != x2.start ) return false;
  if ( x1.stop  != x2.stop  ) return false;
  if ( x1.chrID != x2.chrID ) return false;
  return true;
}



// Convert a line in a BED or BEDgraph file to a chrom_interval.
// Each (non-commented) line in a BED file is tab-delimited.  The first three
// tokens in each line describe the chromosome interval, and the rest give
// details about the interval (which we ignore here).
chrom_interval
chrom_interval_from_BED_line( const string & line )
{
  // Split the line on tabs to get the first three tokens.
  vector<string> tokens;
  boost::split( tokens, line, boost::is_any_of("\t") );
  assert( tokens.size() >= 3 );

  return chrom_interval( tokens[0],
			 boost::lexical_cast<int>( tokens[1] ),
			 boost::lexical_cast<int>( tokens[2] ) );
}
