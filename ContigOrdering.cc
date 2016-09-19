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


// For documentation, see ContigOrdering.h
#include "ContigOrdering.h"
#include "TextFileParsers.h" // TODO: use this here
#include "ChromLinkMatrix.h"
#include "TrueMapping.h"

#include <assert.h>
#include <stdlib.h> // srand48, lrand48
#include <limits.h> // INT_MAX
#include <cmath> // sqrt
#include <iostream>
#include <iomanip> // noboolalpha
#include <fstream>
#include <sstream> // ostringstream
#include <string>
#include <vector>
#include <set>
#include <algorithm> // count, reverse



// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "markov/WDAG.h"




// Set random seed at program initialization.
int set_random_seed() {
  srand48( time(0) );
  return 0;
}

int scratch = set_random_seed();





// Constructor.
ContigOrdering::ContigOrdering( const int N_contigs, const bool all_used )
  : _N_contigs( N_contigs ),
    _N_contigs_used( all_used ? N_contigs : 0 )
{
  assert( N_contigs < INT_MAX ); // we want to make sure N_contigs stays positive and ~N_contigs stays negative

  _contigs_used = vector<bool>( _N_contigs, all_used );

  if ( all_used ) Sort(); // default ordering: in numerical order
  else Clear();
}



// Constructor with a pre-supplied ordering.
ContigOrdering::ContigOrdering( const int N_contigs, const vector<int> & data )
  : _N_contigs( N_contigs ),
    _N_contigs_used( data.size() ),
    _data( data )
{
  assert( N_contigs < INT_MAX ); // we want to make sure N_contigs stays positive and ~N_contigs stays negative
  assert( N_contigs >= _N_contigs_used ); // can't use more contigs than we have!


  // Mark which contigs are in the ordering.  Make sure the ordering contains reasonable contig IDs.
  _contigs_used = vector<bool>( _N_contigs, false );
  for ( int i = 0; i < _N_contigs_used; i++ ) {
    int c = data[i];
    assert( c >= 0 && c < _N_contigs );
    assert( !_contigs_used[c] );
    _contigs_used[c] = true;
  }
}



ContigOrdering::ContigOrdering( const int N_contigs, vector<bool> contigs_used )
  : _N_contigs( N_contigs ),
    _contigs_used( contigs_used )
{
  assert( N_contigs < INT_MAX ); // we want to make sure N_contigs stays positive and ~N_contigs stays negative

  // Find how many contigs are not being used (because they have no data.)
  _N_contigs_used = count( contigs_used.begin(), contigs_used.end(), true );

  Sort(); // default ordering: in order.
}


// Create a sub-ordering containing only the contigs in [start,stop).
ContigOrdering::ContigOrdering( const ContigOrdering & order, const int start, const int stop )
  : _N_contigs( order.N_contigs() ),
    _N_contigs_used( 0 )
{
  assert( start >= 0 );
  assert( start < stop );
  assert( stop <= order.N_contigs_used() );

  _contigs_used.resize( _N_contigs, false );

  // Add only the contigs in the subrange of the input ContigOrdering.
  for ( int i = start; i < stop; i++ )
    AddContig( order.contig_ID(i), -1, order.contig_rc(i) );

}



// ReadFile: Input the file format created by WriteFile().  See WriteFile() for the exact format.
void
ContigOrdering::ReadFile( const string & order_file )
{
  assert( boost::filesystem::is_regular_file( order_file ) );

  static const unsigned LINE_LEN = 1000;
  char line[LINE_LEN];
  vector<string> tokens;

  bool first_line = true;
  bool has_Q = false, has_gaps = false;
  int N_contigs_to_read = 0;
  bool old_version = false;


  ifstream in( order_file.c_str(), ios::in );


  // Read the file line-by-line.
  while ( 1 ) {
    in.getline( line, LINE_LEN );
    assert( strlen(line)+1 < LINE_LEN );
    if ( in.fail() ) break;
    boost::split( tokens, line, boost::is_any_of("\t") );

    // There are two possible file formats.  The new version includes commented lines.  The old version (included for backwards compatibility) does not.
    // Use the first line to figure out which format it is.

    if ( first_line ) {

      old_version = ( line[0] != '#' );

      if ( old_version ) {

	// OLD VERSION:
	// The first line in the file indicates the number of contigs.  Once we know this number, we can recreate this ContigOrdering object.
	// It also indicates, by the number of tokens it has, whether or not there are orientation quality scores stored in this file.
	assert( tokens.size() == 1 || tokens.size() == 2 );
	has_Q = ( tokens.size() == 2 );
	has_gaps = false;

	_N_contigs = boost::lexical_cast<int>( tokens[0] );
	*this = ContigOrdering( _N_contigs, false ); // invoke constructor with no contigs used by default
      }

      first_line = false;
    }

    else {

      // In the old version, all lines after the first indicate a contig that is used.  They may also include a contig's orientation quality score.
      if ( old_version ) {

	if ( has_Q ) assert( tokens.size() == 3 );
	else         assert( tokens.size() == 2 );

	int ID = boost::lexical_cast<int>( tokens[0] );
	bool rc = bool( boost::lexical_cast<int>( tokens[1] ) );
	AddContig( ID, -1, rc );

	if ( has_Q ) _orient_Q.push_back( boost::lexical_cast<double>( tokens[2] ) );
      }

      // In the new format version, all lines after the first are either commented, in which case they're part of the header, or uncommented, in which case
      // they describe a contig in the ordering.
      else {

	bool in_header = ( line[0] == '#' );

	// The only header lines we care about are the ones that describe important variables.  These are also the only header lines with a tab.
	if ( in_header ) {

	  if ( tokens.size() < 4 ) continue; // filter out lines without (at least 3) tab characters

	  if ( tokens[1] == "N_contigs" )      {
	    _N_contigs = boost::lexical_cast<int>( tokens[2] );
	    *this = ContigOrdering( _N_contigs, false );
	  }
	  if ( tokens[1] == "N_contigs_used" ) N_contigs_to_read = boost::lexical_cast<int>( tokens[2] );
	  if ( tokens[1] == "has_Q_scores" ) has_Q    = ( tokens[2] == "1" );
	  if ( tokens[1] == "has_gaps"     ) has_gaps = ( tokens[2] == "1" );
	  assert( has_Q || !has_gaps ); // can't have gaps but not quality scores

	}

	// Parse non-header lines.  They should contain 5 tokens: local contig ID, global contig name, contig orientation, orientation quality, gap size.
	else {

	  assert( tokens.size() == 5 );

	  // Get the contig ID and orientation.  Note that we don't actually care about global contig name.
	  int ID = boost::lexical_cast<int>( tokens[0] );
	  bool rc = bool( boost::lexical_cast<int>( tokens[2] ) );
	  assert( ID < _N_contigs );

	  double orient_Q = has_Q ? boost::lexical_cast<double>( tokens[3] ) : -1;
	  int gap      = has_gaps ? boost::lexical_cast<int>   ( tokens[4] ) : -1;

	  AddContig( ID, -1, rc, orient_Q, gap );
	}


      }


    }

  }

  // Sanity checks on the input data.
  assert( _N_contigs_used == N_contigs_to_read );
  assert( _N_contigs_used <= _N_contigs );
  if ( has_gaps ) assert( _gaps.back() == -1 );
}



// WriteFile: Write files in the ContigOrdering format.  The format consists of a header with commented lines; then one line for each
// contig used in the ContigOrdering, with five columns: local ID, global contig name, orientation (1=rc), orientation quality, gap size.
// If global_IDs and global_contig_names aren't given, the contig name column is filled with '.'s.  Likewise for the quality column if !has_Q_scores().
void
ContigOrdering::WriteFile( const string & order_file, const set<int> & global_IDs, const vector<string> * global_contig_names ) const
{
  ofstream out( order_file.c_str(), ios:: out );

  const bool output_old_version = false;

  if ( output_old_version ) {
    // OLD VERSION: before 2013-07-10

    // On the first line, write the number of contigs.  If there are quality scores, add another token, just labeled 'Q'.
    out << _N_contigs;
    if ( has_Q_scores() ) out << "\tQ";
    out << endl;


    // For each contig, write two tokens: the contig ID and orientation.  If orientation quality scores are available, write them as a third token.
    for ( int i = 0; i < _N_contigs_used; i++ ) {
      out << contig_ID(i) << '\t' << noboolalpha << contig_rc(i);
      if ( has_Q_scores() ) out << '\t' << contig_orient_Q(i);
      out << endl;
    }

  }

  else {
    // NEW VERSION: since 2013-07-10


    // Output the header with basic information.  The values of _N_contigs and _N_contigs_used must be output here because they will be parsed by ReadFile().
    out << noboolalpha;
    out << "# ContigOrdering file.  This file represents a ContigOrdering object, which describes an ordering of contigs within a Lachesis cluster." << endl;
    out << "# See ContigOrdering.h for more documentation." << endl;
    out << "#\n";
    out << "# Important numbers:" << endl;
    out << "#\tN_contigs\t" << _N_contigs << "\t(Number of contigs in the input cluster)" << endl;
    out << "#\tN_contigs_used\t" << _N_contigs_used << "\t(Number of contigs ordered by this ordering)" << endl;
    out << "#\thas_Q_scores\t" << has_Q_scores() << "\t(Boolean: has orientation quality scores?)" << endl;
    out << "#\thas_gaps\t" << has_gaps() << "\t(Boolean: have gap sizes between contigs been estimated?)" << endl;
    out << "#\n";
    out << "# Columns:" << endl;
    out << "#contig_ID(local)\tcontig_name\tcontig_rc\torientation_Q_score\tgap_size_after_contig" << endl;

    // Determine whether or not the global contig information has been supplied.  Without it, we can't write out the global contig names.
    bool has_contig_names = !global_IDs.empty() && global_contig_names != NULL;

    // If there's global ID information, convert it from set to vector form.
    vector<int> global_IDs_vec;
    if ( has_contig_names ) {
      // The following assert will fail if the number of contigs in the group used to make this ContigOrdering does not match the number of contigs in the
      // ChromLinkMatrix where the data was taken from.  You may be able to solve the problem by setting the option OVERWRITE_CLMS=1.
      if ( (int) global_IDs.size() != _N_contigs ) {
	cout << "ERROR: global_IDs.size() == " << global_IDs.size() << ", _N_contigs = " << _N_contigs << ".  Try setting OVERWRITE_CLMS=1." << endl;
	assert( (int) global_IDs.size() == _N_contigs );
      }
      for ( set<int>::const_iterator it = global_IDs.begin(); it != global_IDs.end(); it++ )
	global_IDs_vec.push_back(*it);
    }


    // For each contig, write five tokens: local contig ID, global name ID, orientation, orientation quality score, gap size.
    // If any of these pieces of information are unavailable, write '.' as a placeholder.
    for ( int i = 0; i < _N_contigs_used; i++ ) {
      out << contig_ID(i) // local contig ID
	  << '\t' << ( has_contig_names ? (*global_contig_names)[ global_IDs_vec[ contig_ID(i) ] ] : "." ) // global name ID
	  << '\t' << contig_rc(i); // orientation
      if ( has_Q_scores() ) out << '\t' << _orient_Q[i]; // orientation quality score
      else out << "\t.";
      if ( has_gaps() ) out << '\t' << _gaps[i]; // orientation quality score
      else out << "\t.";
      out << endl;
    }

  }

  out.close();
}





// Add a previously unused contig in position pos (pos = -1: add to end).  Use orientation FW or RC.
// orient_Q and gap both default to -1.  If this ContigOrdering has quality scores and/or gaps, they must be set to something else so they can be loaded in.
void
ContigOrdering::AddContig( const int contig_ID, const int pos, const bool rc, const double orient_Q_score, const int gap )
{
  // Sanity checks.
  assert( contig_ID >= 0 );
  assert( contig_ID < _N_contigs );
  assert( !_contigs_used[contig_ID] );

  if ( orient_Q_score == -1 ) assert( !has_Q_scores() );

  int contig_ID_rc = rc ? ~contig_ID : contig_ID;

  // Insert the contig.
  if ( pos == -1 ) {
    _data.push_back( contig_ID_rc );
    if ( orient_Q_score != -1 ) _orient_Q.push_back( orient_Q_score );
    if ( has_gaps() || gap != -1 ) _gaps.push_back( gap );
  }
  else {
    assert( pos <= _N_contigs_used ); // equality is ok here - just add to end
    _data.insert( _data.begin() + pos, contig_ID_rc );
    if ( orient_Q_score != -1 ) _orient_Q.insert( _orient_Q.begin() + pos, orient_Q_score );
    if ( has_gaps() || gap != -1 ) _gaps.insert( _gaps.begin() + pos, gap );
  }


  // Bookkeeping.
  _contigs_used[contig_ID] = true;
  _N_contigs_used++;
}


// Adds a set of previously used contigs in position pos (pos = -1: add to end).  Always use orientation FW.  Quality scores and gaps can't be added this way.
void
ContigOrdering::AddContigs( vector<int> contig_IDs, const int pos )
{
  // Sanity checks.
  int size = contig_IDs.size();
  assert( size <= N_contigs_unused() );
  for ( int i = 0; i < size; i++ ) {
    assert( contig_IDs[i] < _N_contigs );
    assert( !_contigs_used[ contig_IDs[i] ] );
    _contigs_used[ contig_IDs[i] ] = true;
  }
  assert( !has_Q_scores() && !has_gaps() ); // don't add Q-scores or gaps and then modify the underlying ordering!

  // Bookkeeping.
  _N_contigs_used += size;

  // Insert the contigs.
  if ( pos == -1 )
    _data.insert( _data.end(), contig_IDs.begin(), contig_IDs.end() );
  else {
    assert( pos <= _N_contigs_used );
    _data.insert( _data.begin() + pos, contig_IDs.begin(), contig_IDs.end() );
  }
}


// Remove a previously used contig.  Note that the contig is input by its current position in the ContigOrdering, rather than its ID.
void
ContigOrdering::RemoveContig( const int pos )
{
  // Sanity checks.
  assert( pos >= 0 && pos < _N_contigs_used );
  assert( !has_Q_scores() && !has_gaps() ); // don't add Q-scores or gaps and then modify the underlying ordering!

  // Bookkeeping.
  _contigs_used[ contig_ID(pos) ] = false;
  _N_contigs_used--;

  // Remove the contig from the ordering.
  _data.erase( _data.begin() + pos );
}




// Move a contig from position old_pos to new_pos, without changing any contig orientations.
// Contigs in between old_pos and new_pos get shifted by 1 to make up the difference.
void
ContigOrdering::MoveContig( const int old_pos, const int new_pos )
{
  assert( old_pos >= 0 && old_pos < _N_contigs_used );
  assert( new_pos >= 0 && new_pos < _N_contigs_used );
  assert( !has_Q_scores() && !has_gaps() ); // don't add Q-scores or gaps and then modify the underlying ordering!
  if ( old_pos == new_pos ) return; // no work to do

  // Find the contig that's being moved.
  int contig_ID = _data[old_pos];
  //cout << "Moving a contig (" << _data[old_pos] << ") from " << old_pos << " to " << new_pos << endl;

  // Shift the contigs in between as necessary.
  if ( new_pos > old_pos )
    for ( int i = old_pos; i < new_pos; i++ )
      _data[i] = _data[i+1];
  else
    for ( int i = old_pos; i > new_pos; i-- )
      _data[i] = _data[i-1];

  // Put the moved contig into its new place.
  _data[new_pos] = contig_ID;
}





// Invert: Flip the order of the numbers in the range [start,stop].
void
ContigOrdering::Invert( const int start, const int stop )
{
  assert( start >= 0 );
  assert( start <= stop );
  assert( stop < _N_contigs_used );

  // We accomplish the inversion by swapping pairs of numbers in-place (and also inverting them.)
  int swap1 = start, swap2 = stop;
  while ( swap1 < swap2 ) {
    int swap = _data[swap1];
    _data[swap1] = ~_data[swap2];
    _data[swap2] = ~swap;
    swap1++;
    swap2--;
  }

  if ( swap1 == swap2 ) _data[swap1] = ~_data[swap1];


  // Also reverse the quality scores and gaps, if there are any.
  if ( has_Q_scores() ) {
    vector<double>::iterator it1 = _orient_Q.begin(), it2 = _orient_Q.begin();
    it1 += start;
    it2 += stop;
    reverse( it1, it2 );
  }
  if ( has_gaps() ) {
    vector<int>::iterator it1 = _gaps.begin(), it2 = _gaps.begin();
    it1 += start;
    it2 += stop - 1; // the -1 is necessary because _gaps[i] is the gap between contig i and i+1
    reverse( it1, it2 );
  }
}



// InvertRandom: Apply one or more random inversiona via Invert().
void
ContigOrdering::InvertRandom( const int N )
{
  int start, stop;

  for ( int i = 0; i < N; i++ ) {
    do {
      start = lrand48() % _N_contigs_used;
      stop  = lrand48() % _N_contigs_used;
    } while ( start >= stop ); // require start < stop before proceeding

    Invert( start, stop );
  }
}

// PerturbRandom: Apply one or more random changes, either MoveContig() or Invert().
void
ContigOrdering::PerturbRandom( const int N )
{
  int start, stop;
  int N_contigs_squared_m1 = _N_contigs_used * _N_contigs_used - 1;

  for ( int i = 0; i < N; i++ ) {

    // First, choose a random operation: either MoveContig() or Invert().
    bool invert = lrand48() & 1;

    // Next, choose a random distance over which to apply the operation.
    // The following line of code generates a random distance in the range [1,N_contigs_used) and favors small distances over large ones.
    int dist = int( _N_contigs_used - sqrt( 1 + lrand48() % N_contigs_squared_m1 ) );
    //cout << "dist = " << dist << endl;

    // Next, choose a random starting place for the random perturbation.
    start = lrand48() % ( _N_contigs_used - dist );
    stop = start + dist;

    // Lastly, if this is a MoveContig() (not an Invert) then maybe switch the positions.
    if ( !invert && lrand48() & 1 ) { int swap = start; start = stop; stop = swap; }
    //cout << "Perturbation is a " << ( invert ? "inversion" : "move" ) << " between " << start << " and " << stop << endl;

    // Apply the random perturbation.
    if ( invert ) Invert( start, stop );
    else MoveContig( start, stop );
  }

}




void
ContigOrdering::Clear()
{
  _data.clear();
  _N_contigs_used = 0;
  _contigs_used = vector<bool>( _N_contigs, false );
  _orient_Q.clear();
  _gaps.clear();
}



// Put the contigs (the ones currently used) in ascending order with all fw orientations.
// If the contigs are in reference order - e.g., a non-de novo assembly - this is the "correct" order.
void
ContigOrdering::Sort()
{
  // If there are quality scores or gaps, scrap them, because they're about to lose meaning.
  _orient_Q.clear();
  _gaps.clear();

  _data.clear();

  for ( int i = 0; i < _N_contigs; i++ )
    if ( _contigs_used[i] )
      _data.push_back( i ); // i instead of ~i means fw instead of rc

  assert( (int) _data.size() == _N_contigs_used );
}



// Randomize the order and orientation of all contigs in this ContigOrdering, without changing the set of contigs used.
void
ContigOrdering::Randomize()
{
  // If there are quality scores or gaps, scrap them, because they're about to lose meaning.
  _orient_Q.clear();
  _gaps.clear();

  _data.clear();

  // Vector of which contigs are available to add to the ordering.  Unused contigs are pre-marked as unavailable.
  vector<bool> avail = _contigs_used;

  // For each value i, choose a random integer among the ones that haven't already been chosen.  Then add this integer to the ordering.
  for ( int i = 0; i < _N_contigs_used; i++ ) {
    int x_ID = lrand48() % (_N_contigs_used-i); // x_ID is the index (among not-yet-chosen integers) of the integer to choose
    int x = 0; // x is the integer to choose
    int n_seen = 0; // n_seen is the number of not-yet-chosen integers seen so far
    while ( n_seen < x_ID || !avail[x] ) {
      if ( avail[x] ) n_seen++;
      x++;
    }
    assert( x < _N_contigs );
    avail[x] = false;

    // Choose a random orientation for this contig.
    if ( lrand48() % 2 ) x = ~x;
    _data.push_back(x);
  }

}


// Canonicalize: Flip the entire ordering, if necessary, so that _data[0] < _data[last].
void
ContigOrdering::Canonicalize()
{
  if ( contig_ID(0) > contig_ID( _N_contigs_used-1 ) )
    Invert( 0, _N_contigs_used-1 );
}



// Add all previously unused contigs to the end of the ContigOrdering, in forward order.
void
ContigOrdering::AppendUnusedContigs()
{
  assert( !has_Q_scores() && !has_gaps() ); // don't add Q-scores or gaps and then modify the underlying ordering!

  for ( int contig_ID = 0; contig_ID < _N_contigs; contig_ID++ )
    if ( !_contigs_used[contig_ID] )
      _data.push_back(contig_ID);

  assert( (int) _data.size() == _N_contigs );

  // Mark all contigs as used.
  _contigs_used = vector<bool>( _N_contigs, true );
  _N_contigs_used = _N_contigs;
}


// Add an orientation quality score.
void
ContigOrdering::AddOrientQ( const int pos, const double Q )
{
  assert( pos >= 0 && pos < _N_contigs_used );

  // Create the quality score vector, if necessary.
  if ( _orient_Q.empty() ) _orient_Q.resize( _N_contigs_used, 0 );

  _orient_Q[pos] = Q;
}



// Set an element in the _gaps vector.  If necessary, create the vector (and set all other gaps to -1, indicating no data yet.)
void
ContigOrdering::SetGap( const int pos, const int gap_size )
{
  assert( pos >= 0 && pos + 1 < _N_contigs_used );

  if ( !has_gaps() ) {
    _gaps = vector<int>( _N_contigs_used - 1, -1 );
    _gaps.push_back(-1); // add the 'backstop'
  }

  _gaps[pos] = gap_size;
}



// Set the _gaps vector.
void
ContigOrdering::SetGaps( const vector<int> & gaps )
{
  assert( (int) gaps.size() + 1 == _N_contigs_used );
  _gaps = gaps;
  _gaps.push_back(-1); // add the 'backstop'
}




/* OrientationWDAG: Make a WDAG representing contig orientations in this ContigOrdering.
 *
 * The WDAG has a start node, an end node, and two nodes in between for each contig, in the following configuration:
 *
 *         C1_fw --- C2_fw --- ... --- Cn_fw
 *       /       \ /       \ /     \ /       \
 * start          X         X       X          end
 *       \       / \       / \     / \       /
 *         C1_rc --- C2_rc --- ... --- Cn_rc
 *
 * Each edge between two contigs has a weight equal to the log-likelihood of the set of observed link sizes between the two contigs, in the given orientation.
 * For an ASCII illustration of the four possible contig orientations and the implied link sizes, see ChromLinkMatrix -> LoadLinkMatricesFromSAM().
 *
 *************************************************************************************************************************************************************/
WDAG
ContigOrdering::OrientationWDAG( const ChromLinkMatrix * clm ) const
{

  // Build a WDAG representing the possible paths through this ContigOrdering with different orientations of the contigs.
  // This method for building a WDAG follows HMM:to_WDAG().
  WDAG wdag;

  wdag.Reserve( 2 * _N_contigs_used + 2 );

  // Vectors to keep track of the WDAGNode objects.
  vector<WDAGNode *> contig_fw( _N_contigs_used, NULL );
  vector<WDAGNode *> contig_rc( _N_contigs_used, NULL );
  char edge_name[50];

  // Create the start node.
  WDAGNode * start_node = wdag.AddNode();
  wdag.SetReqStart( start_node );


  // Step through the set of contigs and create two new nodes for each contig.
  for ( int i = 0; i < _N_contigs_used; i++ ) {

    contig_fw[i] = wdag.AddNode();
    contig_rc[i] = wdag.AddNode();

    // For the first contig, connect the fw and rc nodes to the start nodes.
    // The input weights are 0 because there's no a priori reason to favor one orientation over the other.
    if ( i == 0 ) {
      sprintf( edge_name, "S__%d_fw", i ); // "S" = start
      contig_fw[i]->AddEdge( start_node, edge_name, 0 );
      sprintf( edge_name, "S__%d_rc", i ); // "S" = start
      contig_rc[i]->AddEdge( start_node, edge_name, 0 );
    }

    // For all contigs past the first, there are four edges that need to be made between this node and the previos node, corresponding to the four possible
    // combined orientations of the two contigs.
    else {

      for ( int rc1 = 0; rc1 < 2; rc1++ )
	for ( int rc2 = 0; rc2 < 2; rc2++ ) {

	  int contig1 = contig_ID(i-1);
	  int contig2 = contig_ID(i);

	  // Each oriented pair of contigs points to an element in the ChromLinkMatrix, which is a vector<int> giving the distance between the reads in those
	  // two contigs, assuming the contigs are immediately adjacent with the specified orientations.
	  double log_like = clm->ContigOrientLogLikelihood( contig1, rc1, contig2, rc2 );

	  // Make the edge.
	  WDAGNode * node1 = rc1 ? contig_rc[i-1] : contig_fw[i-1];
	  WDAGNode * node2 = rc2 ? contig_rc[i  ] : contig_fw[i  ];
	  sprintf( edge_name, "T__%d_%s__%d_%s", i-1, (rc1 ? "rc" : "fw"), i, (rc2 ? "rc" : "fw") ); // "T" = transition
	  node2->AddEdge( node1, edge_name, log_like );

	}
    }


  }


  // Finally, create the ending node.  This node's input weights are all 0.
  WDAGNode * end_node = wdag.AddNode();
  end_node->AddEdge( _N_contigs_used == 0 ? start_node : contig_fw[_N_contigs_used-1], "F", 0 ); // "F" = finish
  end_node->AddEdge( _N_contigs_used == 0 ? start_node : contig_rc[_N_contigs_used-1], "F", 0 ); // "F" = finish
  wdag.SetReqEnd( end_node );

  assert( wdag.N() == 2 * _N_contigs_used + 2 );


  return wdag;
}




// Print the ContigOrdering in human-readable string format - e.g., "1_fw,2_fw,4_rc,3_rc,5_fw".
// TODO: add an alternate printing format that includes gaps
string
ContigOrdering::as_string() const
{
  ostringstream oss;
  for ( int i = 0; i < _N_contigs_used; i++ ) {
    if ( i > 0 ) oss << ',';
    if ( _data[i] >= 0 ) oss <<  _data[i] << "_fw";
    else                        oss << ~_data[i] << "_rc";
  }
  return oss.str();
}



// Print a full description of this ContigOrdering, including as_string().
void
ContigOrdering::Print( ostream & out ) const
{
  out << "ContigOrdering with " << _N_contigs_used << " contigs used (" << N_contigs_unused() << " unused):\t" << as_string() << endl;
}




// DrawDotplot: Create a visual dotplot of this ordering using QuickDotplot.
// In the dotplot, fw contigs will be red dots and rc contigs will be blue dots.
void
ContigOrdering::DrawDotplot( const string & file ) const
{
  // Open the dotplot file for output.
  ofstream out( file.c_str(), ios::out );


  // Plot the points.  Note that the X axis is contig ID (which may indicate true position of a contig) and the Y axis is the order in this ContigOrdering.
  for ( int i = 0; i < _N_contigs_used; i++ )
    out << contig_ID(i) << '\t' << i << '\t' << ( contig_rc(i)?"rc":"fw" ) << endl;

  out.close();


  // Run the QuickDotplot script to generate a dot plot image, which gets placed at out/<file>.jpg.
  // For details on how this script works, see the script itself.
  string cmd = "QuickDotplot " + file;
  system( cmd.c_str() );
}



// DrawDotplotVsTruth: Use QuickDotplot to create a visual dotplot of this ordering compared to the true ordering of the contigs in this ordering.
void
ContigOrdering::DrawDotplotVsTruth( const set<int> & cluster, const TrueMapping & true_mapping, const string & file ) const
{
  cout << "DrawDotplotVsTruth" << endl;

  // First, figure out the true order and orientation of the contigs in this ordering.
  // Note that we don't bother with contigs that are in this cluster but not in this ordering.
  // TODO: for now, just do starting position, no orientation, or true chrom; do more later
  vector<int> true_pos; // maps _local_ contig ID to true start position on chromosome
  for ( set<int>::const_iterator it = cluster.begin(); it != cluster.end(); it++ ) {
    int start = true_mapping.QTargetStart(*it);
    int stop  = true_mapping.QTargetStop (*it);
    int pos = ( start + stop ) / 2;
    true_pos.push_back( pos );
  }
  PRINT3( true_pos.size(), N_contigs(), N_contigs_used() );

  // Open the dotplot file for output.
  ofstream out( file.c_str(), ios::out );

  // Now step through the contigs in the order they appear here.
  for ( int i = 0; i < N_contigs_used(); i++ ) {
    int ID = contig_ID(i);
    assert( _contigs_used[ID] );
    int start = true_pos[ID];
    out << i << "\t" << start << "\tTHING" << endl;
  }

  out.close();


  // Run the QuickDotplot script to generate a dot plot image, which gets placed at out/<file>.jpg.
  // For details on how this script works, see the script itself.
  string cmd = "QuickDotplot " + file;
  system( cmd.c_str() );
}
