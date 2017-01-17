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
 * ContigOrdering.h
 *
 * A ContigOrdering is a defined ordering of the contigs in a ChromLinkMatrix, possibly including an orientation and relative spacing of the contigs as well.
 *
 * Contigs are stored in a ContigOrdering with local IDs, indicating their indices within the set of contigs in the group.  The ContigOrdering object does not
 * know anything about the contigs in its ordering: not their global IDs, not their lengths, not the linkage between them.
 *
 * ContigOrderings are initially created in the ChromLinkMatrix functions MakeTrunkOrder() and MakeFullOrder().  They are assigned orientations by
 * ChromLinkMatrix::OrientContigs(), and assigned spacings by ChromLinkMatrix::SpaceContigs().  The goal of all these algorithmiv functions is to produce a
 * ContigOrdering that best represents the link data, as measured by the function ChromLinkMatrix::OrderingScore().
 *
 * Note that some of the contigs in a group may be "unused", that is, not in the ordering even though they're in the cluster.  This reflects the reality that
 * we sometimes can't tell where a contig is on a chromosome even if we know it's on that chromosome.
 *
 *
 *
 * Josh Burton
 * December 2012
 *
 *************************************************************************************************************************************************************/


#ifndef _CONTIG_ORDERING__H
#define _CONTIG_ORDERING__H

#include "TrueMapping.h"

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "markov/WDAG.h"



#include <iostream>
#include <vector>
#include <set>
using namespace std;


class ChromLinkMatrix;

// This enum indicates orientations. 0 = false = FW; 1 = true = RC.
enum { FW, RC };



class ContigOrdering
{
 public:

  /* CONSTRUCTORS */

  ContigOrdering( const int N_contigs, const bool all_used = true ); // constructor with all contigs used (or none)
  ContigOrdering( const int N_contigs, const vector<int> & data );
  ContigOrdering( const int N_contigs, const vector<bool> contigs_used );
  ContigOrdering( const ContigOrdering & order, const int start, const int stop ); // create a sub-ordering containing only the contigs in [start,stop)
  ContigOrdering( const string & order_file ) { ReadFile( order_file ); }


  /* FILE I/O */

  // ReadFile, WriteFile: Read and write files in the ContigOrdering format.  The format consists of a header with commented lines; then one line for each
  // contig used in the ContigOrdering, with five columns: local ID, global contig name, orientation (1=rc), orientation quality, gap size.
  // If global_IDs and global_contig_names aren't given, the contig name column is filled with '.'s.  Likewise for the quality column if !has_Q_scores().
  void ReadFile ( const string & order_file );
  void WriteFile( const string & order_file, const set<int> & global_IDs = set<int>(), const vector<string> * global_contig_names = NULL ) const;


  /* MODIFICATIONS TO THE ORDERING OF CONTIGS */

  // Local modifications
  void AddContig( const int contig_ID, const int pos = -1, const bool rc = FW, const double orient_Q_score = -1, const int gap = -1 ); // pos = -1: add to end
  void AddContigs( vector<int> contig_IDs, const int pos = -1 ); // adds a set of previously used contigs (orientation always fw)
  void RemoveContig( const int pos ); // remove a contig (note that input integer is a position, not an ID!)
  void MoveContig( const int old_pos, const int new_pos ); // moves a contig from position old_pos to new_pos; doesn't change orientation
  void Invert( const int start ) { Invert(start,start); } // flip the order of one contig
  void Invert( const int start, const int stop ); // flip the order of the numbers in the range [start,stop]
  void InvertRandom( const int N = 1 ); // apply N random inversions via Invert()
  void PerturbRandom( const int N = 1 ); // apply N random changes: either MoveContig() or Invert()

  // Global modifications
  void Clear(); // un-use all contigs
  void Sort(); // puts the contigs in ascending order with fw orientation
  void Randomize(); // creates a totally random ordering
  void Canonicalize(); // flip the entire ordering, if necessary
  void AppendUnusedContigs(); // add all previously unused contigs to the end of the ContigOrdering

  // Quality scores
  void AddOrientQ( const int pos, const double Q ); // see notes for _orient_Q below

  // Gaps
  void SetGap( const int pos, const int gap_size ); // set an element in the _gaps vector; create the vector if necessary
  void SetGaps( const vector<int> & gaps ); // set the _gaps vector
  void ClearGaps() { _gaps.clear(); }

  /* QUERY FUNCTIONS */

  // Queries to see whether or not this ContigOrdering has had its contigs oriented and/or spaced.
  bool has_Q_scores() const { return !_orient_Q.empty(); } // orientation happens in ChromLinkMatrix::OrientContigs()
  bool has_gaps()     const { return !_gaps.empty(); }     // spacing happens in ChromLinkMatrix::SpaceContigs()

  int N_contigs()        const { return _N_contigs; }
  int N_contigs_used()   const { return _N_contigs_used; }
  int N_contigs_unused() const { return _N_contigs - _N_contigs_used; }
  // The following four functions input integers for position in the ContigOrdering, NOT contig IDs.
  int    contig_ID      ( const int pos ) const { assert(pos>=0 && pos<_N_contigs_used); return ( _data[pos] >= 0 ? _data[pos] : ~_data[pos] ); }
  bool   contig_rc      ( const int pos ) const { assert(pos>=0 && pos<_N_contigs_used); return ( _data[pos] < 0 ); }
  double contig_orient_Q( const int pos ) const { assert(pos>=0 && pos<_N_contigs_used); assert( has_Q_scores() ); return _orient_Q[pos]; }
  int    gap_size       ( const int pos ) const { assert(pos>=0 && pos<_N_contigs_used); if ( !has_gaps() ) return -1; return _gaps[pos]; } // gap size after contig
  // The following function inputs contig IDs, NOT integers for position.
  bool contig_used( const int contig_ID ) const { return _contigs_used.at(contig_ID); }


  // Orientation quality scores.  These must be loaded via AddOrientQC() or via ReadFile().


  // Make a WDAG representing contig orientations in this ContigOrdering.
  WDAG OrientationWDAG( const ChromLinkMatrix * clm ) const; // clm is used just to call ContigOrientLogLikelihood()


  /* OUTPUT FUNCTIONS */

  string as_string() const;
  void Print( ostream & out = cout ) const;
  // DrawDotplot: Use QuickDotplot to create a visual dotplot of this ordering.
  void DrawDotplot( const string & file ) const;
  // DrawDotplotVsTruth: Use QuickDotplot to create a visual dotplot of this ordering compared to the true ordering of the contigs in this ordering.
  void DrawDotplotVsTruth( const set<int> & cluster, const TrueMapping & true_mapping, const string & file ) const;



 private:

  /* DATA */

  int _N_contigs; // number of contigs total
  vector<bool> _contigs_used; // flags indicating which contigs are in _data
  int _N_contigs_used; // equal to size of _data; also number of true values in _contigs_used

  /* MAIN DATA STRUCTURE: a vector representing the positions and orientations of contigs.
   * This vector should contain _N_contigs_used distinct integers in range [_N_contigs,_N_contigs).
   * Forward contigs are represented by positive numbers, while reversed contigs are represented by negative numbers: contig x in rc is given the number ~x. */
  vector<int> _data;

  // A vector representing the estimated gap sizes between contigs.  _gaps[i] represents the gap *after* contig #i.  As a placeholder, _gaps.back() = 0.
  // This vector will be empty until one of SetGap and SetGaps is called (by ChromLinkMatrix::SpaceContigs()); then it will have length _N_contigs_used.
  vector<int> _gaps;


  // Quality scores for orientations.  This vector does NOT parallel the _data vector; it is only used in ContigOrderings in which AddOrientQ() has been
  // called, which should only happen in ContigOrderings that have been fixed.  If you call AddOrientQ() and then any other modification function, the quality
  // scores vector will stop making sense.
  // All orientation quality scores are initially calculated in ChromLinkMatrix::OrientContigs() and loaded into this class via AddOrientQ().
  vector<double> _orient_Q;
};







#endif
