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
 * ClusterVec.h
 *
 * A ClusterVec (subtle pun intended) is a class of type vector< set<int> > that describes a set of clusters of contigs.  It represents the result (possibly
 * a work-in-progress) of the clustering algorithm in the GenomeLinkMatrix class, in which contigs are grouped into clusters that predict what chromosome they
 * are likely to be from
 *
 * The ClusterVec is a low-level object.  It contains no data structures except the vector< set<int> > and the number of contigs.  It only knows the contigs by
 * their IDs and knows nothing about them except what cluster they are in.
 *
 *
 *
 * Josh Burton
 * March 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _CLUSTER_VEC__H
#define _CLUSTER_VEC__H



#include <inttypes.h> // uint64_t
#include <iostream>
#include <vector>
#include <set>
using namespace std;




class ClusterVec : public vector< set<int> >
{
 public:

  /* CONSTRUCTORS */
  ClusterVec() { _N_contigs = 0; }
  ClusterVec( const int N_contigs ) { _N_contigs = N_contigs; }
  ClusterVec( const int N_clusters, const int N_contigs ) { resize( N_clusters, set<int>() ); _N_contigs = N_contigs; }
  ClusterVec( const string & file, const vector<string> * contig_names = NULL ) { ReadFile( file, contig_names ); }
  ClusterVec( const vector<int> & cluster_IDs, const bool remove_singletons = true ); // input vector: contig to cluster ID, or -1 for contigs not in a cluster

  vector<int> cluster_IDs() const; // make an output vector like the input vector: contig to cluster ID, or -1 for contigs not in a cluster

  /* FILE I/O: Read/write clusters using a simple file format. */
  // The list of contig_names is optional.  If contig names are given to WriteFile(), it will be more human-readable.  However, when ReadFile() is called, it
  // must be given the _same_ set of contig names that was used in the WriteFile() call (or no contig names, if none were given to WriteFile().)
  void ReadFile ( const string & file, const vector<string> * contig_names = NULL );
  void WriteFile( const string & file, const vector<string> * contig_names = NULL ) const;

  // Print the IDs in a cluster, tab-delimited, with a newline at the end.
  void PrintCluster( const int i, ostream & out = cout, const vector<string> * contig_names = NULL ) const;

  /* QUERIES */
  int N_contigs() const { return _N_contigs; }
  uint64_t SizeSum() const; // total size of all clusters


  void RemoveSingletons(); // remove all single-element clusters, which have no useful info in them

  /* SORTING METHODS - these methods remove any empty clusters */

  void SortBySmallest();
  void SortByMedian();


 private:
  int _N_contigs;
};







#endif
