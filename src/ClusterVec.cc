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


// For documentation, see ClusterVec.h
#include "ClusterVec.h"
#include "TextFileParsers.h" // TokenizeFile

#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>


// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>






// Constructor: Create a ClusterVec from a vector that indicates, for a set of contigs, which cluster each contig belongs to ( ID<0 means no cluster.)
// The set of cluster IDs in the input vector may not be continuous; the clusters are ordered by the order in which their cluster_IDs first appear.
ClusterVec::ClusterVec( const vector<int> & cluster_IDs, const bool remove_singletons )
{
  assert( !cluster_IDs.empty() );
  clear();

  _N_contigs = cluster_IDs.size();

  map<int,int> vec_to_local; // map of cluster ID in the input vector to cluster ID in this ClusterVec
  int local_ID = 0;

  // Step through the input vector.
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( cluster_IDs[i] < 0 ) continue; // negative cluster ID means this contig isn't in any clusters

    // Determine whether or not this cluster ID has been seen before.
    // If not, create a new cluster.
    map<int,int>::const_iterator it = vec_to_local.find( cluster_IDs[i] );
    if ( it == vec_to_local.end() ) {
      vec_to_local.insert( make_pair( cluster_IDs[i], local_ID++ ) );
      set<int> new_cluster;
      new_cluster.insert(i);
      push_back(new_cluster);
    }
    else
      at( it->second ).insert(i);
  }


  // Remove singleton clusters, which contain no information whatsoever (except in the case of chromosome-spanning contigs with a CEN,
  // which occurs in some yeast assemblies.)
  if ( remove_singletons ) RemoveSingletons();

  // Sort the clusters to make them nicer looking for dotplots.
  SortByMedian();
}



// cluster_IDs: Make an output vector like the input vector: contig to cluster ID, or -1 for contigs not in a cluster.
vector<int>
ClusterVec::cluster_IDs() const
{
  //PRINT2( _N_contigs, SizeSum() );
  vector<int> out( _N_contigs, -1 );

  for ( size_t i = 0; i < size(); i++ )
    for ( set<int>::const_iterator it = at(i).begin(); it != at(i).end(); it++ )
      out[*it] = i;

  return out;
}



// Read a file created with ClusterVec::WriteFile().  contig_names must be the same vector (possibly NULL) that was passed to to WriteFile() call.
void
ClusterVec::ReadFile( const string & file, const vector<string> * contig_names )
{
  clear();

  // If contig names are given, make a lookup table.
  map<string,int> contig_name_to_ID;
  if ( contig_names != NULL )
    for ( size_t i = 0; i < contig_names->size(); i++ )
      contig_name_to_ID.insert( make_pair( (*contig_names)[i], i ) );


  // Parse the whole input file into tokens.
  vector< vector<string> > file_as_tokens;
  TokenizeFile( file, file_as_tokens );
  _N_contigs = 0;

  // Loop over every line of tokens.
  for ( size_t i = 0; i < file_as_tokens.size(); i++ ) {
    const vector<string> & tokens = file_as_tokens[i];

    // Ignore the commented lines in the header, except the one that gives the number of contigs.
    if ( tokens[0][0] == '#' ) {

      // Line: "# N_contigs = 1"
      if ( tokens.size() == 4 && tokens[1] == "N_contigs" )
	_N_contigs = boost::lexical_cast<int> ( tokens[3] );

      continue;
    }

    // For each non-commented line, reconstitute a cluster.
    set<int> cluster;
    for ( size_t j = 0; j < tokens.size(); j++ )
      // If necessary, map the contig name onto the contig ID.
      if ( contig_names == NULL )
	cluster.insert( boost::lexical_cast<int>( tokens[j] ) );
      else {
	if ( contig_name_to_ID.find( tokens[j] ) == contig_name_to_ID.end() ) {
	  cerr << "ERROR: Contig name `" << tokens[j] << "' appears in clusters file " << file << " but not in contig names list.  Take this contig out of the clusters file and try again." << endl;
	  assert(0);
	}
	cluster.insert( contig_name_to_ID.at( tokens[j] ) );
	//cout << "Cluster #" << size() << " gets contig #" << j << ": name = " << tokens[j] << ", global ID = " << contig_name_to_ID.at( tokens[j] ) << endl;
      }
    push_back( cluster );
  }

  if ( _N_contigs == 0 ) _N_contigs = SizeSum(); // for backward compatibility with files that don't have the "N_contigs" line
}





// Spit a ClusterVec out to a file.  Each set<int> becomes one tab-delimited line.
// If contig names are given (optional; could be NULL), output contig names instead of IDs.  This file is more human-readable but can't be read by ReadFile()
// unless ReadFile() is also supplied with the same set of contig names.
void
ClusterVec::WriteFile( const string & file, const vector<string> * contig_names ) const
{
  // Sanity check if contig names are given.
  //if ( contig_names != NULL ) PRINT2( _N_contigs, contig_names->size() );
  if ( contig_names != NULL ) assert( _N_contigs == (int) contig_names->size() );

  ofstream out( file.c_str(), ios::out );

  // Write a header.
  out << "# ClusterVec file - see ClusterVec.h for documentation of this object type" << endl;
  out << "#\n";
  out << "# N_contigs = " << _N_contigs << endl;
  out << "# Number of clusters: " << size() << endl;
  out << "# Number of contigs in clusters: " << SizeSum() << endl;
  out << "#\n";
  out << "# There is one (non-commented) line in this file for each cluster." << endl;
  out << "# Each line lists all the contigs, indicated by their ID in the draft assembly, in that cluster." << endl;
  out << "#\n";

  // Write each cluster out on a separate line, in the form of a tab-delimited series of integers (or strings, if contig_names are given.)
  for ( size_t i = 0; i < size(); i++ )
    PrintCluster( i, out, contig_names );

  out.close();
}



// Print the IDs in a cluster, tab-delimited, with a newline at the end.
// If contig_names is non-NULL, print contig names instead of IDs.
void
ClusterVec::PrintCluster( const int i, ostream & out, const vector<string> * contig_names ) const
{
  assert( i < (int) size() );

  for ( set<int>::const_iterator it = at(i).begin(); it != at(i).end(); it++ ) {
    if ( it != at(i).begin() ) out << '\t';
    if ( contig_names != NULL ) out << contig_names->at(*it);
    else out << *it;
  }

  out << '\n';

}



// Find the total size of all clusters.
uint64_t
ClusterVec::SizeSum() const
{
  uint64_t sum = 0;
  for ( size_t i = 0; i < size(); i++ )
    sum += at(i).size();
  return sum;
}







// Sort the clusters by increasing value of the smallest element in each cluster.
void
ClusterVec::SortBySmallest()
{
  set< set<int> > clusters_set;
  for ( size_t i = 0; i < size(); i++ )
    if ( !at(i).empty() )
      clusters_set.insert( at(i) );

  clear();
  for ( set< set<int> >::const_iterator it = clusters_set.begin(); it != clusters_set.end(); ++it )
    push_back( *it );

}


// Remove singleton clusters, which contain no information whatsoever.
void
ClusterVec::RemoveSingletons()
{
  // Make a mapping of old cluster ID -> new cluster ID.  Singleton clusters will be mapped to -1, indicating they'll be destroyed.
  vector<int> new_IDs( size(), -1 );

  int n_non_singletons = 0;
  for ( size_t i = 0; i < size(); i++ )
    if ( at(i).size() > 1 )
      new_IDs[i] = n_non_singletons++;

  // Move all clusters into their new places.
  for ( int i = 0; i < (int) size(); i++ )
    if ( new_IDs[i] != -1 && new_IDs[i] < i ) // new_IDs[i] <= i for all i, but if they're equal, nothing needs to be done
      at( new_IDs[i] ) = at(i);

  resize( n_non_singletons );
}





// Sort the clusters by increasing value of the median element in each cluster.
void
ClusterVec::SortByMedian()
{
  map<int,int> medians;

  // Find the median of each cluster.  (Technically this is just the (n/2)-th element, which isn't quite the same thing as the median.  But calculating the
  // median would (a) be harder and (b) create the possibility of two clusters with the same median, which would be bad for this map.)
  for ( size_t i = 0; i < size(); i++ ) {
    if ( at(i).empty() ) continue;
    int med_idx = at(i).size() / 2;
    set<int>::const_iterator it = at(i).begin();
    for ( int j = 0; j < med_idx; j++ ) ++it;
    medians.insert( make_pair( *it, i ) );
  }

  // Rebuild the cluster order.
  ClusterVec sorted;
  sorted._N_contigs = _N_contigs;
  for ( map<int,int>::const_iterator it = medians.begin(); it != medians.end(); ++it )
    sorted.push_back( at( it->second ) );

  *this = sorted;
}
