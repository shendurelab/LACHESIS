// For documentation, see ChromLinkMatrix.h
#include "ChromLinkMatrix.h"
#include "ClusterVec.h"
#include "ContigOrdering.h"
#include "TextFileParsers.h" // ParseTabDelimFile
#include "TrueMapping.h"


#include <assert.h>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <iomanip> // setprecision, boolalpha
#include <algorithm> // count, max_element
#include <numeric> // accumulate

// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>


// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "markov/WDAG.h"
#include "gtools/SAMStepper.h" // SAMStepper, NTargetsInSAM

// Modules in <samtools-dir> (must add -I~/include)
#include "sam.h"




// most_distant_node: Graph theory helper function for find_longest_path.
// Given a graph (in the form of an adjacency list) and a node A, find the node B that is most distant from A (in the same component as A.)
// Return the node B; also returns a path from A to B, inclusive (so path.size() is one more than the distance from A to B.)
int
most_distant_node( const vector< vector<int> > & adj_list, const int A, vector<int> & path_AB )
{
  assert( A < (int) adj_list.size() );
  
  bool verbose = false;
  
  // Make a vector to keep track of each node's distance from A.
  vector<int> dist_from_A( adj_list.size(), -1 );
  vector< vector<int> > path_from_A( adj_list.size(), vector<int>(0) );
  dist_from_A[A] = 0;
  path_from_A[A].push_back(A);
  
  // Initially, set B = A.
  int B = A;
  int dist_AB = 0;
  path_AB.resize( 1, A );
  
  // Set up a queue of "live" nodes with which to step through the graph, starting at A.
  queue<int> live_nodes;
  live_nodes.push(A);
  
  // Step through the graph (or at least the component containing A) and find each node's distance from A.
  while ( !live_nodes.empty() ) {
    int node = live_nodes.front();
    live_nodes.pop();
    int dist = dist_from_A[node];
    vector<int> & path = path_from_A[node];
    if ( verbose ) cout << "node = " << node << "; dist = " << dist << "; path has length " << path.size() << "; adj_list = ";
    
    // Check the distance of this node.
    if ( dist > dist_AB ) {
      path_AB = path;
      dist_AB = dist;
      B = node;
    }
    
    // Add this node's neighbors to the queue.
    for ( size_t i = 0; i < adj_list[node].size(); i++ ) {
      int node2 = adj_list[node][i];
      if ( verbose ) cout << "\t" << node2;
      if ( dist_from_A[node2] == -1 ) { // e.g., if this node hasn't been seen yet
	dist_from_A[node2] = dist + 1;
	path_from_A[node2] = path;
	path_from_A[node2].push_back( node2 );
	live_nodes.push(node2);
      }
    }
    if ( verbose ) cout << endl;
  }
  
  if ( verbose ) {
    cout << "MAXIMUM DISTANCE IS " << dist_AB << " AT NODE " << B << endl;
    cout << "PATH:";
    for ( size_t i = 0; i < path_AB.size(); i++ ) cout << "\t" << path_AB[i];
    cout << endl;
  }
  
  return B;
}



// find_longest_path: Graph theory function to find the longest path in a tree (given in the form of an adjacency list.)
// This is done via a two-step algorithm, as outlined here: http://stackoverflow.com/questions/13187404/linear-time-algorithm-for-longest-path-in-tree
// Runtime is O(n) for trees; longer for arbitrary graphs.  This function has only been tested for spanning trees, which have only one non-trivial component.
vector<int>
find_longest_path( const vector< vector<int> > & adj_list )
{
  vector<int> path(0);
  srand48( time(0) );
  
  // If there are no adjacencies left in this graph, return an empty vector.
  int size = adj_list.size();
  bool has_adj = false;
  for ( int i = 0; i < size; i++ )
    if ( !adj_list[i].empty() ) has_adj = true;
  if ( !has_adj ) return path;
  
  
  // Longest Path Step 1: Pick a random node A and find the node B that is furthest from A.
  // We do this repeatedly with different nodes A (of degree at least 1), just to make sure we have the right answer.
  // This always finds the largest path on one of the non-trivial components of the graph.  Spanning trees have only one non-trivial component, although they
  // may have components of just one node, representing a contig with no data.
  bool stochastic = false;
  int A = 0;
  while ( adj_list[A].empty() && A < size ) A++;
  
  if ( stochastic ) {
    int A = lrand48() % size;
    while ( adj_list[A].empty() && A < size ) A = lrand48() % size;
  }
  
  assert( A < size );
  
  int B = most_distant_node( adj_list, A, path );
  
  // Longest Path Step 2: Find the node C that is furthest from B.  The path from B to C is guaranteed to be the longest path in the tree.
  most_distant_node( adj_list, B, path );
  
  // Canonicalize the path order, just to make this function deterministic (in the case of a spanning tree with only one component.)
  if ( path[0] > path[path.size()-1] ) reverse( path.begin(), path.end() );
  
  return path;
}




// AssertFilesExist: Check for the existence of a set of filenames.  It's a good idea to run this on a set of SAM files before loading them in, lest we spend
// hours loading in six SAM files only to fail an assert when the seventh file fails to exist.
void
AssertFilesExist( const vector<string> & files )
{
  for ( size_t i = 0; i < files.size(); i++ )
    assert( boost::filesystem::is_regular_file( files[i] ) );
}




// Null constructor.
ChromLinkMatrix::ChromLinkMatrix()
{
  _matrix_init = false;
  _N_contigs = -1;
  _contig_size = -1;
  _contig_lengths.clear();
  _longest_contig = -1;
  _contig_RE_sites.clear();
  _most_contig_REs = -1;
  _repeat_factors.clear();
  _tree.clear();
  _SAM_files.clear();
  _CP_score_dist = 1e7;
}



// Create an empty non-de novo ChromLinkMatrix with a contig size and chromosome length.  This is designed to be used in calls to LoadNonDeNovoCLMsFromSAM().
ChromLinkMatrix::ChromLinkMatrix( const string & species, const int contig_size, const int chrom_len )
  : _species( species )
{
  assert( contig_size > 0 );
  
  _matrix_init = false;
  _N_contigs = ceil( double(chrom_len) / contig_size );
  _contig_size = contig_size; // this implies !DeNovo()
  _contig_lengths.resize( _N_contigs, 0 ); // for non-de novo CLMs, this remains at all 0's, and _longest_contig remains at -1
  _longest_contig = -1;
  _contig_RE_sites.resize( _N_contigs, 0 );
  _most_contig_REs = -1;
  _repeat_factors.clear();
  _tree.clear();
  _SAM_files.clear();
  _CP_score_dist = 1e7;
  
  int N_bins = 2 * _N_contigs;
  cout << Time() << ": Creating a new ChromLinkMatrix for a chromosome with " << _N_contigs << " contigs of size " << _contig_size << " (matrix size = " << N_bins << "x" << N_bins << ")" << endl;
  
  InitMatrix();
}






// Create an empty ChromLinkMatrix with a specific number of contigs.  This is designed to be used in calls to LoadDeNovoCLMsFromSAM().
// Note that this CLM can't know how long its contigs are because no SAM file has been loaded in yet.
ChromLinkMatrix::ChromLinkMatrix( const string & species, const int cluster_size )
  : _species( species ),
    _N_contigs( cluster_size )
{
  assert( cluster_size != 0 );
  
  _matrix_init = false;
  _contig_size = 0; // this implies DeNovo()
  _longest_contig = -1;
  _tree.clear();
  _SAM_files.clear();
  _CP_score_dist = 1e7;
  
  cout << Time() << ": Creating a new ChromLinkMatrix for a cluster with " << _N_contigs << " contigs (matrix size = " << N_bins() << "x" << N_bins() << ")" << endl;
  
  InitMatrix();
  
}



// Constructor: Load a de novo ChromLinkMatrix from a set of SAM files, and only include contigs with IDs in cluster #cluster_ID.
ChromLinkMatrix::ChromLinkMatrix( const string & species, const vector<string> & SAM_files, const string & RE_sites_file, const ClusterVec & clusters, const int cluster_ID )
  : _species( species )
{
  assert( !SAM_files.empty() );
  assert( !clusters.empty() );
  
  AssertFilesExist( SAM_files );
  
  _matrix_init = false;
  _N_contigs = clusters[cluster_ID].size();
  _contig_size = 0; // this implies DeNovo()
  _longest_contig = -1;
  _tree.clear();
  _SAM_files = SAM_files;
  _CP_score_dist = 1e7;
  
  InitMatrix();
  
  LoadFromSAMDeNovo( SAM_files, RE_sites_file, clusters, cluster_ID );
  return;
}




// Destructor.
ChromLinkMatrix::~ChromLinkMatrix()
{
  if ( _matrix_init ) FreeMatrix();
}








// ReadFile: Read the data from file CLM_file into this ChromLinkMatrix.
// Overwrite any existing data.
// The file CLM_file should have been created by a previous call to ChromLinkMatrix::WriteFile() with heatmap=false.
void
ChromLinkMatrix::ReadFile( const string & CLM_file )
{
  cout << Time() << ": ChromLinkMatrix::ReadFile   <-  " << CLM_file << "\t" << flush;
  assert( boost::filesystem::is_regular_file( CLM_file ) );
  
  // Set initial values that MUST be overwritten later.
  _contig_size = -1;
  _SAM_files.clear();
  _matrix_init = false;
  
  static const unsigned LINE_LEN = 10000000;
  char line[LINE_LEN];
  vector<string> tokens;
  
  string contig_lens_file = "";
  string RE_sites_file = "";
  bool seen_data = false;
  
  ifstream in( CLM_file.c_str(), ios::in );
  
  // Read the file line-by-line.
  
  while ( 1 ) {
    in.getline( line, LINE_LEN );
    assert( strlen(line)+1 < LINE_LEN );
    if ( in.fail() ) break;
    
    // Examine the header lines carefully.  These header lines should have been generated in ChromLinkMatrix::WriteFile.
    if ( line[0] == '#' ) {
      boost::split( tokens, line, boost::is_any_of(" ") );
      assert( tokens.size() > 3 );
      
      if ( tokens[1] == "Species" ) // line: "# species = human"
	_species = tokens[3];
      else if ( tokens[1] == "N_contigs" ) { // line: "# N_contigs = 1"
	if ( _matrix_init ) assert( _N_contigs == boost::lexical_cast<int>( tokens[3] ) );
	else {
	  _N_contigs = boost::lexical_cast<int>( tokens[3] );
	  InitMatrix();
	}
      }
      else if ( tokens[1] == "contig_size" ) // line: "# contig_size = 1"
	_contig_size = boost::lexical_cast<int>( tokens[3] );
      else if ( tokens[1] == "contig_lens_file" ) // line: "# contig_lens_file = <filename>" or "# contig_lens_file = ."
	contig_lens_file = tokens[3];
      else if ( tokens[1] == "RE_sites_file" ) // line: "# RE_sites_file = <filename>"
	RE_sites_file = tokens[3];
      else if ( tokens[0] == "heatmap" ) // line: "# heatmap = false"
	assert( tokens[3] == "false" );
      else if ( tokens[1] == "SAM" ) { // line: "# SAM files used in generating this dataset: test.sam"
	for ( size_t i = 8; i < tokens.size(); i++ )
	  _SAM_files.push_back( tokens[i] );
	AssertFilesExist( _SAM_files );
      }
    }
    
    else if ( line[0] == 'X' ) continue; // skip the header line of the matrix itself
    
    // If this is not a header line, put data in the matrix.
    // The data may be sparse (i.e., not every X,Y pair appearing) and that's ok.
    else {
      seen_data = true;
      boost::split( tokens, line, boost::is_any_of("\t") );
      int X = boost::lexical_cast<int>( tokens[0] );
      int Y = boost::lexical_cast<int>( tokens[1] );
      int Z = boost::lexical_cast<int>( tokens[2] );
      if ( Z == 0 ) continue;
      assert( (int) tokens.size() == 3 + Z );
      for ( int i = 0; i < Z; i++ )
	_matrix[X][Y].push_back( boost::lexical_cast<int>( tokens[3+i] ) );
    }
  }
  
  in.close();
  
  cout << "\tN contigs = " << _N_contigs << endl;
  
  if ( _N_contigs > 1 && !seen_data )
    cerr << "WARNING: ChromLinkMatrix::ReadFile: CLM file '" << CLM_file << "' has multiple contigs but no link data" << endl;
  //PRINT2( seen_data, has_links() );
  assert( seen_data == has_links() );
  
  assert( _contig_size != -1 );
  assert( !_SAM_files.empty() );
  assert( _matrix_init );
  
  
  // If this is a de novo CLM, read contig lengths from the contig_lens_file and contig RE sites from the contig_RE_sites file.
  if ( DeNovo() ) { // this is equivalent to _contig_size == 0
    assert( contig_lens_file != "" && contig_lens_file != "." );
    _contig_lengths = ParseTabDelimFile<int>( contig_lens_file, 0 );
    assert( _N_contigs == (int) _contig_lengths.size() ); // if this fails, the contig_lens file has the wrong number of lines
    FindLongestContig();
    
    LoadRESitesFile( RE_sites_file );
  }
  
  
  // If this isn't a de novo CLM, the contig_lens_file should have been marked as ".".
  else {
    assert( contig_lens_file == "." );
    _longest_contig = -1;
  }
  
  
  CalculateRepeatFactors();
}






// WriteFile: Write the data in this ChromLinkMatrix to file CLM_file.
// The output format is a long tall table with (4*_N_contigs^2) rows and three (or more) columns: X, Y, Z, [data].
// X and Y are bin IDs (bin ID = 2 * contig ID + contig orientation).  Z is an integer representing the amount of data in bin [X,Y].
// If heatmap = false (default), then following Z is actually a tab-separated list of all the integers in bin [X,Y].
// This format is designed for easy input to ChromLinkMatrix::ReadFile() (if heatmap = false), or to R (if heatmap = true).
// If this is a de novo CLM, also write the contig lengths to an auxiliary file.
void
ChromLinkMatrix::WriteFile( const string & CLM_file, const bool heatmap ) const
{
  cout << Time() << ": ChromLinkMatrix::WriteFile(" << ( heatmap ? "heatmap" : "CLM" ) << ") -> " << CLM_file << endl;
    
  
  // If this is a de novo CLM, set filenames for the auxiliary contig lengths and contig RE sites file.
  string contig_lens_file = DeNovo() ? CLM_file + ".lens" : ".";
  string contig_RE_sites_file = DeNovo() ? CLM_file + ".RE_sites" : ".";
  
  bool seen_data = false;
  ofstream out( CLM_file.c_str(), ios::out );
  
  // Print a header to file.  The header is for easier human reading and also contains numbers used by ChromLinkMatrix::ReadFile().
  out << "# ChromLinkMatrix file - see ChromLinkMatrix.h for documentation of this object type" << endl;
  out << "# Species = " << _species << endl;
  out << "# De novo CLM? " << boolalpha << DeNovo() << endl;
  out << "# N_contigs = " << _N_contigs << endl;
  out << "# contig_size = " << _contig_size << " (ignored if a contig_lens_file is supplied)" << endl;
  out << "# contig_lens_file = " << contig_lens_file << endl;
  out << "# RE_sites_file = " << contig_RE_sites_file << endl;
  out << "# heatmap = " << boolalpha << heatmap << endl;
  out << "# SAM files used in generating this dataset:";
  for ( size_t i = 0; i < _SAM_files.size(); i++ )
    out << " " << _SAM_files[i];
  out << endl;
  
  // Print the table to file.
  out << "X\tY\tZ" << ( heatmap ? "" : "\tlink_lengths" ) << endl;
  
  for ( int X = 0; X < 2*_N_contigs; X++ )
    for ( int Y = 0; Y < 2*_N_contigs; Y++ ) {
      
      vector<int> & Z = _matrix[X][Y];
      //PRINT3( X, Y, Z.size() );
      if ( Z.empty() ) continue; // this makes the matrix sparse
      
      // If writing a heatmap, write just the vector size, and only once for each pair of contigs (i.e., don't differentiate by orientation.)
      if ( heatmap ) {
	if ( X%2 || Y%2 ) continue;
	seen_data = true;
	out << X/2 << '\t' << Y/2 << '\t' << Z.size() << endl;
      }
      // If writing a CLM output file, write the entire vector, and write for all four contig orientations.
      else {
	seen_data = true;
	out << X << '\t' << Y << '\t' << Z.size();
	for ( size_t i = 0; i < Z.size(); i++ )
	  out << '\t' << Z[i];
	out << endl;
      }
    }
  
  out.close();
  
  if ( _N_contigs > 1 && !seen_data )
    cerr << "WARNING: ChromLinkMatrix::ReadFile: CLM file '" << CLM_file << "' has multiple contigs but no link data" << endl;
  //PRINT2( seen_data, has_links() );
  assert( seen_data == has_links() );
  
  // If this is a de novo CLM, also write the contig lengths and RE sites to auxiliary files.
  if ( DeNovo() ) {
    
    ofstream out2( contig_lens_file.c_str(), ios::out );
    for ( int i = 0; i < _N_contigs; i++ )
      out2 << _contig_lengths[i] << endl;
    out2.close();
    
    ofstream out3( contig_RE_sites_file.c_str(), ios::out );
    for ( int i = 0; i < _N_contigs; i++ )
      out3 << ( _contig_RE_sites[i] - 1 ) << endl; // subtract 1 to make up for the 1 added in LoadRESites()
    out3.close();
    
  }
}





// DrawHeatmap: Call WriteFile("heatmap.txt"), then run the R script "heatmap.R", which uses R and ggplot2 to make a heatmap image of this ChromLinkMatrix.
void
ChromLinkMatrix::DrawHeatmap( const string & heatmap_file ) const
{
  WriteFile( "heatmap.txt", true ); // true means to write the file in heatmap format
  
  // The R script "heatmap.R" is hardwired to take "heatmap.txt" as input and write to ~/public_html/heatmap.jpg.
  // For details on how this script works, see the script itself.
  cout << Time() << ": Plotting a heatmap at ~/public_html/" << heatmap_file << endl;
  system( "heatmap.R" );
  
  // Copy the file heatmap.jpg into the place desired.
  if ( heatmap_file != "" && heatmap_file != "heatmap.jpg" )
    system( ( "cp ~/public_html/heatmap.jpg ~/public_html/" + heatmap_file ).c_str() );
  
  cout << Time() << ": Done plotting!" << endl;
}




// LoadFromNonDeNovo: Fill this de novo ChromLinkMatrix with data from one or more SAM/BAM files, but only for the contigs in cluster #cluster_ID.
// This is a wrapper to the function LoadDeNovoCLMsFromSAM, which loads data into one or more SAM files.
void
ChromLinkMatrix::LoadFromSAMDeNovo( const string & SAM_file, const string & RE_sites_file, const ClusterVec & clusters, const int cluster_ID )
{
  // Make a list of ChromLinkMatrix pointers in which the only non-null pointer is to this ChromLinkMatrix.  This is needed for LoadDeNovoCLMsFromSAM.
  vector<ChromLinkMatrix *> CLMs( clusters.size(), NULL );
  CLMs[cluster_ID] = this;
  
  LoadDeNovoCLMsFromSAM( SAM_file, RE_sites_file, clusters, CLMs );
}



// This multi-SAM-file function is a wrapper for the one-SAM-file function.
void
ChromLinkMatrix::LoadFromSAMDeNovo( const vector<string> & SAM_files, const string & RE_sites_file, const ClusterVec & clusters, const int cluster_ID )
{
  AssertFilesExist( SAM_files );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadFromSAMDeNovo( SAM_files[i], RE_sites_file, clusters, cluster_ID );
}






// LoadFromSAMNonDeNovo: Fill this non-de novo ChromLinkMatrix with data from one or more SAM/BAM files describing Hi-C reads aligned to a reference genome.
// The chromosome ID is specified by the chromosome order in the reference fasta file.
// This is a wrapper to the function LoadNonDeNovoCLMsFromSAM, which loads data into one or more SAM files.
void
ChromLinkMatrix::LoadFromSAMNonDeNovo( const string & SAM_file, const int chrom_ID )
{
  // Look at the SAM file's header to determine the number of chromosomes in the reference genome.
  int N_chroms = NTargetsInSAM( SAM_file );
  assert( chrom_ID >= 0 && chrom_ID < N_chroms );
  
  // Make a list of ChromLinkMatrix pointers in which the only non-null pointer is to this ChromLinkMatrix.  This is needed for LoadNonDeNovoCLMsFromSAM.
  vector<ChromLinkMatrix *> CLMs( N_chroms, NULL );
  CLMs[chrom_ID] = this;
  
  LoadNonDeNovoCLMsFromSAM( SAM_file, CLMs );
}




// This multi-SAM-file function is a wrapper for the one-SAM-file function.
void
ChromLinkMatrix::LoadFromSAMNonDeNovo( const vector<string> & SAM_files, const int chrom_ID )
{
  AssertFilesExist( SAM_files );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadFromSAMNonDeNovo( SAM_files[i], chrom_ID );
}






bool
ChromLinkMatrix::has_links() const
{
  for ( int i = 0; i < _N_contigs; i++ )
    for ( int j = i+1; j < _N_contigs; j++ )
      if ( !_matrix[ 2*i ][ 2*j ].empty() )
	return true;
  
  return false;
}



// EmptyRows: Return a vector indicating which contigs in the ChromLinkMatrix have data at all.  Contigs in centromeres will end up as false.
// If flag_adjacent = true, mark contigs adjacent to data-free contigs as unused as well.
vector<bool>
ChromLinkMatrix::ContigsUsed( const bool flag_adjacent ) const
{
  if ( _N_contigs == 1 ) return vector<bool>( 1, true ); // handle edge case
  
  vector<bool> used( _N_contigs, true );
  
  
  // Loop over all bins in a row.
  for ( int i = 0; i < _N_contigs; i++ ) {
    bool has_data = false;
    
    for ( int j = 0; j < _N_contigs; j++ )
      if ( !_matrix[ 2*min(i,j) ][ 2*max(i,j) ].empty() ) {
	has_data = true;
	break;
      }
    
    // If there's no data in this row, flag it (and maybe its neighbors).
    if ( !has_data ) {
      used[i] = false;
      if ( flag_adjacent ) {
	if ( i > 0 )            used[i-1] = false;
	if ( i+1 < _N_contigs ) used[i+1] = false;
      }
    }
  }
  
  return used;
}




// ContigOrientLogLikelihood: Return the log-likelihood of observing two contigs in a given orientation, as defined by the links between the contigs.
// Note that the result is not normalized to the total number of links.  Hence if two contigs have a lot of links between them, the LogLikelihood will be
// very small - for ALL possible orientations of these contigs.
double
ChromLinkMatrix::ContigOrientLogLikelihood( const int c1, const bool rc1, const int c2, const bool rc2 ) const
{
  double log_like = 0;
  
  // Find the vector of contig distances that represents this pair of contigs with this orientation.
  // For an ASCII illustration of these orientations, see AddLinkToMatrix().
  const vector<int> & dists = _matrix[ 2*c1 + (rc1?1:0) ][ 2*c2 + (rc2?1:0) ];
  
  // Each link of distance x makes a contribution of 1/x to the likelihood; hence -ln(x) to the log-likelihood.
  int N_links = dists.size();
  for ( int j = 0; j < N_links; j++ )
    log_like -= log( double(dists[j]) );
  
  return log_like;
}









// OrderingScore: Find the "score" of this ContigOrdering, indicating how well it matches up with the Hi-C link data.
// Specifically: This ContigOrdering implies a certain distribution of Hi-C link sizes.  The ContigOrdering's OrderingScore measures how well this
// implied distribution matches the theoretical distribution, in which the contact probability at a distance x drops off as 1/x.
// This function simply calculates the sum of 1/x over all links.
// (This theoretical distribution comes from Figure 4 in the original 2009 Hi-C paper: http://www.sciencemag.org/content/326/5950/289.full#F4).
// This function is computationally very expensive!
// If oriented = true (default), takes orientation into account.  oriented=false is faster.  Scores with and without orientation can't be directly compared.
// If range_start and range_stop are given, the score is only calculated between contig pairs for which at least one contig falls within the range
// [range_start,range_stop).  When testing the effect of adding a single contig or a shred of contigs, use this: it's O(N) instead of O(N^2)!
// If range_start is given, but range_stop = -1, the score is only calculated between contig pairs on either side of range_start.
double
ChromLinkMatrix::OrderingScore( const ContigOrdering & order, const bool oriented, const int range_start, const int range_stop ) const
{
  assert( order.N_contigs() == _N_contigs ); // sanity check
  
  double score = 0;
  
  
  // Loop over all distinct values of contig1 and contig2, such that contig1 < contig2.
  for ( int i1 = 0; i1+1 < order.N_contigs_used(); i1++ ) {
    int contig1 = order.contig_ID(i1);
    int rc1     = order.contig_rc(i1);
    
    // Keep track of the distance between contig1 and contig2 (which is occupied by other intervening contigs.)
    int contig_dist = 0;
    
    for ( int i2 = i1+1; i2 < order.N_contigs_used(); i2++ ) {
      if ( range_start != -1 && range_stop == -1 )
	if ( i1 >= range_start || i2 < range_start ) continue;
      if ( range_stop != -1 )
	if ( i1 >= range_stop || i2 < range_start ) continue;
      int contig2 = order.contig_ID(i2);
      int rc2     = order.contig_rc(i2);
      
      // Each oriented pair of contigs points to an element in the ChromLinkMatrix, which is a vector<int> giving the distance between the reads in those
      // two contigs, assuming the contigs are immediately adjacent with the specified orientations.
      // For an ASCII illustration of these distances, see AddLinkToMatrix().
      const vector<int> & dists = _matrix[ 2*contig1+rc1 ] [ 2*contig2+rc2 ];
      
      // If we take orientation into account, we have to do a bunch of computation on the individual links.
      if ( oriented ) {
	
	// Adjust the distances to account for the space between contig1 and contig2 in this ContigOrdering.
	// Each link of distance x makes a contribution to the score that is equal to 1/x.
	for ( size_t i = 0; i < dists.size(); i++ ) {
	  if ( dists[i] == 0 ) PRINT4( i1, i2, i, dists[i] );
	  score += 1.0 / double( dists[i] + contig_dist );
	}
      }
      
      contig_dist += ( _contig_size != 0 ? _contig_size : _contig_lengths[contig2] );
      if ( contig_dist > _CP_score_dist ) break;
      
      if ( !oriented ) {
	// Just count the number of links between the two contigs.
	score += dists.size() / double( contig_dist );
      }
    }
  }
  
  assert( !isnan(score) );
  return score;
}



// EnrichmentScore: Find the "enrichment" for this ContigOrdering: the degree to which it causes the Hi-C links between contigs to be between close contigs.
// It's analogous to the concentration of signal along the main diagonal of the heatmap.
// To calculate enrichment, the OrderingScore is divided by the "null score", which is the value that would result from this ContigOrdering under the
// condition of constant link density.
double
ChromLinkMatrix::EnrichmentScore( const ContigOrdering & order ) const
{
  assert( order.N_contigs() == _N_contigs ); // sanity check
  
  if ( order.N_contigs() <= 2 || order.N_contigs_used() <= 1 ) return 0; // handle degenerate case
  
  
  // Find the total number of links and the total (squared) distance in which the links occur.  Use these numbers to calculate a link density.
  // Also calculate the null score, which depends on the set of contigs used in the ContigOrdering, but not on their order or orientation.
  int64_t N_links = 0, total_len = 0, total_len_sq = 0;
  double null_score = 0;
  
  for ( int i1 = 0; i1 < order.N_contigs_used(); i1++ ) {
    int contig1 = order.contig_ID(i1);
    
    total_len += _contig_lengths[contig1];
    
    // Keep track of the distance between contig1 and contig2 (which is occupied by other intervening contigs.)
    int contig_dist = 0;
    
    for ( int i2 = i1+1; i2 < order.N_contigs_used(); i2++ ) {
      int contig2 = order.contig_ID(i2);
      
      // If the distance is too far, skip this contig pair.
      contig_dist += ( _contig_size != 0 ? _contig_size : _contig_lengths[contig2] );
      if ( contig_dist > _CP_score_dist ) break;
      
      // Add to the contig length tallies.
      N_links += _matrix[ 2*contig1 ] [ 2*contig2 ].size();
      int64_t len_sq = (int64_t) _contig_lengths[contig1] * (int64_t) _contig_lengths[contig2];
      total_len_sq += len_sq;
      
      // Add to the null score.  Note that at the moment, the null score is purely a measure of the relative lengths and positions of contigs.
      null_score += len_sq / contig_dist;
      if ( isnan( null_score ) ) PRINT7( i1, i2, contig1, contig2, null_score, len_sq, contig_dist ); // trying to track down a weird bug
    }
  }
  
  // Calculate the average link density per bp of sequence.  Multiply the null score by this to convert it from a measure of length to a density of links.
  double link_density = double( N_links ) / total_len_sq;
  null_score *= link_density;
  assert( null_score != 0 );
  assert( !isnan( null_score ) );
  
  
  // Lastly, find the OrderingScore and divide it by the null score.
  double score = OrderingScore(order);
  
  //PRINT5( score, null_score, N_links, total_len_sq, link_density );
  return score / null_score;
}






// PrefilterLinks: Find contig pairs in which the distribution of Hi-C link positions on the contigs suggest long-range rather than short-range contacts.
// Flag these pairs of contigs and filter them out.
// The principle is that if two *small* regions have a very high frequency of links between them, these links probably represent either technical noise or a
// long-range chromatin interaction rather than the usual signal of proximity.  If two contigs have a lot of links between them but the links are all localized
// in this fashion, we don't want to mistake them for close contigs.  Contigs that are truly close on the genome will have a lot of links that are spread out
// all across the contigs, and those links will not be flagged here.
// The TrueMapping is included for evaluation.
void
ChromLinkMatrix::PrefilterLinks( const set<int> & cluster, const TrueMapping * mapping )
{
  assert(0); // this function was test code and is now deprecated
  
  cout << Time() << ": PrefilterLinks" << endl;
  cout << "THING\tX\tY\tZ\n";
  
  // Make a map of local contig ID to true position on the chromosome.
  map<int,int> true_starts;
  set<int>::const_iterator it = cluster.begin();
  for ( int i = 0; i < _N_contigs; i++ ) {
    int start = mapping->QTargetStart( *it );
    if ( _contig_lengths[i] < 50000 ) continue; // only include large contigs
    true_starts[start] = i;
    //PRINT3( *it, mapping->QTargetName(*it), start );
    it++;
  }
  vector<int> true_pos( _N_contigs, -1 );
  int i = 0;
  for ( map<int,int>::const_iterator it = true_starts.begin(); it != true_starts.end(); it++ ) {
    true_pos[i++] = it->second;
    //PRINT( it->second );
  }
  
  
  // Loop over all contig pairs.
  for ( int contig1 = 0; contig1 < _N_contigs; contig1++ ) {
    for ( int contig2 = contig1+1; contig2 < _N_contigs; contig2++ ) {
      const vector<int> & links = _matrix[ 2*contig1 ][ 2*contig2 ];
      int N_links = links.size();
      if ( N_links == 0 ) continue; // no links, so nothing to do
      
      int len1 = _contig_lengths[contig1], len2 = _contig_lengths[contig2];
      
      // The possible range of link distances between these contigs is [0,len1+len2).  If the links are evenly distributed on the contigs, we expect the
      // distribution of link sizes to be roughly trapezoidal, as follows:
      //
      //            |  L1      L2
      //            |  ._______.
      // frequency  |  /       \     ^
      //  of links  | /		\    | expected_link_freq
      //            |/           \   v
      //            +----------------
      //            0           L1+L2
      //
      // We want to measure the extent to which the distribution of link distances deviates from this.
      double expected_link_freq = double( N_links ) / len2;
      
      // Make bins for possible link distances.
      int bin_size = 10000;
      int N_bins = (len1+len2)/bin_size + 1;
      
      // Determine how many links actually fall into each bin.
      vector<int> N_links_in_bin( N_bins, 0 );
      for ( int i = 0; i < N_links; i++ )
	N_links_in_bin[ links[i] / bin_size ]++;
      
	
      
      double esum = 0;
      double var = 0;
      
      for ( int i = 0; i < N_bins; i++ ) {
	int bin_start = i*bin_size, bin_stop = (i+1)*bin_size;
	double expected_N_links_in_bin = expected_link_freq * bin_size;
	if ( bin_start < len1 ) {
	  double mid = ( bin_start + min( bin_stop, len1 ) ) / 2.0;
	  expected_N_links_in_bin *= ( mid / len1 );
	}
	else if ( bin_stop > len2 ) {
	  double mid = ( bin_stop + max( bin_start, len2 ) ) / 2.0;
	  expected_N_links_in_bin *= ( ( len1+len2 - mid ) / len1 );
	}
	//PRINT2( expected_N_links_in_bin, N_links_in_bin[i] );
	double diff = expected_N_links_in_bin - N_links_in_bin[i];
	var += diff * diff;
	esum += expected_N_links_in_bin;
      }
      
      var /= N_links;
      
      //PRINT6( contig1, contig2, N_links, N_bins, esum, var );
      int pos1 = true_pos[contig1];
      int pos2 = true_pos[contig2];
      if ( pos1 != -1 && pos2 != -1 ) { // skip cases of missing data in the true_pos map (because of short contigs or because it's not properly one-to-one)
	cout << "THING\t" << pos1 << "\t" << pos2 << "\t" << var << endl;
	cout << "THING\t" << pos2 << "\t" << pos1 << "\t" << var << endl;
      }
      
      /*
      int bin_counts[100];
      for ( size_t i = 0; i < 100; i++ ) bin_counts[i] = 0;
      
      // Find the number of links landing in each sub-bin.  We will calculate the variance of this distribution and use it as a statistic of non-randomness.
      for ( int i = 0; i < N_links; i++ )
	
      
      
      // We've gotten all the bins.  Now calculate the variance of the bin contents.
      double avg_N_links = N_links / 100.0;
      double var = 0;
      
      for ( int i = 0; i < 100; i++ ) {
	//PRINT2( i, bin_counts[i] );
	double diff = bin_counts[i] - avg_N_links;
	var += diff * diff;
      }
      
      var /= N_links;
      
      PRINT4( contig1, contig2, N_links, var );
      */
    }
  }
  
  exit(0);
  
}



ContigOrdering
ChromLinkMatrix::MakeTrunkOrder( const int min_N_REs ) const
{
  cout << Time() << ": MakeTrunkOrder!" << endl;
  assert ( !_contig_RE_sites.empty() );
  
  // Handle the trivial case, where there are fewer than two contigs or there are no links between the contigs.
  if ( !has_links() ) return ContigOrdering( _N_contigs, false );
  
  // First, find the minimum spanning tree (MST) that spans a graph version of this ChromLinkMatrix.  The minimum spanning tree is a way to connect all the
  // nodes (contigs) such that the total weight of all the connections (measured as 1 / the number of links in support of a connection) is at a minimum.
  _tree = FindSpanningTree( min_N_REs );
  
  
  // SmoothThornsInTree: Remove from the tree as many "thorns" (i.e., single-vertex spurs from the main trunk) as possible.
  SmoothThornsInTree(_tree);
  
  // Optional output: Make a dotplot of the tree.
  if (0) {
    string file = "dotplot.tree.txt";
    cout << Time() << ": Drawing a dotplot of this spanning tree at ~/public_html/" << file << endl;
    ofstream out( file.c_str(), ios::out );
    
    for ( int i = 0; i < _N_contigs; i++ )
      for (size_t j = 0; j != _tree[i].size(); ++j)
	if ( (int)_tree[i][j] > i )
	  out << i << "\t" << _tree[i][j] << "\n";
    out.close();
    
    // Run the QuickDotplot script to generate a dot plot image, which gets placed at ~/public_html/<file>.jpg.
    string cmd = "QuickDotplot " + file;
    std::system( cmd.c_str() );
  }
  
  
  // Optional output: Make a graph image for this spanning tree.
  // NOT RECOMMENDED for graphs over 500 vertices due to runtime problems and unreadability of output.
  if (0) PlotTree( _tree, "tree.png" );
  
  // The MST is a tree, not a linear ordering (or ContigOrdering).  Convert it to a ContigOrdering by finding the longest path.
  ContigOrdering trunk = TreeTrunk( _tree, true );
  OrientContigs( trunk );
  
  //trunk.DrawDotplot( "order_dotplot.1.trunk.txt" );
  
  return trunk;
}




ContigOrdering
ChromLinkMatrix::MakeFullOrder( const int min_N_REs, const bool use_CP_score ) const
{
  cout << Time() << ": MakeFullOrder!" << endl;
  
  // Handle the trivial case, where there are fewer than two contigs or there are no links between the contigs.
  if ( !has_links() ) return ContigOrdering( _N_contigs, false );
  
  
  assert( !_tree.empty() ); // if this fails, you need to run MakeTrunkOrder() first
  
  
  // Reinsert the pruned contigs into the ordering.
  ContigOrdering ordering = ReinsertShreds( _tree, min_N_REs, use_CP_score );
  
  
  // Determine the proper orientation of contigs in this ContigOrdering.
  OrientContigs( ordering );
  ordering.Print();
  
  
  
  // Find the score for the "true" ordering.  (This is worthless except in simulated assemblies in which the contigs are in the proper order.)
  //double true_score = OrderingScore( ContigOrdering( _N_contigs, true ) );
  cout << setprecision(8);
  
  // Report on the size of this ContigOrdering.  Also report on its "enrichment", which describes the degree to which the Hi-C links between
  // the placed contigs are in close contigs - i.e., the amount of signal along the main diagonal of the heatmap.
  cout << "Full ordering:\t";
  ReportOrderingSize(ordering);
  cout << "..... ENRICHMENT SCORE = " << EnrichmentScore(ordering) << endl;
  //ordering.DrawDotplot("order_dotplot.2.shreds+oriented.txt");
  
  ContigOrdering true_ordering = ordering;
  true_ordering.Sort();
  cout << "TRUTH ENRICHMENT SCORE = " << EnrichmentScore( true_ordering ) << endl;
  
  // Try out all possible inversions.  Use ordering scores (non-normalized) instead of enrichment scores, to save CPU time.
  /*
  double base_score = OrderingScore(ordering);
  double true_score = OrderingScore(true_ordering);
  PRINT2( base_score, true_score );
  int N_contigs_used = ordering.N_contigs_used();
  for ( int i = 0; i < N_contigs_used; i++ )
    for ( int j = i+1; j < N_contigs_used; j++ ) {
      ContigOrdering new_ordering = ordering;
      new_ordering.Invert(i,j-1);
      double score = OrderingScore(new_ordering);
      if ( score > base_score )
	cout << "BETTER:\ti = " << i << "\tj = " << j << "\tscore = " << score << endl;
    }
  */
  
  
  return ordering;
}



// SpaceContigs: No, not contigs in space.
// Given an ordering of the contigs in this ChromLinkMatrix, estimate the spacing between them.  Fill the variable _gaps in the ContigOrdering.
// Once a ContigPermutation has had its gaps estimated by SpaceContigs, it can be converted into a super-scaffold.
// Method: For each pair of adjacent contigs, consider the set of link distances between them.  (Note that their orientations must have already been determined
// by OrientContigs.)  Determine what number to add to these link distances to make them most concordant with the expectations of the LinkSizeDistribution.
void
ChromLinkMatrix::SpaceContigs( ContigOrdering & order, const LinkSizeDistribution & link_size_distribution ) const
{
  cout << Time() << ": SpaceContigs" << endl;
  
  // Sanity checks.
  assert( _SAM_files == link_size_distribution.SAM_files() );
  assert( order.N_contigs() == _N_contigs );
  assert( order.has_Q_scores() ); // can't space contigs if they're not already oriented
  
  vector<int> gaps;
  
  // Loop over all the pairs of adjacent contigs.
  for ( int i = 0; i+1 < order.N_contigs_used(); i++ ) {
    //if ( i != 4 ) continue; // TEMP
    
    // Find the vector of link lengths associated with this adjacency.  This vector describes the distance that the links would have if the contigs were
    // immediately adjacent.  For an ASCII illustration of these distances, see AddLinkToMatrix().
    int contig1 = order.contig_ID(i);
    int rc1     = order.contig_rc(i);
    int contig2 = order.contig_ID(i+1);
    int rc2     = order.contig_rc(i+1);
    if ( contig1 > contig2 ) { int swap = contig1; contig1 = contig2; contig2 = swap; swap = !rc1; rc1 = !rc2; rc2 = swap; }
    PRINT4( contig1, contig2, rc1, rc2 );
    const vector<int> & dists = _matrix[ 2*contig1+rc1 ] [ 2*contig2+rc2 ];
    
    
    // Assuming these contigs are separated at a distance D, the actual size of all the Hi-C links is D higher than the reported number.
    // Determine the value of D that makes this set of links most concordant with the expectations of the LinkSizeDistribution.
    int D = link_size_distribution.FindDistanceBetweenLinks( _contig_lengths[contig1], _contig_lengths[contig2], dists );
    
    double RF = _repeat_factors[contig1] * _repeat_factors[contig2];
    cout << Time() << ": Result for gap #" << i << ": " << D << '\t' << RF << endl;
    //D /= sqrt(RF); // TEMP
    gaps.push_back(D);
  }
  
  order.SetGaps( gaps );
}






// Find a spanning tree that spans a graph version of this ChromLinkMatrix.  The output is a vector< vector<int> > with the tree as an adjacency list.
// This goes a long way toward finding the proper contig order (but not orientation).
// This is actually an iterative process, with the goal of finding a spanning tree with an optimally long longest_path: First we make a graph, then we find the
// MST (minimum spanning tree) of the graph, then we use the longest path in this MST to re-create the graph and find a new spanning tree.
// The tree from the first iteration is an MST, but subsequent trees are technically merely spanning trees.
// The code here follows the example given in [boost]/libs/graph/example/prim-example.cpp.
vector< vector<int> >
ChromLinkMatrix::FindSpanningTree( const int min_N_REs ) const
{
  // If this is a small cluster, return a trivial result, rather than failing an assert later.
  assert( _N_contigs > 1 );
  if ( _N_contigs == 2 ) {
    vector< vector<int> > best_tree( _N_contigs );
    for ( int i = 0; i+1 < _N_contigs; i++ ) {
      best_tree[i  ].push_back(i+1);
      best_tree[i+1].push_back(i);
    }
    return best_tree;
  }
  
  
  using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS, property<vertex_distance_t, int>, property < edge_weight_t, double > > Graph;
  typedef pair<int, int> E;
  
  bool use_RE = !_contig_RE_sites.empty();
  const vector<int> & lens = use_RE ? _contig_RE_sites : _contig_lengths;
  
  static const int MAX_EDGE_WEIGHT_PENALTY = 10; // HEUR: setting this higher will make the search for a long path more aggressive
  
  
  // Blank out all small contigs to force them to be left out of the ordering.
  // Reduce the min contig length as necessary so that the cluster has at least 3 long contigs
  // If this isn't done, there will be too few real edges in the graph and prim_minimum_spanning_tree will get confused.
  int N_long_contigs = 0;
  for ( int i = 0; i < _N_contigs; i++ )
    if ( lens[i] >= min_N_REs ) N_long_contigs++;
  
  int min_len_reduced = min_N_REs;
  assert( _N_contigs >= 3 );
  while ( N_long_contigs < 3 && min_len_reduced > 0 ) {
    min_len_reduced -= (use_RE ? 1 : 1000);
    N_long_contigs = 0;
    for ( int i = 0; i < _N_contigs; i++ )
      if ( lens[i] >= min_len_reduced ) N_long_contigs++;
  }
  PRINT2( min_len_reduced, N_long_contigs );
  
  
  vector< vector<int> > best_tree; // this will contain the output
  
  // Set up the edges of the graph.
  int N_edges = _N_contigs * (_N_contigs-1) / 2;
  E * edges = new E[N_edges];
  int max_edge_ID = 0;
  for ( int i = 0; i < _N_contigs; i++ )
    for ( int j = i+1; j < _N_contigs; j++ )
      edges[max_edge_ID++] = make_pair(i,j);
  assert( max_edge_ID == N_edges );
  
  // Allocate memory for the edge weigths.
  double * weights = new double[N_edges];
  
  vector<int> trunk; // this will be repeatedly updated
  int edge_weight_penalty = 1;
  size_t trunk_size = 0;
  
  // In each iteration of the loop, we do the following 3 things:
  // 1. Build a graph.  If there was a previous iteration, use the longest path from that iteration to inform the graph building for this iteration.
  // 2. Run Prim's algorithm to find the MST on this graph.
  // 3. Find the "trunk", the longest path in this MST.
  for ( size_t iteration = 0; iteration < 100; iteration++ ) {
    
    // Set up the edge weights to be infinity (i.e., unusable) by default.  The prim_minimum_spanning_tree algorithm can't handle non-connected graphs.  Hence,
    // if two vertices in the graph should have no edge between them, they'll have an edge with weight infinity.
    for ( int i = 0; i < N_edges; i++ )
      weights[i] = INFINITY;
    
    max_edge_ID = 0;
    vector<bool> isolated( _N_contigs, true );
    
    // 1. Build a graph.  If there was a previous iteration, use the longest path from that iteration to inform the graph building for this iteration.
    // Define the vertices, edges, and edge weights.  Each contig becomes a vertex, and the Hi-C links between contigs are represented as edges between them.
    
    // Loop over all pairs of contigs (vertices).
    for ( int i = 0; i < _N_contigs; i++ ) {
      vector<int>::const_iterator i_it = find( trunk.begin(), trunk.end(), i );
      
      for ( int j = i+1; j < _N_contigs; j++ ) {
	max_edge_ID++;
	
	if ( lens[i] < min_len_reduced || lens[j] < min_len_reduced ) continue; // leave out short contigs
	//PRINT7( iteration, i, j, link_weight, lens[i], lens[j], min_len_reduced );
	
	double link_weight = LinkDensity(i,j);
	//PRINT3( i, j, link_weight );
	if ( link_weight == 0 ) continue; // no links between these two contigs
	double edge_weight = 1.0 / link_weight; // better-linked contigs should have a *smaller* edge weight in the graph
	
	// Use the already-observed longest path to guide further development of the tree.  The idea is to use the longest path as a scaffold, and fit into it
	// the contigs that were initially left out of it.  So, for the purposes of extending the longest path, prevent non-adjacent contigs in the longest
	// path from moving directly to each other, and penalize adjacent ones.
	vector<int>::const_iterator j_it = find( trunk.begin(), trunk.end(), j );
	if ( i_it != trunk.end() && j_it != trunk.end() ) {
	  if ( i_it + 1 != j_it && i_it - 1 != j_it ) continue; // these two contigs aren't adjacent in the longest path
	  else                                        edge_weight *= edge_weight_penalty;
	}
	
	// Keep track of which vertices are really isolated so we can remove them from the graph later.
	isolated[i] = false;
	isolated[j] = false;
	
	//cout << "FindSpanningTree iteration #" << iteration << ": Adding edge #" << max_edge_ID << " between " << i << " and " << j << " with weight " << edge_weight << endl;
	weights[max_edge_ID-1] = edge_weight; // subtract 1 because max_edge_weight has already been incremented in this loop
      }
      
    }
    assert( max_edge_ID == N_edges );
    
    // If all edges are isolated, there is effectively no graph - i.e., no sufficiently long contigs have any links to other sufficiently long contigs.
    // In this case, let's just forget about trying to make an ordering and just return an empty one.
    if ( count( isolated.begin(), isolated.end(), false ) == 0 ) {
      best_tree.resize( _N_contigs, vector<int>() );
      return best_tree;
    }
    
    
    if ( _N_contigs > 10000 ) cout << "WARNING: ChromLinkMatrix contains a large number of contigs (" << _N_contigs << ").  About to create a Graph object, which will take a lot of memory.  May run out of memory and throw a bad_alloc error." << endl;
    
    // Build the Graph object.
    Graph g( edges, edges + max_edge_ID, weights, _N_contigs );
    vector < graph_traits<Graph>::vertex_descriptor > MST_adjs( _N_contigs );
    
    // 2. Run Prim's algorithm to find the MST on this graph.
    cout << Time() << ": prim_minimum_spanning_tree (iteration #" << iteration << ", previous trunk_size = " << trunk_size << ")" << endl;
    prim_minimum_spanning_tree(g, &MST_adjs[0]);
    //cout << Time() << ": prim_minimum_spanning_tree done!" << endl;
    
    for ( int i = 0; i < _N_contigs; i++ )
      if ( lens[i] < min_len_reduced || lens[ MST_adjs[i] ] < min_len_reduced )
	MST_adjs[i] = i;
    
    // Convert the tree from the output format of prim_minimum_spanning_tree into adjacency-list format.
    vector< vector<int> > MST( _N_contigs, vector<int>( 0 ) );
    for ( int i = 0; i < _N_contigs; i++ ) {
      int j = MST_adjs[i];
      
      // Don't make connections for isolated contigs.
      if ( i == j || isolated[i] || isolated[j] ) continue;
      
      MST[i].push_back(j);
      MST[j].push_back(i);
    }
    
    
    
    // 3. Find the "trunk", the longest path in this tree.
    trunk = find_longest_path( MST );
    
    // If this tree's path length is an improvement over the last iteration, great!  Save the result.
    if ( trunk.size() > trunk_size ) {
      trunk_size = trunk.size();
      best_tree = MST;
    }
    
    // Otherwise, get more aggressive by increasing the edge weight penalty.  If it's at max, stop iterating.
    else {
      edge_weight_penalty++;
      if ( edge_weight_penalty >= MAX_EDGE_WEIGHT_PENALTY ) break;
      continue;
    }
    //PRINT( trunk_size );
    
    
    // Optional output: Make a dot file and then an image of this complete graph.
    // NOT RECOMMENDED for graphs over 30 vertices, and unlikely to be useful for graphs above 10.
    if (0) {
      ofstream out( "graph.dot", ios::out );
      boost::write_graphviz( out, g, make_label_writer( get(vertex_index,g) ), make_label_writer( get(edge_weight,g) ) );
      out.close();
      std::system ( "dot -Tpng graph.dot > ~/public_html/graph.png" );
      cout << "Drew an image of this graph at ~/public_html/graph.png" << endl;
    }
    
  }
  
  
  
  // Cleanup.
  delete[] edges;
  delete[] weights;
  
  cout << Time() << ": FindSpanningTree: Trunk length = " << trunk_size << endl;
  assert( !best_tree.empty() );
  
  return best_tree;
}






// Find all structures of this type: (a "thorn node" A with one neighbor, B, which has exactly two other neighbors C and D, neither of which is a leaf node)
//
// --- C --- B --- D ----
//           |
//           |
//           A
//
// "Smooth out" these structures, if possible, in one of these two ways:
//
// --- C     B --- D ----                 --- C --- B     D ---
//       \_  |                  OR                  |  _/
//         \ |                                      | /
//           A                                      A
//
void
ChromLinkMatrix::SmoothThornsInTree( vector< vector<int> > & tree ) const
{
  cout << Time() << ": SmoothThornsInTree...\t";
  assert( (int) tree.size() == _N_contigs );
  
  int N_smoothed = 0;
  const double aggression = 3.0; // HEUR: parameter; 1.0 is most aggressive
  
  // First, look for thorns.  Loop over all nodes A that might be thorns.
  // NOTE: Because we modify the tree in place within each iteration, the overall result could be affected by the contig ordering.
  
  for ( int A = 0; A < _N_contigs; A++ ) {
    if ( tree[A].size() != 1 ) continue;
    
    // Find the nodes B,C,D.  Verify that the structure in this area is appropriate for a thorn.
    int B = tree[A][0];
    if ( tree[B].size() != 3 ) continue;
    int C = -1, D = -1;
    for ( int i = 0; i < 3; i++ ) {
      if ( tree[B][i] == A ) continue;
      ( C == -1 ? C : D ) = tree[B][i];
    }
    if ( tree[C].size() == 1 || tree[D].size() == 1 ) continue;
    
    
    // Consider the net change that would happen from smoothing this thorn in either of the two possible ways (designated "AC" and "AD".)
    double N_links_AC = LinkDensity(A,C) - LinkDensity(B,C);
    double N_links_AD = LinkDensity(A,D) - LinkDensity(B,D);
    
    // Pick which change, if either, to apply.  0 = neither; 1 = AC; 2 = AD
    int choice = 0;
    
    // Most of the time, the change will be negative (since the graph was created with the thorn).  If a change is non-negative, go with it.
    if ( N_links_AC >= 0 || N_links_AD >= 0 )
      choice = ( N_links_AC > N_links_AD ? 1 : 2 );
    
    // If both changes would have a net negative effect, take one of them if its negative effect is much less than the other's.
    // If the fold difference in their effects is less than <aggression>, leave the thorn as is.
    else {
      double change_ratio = N_links_AC / N_links_AD;
      if      ( change_ratio < 1/aggression ) choice = 1;
      else if ( change_ratio >   aggression ) choice = 2;
    }
    
    //cout << "THORNY SITUATION:\t(C,B,D) = (" << C << "," << B << "," << D << ")\tA = " << A << "\tN_links_AC = " << N_links_AC << "\tN_links_AD = " << N_links_AD << "\tCHANGE RATIO = " << (N_links_AC/N_links_AD) << "\tCHOICE = " << choice << endl;
    
    
    // Apply the change!
    if ( choice != 0 ) {
      N_smoothed++;
      
      int CorD = choice == 1 ? C : D;
      tree[B]   .erase( find( tree[B]   .begin(), tree[B]   .end(), CorD ) );
      tree[CorD].erase( find( tree[CorD].begin(), tree[CorD].end(), B ) );
      tree[A].push_back(CorD);
      tree[CorD].push_back(A);
    }
  }
  
  
  cout << "N thorns smoothed: " << N_smoothed << endl;
}






// TreeTrunk: Find the "trunk", the longest path within the tree.  This trunk ContigOrdering is used as a scaffold, and contigs not in the trunk will be
// added to it later.
// This is done via a two-step algorithm, as outlined here: http://stackoverflow.com/questions/13187404/linear-time-algorithm-for-longest-path-in-tree
// This function is prepatory to ReinsertShreds().
ContigOrdering
ChromLinkMatrix::TreeTrunk( const vector< vector<int> > & tree, const bool verbose ) const
{
  assert( (int) tree.size() == _N_contigs );
  vector<int> longest_path = find_longest_path( tree );
  if ( _N_contigs == 1 ) longest_path.push_back(0); // handle the trivial case
  ContigOrdering trunk( _N_contigs, longest_path );
  
  if ( verbose ) {
    cout << "Trunk:\t";
    ReportOrderingSize( trunk );
  }
  
  return trunk;
}








// Convert the spanning tree - which is a graph-theory tree of the nodes (contigs) - to a linear ordering of the contigs.
// Method:
// 1. Identify the "trunk", the longest path within the tree.  This should hopefully comprise most of the nodes.  (TreeTrunk())
// 2. "Prune" the tree down to the longest path, by creating a ContigOrdering consisting of just the nodes in the longest path.
//    This will leave out some "pruned" contigs. (TreeTrunk())
// 3. Repeat this pruning and shred the tree into a set of linear paths, by repeatedly finding the longest path and then pruning that path.
// 4. Reinsert the linear "shreds" into the ContigOrdering, each time simply finding the optimal point at which to insert them.
ContigOrdering
ChromLinkMatrix::ReinsertShreds( vector< vector<int> > & tree, const int min_N_REs, const bool use_CP_score ) const
{
  bool verbose = false;
  assert ( !_contig_RE_sites.empty() );
  
  cout << Time() << ": ReinsertShreds with use_CP_score = " << boolalpha << use_CP_score << ", CP_score_dist = " << _CP_score_dist << endl;
  
  assert( _N_contigs == (int) tree.size() );
  
  // Handle a degenerate case that sometimes produces odd results.
  if ( _N_contigs == 1 ) return ContigOrdering( 1, true );
  
  
  vector<int> longest_path = find_longest_path( tree ); // this is redundant but oh well
  ContigOrdering order( _N_contigs, longest_path );
  
  
  // 3. Repeat this pruning and shred the tree into a set of linear paths, by repeatedly finding the longest path and then pruning that path.
  vector< vector<int> > shreds;
  int n_nodes_left = _N_contigs;
  vector<bool> nodes_left( _N_contigs, true );
  
  while ( 1 ) {
    
    // Excise the longest path from the graph, along with any links to other nodes.
    for ( size_t i = 0; i < longest_path.size(); ++i ) {
      int node = longest_path[i];
      
      assert( nodes_left[node] );
      nodes_left[node] = false;
      n_nodes_left--;
      
      for ( size_t j = 0; j < tree[node].size(); ++j ) {
	vector<int> & j_adjs = tree[ tree[node][j] ];
	j_adjs.erase( remove( j_adjs.begin(), j_adjs.end(), node ), j_adjs.end() );
      }
      tree[node].clear();
    }
    
    
    // Find one of the longest remaining paths in the reduced graph.
    // Technically the longest_path function is guaranteed to find one of the longest paths (in case of a tie) in one of the components of the graph.
    // If longest_path returns an empty path, there are no non-trivial paths left, and we're done.
    longest_path = find_longest_path( tree );
    if ( longest_path.empty() ) break;
    shreds.push_back( longest_path );
    //cout << "Next shred:";
    //for ( size_t i = 0; i < longest_path.size(); i++ )
    //  cout << '\t' << longest_path[i];
    //cout << endl;
  }
  
  
  
  // Examine the singleton shreds.  Every contig with no data should be in a singleton shred, but there may also be other singleton shreds.
  vector<bool> contigs_used = ContigsUsed(false);
  for ( int i = 0; i < _N_contigs; i++ ) {
    if ( !contigs_used[i] ) assert( nodes_left[i] );
    if ( nodes_left[i] && contigs_used[i] ) shreds.push_back( vector<int>( 1, i ) );
  }
  
  cout << "Singleton shreds (" << n_nodes_left << "):";
  for ( size_t i = 0; i < nodes_left.size(); i++ )
    if ( nodes_left[i] ) {
      if ( contigs_used[i] ) cout << '\t' << i;
      else cout << "\t(" << i << ")";
    }
  cout << endl;
  
  
  // 4. Reinsert the linear "shreds" into the ContigOrdering, each time simply finding the optimal point at which to insert them.
  // Note that we are reinserting them in decreasing order of their size.
  cout << Time() << ": Reinserting " << shreds.size() << " shreds" << endl;
  
  for ( size_t i = 0; i < shreds.size(); i++ ) {
    vector<int> shred = shreds[i];
    int shred_start = shred[0];
    int shred_end   = shred.back();
    if ( verbose ) { cout << Time() << ": Adding shred #" << i << ":"; for ( size_t j=0; j < shred.size(); j++ ) cout << '\t' << shred[j]; cout << endl; }
    
    
    int shred_len = 0;
    for ( size_t j = 0; j < shred.size(); j++ )
      shred_len += ( _contig_RE_sites.empty() ? _contig_lengths[ shred[j] ] : _contig_RE_sites[ shred[j] ] );
    if ( verbose ) cout << "Shred length: " << shred_len << '\t';
    if ( shred_len > min_N_REs ) {
      if ( verbose ) cout << "PASS!\n";
    }
    else {
      if ( verbose ) cout << "FAIL!\n";
      continue;
    }
    
    
    // Consider all possible positions and orientations in which to insert this shred.
    // Find the position with the most data (immediate links) in support of it.
    double best_N_links = 0;
    double best_score = 0;
    int best_j = -1;
    bool best_rc = false;
    
    for ( int j = 0; j <= order.N_contigs_used(); j++ ) {
      
      for ( int rc = 0; rc < 2; rc++ ) {
	if ( rc && shred_start == shred_end ) continue; // no need to reverse singletons
	if ( verbose ) cout << "Adding at position " << j << " with " << ( rc ? "RC" : "FW" ) << " orientation" << endl;
	
	
	if ( use_CP_score ) {
	  ContigOrdering order_plus = order;
	  order_plus.AddContigs( shred, j );
	  if ( rc ) order_plus.Invert( j, j + shred.size() - 1 );
	  double score;
	  score = OrderingScore( order_plus, true );
	  if ( score > best_score ) {
	    best_score = score;
	    best_j = j;
	    best_rc = bool(rc);
	  }
	}
	
	else {
	  
	  // Count the NET number of Hi-C links in the adjacencies created by inserting this shred here.
	  // Note the edge cases: Inserting the shred at the very beginning or ending only creates one adjacency, but doesn't destroy any adjacencies.
	  double N_links = 0;
	  int cID_prev  = j-1 < 0                       ? -1 : order.contig_ID(j-1);
	  int cID_next  = j   >= order.N_contigs_used() ? -1 : order.contig_ID(j);
	  if (rc) { int swap = cID_prev; cID_prev  = cID_next; cID_next  = swap; }
	  
	  if ( cID_prev  != -1 )                    N_links += LinkDensity( cID_prev, shred_start );
	  if ( cID_next  != -1 )                    N_links += LinkDensity( cID_next, shred_end );
	  if ( cID_prev  != -1 && cID_next  != -1 ) N_links -= LinkDensity( cID_prev, cID_next );
	  
	  
	  if ( N_links > best_N_links ) {
	    best_N_links = N_links;
	    best_j = j;
	    best_rc = bool(rc);
	  }
	  
	  if ( verbose ) cout << "Total (normalized) number of net links supporting this mofo: " << N_links << endl;
	}
      }
    }
    
    
    if ( verbose ) cout << "Best stuff: N_links = " << ( use_CP_score ? best_score : best_N_links ) << " from j=" << best_j << ", rc=" << int(best_rc) << endl;
    
    if ( best_rc ) reverse( shred.begin(), shred.end() );
    order.AddContigs( shred, best_j );
  }
  
  
  
  return order;
}



// OrientContigs: Determine the proper orientation for the contigs in this ContigOrdering (without rearranging the order) and flip contigs accordingly.
// Method: Build a WDAG (Weighed Directed Acyclic Graph) representing all possible ways to orient contigs, then find the best path through the WDAG.
void
ChromLinkMatrix::OrientContigs( ContigOrdering & order ) const
{
  cout << Time() << ": OrientContigs" << endl;
  
  // Build a WDAG (Weighed Directed Acyclic Graph) representing all possible ways to orient contigs in this ContigOrdering
  WDAG wdag = order.OrientationWDAG( this );
  
  // Compute the highest-weight path on this WDAG.
  wdag.FindBestPath();
  //wdag.WriteToFile("wdag.txt");
  vector<int> best_node_IDs = wdag.BestNodeIDs();
  assert( (int) best_node_IDs.size() == order.N_contigs_used() + 2 ); // path includes start,end nodes
  
  // Use the node IDs of the highest-weight path to determine which contigs should be flipped.
  // Due to the node ID numbering, an even node ID indicates that a contig should be flipped.
  for ( int i = 1; i+1 < (int) best_node_IDs.size(); i++ ) {
    int node_ID = best_node_IDs[i];
    assert( node_ID == 2*i-1 || node_ID == 2*i );
    if ( node_ID == 2*i ) order.Invert(i-1);
  }
  
  
  // Calculate quality scores for each contig's orientation.  The quality score is defined as the relative likelihood that the contig belongs in its chosen
  // orientation rather than the opposite, given its neighbors' orientations and the links it shares with them.
  for ( int i = 0; i < order.N_contigs_used(); i++ ) {
    
    int  ID = order.contig_ID(i);
    bool rc = order.contig_rc(i);
    
    // Calculate the log-likelihood of the data, given the orientation as we see it, and the log-likelihood of the data, if this contig were flipped.
    // We must be careful to handle edge cases.
    double LL = 0, LL_alt = 0;
    if ( i > 0 ) {
      LL     += ContigOrientLogLikelihood( order.contig_ID(i-1), order.contig_rc(i-1), ID,  rc );
      LL_alt += ContigOrientLogLikelihood( order.contig_ID(i-1), order.contig_rc(i-1), ID, !rc );
    }
    if ( i+1 < order.N_contigs_used() ) {
      LL     += ContigOrientLogLikelihood( ID,  rc, order.contig_ID(i+1), order.contig_rc(i+1) );
      LL_alt += ContigOrientLogLikelihood( ID, !rc, order.contig_ID(i+1), order.contig_rc(i+1) );
    }
    
    // The difference between log-likelihoods describes how much statistical confidence is behind our call of this contig's orientation.
    double diff = LL - LL_alt;
    assert( diff >= 0 ); // if this fails, the WDAG hasn't done its job properly - there's a bug somewhere
    
    //double ratio = LL_alt / LL;
    //PRINT5( i, LL, LL_alt, diff, ratio );
    
    // Load the quality score into the ContigOrdering.
    order.AddOrientQ( i, diff );
    
  }
  
}






// Initialize the matrix of data and allocate memory for it.
void
ChromLinkMatrix::InitMatrix()
{
  assert( !_matrix_init );
  assert( _N_contigs > 0 );
  
  int N_bins = 2 * _N_contigs; // each contig gets 2 bins: one for each of its possible orientations
  
  // Initialize the matrix of data.  If the matrices are big, the machine may run out of memory, so catch that error.
  try {
    _matrix = new vector<int> * [N_bins];
    for ( int i = 0; i < N_bins; i++ ) {
      _matrix[i] = new vector<int>[N_bins];
      for ( int j = 0; j < N_bins; j++ )
	_matrix[i][j].clear();
    }
  }
  catch ( std::bad_alloc & ba ) {
    cerr << "ERROR: Sorry, there's not enough memory to allocate for this ChromLinkMatrix!  Try running on a machine with more RAM,\nbad_alloc error message: " << ba.what() << endl;
    exit(1);
  }
  
  
  _matrix_init = true;
}




 

// De-allocate memory for the matrix of data.
void
ChromLinkMatrix::FreeMatrix()
{
  assert( _matrix_init );
  _matrix_init = false;
  
  for ( int i = 0; i < 2 * _N_contigs; i++ )
    delete[] _matrix[i];
  delete[] _matrix;
}




// LoadRESitesFile: Fill _contig_RE_sites.
void
ChromLinkMatrix::LoadRESitesFile( const string & RE_sites_file )
{
  cout << Time() << ": Loading contig RE lengths for use in normalization <-\t" << RE_sites_file << endl;
  assert ( DeNovo() );
  
  if ( RE_sites_file == "." ) {
    cerr << "WARNING: Tried to load an RE_sites_file with the name `.'.  Is this taken from a CLM file with no RE_sites_file listed?  Can't normalize to RE sites." << endl;
    return;
  }
  
  // Get token 1 (zero-indexed) from each line of the lengths file.
  _contig_RE_sites = ParseTabDelimFile<int>( RE_sites_file, 0 );
  
  // If this assert fails, the input RE_sites_file is inconsistent with the original dataset (wrong number of contigs).
  assert( (int) _contig_RE_sites.size() == _N_contigs );
  
  // Add 1 to all length to prevent dividing by 0 (some short assembly contigs have no restriction sites)
  for ( int i = 0; i < _N_contigs; i++ )
    _contig_RE_sites[i]++;
  
  // Find the most REs in one contig.  Note that this occurs after adding 1 to lengths.
  _most_contig_REs = *( max_element( _contig_RE_sites.begin(), _contig_RE_sites.end() ) );
}







// AddLinkToMatrix: Add a individual Hi-C link to the matrix.  This function is only used when loading data from SAM files.
void
ChromLinkMatrix::AddLinkToMatrix( const int contig1, const int contig2, const int read1_dist1, const int read1_dist2, const int read2_dist1, const int read2_dist2 )
{
  // Set a minimum offset.  This prevents distances from being equal to 0, which is good because that leads to division-by-0 errors in OrderingScore().
  // To be pedantically accurate, we would want to set OFFSET to twice the read length, but the read length currently isn't tracked.
  static const int OFFSET = 100;
  
  
  /* In this ASCII illustration, read1_dist1 = 8, read1_dist2 = 2, read2_dist1 = 1, read2_dist2 = 9.
   *
   * fw, fw:          R----R               fw, rc:          R------------R
   *          ==========> ==========>               ==========> <==========
   *            contig1     contig2                   contig1     2gitnoc
   *
   *
   * rc, fw:    R----------R               rc, rc:    R------------------R
   *          <========== ==========>               <========== <==========
   *            1gitnoc     contig2                   1gitnoc     2gitnoc
   *
   * Find the distance between these reads, assuming the contigs are immediately adjacent, under any of their four possible orientations:
   *
   *************************************************************************************************/
  int dist_fw_fw = read1_dist2 + read2_dist1 + OFFSET;
  int dist_fw_rc = read1_dist2 + read2_dist2 + OFFSET;
  int dist_rc_fw = read1_dist1 + read2_dist1 + OFFSET;
  int dist_rc_rc = read1_dist1 + read2_dist2 + OFFSET;
  
  
  // Record all four distances in the 2-D matrix object for this ChromLinkMatrix.
  // The conversion from contig ID to bin ID is: bin ID = 2 * contig ID + (1 if rc)
  // Note that the matrix is filled in a symmetric fashion, and should always remain symmetric.
  _matrix[2*contig1  ][2*contig2  ].push_back( dist_fw_fw );
  _matrix[2*contig1  ][2*contig2+1].push_back( dist_fw_rc );
  _matrix[2*contig1+1][2*contig2  ].push_back( dist_rc_fw );
  _matrix[2*contig1+1][2*contig2+1].push_back( dist_rc_rc );
  _matrix[2*contig2+1][2*contig1+1].push_back( dist_fw_fw );
  _matrix[2*contig2  ][2*contig1+1].push_back( dist_fw_rc );
  _matrix[2*contig2+1][2*contig1  ].push_back( dist_rc_fw );
  _matrix[2*contig2  ][2*contig1  ].push_back( dist_rc_rc );
  
}



// CalculateRepeatFactors: Fill the _repeat_factors vector.  This vector contains the 'factor', or multiplicity, of each contig: the ratio by which the total
// density of links in each contig exceeds the density expected by chance.  Used in normalization.  Run this after all the links are loaded in.
void
ChromLinkMatrix::CalculateRepeatFactors()
{
  return; // TEMP: not used at the moment, buggy in fragScaff cases, screw it
  cout << Time() << ": CalculateRepeatFactors for normalization" << endl;
  assert( DeNovo() );
  
  _repeat_factors.clear();
  _repeat_factors.resize( _N_contigs );
    
  // Find the number of Hi-C links on each contig.  This is as simple as adding rows/columns in the matrix.  Also calculate the total sum.
  int64_t N_links_total = 0;
  vector<int64_t> N_links( _N_contigs, 0 );
  
  for ( int i = 0; i < _N_contigs; i++ )
    assert( _contig_RE_sites[i] > 0 ); // not sure why this would fail
  
  for ( int i = 0; i < _N_contigs; i++ )
    for ( int j = 0; j < _N_contigs; j++ )
      if ( !_matrix[2*i][2*j].empty() ) {
	// TODO: there's a floating point exception in here somewhere
	int64_t N_links_norm = _matrix[2*i][2*j].size();
	N_links_norm = N_links_norm * _most_contig_REs / _contig_RE_sites[i] * _most_contig_REs / _contig_RE_sites[j]; // normalize
	//N_links_norm = N_links_norm * _longest_contig / _contig_lengths[i] * _longest_contig / _contig_lengths[j]; // normalize
	N_links_total += N_links_norm;
	N_links[i]    += N_links_norm;
	N_links[j]    += N_links_norm;
      }
  
  // Calculate the 'average' number of links.
  double N_links_avg = 2.0 * N_links_total / _N_contigs;
  
  for ( int i = 0; i < _N_contigs; i++ )
    _repeat_factors[i] = N_links[i] / N_links_avg;
  
  
  // TEMP: Write to file so that Reporter::EvalGapSizes can access these.
  cout << Time() << ": Writing to RFs.txt" << endl;
  ofstream out( "RFs.txt", ios::out );
  for ( int i = 0; i < _N_contigs; i++ )
    out << _repeat_factors[i] << endl;
  out.close();
}


// LinkDensity: Return the number of Hi-C links connecting these two contigs (normalized to contig lengths, if this is a de novo CLM.)
// This function assumes, and does not check, that contig1,contig2 are in the range [0,_N_contigs).
double
ChromLinkMatrix::LinkDensity( const int contig1, const int contig2 ) const
{
  double N_links = _matrix[ 2*contig1 ][ 2*contig2 ].size();
  if ( !DeNovo() ) return N_links;
  
  //N_links /= ( _repeat_factors[contig1] * _repeat_factors[contig2] );
  //return N_links; // TEMP: should already be effectively normalized by contig length/RE sites thru repeat_factor
  
  // Normalize by contig length, 
  const bool use_REs = !_contig_RE_sites.empty();
  if ( use_REs )
    return N_links * _most_contig_REs / _contig_RE_sites[contig1] * _most_contig_REs / _contig_RE_sites[contig2];
  else
    return N_links * _longest_contig / _contig_lengths[contig1] * _longest_contig / _contig_lengths[contig2];
}






// PlotTree: Use grpahviz to print a spanning tree to a graph image at ~/public_html/<filename>.
// Note that graphviz is very slow for large graphs, and the output images themselves are sometimes so large as to cause memory problems.
// Hence this is NOT RECOMMENDED for graphs with over 500 vertices.
void
ChromLinkMatrix::PlotTree( const vector< vector<int> > & tree, const string & filename ) const
{
  assert( (int) tree.size() == _N_contigs );
  cout << Time() << ": PlotTree: Drawing an image of this spanning tree at ~/public_html/" << filename << endl;
  
  
  if ( _N_contigs >= 500 ) cerr << "WARNING: Called ChromLinkMatrix::PlotTree on a large tree of size " << _N_contigs << "; this may take a while, and the output file " << filename << " will take a long time to open" << endl;
  if ( filename.substr( filename.size() - 4, 4 ) != ".png" ) cerr << "WARNING: Called ChromLinkMatrix::PlotTree to make a file called " << filename << "; filename should end in `.png'" << endl;
       
  // Make a dot file and then an image of the spanning tree.
  ofstream out( ( filename + ".dot" ).c_str(), ios::out );
  out << "graph G {\n";
  // Only list edges, not vertices.  This has the effect of leaving out vertices with no edges.
  for ( int i = 0; i != _N_contigs; ++i)
    for ( size_t j = 0; j != tree[i].size(); ++j)
      if ( tree[i][j] > i )
	out << i << "--" << tree[i][j] << ";\n";
  out << "}\n";
  out.close();
  
  string cmd = "dot -Tpng " + filename + ".dot > ~/public_html/" + filename;
  system( cmd.c_str() );
}



// ReportOrderingSize: Report about the number and length of contigs in this ContigOrdering, compared to the total set of contigs in this CLM.
void
ChromLinkMatrix::ReportOrderingSize( const ContigOrdering & order ) const
{
  int64_t order_length = 0, total_contig_length;
  
  if ( DeNovo() ) {
    total_contig_length = std::accumulate( _contig_lengths.begin(), _contig_lengths.end(), int64_t(0) );
    for ( int i = 0; i < order.N_contigs_used(); i++ )
      order_length += _contig_lengths[ order.contig_ID(i) ];
  }
  else {
    total_contig_length = _contig_size * order.N_contigs();
    order_length        = _contig_size * order.N_contigs_used();
  }
  
  double pct_used   = 100.0 * order.N_contigs_used() / _N_contigs;
  double pct_length = 100.0 * order_length / total_contig_length;
  
  cout << "Order contains " << order.N_contigs_used() << " of " << _N_contigs << " contigs (" << pct_used << "%), with a total length of " << order_length << " of " << total_contig_length << " (" << pct_length << "%)." << endl;
}








// LoadDeNovoCLMsFromSAM: Import one or more SAM/BAM files and create a set of de novo ChromLinkMatrices corresponding to a ClusterVec (i.e., a set of contig
// clusters that have been derived by Lachesis' clustering algorithm; see GenomeLinkMatrix).
// As many or as few of the ChromLinkMatrix pointers can be non-NULL.  Whichever ones are non-NULL will be assumed to represent ChromLinkMatrices for the
// cluster ID equal to their index in the vector, and will be filled accordingly.
// If you are creating a set of ChromLinkMatrices for each chromosome, this is much faster than calling LoadFromSAMDeNovo individually for each ChromLinkMatrix
// object because it only reads through the SAM file(s) once.
void
LoadDeNovoCLMsFromSAM( const string & SAM_file, const string & RE_sites_file, const ClusterVec & clusters, vector<ChromLinkMatrix *> CLMs )
{
  assert( CLMs.size() == clusters.size() );
  int N_clusters = clusters.size();
  
  bool verbose = true;
  
  
  // Find the lengths of all of the de novo contigs.  This tells us how many de novo contigs there are in the total dataset.
  vector<int> contig_lengths_orig  = TargetLengths( SAM_file );
  vector<int> contig_RE_sites_orig = ParseTabDelimFile<int>( RE_sites_file, 1 );
  int N_contigs_total = contig_lengths_orig.size();
  
  // For each ChromLinkMatrix that we're actually creating, make a lookup table of the lengths and number of RE sites in the contigs in this cluster.
  vector<int> used;
  for ( int i = 0; i < N_clusters; i++ ) {
    if ( CLMs[i] == NULL ) continue;
    used.push_back(i);
    
    CLMs[i]->_contig_lengths .clear();
    CLMs[i]->_contig_RE_sites.clear();
    for ( set<int>::const_iterator it = clusters[i].begin(); it != clusters[i].end(); ++it ) {
      //cout << "cluster #" << i << ", contig global#" << *it << ": length = " << contig_lengths_orig[*it] << endl;
      CLMs[i]->_contig_lengths .push_back( contig_lengths_orig [*it] );
      CLMs[i]->_contig_RE_sites.push_back( contig_RE_sites_orig[*it] );
    }
    CLMs[i]->FindLongestContig();
    CLMs[i]->_SAM_files.push_back( SAM_file ); // keep track of which SAM files were used to generate this matrix
  }
  
  cout << Time() << ": Filling "
       << ( used.size() == 1 ? "cluster " + boost::lexical_cast<string>( used[0] ) : boost::lexical_cast<string>( used.size() ) + " clusters" )
       << " with Hi-C data from SAM file " << SAM_file
       << ( verbose ? "\t(dot = 1M alignments)" : "" ) << endl;
  
  
  // Convert the clusters vector into two lookup tables, which map contig IDs to cluster IDs, and also onto "local" contig indices for the cluster.
  // For example, if there are 8 contigs in two clusters: {0,2,5,6} and {1,3,7}, then the cluster_IDs table will look like this: [ 0, 1, 0, 1, -1, 0, 0, 1 ].
  // The local_cIDs lookup table will look like this for cluster #0: [0, -1, 1, -1, -1, 2, 3, -1] and like this for cluster #1: [-1, 0, -1, 1, -1, -1, -1, 2].
  vector<int> cluster_IDs( N_contigs_total, -1 );
  vector<int>  local_cIDs( N_contigs_total, -1 );
  vector<int> ID( N_clusters, 0 );
  for ( int i = 0; i < N_clusters; i++ ) {
    for ( set<int>::const_iterator it = clusters[i].begin(); it != clusters[i].end(); ++it ) {
      cluster_IDs[*it] = i;
      local_cIDs [*it] = ID[i]++;
    }
  }
  
  
  /*
  for ( int i = 0; i < N_contigs_total; i++ )
    PRINT3( i, cluster_IDs[i], local_cIDs[i] );
  for ( int i = 0; i < N_clusters; i++ )
    PRINT2( i, ID[i] );
  for ( int i = 0; i < N_clusters; i++ )
    PRINT( CLMs[i] );
  */
  
  
  
  int N_pairs_used = 0;
  
  // Set up a SAMStepper object to read in the alignments.
  SAMStepper stepper(SAM_file);
  stepper.FilterAlignedPairs(); // Only look at read pairs where both reads aligned to the assembly.
  
  // Loop over all pairs of alignments in the SAM file.
  // Note that the next_pair() function assumes that all reads in a SAM file are paired, and the two reads in a pair occur in consecutive order.
  for ( pair< bam1_t *, bam1_t *> aligns = stepper.next_pair(); aligns.first != NULL; aligns = stepper.next_pair() ) {
    
    if ( verbose && stepper.N_aligns_read() % 1000000 == 0 ) cout << "." << flush;
    
    const bam1_core_t & c1 = aligns.first->core;
    const bam1_core_t & c2 = aligns.second->core;
    
    // Ignore reads with mapping quality 0.  (In the GSM862723 dataset this is roughly 4% of reads.)
    if ( c1.qual == 0 || c2.qual == 0 ) continue;
    
    // Sanity checks to make sure the read pairs appear as they should in the SAM file.
    //cout << stepper.as_SAM_line(aligns.first) << endl << stepper.as_SAM_line(aligns.second) << endl << endl;
    assert( c1.tid == c2.mtid );
    assert( c2.tid == c1.mtid );
    assert( c1.pos == c2.mpos );
    assert( c2.pos == c1.mpos );
    assert( c1.tid < N_contigs_total );
    assert( c2.tid < N_contigs_total );
    
    // Find what contigs these reads map to, and which clusters they belong to.
    // We only care about this link if both reads align, to contigs in the same cluster, and that cluster is one of the non-NULL ones that we're building.
    int cluster  = cluster_IDs[c1.tid];
    int cluster2 = cluster_IDs[c2.tid];
    if ( cluster  == -1 ) continue;
    if ( cluster  != cluster2 ) continue;
    if ( CLMs[cluster] == NULL ) continue;
    
    // If the two reads align to the exact same contig, the link isn't informative, so skip it.
    if ( c1.tid == c2.tid ) continue;
    
    // For each read, find the distance to either end of its contig.
    int read1_dist1 = c1.pos;
    int read2_dist1 = c2.pos;
    int read1_dist2 = contig_lengths_orig[c1.tid] - c1.pos;
    int read2_dist2 = contig_lengths_orig[c2.tid] - c2.pos;
    assert( read1_dist2 >= 0 );
    assert( read2_dist2 >= 0 );
    
    // Add the link to the matrix object.
    //cout << "Adding link to cluster #" << cluster << endl;
    CLMs[cluster]->AddLinkToMatrix( local_cIDs[c1.tid], local_cIDs[c2.tid], read1_dist1, read1_dist2, read2_dist1, read2_dist2 );
    
    N_pairs_used++;
  }
  
  if ( verbose ) cout << endl;
  
  
  for ( int i = 0; i < N_clusters; i++ )
    CLMs[i]->CalculateRepeatFactors();
}




// LoadNonDeNovoCLMsFromSAM: Import one or more SAM/BAM files and fill a set of non-de novo ChromLinkMatrices corresponding to each chromosome.
// As many or as few of the ChromLinkMatrix pointers can be non-NULL.  Whichever ones are non-NULL will be assumed to represent LinkMatrices with chrID equal
// to their index in the vector, and will be filled accordingly.
// If you are creating a set of LinkMatrices for each chromosome, this is much faster than calling LoadFromSAMNonDeNovo individually for each ChromLinkMatrix
// object because it only reads through the SAM file(s) once.
void
LoadNonDeNovoCLMsFromSAM( const string & SAM_file, vector<ChromLinkMatrix *> CLMs )
{
  const bool verbose = true; // verbosely read in the file
  
  cout << Time() << ": Reading Hi-C data from SAM file " << SAM_file << (verbose ? "\t(dot = 1M alignments)" : "" ) << endl;
  assert( boost::filesystem::is_regular_file( SAM_file ) );
  
  int N_chroms = NTargetsInSAM( SAM_file );
  assert( N_chroms == (int) CLMs.size() ); // if this fails, the CLMs vector doesn't have a number of elements equal to the number of chromosomes
  
  // Gather up stats on each ChromLinkMatrix object for local use.
  bool has_data = false;
  vector<int> N_contigs   ( N_chroms, -1 );
  vector<int> contig_sizes( N_chroms, -1 );
  for ( int chrID = 0; chrID < N_chroms; chrID++ )
    if ( CLMs[chrID] ) {
      N_contigs   [chrID] = CLMs[chrID]->_N_contigs;
      contig_sizes[chrID] = CLMs[chrID]->_contig_size;
      CLMs[chrID]->_SAM_files.push_back( SAM_file ); // keep track of which SAM files were used to generate this matrix
      
      has_data = true;
    }
  
  
  if (!has_data) {
    cout << "WARNING: Called LoadNonDeNovoCLMsFromSAM with no non-NULL LinkMatrix objects.  Nothing will happen.  Skipping." << endl;
    return;
  }
  
  
  
  
  int N_pairs_used = 0;
  
  // Set up a SAMStepper object to read in the alignments.
  SAMStepper stepper(SAM_file);
  stepper.FilterAlignedPairs(); // Only look at read pairs where both reads aligned to the reference.
  
  // Loop over all pairs of alignments in the SAM file.
  // Note that the next_pair() function assumes that all reads in a SAM file are paired, and the two reads in a pair occur in consecutive order.
  for ( pair< bam1_t *, bam1_t *> aligns = stepper.next_pair(); aligns.first != NULL; aligns = stepper.next_pair() ) {
    
    if ( verbose && stepper.N_aligns_read() % 1000000 == 0 ) cout << "." << flush;
    
    const bam1_core_t & c1 = aligns.first->core;
    const bam1_core_t & c2 = aligns.second->core;
    
    // Sanity checks to make sure the read pairs appear as they should in the SAM file.
    assert( c1.tid == c2.mtid );
    assert( c2.tid == c1.mtid );
    assert( c1.pos == c2.mpos );
    assert( c2.pos == c1.mpos );
    
    // Ignore the oddly large number of cases where both reads map to the same location.  (Maybe these are cases in which really only one read mapped?  Dunno.)
    if ( c1.pos == c2.pos ) continue;
    // Only select one half of each read pair.  This breaks symmetry and therefore has a small stochastic effect in the odd situations of unmatched read pairs.
    if ( c1.pos >  c2.pos ) continue;
    
    // Ignore reads with mapping quality 0.  (In the GSM862723 dataset this is roughly 4% of reads.)
    if ( c1.qual == 0 || c2.qual == 0 ) continue;
    
    if ( c1.tid != c2.tid ) continue; // skip interchromosomal links
    if ( c1.tid >= N_chroms ) continue; // skip links on non-canonical chromosomes
    if ( !CLMs[c1.tid] ) continue; // skip links on chromosomes for which we have no ChromLinkMatrix object
    
    // Find the contig IDs for these reads.  This determines in what bins the distance between them will be stored.
    int contig_size = contig_sizes[c1.tid];
    int contig1 = c1.pos / contig_size;
    int contig2 = c2.pos / contig_size;
    
    // Ignore cases in which both reads map to the same contig, because these are uninformative for scaffolding.
    if ( contig1 == contig2 ) continue;
    
    assert( contig1 < N_contigs[c1.tid] );
    assert( contig2 < N_contigs[c1.tid] );
    
    // For each read, find the distance to either end of its contig.
    int read1_dist1 = c1.pos - contig1 * contig_size;
    int read2_dist1 = c2.pos - contig2 * contig_size;
    int read1_dist2 = (contig1+1) * contig_size - c1.pos;
    int read2_dist2 = (contig2+1) * contig_size - c2.pos;
    
    // Add the link to the matrix object of the appropriate ChromLinkMatrix.
    CLMs[c1.tid]->AddLinkToMatrix( contig1, contig2, read1_dist1, read1_dist2, read2_dist1, read2_dist2 );
    
    N_pairs_used++;
    
  }
  
  if ( verbose ) cout << endl;
  
  cout << Time() << ": Done with " << SAM_file << "!  N aligns/pairs read: " << stepper.N_aligns_read() << "/" << stepper.N_pairs_read() << "; N pairs used: " << N_pairs_used << endl;
  
  for ( int i = 0; i < N_chroms; i++ )
    CLMs[i]->CalculateRepeatFactors();
}





// These multi-SAM-file functions are wrappers for their respective one-SAM-file functions.
void
LoadDeNovoCLMsFromSAM( const vector<string> & SAM_files, const string & RE_sites_file, const ClusterVec & clusters, vector<ChromLinkMatrix *> CLMs )
{
  AssertFilesExist( SAM_files );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadDeNovoCLMsFromSAM( SAM_files[i], RE_sites_file, clusters, CLMs );
}


void
LoadNonDeNovoCLMsFromSAM( const vector<string> & SAM_files, vector<ChromLinkMatrix *> CLMs )
{
  AssertFilesExist( SAM_files );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadNonDeNovoCLMsFromSAM( SAM_files[i], CLMs );
}

