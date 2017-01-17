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


// For documentation, see GenomeLinkMatrix.h
#include "GenomeLinkMatrix.h"
#include "TrueMapping.h"
#include "TextFileParsers.h" // ParseTabDelimFile

#include <sys/time.h> // struct timeval, gettimeofday
#include <assert.h>
#include <set>
#include <map> // map, multimap
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip> // boolalpha
#include <algorithm> // max_element
#include <numeric> // accumulate

// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp> // mapped_matrix, compressed_matrix


// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/HumanGenome.h"
#include "gtools/SAMStepper.h" // NTargetsInSAM(), TargetLengths()




// Options for reporting validation & dotplots.
const bool exclude_noninformative_contigs_from_dotplot = true;
//const bool exclude_low_Q_scores = false; // don't use contigs with insufficiently high quality scores of their mapping to reference
const bool cluster_noninformative_into_singletons = true; // allow singleton clusters to become non-singletons by adding noninformative contigs to them (TEMP: do this for yeasts only)


// Print a set in the following format: "[1-3,5,7]".
string
PrintSet( const set<int> & s ) {

  if ( s.empty() ) return "[]";

  // Beginning.
  ostringstream oss;
  oss << '[' << *(s.begin());

  bool start_streak = false;
  bool in_streak = false;
  int prev = -1;

  // Middle loop.
  for ( set<int>::const_iterator it = s.begin(); it != s.end(); it++ ) {

    // Look for streaks of consecutive numbers.
    if ( it != s.begin() ) {

      if ( *it == prev + 1 ) {
	start_streak = !in_streak;
	in_streak = true;
      }
      else {

	// End streak, if there is one.
	if ( in_streak && !start_streak ) oss << '-' << prev;
	start_streak = false;
	in_streak = false;

	oss << ',' << *it;
      }
    }

    prev = *it;
  }


  // End.
  if ( in_streak ) oss << '-' << prev;

  oss << ']';
  return oss.str();
}



int
median_of_set( const set<int> & s )
{
  set<int>::const_iterator it = s.begin();
  for ( size_t i = 0; i < s.size() / 2; i++ ) ++it;
  return *it;
}





// Make a non-de novo GenomeLinkMatrix with no Hi-C link data.  This isn't good for much other than calling NonDeNovoTrueMapping().
// In fact, many other function calls will fail on this GenomeLinkMatrix because _SAM_files is empty.
GenomeLinkMatrix::GenomeLinkMatrix( const string & species, const int bin_size )
  : _species( species ),
    _bin_size( bin_size )
{
  assert( bin_size > 0 );
  _N_bins = 0;
  _SAM_files.clear();
}




// Load a non-de novo GenomeLinkMatrix with a set of SAM files representing human-genome alignments.  Split the GLM into bins of size bin_size.
GenomeLinkMatrix::GenomeLinkMatrix( const vector<string> & SAM_files, const int bin_size )
{
  assert( !SAM_files.empty() );
  assert( bin_size > 0 );
  _bin_size = bin_size; // setting a non-zero bin_size indicates that this is a non-de novo GLM
  _species = "human"; // these SAM files must be human


  cout << "Creating a new GenomeLinkMatrix with bin_size = " << bin_size << " and SAM_files =";
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    cout << "  " << SAM_files[i];
  cout << endl;



  // Determine the number of bins in this GenomeLinkMatrix, by splitting up every chromosome into bins of size _bin_size.
  vector<int> bins_per_chrom = BinsPerChromInHumanGenome();
  _N_bins = std::accumulate( bins_per_chrom.begin(), bins_per_chrom.end(), 0 );


  cout << "Number of bins in entire genome = " << _N_bins << endl;

  // Set the original order to identity (i.e., no rearrangement until ReorderContigsByRef() gets called.)
  _contig_orig_order.clear();
  for ( int i = 0; i < _N_bins; i++ )
    _contig_orig_order.push_back(i);

  // Set a bunch of stuff that remains empty for non-de novo GLMs.
  _RE_sites_file = "";
  _contig_lengths .clear();
  _contig_RE_sites.clear();
  _contig_skip = vector<bool>(_N_bins,false);

  // Create an empty matrix for these bins.
  InitMatrix();


  // Fill the matrix with data from the SAM files.
  LoadFromSAMNonDeNovo( SAM_files );
}





// Load a de novo GenomeLinkMatrix with a set of SAM files representing alignments to contigs.
// Also optionally load a list of the number of restriction sites per contig.  If a file is given, contigs' RE lengths will be used for normalization, instead
// of their lengths in bp.
GenomeLinkMatrix::GenomeLinkMatrix( const string & species, const vector<string> & SAM_files, const string & RE_sites_file )
{
  assert( !SAM_files.empty() );
  _bin_size = 0; // setting this indicates at this a de novo GLM
  _species = species;
  _RE_sites_file = RE_sites_file;


  cout << "Creating a new GenomeLinkMatrix for an assembly in progress, with SAM_files =";
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    cout << "  " << SAM_files[i];
  cout << endl;

  for ( size_t i = 0; i < SAM_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( SAM_files[i] ) );


  // Determine the number of bins in this GenomeLinkMatrix by reading the SAM file headers.  Each target contig is a bin.
  _N_bins = NTargetsInSAM( SAM_files[0] );

  // Set the original order to identity (i.e., no rearrangement until ReorderContigsByRef() gets called.)
  _contig_orig_order.clear();
  for ( int i = 0; i < _N_bins; i++ )
    _contig_orig_order.push_back(i);

  // Set the contig lengths in accordance with the SAM files.
  _contig_lengths = TargetLengths( SAM_files[0] );



  // Set the number of restriction sites per contig in accordance with the input file, if one was given.
  _contig_RE_sites.clear();
  if ( RE_sites_file != "" ) LoadRESitesFile( RE_sites_file );


  _contig_skip = vector<bool>(_N_bins,false);

  cout << "Number of contigs = " << _N_bins << endl;

  // Create an empty matrix for these bins.
  InitMatrix();


  // Fill the matrix with data from the SAM files.
  LoadFromSAMDeNovo( SAM_files );
}








static const unsigned LINE_LEN = 50000;




// ReadFile: Read the data from file GLM_file into this GenomeLinkMatrix.
// The file GLM_file should have been created by a previous call to GenomeLinkMatrix::WriteFile(), and needs to have a commented header as defined there.
void
GenomeLinkMatrix::ReadFile( const string & GLM_file )
{
  cout << "GenomeLinkMatrix::ReadFile   <-  " << GLM_file << endl;
  assert( boost::filesystem::is_regular_file( GLM_file ) );


  // Set initial values that MUST be overwritten later.
  _N_bins = -1;
  _bin_size = -1;
  _species = "";

  char line[LINE_LEN];
  vector<string> tokens;

  ifstream in( GLM_file.c_str(), ios::in );
  bool first = true;

  // Read the file line-by-line.

  while ( 1 ) {
    in.getline( line, LINE_LEN );
    assert( strlen(line)+1 < LINE_LEN );
    if ( in.fail() ) break;

    // If this is a header line, look for useful data in it.  These header lines were generated in GenomeLinkMatrix::WriteFile, below.
    if ( line[0] == '#' ) {
      boost::split( tokens, line, boost::is_any_of(" ") );
      assert( tokens.size() > 3 );

      if ( tokens[1] == "Species" ) // line: "# species = human"
	_species = tokens[3];

      else if ( tokens[1] == "N_bins" ) { // line: "# N_bins = 1"
	_N_bins = boost::lexical_cast<int>( tokens[3] );
	cout << "N_bins = " << _N_bins << endl;
	InitMatrix();
      }

      else if ( tokens[1] == "bin_size" ) // line: "# bin_size = 1"
	_bin_size = boost::lexical_cast<int>( tokens[3] );

      else if ( tokens[1] == "RE_sites_file" ) { // line: "# RE_sites_file = <filename>"
	if ( !boost::filesystem::is_regular_file( tokens[3] ) ) {
	  cout << "ERROR: Trying to load RE sites file that doesn't seem to exist.  You need to create the following file: " << tokens[3] << endl;
	  cout << "Try running `CountMotifsInFasta.pl <ref fasta> <RE motif>`.\n";
	}
	LoadRESitesFile( tokens[3] );
      }

      else if ( tokens[1] == "SAM" ) // line: "# SAM files used in generating this dataset: test.sam"
	for ( size_t i = 8; i < tokens.size(); i++ )
	  _SAM_files.push_back( tokens[i] );
    }

    else if ( line[0] == 'X' ) continue; // skip the header line of the matrix itself

    // If this is not a header line, put data in the matrix.
    // The data may be sparse, and that's ok - the matrix will just contain 0s.  However, there may not be data on the diagonal.
    else {
      if ( first ) { first = false; cout << "Loading matrix data..." << endl; }

      boost::split( tokens, line, boost::is_any_of("\t") );
      int X = boost::lexical_cast<int>( tokens[0] );
      int Y = boost::lexical_cast<int>( tokens[1] );
      int Z = boost::lexical_cast<int>( tokens[2] );
      //assert( X != Y ); // comment out to save runtime
      _matrix(X,Y) = Z;
    }
  }

  in.close();



  assert( _N_bins != -1 );
  assert( _bin_size != -1 );
  assert( _species != "" );
  assert( !_SAM_files.empty() );


  // Set the original order to identity (i.e., no rearrangement until ReorderContigsByRef() gets called.)
  _contig_orig_order.clear();
  for ( int i = 0; i < _N_bins; i++ )
    _contig_orig_order.push_back(i);

  _contig_skip = vector<bool>( _N_bins, false );

  // If this is a de novo GLM, set the contig lengths in accordance with the SAM files.
  if ( _bin_size == 0 )
    _contig_lengths = TargetLengths( _SAM_files[0] );
}






// WriteFile: Write the data in this GenomeLinkMatrix to file GLM_file.
// The output format is a long tall table with (_N_bins^2) rows and three columns: X, Y, Z.  X and Y are bin IDs; Z is the value in the bin.
// It's printed as a "sparse matrix", so that lines with Z=0 are not printed.
// This format is designed for easy input to GenomeLinkMatrix::ReadFile(), or to R and Perl.  (NOTE: R and Perl may not like the new sparse-matrix format.)
void
GenomeLinkMatrix::WriteFile( const string & GLM_file ) const
{
  cout << "GenomeLinkMatrix::WriteFile  ->  " << GLM_file << endl;

  ofstream out( GLM_file.c_str(), ios::out );

  // Print a header to file.  The header is for easier human reading and also contains numbers used by GenomeLinkMatrix::ReadFile().
  out << "# GenomeLinkMatrix file - see GenomeLinkMatrix.h for documentation of this object type" << endl;
  out << "# Species = " << _species << endl;
  out << "# N_bins = " << _N_bins << endl;
  out << "# bin_size = " << _bin_size << endl;
  out << "# RE_sites_file = " << ( _RE_sites_file != "" ? _RE_sites_file : "." ) << endl;
  out << "# SAM files used in generating this dataset:";
  for ( size_t i = 0; i < _SAM_files.size(); i++ )
    out << " " << _SAM_files[i];
  out << endl;

  // Print the table to file.
  out << "X\tY\tZ" << endl;
  for ( int x = 0; x < _N_bins; x++ )
    for ( int y = 0; y < _N_bins; y++ ) {

      // Skip the main diagonal.  We don't need this data in our histogram because we're interested in links *between* bins, not within a bin.
      if ( x == y ) continue;

      // Also skip empty bins.  This makes the matrix sparse.
      if ( _matrix(x,y) == 0 ) continue;
      out << x << '\t' << y << '\t' << _matrix(x,y) << endl;
    }

  out.close();
}





// LoadForSAMDeNovo: A wrapper for LoadFromSAM for de novo GLMs.
void
GenomeLinkMatrix::LoadFromSAMDeNovo( const string & SAM_file )
{
  assert( DeNovo() ); // this indicates that the input contigs aren't split into bins
  assert( _N_bins > 0 );
  LoadFromSAM( SAM_file, vector<int>( _N_bins, 1 ) );
}



// This multi-SAM-file function is a wrapper for the one-SAM-file function.
void
GenomeLinkMatrix::LoadFromSAMDeNovo( const vector<string> & SAM_files )
{
  // Check the existence of *all* of the files, before taking the time to load in *any* of the files.
  assert( !SAM_files.empty() );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( SAM_files[i] ) );

  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadFromSAMDeNovo( SAM_files[i] );
}




// LoadFromSAMNonDeNovo: A wrapper for LoadFromSAM for non-de novo GLMs.
// Place the Hi-C links in bins whose numbering is determined by BinsPerChromInHumanGenome().  Each chromosome thus gets split up into several contigs of size
// bin_size, with the intent of putting those contigs back together again later.
void
GenomeLinkMatrix::LoadFromSAMNonDeNovo( const string & SAM_file )
{
  assert( _species == "human" );
  assert( !DeNovo() );
  LoadFromSAM( SAM_file, BinsPerChromInHumanGenome() );
}



// This multi-SAM-file function is a wrapper for the one-SAM-file function.
void
GenomeLinkMatrix::LoadFromSAMNonDeNovo( const vector<string> & SAM_files )
{
  // Check the existence of *all* of the files, before taking the time to load in *any* of the files.
  assert( !SAM_files.empty() );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( SAM_files[i] ) );

  for ( size_t i = 0; i < SAM_files.size(); i++ )
    LoadFromSAMNonDeNovo( SAM_files[i] );
}







// NormalizeToDeNovoContigLengths: In a de novo GLM, adjust the matrix data to account for the fact that some bins are smaller.
void
GenomeLinkMatrix::NormalizeToDeNovoContigLengths( const bool use_RE_sites )
{
  cout << "NormalizeToDeNovoContigLengths" << ( use_RE_sites ? " (using RE sites)" : " (using lengths in bp)" ) << endl;

  assert( DeNovo() ); // don't use on binned-human-chromosome data

  // This function should not be called more than once!
  assert( !_normalized );
  _normalized = true;


  // Determine whether to use contig lengths, or number of restriction sites per contig.
  if ( use_RE_sites ) assert( !_contig_RE_sites.empty() );
  const vector<int> & lens = use_RE_sites ? _contig_RE_sites : _contig_lengths;


  // Find the longest contig.  This will be used as a denominator in normalization.
  int64_t longest_contig = *( max_element( lens.begin(), lens.end() ) );
  int64_t longest_squared = longest_contig * longest_contig;


  // Now, normalize!  The math is performed in a single step for each bin, to minimize the effect of rounding error.
  for ( int i = 0; i < _N_bins; i++ )
    for ( int j = 0; j < _N_bins; j++ )
      if ( _matrix(i,j) != 0 ) {
	//cout << "DEBUG:\ti = " << i << "\tj = " << j << "\tmatrix = " << _matrix(i,j) << "\tlongest_squared = " << longest_squared << "\tlen[i] = " << lens[i] << "\tlen[j] = " << lens[j] << "\t-> RESULT: ";
	_matrix(i,j) = _matrix(i,j) * longest_squared / double( (int64_t) lens[i] * (int64_t) lens[j] );
	//cout << _matrix(i,j) << endl;
      }
}






// ReorderContigsByRef: In a de novo GLM, re-order the contigs in _matrix (and _contig_lengths, _contig_RE_sites, _contig_skip) in accordance with this
// TrueMapping.
// Also rearrange the TrueMapping itself so that indexing in it will continue to work.
// This modifies _matrix but does not modify _clusters, so it should be called after all calls to Load...() and Normalize...() but before Cluster().
void
GenomeLinkMatrix::ReorderContigsByRef( TrueMapping & true_mapping )
{
  cout << "ReorderContigsByRef" << endl;

  // Get the full ordering of the contigs on reference.
  vector<int> contig_order = true_mapping.QueriesToGenomeOrder();
  if ( (int) contig_order.size() != _N_bins ) cout << "contig_order.size() = " << contig_order.size() << "\tTrueMapping.NQueries() = " << true_mapping.NQueries() << "\tN_bins = " << _N_bins << endl;
  assert( (int) contig_order.size() == _N_bins ); // if this fails, you're using the wrong TrueMapping object for this dataset

  // Update the TrueMapping.
  true_mapping.ReorderQueries( contig_order );



  // Convert the link matrix to the new ordering.  This can be time-consuming for large matrices.
  // The use of iterators here follows the example at: http://www.guwi17.de/ublas/matrix_sparse_usage.html#Q1
  boost::numeric::ublas::mapped_matrix<int64_t> new_matrix( _N_bins, _N_bins );

  boost::numeric::ublas::compressed_matrix<int64_t>::const_iterator1 it1;
  boost::numeric::ublas::compressed_matrix<int64_t>::const_iterator2 it2;

  for ( it1 = _matrix.begin1(); it1 != _matrix.end1(); ++it1 )
    for ( it2 = it1.begin(); it2 != it1.end(); ++it2 )
      new_matrix( contig_order[ it2.index1() ], contig_order[ it2.index2() ] ) = *it2;

  _matrix = new_matrix;

  // Convert the _contig_lengths, _contig_RE_sites, and _contig_skip vectors to the new ordering.
  vector<int>  old_contig_lengths  = _contig_lengths;
  vector<int>  old_contig_RE_sites = _contig_RE_sites;
  vector<bool> old_contig_skip     = _contig_skip;
  for ( int i = 0; i < _N_bins; i++ ) {
    _contig_orig_order[ contig_order[i] ] = i;
    _contig_lengths   [ contig_order[i] ] = old_contig_lengths[i];
    _contig_RE_sites  [ contig_order[i] ] = old_contig_RE_sites[i];
    _contig_skip      [ contig_order[i] ] = old_contig_skip[i];
  }

}



// SkipShortContigs: Skip contigs of length below min_len.
void
GenomeLinkMatrix::SkipShortContigs( const int & min_len )
{
  cout << "SkipShortContigs with min_len = " << min_len << endl;

  int N_short = 0;
  int64_t short_len = 0;

  for ( int i = 0; i < _N_bins; i++ )
    if ( _contig_lengths[i] < min_len ) {
      N_short++;
      short_len += _contig_lengths[i];
      _contig_skip[i] = true;
    }

  double avg_len = N_short == 0 ? 0 : double(short_len) / N_short;

  // The number of contigs reported as small includes contigs that may have already been marked for skipping, e.g., by SkipRepeats().
  cout << "Marked " << N_short << " contigs (avg length " << avg_len << ") as too short to inform clustering." << endl;
}


// SkipContigsWithFewREs: Skip contigs with fewer than min_N_REs restriction sites.  This can only be done in a de novo GLM in which RE_sites_file != "".
void
GenomeLinkMatrix::SkipContigsWithFewREs( const int & min_N_REs )
{
  cout << "SkipContigsWithFewREs with min_N_REs = " << min_N_REs << endl;
  assert( !_contig_RE_sites.empty() ); // this function can only be used if RE site lengths have been loaded in

  int N_short = 0;
  int64_t short_len = 0;
  int short_N_REs = 0;

  for ( int i = 0; i < _N_bins; i++ )
    if ( _contig_RE_sites[i] < min_N_REs ) {
      N_short++;
      short_len += _contig_lengths[i];
      short_N_REs += _contig_RE_sites[i];
      _contig_skip[i] = true;
    }

  double avg_len   = N_short == 0 ? 0 : double(short_len)   / N_short;
  double avg_N_REs = N_short == 0 ? 0 : double(short_N_REs) / N_short;

  // The number of contigs reported as small includes contigs that may have already been marked for skipping, e.g., by SkipRepeats().
  cout << "Marked " << N_short << " contigs (avg len " << avg_len << ", avg number of RE sites " << avg_N_REs << ") as having too few RE sites to inform clustering (CLUSTER_MIN_RE_SITES = " << min_N_REs << ")." << endl;
}




// SkipRepeats: Skip contigs suspected to be from repetitive regions in the genome.  Contigs are repetitive if they have
// at least <repeat_multiplicity> times as many links as the average contig.  This should be run AFTER normalizing for contig length.
void
GenomeLinkMatrix::SkipRepeats( const double & repeat_multiplicity, const bool flip )
{
  cout << "SkipRepeats with repeat_multiplicity = " << repeat_multiplicity << ", flip = " << noboolalpha << flip << endl;
  bool verbose = false;

  // Find the number of Hi-C links on each contig.  This is as simple as adding rows/columns in the matrix.  Also calculate the total sum.
  // The use of iterators here follows the example at: http://www.guwi17.de/ublas/matrix_sparse_usage.html#Q1
  int64_t N_links_total = 0;
  vector<int64_t> N_links( _N_bins, 0 );

  boost::numeric::ublas::compressed_matrix<int64_t>::const_iterator1 it1;
  boost::numeric::ublas::compressed_matrix<int64_t>::const_iterator2 it2;
  for ( it1 = _matrix.begin1(); it1 != _matrix.end1(); ++it1 )
    for ( it2 = it1.begin(); it2 != it1.end(); ++it2 ) {
      N_links_total           += _matrix( it2.index1(), it2.index2() );
      N_links[ it2.index1() ] += _matrix( it2.index1(), it2.index2() );
      N_links[ it2.index2() ] += _matrix( it2.index1(), it2.index2() );
    }


  // Calculate how many links should be used as a threshold in determining whether a contig is repetitive.
  double N_links_avg = 2.0 * N_links_total / _N_bins;
  if ( verbose ) {
    cout << "N_links_total = " << N_links_total << endl;
    cout << "N_links_avg = " << N_links_avg << endl;
    cout << "repeat_multiplicity = " << repeat_multiplicity << endl;
  }
  assert( N_links_total >= 0 ); // if this fails, there's been an overflow problem

  int N_repetitive = 0;
  int64_t repetitive_len = 0;



  // Find the repetitiveness factor for each contig: the number of links it contains, divided by average.
  for ( int i = 0; i < _N_bins; i++ ) {
    double factor = N_links[i] / N_links_avg;

    // Adjust all link densities by their repetitiveness factors.  This mitigates the effect of mappability and repeat-mediated mapping variation.
    for ( int j = 0; j < _N_bins; j++ )
      if ( _matrix(i,j) != 0 )
	_matrix(i,j) /= factor;

    // Regions that are too repetitive can't be trusted, so skip them entirely.
    if ( factor >= repeat_multiplicity ) {
      if ( verbose ) cout << "CONTIG #" << i << " has " << factor << " x the average number of Hi-C links -> MARKED AS REPETITIVE!" << endl;
      N_repetitive++;
      repetitive_len += ( DeNovo() ? _contig_lengths[i] : _bin_size );
      _contig_skip[i] = true;
    }
  }

  double avg_len = N_repetitive == 0 ? 0 : double(repetitive_len) / N_repetitive;

  // The number of contigs reported as repetitive includes contigs that may have already been marked for skipping, e.g., by SkipShortContigs().
  cout << "Marked " << N_repetitive << " contigs (avg length " << avg_len << ") as too repetitive to inform clustering (CLUSTER_MAX_LINK_DENSITY = " << repeat_multiplicity << ")." << endl;
}






// AHClustering: Apply a greedy agglomerative hierarchical clustering algorithm to cluster the contigs into scaffolds.
// The distance metric between clusters is "average linkage", as described here: http://www2.statistics.com/resources/glossary/a/avglnkg.php
// CEN_contigs, if not an empty vector, lists a set of contig IDs (0-indexed) for contigs containing centromeres.  These contigs will NOT be merged.
// Sets: _clusters (via SetClusters())
void
GenomeLinkMatrix::AHClustering( const int N_CLUSTERS_MIN, const vector<int> & CEN_contigs, const double MIN_AVG_LINKAGE, const double NONINFORMATIVE_RATIO, const bool DRAW_DOTPLOT, const TrueMapping * true_mapping )
{
  // Determine the number of non-skipped contigs (contigs not marked either short or reptitive by the Skip() functions.)  If there are none, throw an error.
  int N_non_skipped = count( _contig_skip.begin(), _contig_skip.end(), false );
  if ( N_non_skipped == 0 ) {
    cerr << "ERROR: There are no informative contigs for clustering.  Every contig is either marked as short (fewer than CLUSTER_MIN_RE_SITES sites) or repetitive (link density more than CLUSTER_MAX_LINK_DENSITY.)  We can't cluster this!" << endl;
  }
  assert ( N_non_skipped > 0 );

  cout << "AHClustering!  (N informative contigs = " << N_non_skipped << ", N_CLUSTERS_MIN=" << N_CLUSTERS_MIN << ", MIN_AVG_LINKAGE=" << MIN_AVG_LINKAGE << ", NONINFORMATIVE_RATIO=" << NONINFORMATIVE_RATIO << ")" << endl;

  // Check that the CEN_contigs vector makes sense, if it's non-empty.
  // All values must be in the range of contig IDs, and there can't be more values than clusters.
  int N_CEN_contigs = CEN_contigs.size(); // may be 0, in which case no CEN-based filtering happens, and that's ok

  assert( N_CEN_contigs <= N_CLUSTERS_MIN );
  for ( int i = 0; i < N_CEN_contigs; i++ ) {
    assert( CEN_contigs[i] >= 0 ); // this was checked when CLUSTER_CONTIGS_WITH_CENS was input
    if ( CEN_contigs[i] >= _N_bins ) { // this wasn't, so report informatively on it now
      cout << "ERROR: Contig #" << CEN_contigs[i] << " was marked as centromeric (CLUSTER_CONTIGS_WITH_CENS).  But there are only " << _N_bins << " contigs in the draft assembly.  What gives?" << endl;
      exit(1);
    }
    cout << "\tCLUSTER_CONTIGS_WITH_CENS: Contig marked with a CEN, will not be merged with other CEN contigs:\t" << CEN_contigs[i] << "\tlength = " << _contig_lengths[ CEN_contigs[i] ] << " bp" << endl;
  }




  set<int>::const_iterator set_it1, set_it2;

  // HEUR: Range for number of clusters that will get printed out to SKY dotplots at out/dotplot.SKY.<n>.jpg
  static const int N_CLUSTERS_MAX = DRAW_DOTPLOT ? 1.25 * N_CLUSTERS_MIN : N_CLUSTERS_MIN;




  const int PRUNE_RATE = 10; // prune after each time we do 1/PRUNE_RATE of the total remaining number of merges; this only affects runtime

  // These objects will contain intermediate products of the algorithm.
  vector<bool> cluster_exists( 2 * _N_bins, false );
  vector<int> cluster_size( 2 * _N_bins, 0 );
  vector<int> bin_to_clusterID( _N_bins, -1 );

  // Put each contig in its own distinct cluster (except contigs marked with contig_skip, which get left out for now).
  // This is the starting point of agglomerative clustering.
  int N_contigs_skipped = 0;
  for ( int i = 0; i < _N_bins; i++ ) {
    if ( _contig_skip[i] ) { N_contigs_skipped++; continue; }
    cluster_exists[i] = true;
    cluster_size[i] = 1;
    bin_to_clusterID[i] = i;
  }
  assert( N_contigs_skipped == _N_bins - N_non_skipped );



  // Calculate all possible "merge scores" for all pairs of clusters.  Make a list sorted by distance.
  // The initial "merge score" values for the initial (one-bin) clusters is simply the amount of link data between each pair of bins.
  // Also make a matrix to keep track of the merge scores for each pair of clusters.  This is the same data as merge_score_map, but in a different form.
  cout << "Creating a 'merge score map'..." << endl;
  multimap< double, pair<int64_t,int>, greater<double> > merge_score_map;
  for ( int i = 0; i < _N_bins; i++ )
    if ( !_contig_skip[i] )
      for ( int j = i+1; j < _N_bins; j++ )
	if ( !_contig_skip[j] )
	  if ( _matrix(i,j) > MIN_AVG_LINKAGE )
	    merge_score_map.insert( make_pair( _matrix(i,j), make_pair(i,j) ) );

  // TEMP
  vector< pair<int,int> > merges_to_do;
  /*
  merges_to_do.push_back( make_pair(17,57) ); // 57 = merger of 0,5
  merges_to_do.push_back( make_pair(0,5) );
    Merges to do for NRRL8004
  merges_to_do.push_back( make_pair(2,6) );
  merges_to_do.push_back( make_pair(3,10) );
  merges_to_do.push_back( make_pair(4,7) );
  merges_to_do.push_back( make_pair(9,47) ); // 47 = merger of 5,8
  merges_to_do.push_back( make_pair(5,8) );
  */


  int N_merges = 0;
  int N_merges_since_prune = 0;
  int N_non_singleton_clusters = 0;



  // Hierarchical clustering!  Repeatedly perform the following three steps:
  // 1. Find the pair of clusters with the highest "merge score".
  // 2. Merge this pair and record the merging.
  // 3. Update the multimap<> of merge scores to reflect the merging.

  while ( 1 ) {
    if ( merge_score_map.empty() )  { cout << "empty merge score map. weird. maybe everything is over-clustered." << endl; break; }

    // 1. Find the pair of clusters with the highest "merge score".
    // This is easy - it's simply the first item in the map for which two clusters actually exist.
    multimap< double, pair<int64_t,int>, greater<double> >::iterator it;

    for ( it = merge_score_map.begin(); it != merge_score_map.end(); ++it )
      if ( cluster_exists[ it->second.first ] && cluster_exists[ it->second.second ] ) {


	// If the CEN_contigs vector is set, it indicates contigs known to contain yeast centromeres.  Don't allow a merge that combines multiple such contigs.
	int N_CENs_in_clusters = 0;
	for ( int i = 0; i < N_CEN_contigs; i++ ) {
	  int cID = bin_to_clusterID[ CEN_contigs[i] ];
	  if ( cID == it->second.first || cID == it->second.second ) N_CENs_in_clusters++;
	}

	if ( N_CENs_in_clusters > 1 ) {
	  cout << "Disallowing a merge because it would put " << N_CENs_in_clusters << " CEN contigs (#" << it->second.first << ",#" << it->second.second << ") into the same cluster" << endl;
	  continue;
	}


	break;
      }

    if ( it == merge_score_map.end() ) {
      cout << "No more merges to do (because MIN_AVG_LINKAGE = " << MIN_AVG_LINKAGE << "); so clustering is done after " << N_merges << " merges" << endl;
      break;
    }

    //cout << "CHOOSING MERGE: " << it->first << "\t" << it->second.first << "," << it->second.second << endl;
    double best_linkage = it->first;
    assert( best_linkage > 0 );
    int best_i = it->second.first;
    int best_j = it->second.second;
    if ( best_i > best_j ) { int swap = best_i; best_i = best_j; best_j = swap; }

    // TEMP: assign certain merges-to-do
    if ( !merges_to_do.empty() ) {
      best_i = merges_to_do.back().first;
      best_j = merges_to_do.back().second;
      merges_to_do.pop_back();
      best_linkage = 0;
      cout << "FORCED MERGE BY USER REQUEST: " << best_i << "," << best_j << endl;
    }


    // 2. Merge this pair and record the merging.
    int new_cluster_ID = _N_bins + N_merges;
    N_merges++;
    N_merges_since_prune++;
    //int N_merges_remaining = _N_bins - N_merges - N_contigs_skipped;
    //cout << "MERGE #" << N_merges << ": Best linkage has value = " << best_linkage << "\tbetween clusters " << PrintSet( _clusters[best_i] ) << " and " << PrintSet( _clusters[best_j] ) << endl;



    // Merge together clusters best_i, best_j into a new cluster.
    // We do this instead of merging one cluster into another so that we don't have to remove all the existing entries in the merge score map that correspond
    // to the existing clusters; we can just mark these clusters as unused so that we can skip over those entries.
    assert( cluster_exists[best_i] );
    assert( cluster_exists[best_j] );
    cluster_exists[best_i] = false;
    cluster_exists[best_j] = false;
    assert( !cluster_exists[new_cluster_ID] );
    cluster_exists[new_cluster_ID] = true;
    cluster_size[new_cluster_ID] = cluster_size[best_i] + cluster_size[best_j];
    if ( best_i < _N_bins ) N_non_singleton_clusters++;
    if ( best_j < _N_bins ) N_non_singleton_clusters++;
    N_non_singleton_clusters--;


    set<int> new_cluster;
    for ( int bin = 0; bin < _N_bins; bin++ )
      if ( bin_to_clusterID[bin] == best_i || bin_to_clusterID[bin] == best_j ) {
	bin_to_clusterID[bin] = new_cluster_ID;
	new_cluster.insert(bin);
      }


    // 3. Update the multimap<> of merge scores to reflect the merging.  This has two sub-parts, 3a and 3b.

    // 3a. Find all merge scores involving clusters that have been used, and remove them.
    // We only do this rarely because it's very time-consuming.
    if ( N_merges_since_prune > ( _N_bins - N_merges ) / PRUNE_RATE ) {
      N_merges_since_prune = 0;

      for ( it = merge_score_map.begin(); it != merge_score_map.end(); it++ ) {
	while ( it != merge_score_map.end() &&
		( !cluster_exists[ it->second.first  ] ||
		  !cluster_exists[ it->second.second ] ) )
	  //cout << "MERGE_SCORE_MAP\tERASING:\t" << it->first << "\t" << it->second.first << "," << it->second.second << endl;
	  merge_score_map.erase(it++);

	if ( it == merge_score_map.end() ) break;
      }

    }


    // 3b. Calculate new score entries for the new cluster - that is, the average linkage from this cluster to each other cluster.
    vector<int64_t> total_linkage_by_cluster( 2 * _N_bins, 0 );
    //set<int, greater<int> > linkages_by_cluster;

    for ( int i = 0; i < _N_bins; i++ ) {
      int cluster_ID = bin_to_clusterID[i];
      if ( cluster_ID == new_cluster_ID ) continue; // no need to calculate linkages within a cluster
      if ( cluster_ID == -1 ) continue; // this happens if _contig_skip[i]

      // This loop is a runtime bottleneck, but attempts to fix it haven't worked...
      for ( set<int>::const_iterator it = new_cluster.begin(); it != new_cluster.end(); ++it ) {
	total_linkage_by_cluster[cluster_ID] += _matrix( i, *it );
	//linkages_by_cluster.insert( _matrix( i, *it ) );
      }
    }


    for ( int i = 0; i < 2 * _N_bins; i++ )
      if ( total_linkage_by_cluster[i] > 0 ) {
	assert( cluster_exists[i] );
	double avg_linkage = double( total_linkage_by_cluster[i] ) / cluster_size[i] / cluster_size[new_cluster_ID];
	//cout << "ADDING LINK TO MERGE_SCORE_MAP: \t" << min(i,new_cluster_ID) << ',' << max(i,new_cluster_ID) << "\tScore = " << avg_linkage << endl;

	if ( avg_linkage < MIN_AVG_LINKAGE ) continue;
	//PRINT2( avg_linkage, THING );
	//avg_linkage += THING;
	//avg_linkage /= 2;
	merge_score_map.insert( make_pair( avg_linkage, make_pair( min(i,new_cluster_ID), max(i,new_cluster_ID) ) ) );
      }

    //PRINT4( N_merges,  N_non_singleton_clusters, N_CLUSTERS_MAX, N_CLUSTERS_MIN );

    // If the number of clusters remaining is sufficiently small, analyze the results.
    if ( N_merges > N_non_skipped / 2 && N_non_singleton_clusters <= N_CLUSTERS_MAX ) {
      SetClusters( bin_to_clusterID, NONINFORMATIVE_RATIO );

      // Report on contigs that appear to be mis-clustered (i.e., they belong to a different chromosome than the plurality of other contigs in their cluster.)
      if ( true_mapping != NULL ) ReportMisclusteredContigs( *true_mapping );

      // Make a dotplot.
      if ( true_mapping != NULL && DRAW_DOTPLOT ) {
	DrawClusterDotplot( *true_mapping );
	string cmd = "mv dotplot.txt dotplot." + boost::lexical_cast<string>(N_non_singleton_clusters) + ".txt";
	system( cmd.c_str() );
	cmd = "mv out/dotplot.SKY.jpg out/dotplot.SKY." + boost::lexical_cast<string>(N_non_singleton_clusters) + ".jpg";
	system( cmd.c_str() );
      }

      // If the number of clusters remaining is SUPER small, we're done.
      //if ( N_merges_remaining == N_CLUSTERS_MIN ) {
      if ( N_non_singleton_clusters == N_CLUSTERS_MIN ) {
	cout << N_merges << " merges made so far; this leaves " << N_CLUSTERS_MIN << " clusters, and so we're done!" << endl;
	break;
      }

    }

    cout << "Merge #" << N_merges << ": Clusters\t#" << best_i << "," << best_j << "\t-> " << new_cluster_ID << "\tLinkage = " << best_linkage << endl;

  }


  // We've broken out of the clustering loop, so we're done, but set clusters one last time.
  SetClusters( bin_to_clusterID, NONINFORMATIVE_RATIO );

}







// Remove from the clusters all contigs whose alignments to reference are sketchy.  "Quality" here refers to the quality of the contig alignment on reference.
// This is a bit of cheating because it implies that we already know the alignments to a reference.  Unlike clustering, it REQUIRES a TrueMapping object.
// Also note that this happens AFTER clustering, unlike the Skip...() functions.  It doesn't affect the clustering; it just prunes the result.
void
GenomeLinkMatrix::ExcludeLowQualityContigs( const TrueMapping & true_mapping )
{
  assert( !_clusters.empty() );
  cout << "ExcludeLowQualityContigs (this is kind of cheating and it must be the LAST modification step on the clusters)" << endl;

  // Create a new ClusterVec consisting only of contigs that pass the quality threshold.
  ClusterVec new_clusters( _clusters.size(), _clusters.N_contigs() );

  for ( size_t i = 0; i < _clusters.size(); i++ )
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      if ( true_mapping.QHighQuality(*it) )
	new_clusters[i].insert( *it );

  _clusters = new_clusters;
}



// MoveContigsInClusters: A heuristic tool to improve _clusters.
// For each contig, try to place it in all other possible clusters, and see if that improves things.  Method:
// STEP 1. Loop over all contigs in clusters.
// STEP 2. Find the number of links between this contig and all other contigs in it cluster, and between this contig and all contigs in other clusters.
// STEP 3. Find which cluster has the most average links to this contig.  If it's not the cluster the contig is already in, move it to that cluster.
// STEP 4: Repeat Steps 1-3 until no more changes are made.
// To move a contig (Step 3), its number of links must be at least <annealing_factor> times as much as the links for the cluster the contig is already in.
// Hence a low annealing_factor (close to 1) is more permissive; a higher annealing_factor is less so.
void
GenomeLinkMatrix::MoveContigsInClusters( const double annealing_factor )
{
  cout << "MoveContigsInClusters with annealing_factor = " << annealing_factor << endl;
  assert( annealing_factor >= 1 );

  // Pre-processing: Make a lookup table of contig ID to cluster ID, and a table of the total contig length of each cluster.
  int N_clusters = _clusters.size();

  vector<int> bin_to_clusterID( _N_bins, -1 );
  vector<int> cluster_length( N_clusters, 0 );
  for ( int i = 0; i < N_clusters; i++ )
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it ) {
      bin_to_clusterID[ *it ] = i;
      cluster_length[i] += _contig_lengths[*it];
    }




  // STEP 4: Repeat Steps 1-3 until no more changes are made.  (In practice, only iterate 20 times, to avoid infinite loops.)
  for ( int x = 0; x < 20; x++ ) {

    int N_changes = 0;

    // STEP 1. Loop over all contigs in clusters.
    for ( int i = 0; i < _N_bins; i++ ) {
      int cluster_i = bin_to_clusterID[i];
      if ( cluster_i == -1 ) continue;

      // STEP 2. Find the number of links between this contig and all other contigs in its cluster, and between this contig and all contigs in other clusters.
      //vector<int> total_N_links( N_clusters, 0 );
      //vector<int>   max_N_links( N_clusters, 0 );
      vector<int64_t> total_N_links_by_len( N_clusters, 0 );
      //vector< set<int> > N_links_set( N_clusters );

      // Loop over all contigs and find the total number of links connecting this contig (i) to all other contigs (j).
      for ( int j = 0; j < _N_bins; j++ ) {
	if ( i == j ) continue; // don't do self-links
	int cluster_j = bin_to_clusterID[j];
	if ( cluster_j == -1 ) continue;

	if ( _matrix(i,j) != 0 ) {
	  //total_N_links       [ cluster_j ] += _matrix(i,j);
	  //max_N_links         [ cluster_j ] = max( max_N_links[cluster_j], _matrix(i,j) );
	  total_N_links_by_len[ cluster_j ] += _matrix(i,j) * _contig_lengths[j];
	  //N_links_set[cluster_j].insert( _matrix(i,j) );
	}
      }

      // STEP 3. Find which cluster has the most average links to this contig.  If it's not the cluster the contig is already in, move it to that cluster.

      // Find the average number of links between this contig and all other contigs in its cluster.  Use this as a starting point.
      int cluster_i_size = _clusters[cluster_i].size() - 1; // subtract 1 to discount this contig in normalization
      int best_cluster = cluster_i;
      double best_cluster_avg_N_links = double( total_N_links_by_len[cluster_i] ) / cluster_i_size;
      //best_cluster_avg_N_links = median_of_set( N_links_set[cluster_i] );
      //best_cluster_avg_N_links = max_N_links[cluster_i];
      best_cluster_avg_N_links *= annealing_factor;
      //PRINT2( i, max_N_links[i] );

      // Normalize the link numbers to find the average distance for each cluster.
      for ( int j = 0; j < N_clusters; j++ ) {

	if ( j == cluster_i ) continue;
	double avg_N_links = double( total_N_links_by_len[j] ) / ( _clusters[j].size() );
	//avg_N_links = median_of_set( N_links_set[j] );
	//avg_N_links = max_N_links[j];

	if ( avg_N_links > best_cluster_avg_N_links ) {
	  //if ( max_N_links[j] > best_cluster_max_N_links ) {
	  best_cluster_avg_N_links = avg_N_links;
	  best_cluster = j;
	}
	//PRINT4( i, j, total_N_links_by_len[j], avg_N_links );
      }

      assert( best_cluster != -1 ); // if this fails, the contig has no links to any contigs in any clusters

      // If this contig is in the right cluster already, we're done.
      if ( cluster_i == best_cluster ) continue;

      N_changes++;

      // Move this contig to the other cluster, and update data structures.
      bin_to_clusterID[i] = best_cluster;
      cluster_length[cluster_i]    -= _contig_lengths[i];
      cluster_length[best_cluster] += _contig_lengths[i];
      _clusters[cluster_i]   .erase( i );
      _clusters[best_cluster].insert( i );


      //PRINT4( i, cluster_i, best_cluster, best_cluster_avg_N_links );


    }

    // STEP 4: Repeat Steps 1-3 until no more changes are made.
    cout << "ITERATION #" << x << ": N changes made: " << N_changes << endl;
    if ( N_changes == 0 ) break;

  }

}








// Create a TrueMapping object for a non-de-novo GenomeLinkMatrix.
TrueMapping
GenomeLinkMatrix::NonDeNovoTrueMapping() const
{
  assert( _species == "human" ); // in the future, load chromosome lengths for non-human species from a FastaSize file

  // Get chromosome lengths.  For human, adjust lengths to account for centromeres.
  const vector<string> chroms = HumanGenome_chroms();
  map<string,int> chrom_lens = HumanGenome_chrom_lengths();
  for ( map<string,int>::iterator it = chrom_lens.begin(); it != chrom_lens.end(); it++ )
    it->second += HumanGenome_centro_size;

  return TrueMapping( _species, _bin_size, chroms, chrom_lens );
}



// Some validation exercises for clustering.
// If true_mapping != NULL (i.e., assembling with a reference) then there will be reference-based validation; otherwise there's just a summary of results.
void
GenomeLinkMatrix::ValidateClusters( const TrueMapping * true_mapping, const bool draw_dotplot ) const
{
  int N_clusters = _clusters.size();

  /* REPORTING ON CLUSTER LENGTHS */

  // Count up the total contig length in each cluster.
  int N_singletons = 0;
  int64_t singleton_len = 0;
  vector<int64_t> cluster_len( N_clusters, 0 );
  for ( int i = 0; i < N_clusters; i++ ) {

    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      cluster_len[i] += _contig_lengths[*it];

    // Also count the number and length of singleton clusters.
    // Singleton clusters tend to come from noisy data and contain badly behaving contigs (e.g., from heterochromatin) rather than real chromosomes.
    if ( _clusters[i].size() == 1 ) {
      N_singletons++;
      singleton_len += _contig_lengths[ *( _clusters[i].begin() ) ];
      cout << "Cluster #" << i << " has 1 contig (singleton) of length\t" << cluster_len[i] << "." << endl;
    }
    else
      cout << "Cluster #" << i << " has\t" << _clusters[i].size() << " contigs, with total length\t" << cluster_len[i] << "." << endl;
  }

  PRINT2( N_clusters, N_singletons );
  cout << "Total length in singleton clusters = " << singleton_len << endl;

  int64_t total_cluster_len = std::accumulate(     cluster_len.begin(),     cluster_len.end(), int64_t(0) );
  int64_t total_contigs_len = std::accumulate( _contig_lengths.begin(), _contig_lengths.end(), int64_t(0) );
  double pct_clustered = 100.0 * total_cluster_len / total_contigs_len;
  cout << "Total length of all clustered contigs = " << total_cluster_len << " out of " << total_contigs_len << " (" << pct_clustered << "%)" << endl;



  /* REPORTING ON CONCENTRATION OF LINKS WITHIN CLUSTERS */

  // Make a lookup table of contig ID to cluster ID.  Unclustered contigs get -1.
  vector<int> cluster_ID( _N_bins, -1 );
  for ( int i = 0; i < N_clusters; i++ )
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      cluster_ID[*it] = i;


  // Add up the link densities for all links between contigs in the same cluster, and for all links between contigs in different clusters.
  int64_t N_links_same_cluster = 0, N_links_diff_cluster = 0;
  int64_t total_len_same_cluster = 0, total_len_diff_cluster = 0;
  int N_bitshifts = 20; // number of base-2 bits to round by, to avoid overflow: start at 20 (e.g., 2^20 ~= 1M)

  for ( int i = 0; i < _N_bins; i++ )
    for ( int j = 0; j < _N_bins; j++ ) {
      if ( cluster_ID[i] == -1 || cluster_ID[j] == -1 ) continue; // both contigs must be clustered
      bool same_cluster = ( cluster_ID[i] == cluster_ID[j] );

      // To weight the link densities properly, it's necessary to find the product of the contig lengths.
      int64_t len_product = int64_t( _contig_lengths[i] ) * int64_t( _contig_lengths[j] ) >> N_bitshifts; // divide by 2^20 to avoid overflow

      assert( ( same_cluster ? N_links_same_cluster : N_links_diff_cluster ) >= 0 );

      ( same_cluster ? total_len_same_cluster : total_len_diff_cluster ) += len_product;
      ( same_cluster ?   N_links_same_cluster :   N_links_diff_cluster ) += len_product * _matrix(i,j);

    }

  double density_same_cluster = double( N_links_same_cluster ) / total_len_same_cluster;
  double density_diff_cluster = double( N_links_diff_cluster ) / total_len_diff_cluster;

  cout << "Average (normalized) intra-cluster link density: " << density_same_cluster << endl;
  cout << "Average (normalized) inter-cluster link density: " << density_diff_cluster << endl;
  cout << "\t-> GLOBAL INTRA-CLUSTER ENRICHMENT: " << density_same_cluster / density_diff_cluster << endl;
  //cout << "\t-> NORMALIZED TO NUMBER OF CLUSTERS: " << ( density_same_cluster / density_diff_cluster / N_clusters ) << endl;

  vector<int> bin_to_clusterID = _clusters.cluster_IDs();

  cout << "Done clustering!" << endl;

  multimap< double, int, greater<double> > cluster_linkages;

  if (0)
  // Find the individual amount of intra-cluster enrichment for each contig.  This is akin to a quality score for each individual contig.
  for ( int i = 0; i < _N_bins; i++ ) {
    int cluster = bin_to_clusterID[i];
    if ( cluster == -1 ) continue; // don't calculate quality scores for non-clustered contigs

    // Find the average linkage between this contig and all other contigs in its cluster.  Also find the average linkage with all other contigs.
    //FindClusterLinkages( i, cluster_linkages );
    //if ( cluster_linkages.empty() ) { cout << "CONTIG #" << i << ": NO LINKS" << endl; continue; }

    int64_t N_in_cluster = 0, N_out_of_cluster = 0;
    double linkage_in_cluster = 0, linkage_out_of_cluster = 0;
    for ( int j = 0; j < _N_bins; j++ ) {
      if ( i == j ) continue;
      ( bin_to_clusterID[j] == cluster ?       N_in_cluster :       N_out_of_cluster ) ++;
      ( bin_to_clusterID[j] == cluster ? linkage_in_cluster : linkage_out_of_cluster ) += _matrix( i, j );
    }

    double ratio = ( linkage_in_cluster / N_in_cluster ) / ( linkage_out_of_cluster / N_out_of_cluster );
    if ( N_out_of_cluster == 0 ) ratio = INT_MAX; // no linkage at all out of cluster: awesome! great job!
    cout << "CONTIG #" << i << "\tENRICHMENT RATIO: " << ratio << endl;
  }




  // Everything past here is reference-based validation.
  if ( true_mapping == NULL ) return;

  if ( true_mapping->species() != _species ) PRINT2( true_mapping->species(), _species );
  assert( true_mapping->species() == _species );


  // Report on contigs that appear to be mis-clustered (i.e., they belong to a different chromosome than the plurality of other contigs in their cluster.)
  ReportMisclusteredContigs( *true_mapping );


  // Draw a dotplot!
  PRINT2( draw_dotplot, true_mapping );
  if ( draw_dotplot ) DrawClusterDotplot( *true_mapping );
}




// Report on contigs that appear to be mis-clustered (i.e., they belong to a different chromosome than the plurality of other contigs in their cluster.)
void
GenomeLinkMatrix::ReportMisclusteredContigs( const TrueMapping & true_mapping ) const
{
  int N_chroms = true_mapping.NTargets();

  int N_non_singleton_clusters = 0;
  int N_aligned_total = 0;
  int N_misclustered = 0;
  int N_informative_total = 0;
  int N_informative_misclustered = 0;
  int64_t len_aligned_total = 0;
  int64_t len_misclustered = 0;
  int64_t len_informative_total = 0;
  int64_t len_informative_misclustered = 0;

  vector<string> cluster_chrom_names( _clusters.size(), "" );

  for ( size_t i = 0; i < _clusters.size(); i++ ) {
    if ( _clusters[i].empty() ) continue;
    if ( _clusters[i].size() > 1 ) N_non_singleton_clusters++;

    // Tallies of which chromosomes contain the contigs in this cluster.
    vector<int> aligns( N_chroms, 0 );
    vector<int> aligns_informative( N_chroms, 0 );
    vector<int64_t> align_len( N_chroms, 0 );
    vector<int64_t> aligns_informative_len( N_chroms, 0 );

    // First, find the chromosome (if any) to which each contig in this cluster aligns.
    int N_aligned = 0;
    int N_informative = 0;
    int64_t len_aligned = 0;
    int64_t len_informative = 0;
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it ) {
      int chrID = true_mapping.QTargetID(*it);
      //cout << "Cluster #" << i << "\tcontains contig " << *it << " which is on chrom " << chrID << endl;
      if ( chrID != -1 ) {
	int len = DeNovo() ? _contig_lengths[*it] : _bin_size;
	N_aligned++;
	len_aligned += len;
	aligns[chrID]++;
	align_len[chrID] += len;
	if ( !_contig_skip[*it] ) { N_informative++; len_informative += len; aligns_informative[chrID]++; aligns_informative_len[chrID] += len; }
      }
    }

    // Now find which chromosome that contains the plurality of these alignments (measured by length, not number of contigs.)
    int N_aligned_correct = 0;
    int N_informative_correct = 0;
    int64_t len_aligned_correct = 0;
    int64_t len_informative_correct = 0;
    for ( size_t j = 0; j < aligns.size(); j++ )
      if ( len_aligned_correct < align_len[j] ) {
	cluster_chrom_names[i] = true_mapping.TargetName(j);
        N_aligned_correct   = aligns[j];
        len_aligned_correct = align_len[j];
	N_informative_correct = aligns_informative[j];
	len_informative_correct = aligns_informative_len[j];
      }

    //cout << "Median chromosome is " << median_chrID << " with " << N_aligned_correct << " aligns!" << endl;

    // Add to global tallies.
    N_aligned_total += N_aligned;
    N_misclustered  += N_aligned - N_aligned_correct;
    N_informative_total += N_informative;
    N_informative_misclustered += N_informative - N_informative_correct;
    len_aligned_total += len_aligned;
    len_misclustered  += len_aligned - len_aligned_correct;
    len_informative_total += len_informative;
    len_informative_misclustered += len_informative - len_informative_correct;
    //cout << "N_misclustered = " << N_misclustered << endl;
  }




  // Report to stdout!
  cout << "REPORT" << endl;
  cout << "N clusters = " << _clusters.size() << endl;
  cout << "N non-singleton clusters = " << N_non_singleton_clusters << endl;
  cout << "Among the " << N_aligned_total << " contigs aligned to canonical chromosomes..." << endl;
  cout << "Number of contigs mis-clustered = " <<   N_misclustered << " of " <<   N_aligned_total << " (" << ( 100.0 *   N_misclustered ) /   N_aligned_total << "%)" << endl;
  cout << "Length of contigs mis-clustered = " << len_misclustered << " of " << len_aligned_total << " (" << ( 100.0 * len_misclustered ) / len_aligned_total << "%)" << endl;
  cout << "Number of informative contigs (used in clustering) mis-clustered = " <<   N_informative_misclustered << " of " <<   N_informative_total << " (" << ( 100.0 *   N_informative_misclustered ) /   N_informative_total << "%)" << endl;
  cout << "Length of informative contigs (used in clustering) mis-clustered = " << len_informative_misclustered << " of " << len_informative_total << " (" << ( 100.0 * len_informative_misclustered ) / len_informative_total << "%)" << endl;
}




// DrawClusterDotplot: Create a visual dotplot of a set of clusters, using QuickDotplot.
// Input a mapping of contig IDs to true chromosome IDs.
void
GenomeLinkMatrix::DrawClusterDotplot( const TrueMapping & true_mapping ) const
{
  {
    string outfile = "dotplot.jpg";
    if ( _species == "human" ) outfile = "dotplot.SKY.jpg";
    cout << "DrawClusterDotplot -> out/" << outfile << endl;
  }


  // Set minimum thresholds for the size of a cluster to plot, and the size of contigs to plot (if de novo).
  // Default values of 1,1: plot everything.
  const unsigned MIN_CLUSTER_SIZE = 1; // in sim datasets, set this to 2 to exclude singletons in centromeres; in de novo it can be 1 (but doesn't have to be)
  const int MIN_CONTIG_LEN = 1;

  int y = 0;

  // Make a local copy of the ClusterVec and sort its clusters so the dotplot looks nicer.
  ClusterVec clusters = _clusters;
  clusters.SortByMedian();


  // Write to the output file.
  string file = "dotplot.txt";
  ofstream out( file.c_str(), ios::out );

  for ( size_t i = 0; i < clusters.size(); i++ ) {
    if ( clusters[i].size() < MIN_CLUSTER_SIZE ) continue; // don't plot small clusters

    for ( set<int>::const_iterator it = clusters[i].begin(); it != clusters[i].end(); ++it ) {
      int chrom = true_mapping.QTargetID(*it);
      if ( chrom == -1 ) continue; // this indicates no data or a non-canonical chromosome

      if ( DeNovo() && _contig_lengths[*it] < MIN_CONTIG_LEN ) continue;

      // Don't plot contigs that weren't used in clustering.
      if ( exclude_noninformative_contigs_from_dotplot && DeNovo() && _contig_skip[*it] ) continue;

      // Plot (x,y) points for the sake of prettiness.  The points are ordered by y so that small clusters can be filtered out.
      // For humans, map onto chromosome names including X,Y.  This allows the QuickDotplot script to produce nicer-looking chromosome colors.
      out << *it << '\t' << y << '\t' << true_mapping.QTargetName(*it) << endl;
    }

    y++;
  }

  out.close();


  // Run the QuickDotplot script to generate a dot plot image, which gets placed at out/<file>.jpg.
  // For details on how this script works, see the script itself.
  string cmd = "QuickDotplot " + file;
  if ( _species == "human" ) cmd = "QuickDotplot.SKY.R " + file;
  system( cmd.c_str() );
}



// Return the ClusterVec describing contig clusters.  This is non-trivial because if ReorderContigsByRef has been run, we must de-re-order the contigs.
ClusterVec
GenomeLinkMatrix::GetClusters() const
{
  ClusterVec clusters2( _clusters.size(), _clusters.N_contigs() );

  for ( size_t i = 0; i < _clusters.size(); i++ )
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      clusters2[i].insert( _contig_orig_order[*it] );

  return clusters2;
}



// DrawHeatmap: Call WriteFile("heatmap.txt"), then run the R script "heatmap.R", which uses R and ggplot2 to make a heatmap image of this GenomeLinkMatrix.
void
GenomeLinkMatrix::DrawHeatmap( const string & heatmap_file ) const
{
  WriteFile("heatmap.txt");

  // The R script "heatmap.R" is hardwired to take "heatmap.txt" as input and write to out/heatmap.jpg.
  // For details on how this script works, see the script itself.
  cout << "Plotting a heatmap at out/" << heatmap_file << endl;
  system("./heatmap.R");

  // Copy the file heatmap.jpg into the place desired.
  if ( heatmap_file != "" && heatmap_file != "heatmap.jpg" )
    system( ( "cp out/heatmap.jpg out/" + heatmap_file ).c_str() );

  cout << "Done plotting!" << endl;
}












// Initialize the matrix of data and allocate memory for it.
void
GenomeLinkMatrix::InitMatrix()
{
  assert( _N_bins > 0 );
  _matrix.resize( _N_bins, _N_bins, 0 );
  _normalized = false;
}




// LoadRESitesFile: Fill _RE_sites_file, _contig_RE_sites.
void
GenomeLinkMatrix::LoadRESitesFile( const string & RE_sites_file )
{
  cout << "Loading contig RE lengths for use in normalization <-\t" << RE_sites_file << endl;
  assert ( DeNovo() );

  if ( RE_sites_file == "." ) {
    cerr << "WARNING: Tried to load an RE_sites_file with the name `.'.  Is this taken from a GLM file with no RE_sites_file listed?  Can't normalize to RE sites." << endl;
    _RE_sites_file = "";
    return;
  }

  _RE_sites_file = RE_sites_file;

  // Get token 1 (zero-indexed) from each line of the lengths file.
  _contig_RE_sites = ParseTabDelimFile<int>( RE_sites_file, 1 );

  // If this assert fails, the input RE_sites_file is inconsistent with the original dataset (wrong number of contigs).
  if( (int) _contig_RE_sites.size() != _N_bins ) {
    cout << "ERROR: Contig RE sites file (" << RE_sites_file << ") seems to imply a different number of contigs (" << _contig_RE_sites.size()
	 << ") than the SAM file[s] do (" << _N_bins << ")." << endl;
    assert(0);
  }

  // Add 1 to all length to prevent dividing by 0 (some short assembly contigs have no restriction sites)
  for ( int i = 0; i < _N_bins; i++ )
    _contig_RE_sites[i]++;
}






// LoadFromSAM: Fill this GenomeLinkMatrix with data from one or more SAM/BAM files.
// Note that this is NOT the same function as ChromLinkMatrix::LoadFromSAM() because the two objects store Hi-C data differently.
// DO NOT CALL THIS FUNCTION DIRECTLY - instead call the wrappers LoadFromSAMDeNovo() or LoadFromSAMNonDeNovo(), which fill bins_per_contig.
void
GenomeLinkMatrix::LoadFromSAM( const string & SAM_file, const vector<int> & bins_per_contig )
{

  bool verbose = true;

  cout << "Reading Hi-C data (" << _species << ") from SAM file " << SAM_file << (verbose ? "\t(dot = 1M alignments)" : "" ) << endl;
  assert( boost::filesystem::is_regular_file( SAM_file ) );
  _SAM_files.push_back( SAM_file );

  // Using the number of bins in each contig, find the offset of each contig in the vector of bins.  This is just an indexing exercise.
  int contig_start = 0;
  vector<int> contig_offsets;
  int N_chroms = bins_per_contig.size();
  for ( int i = 0; i < N_chroms; i++ ) {
    contig_offsets.push_back( contig_start );
    contig_start += bins_per_contig[i];
  }
  assert( contig_start == _N_bins ); // this will fail if _N_bins is not the sum of bins_per_contig.

  vector<int> contig_lens = TargetLengths( SAM_file );


  // Set up a mapped matrix to handle the data as it comes in.  Later this will be converted to a compressed_matrix object for referencing.
  boost::numeric::ublas::mapped_matrix<double> mapped_matrix( _N_bins, _N_bins );



  // Set up a SAMStepper object to read in the alignments.
  SAMStepper stepper(SAM_file);
  stepper.FilterAlignedPairs(); // Only look at read pairs where both reads aligned to the assembly.

  // Loop over all alignments.
  for ( bam1_t * align = stepper.next_read(); align != NULL; align = stepper.next_read() ) {
    const bam1_core_t & c = align->core;

    //cout << "ALIGNMENT POSITIONS: " << c.tid << "." << c.pos << "\t" << c.mtid << "." << c.mpos << endl;

    if ( verbose && stepper.N_aligns_read() % 1000000 == 0 ) cout << "." << flush;

    // Ignore the oddly large number of cases where both reads map to the same location.  (Maybe these are cases in which really only one read mapped?  Dunno.)
    if ( c.pos == c.mpos ) continue;

    // Ignore reads whose partner is not aligned.
    if ( c.mtid == -1 ) continue;

    // Ignore reads with mapping quality 0.  (In the GSM862723 dataset this is roughly 4% of reads.)
    if ( c.qual == 0 ) continue;

    // Skip links involving non-canonical chromosomes.
    if ( c.tid >= N_chroms || c.mtid >= N_chroms ) continue;

    // Find the bin ID of each read.  If the contigs are not being divided, then the contig IDs correspond with bin IDs, and this is easy.
    int bin1 = contig_offsets[c. tid];
    int bin2 = contig_offsets[c.mtid];

    // For non-de novo GLMs, adjust the bin IDs by finding the bin ID of each read's position *within its chromosome*, and then applying the chromosomal offset
    // to the bin IDs.  This makes each genomic position (within a range of bin_size) into a unique bin ID.
    if ( !DeNovo() ) {
      bin1 += c. pos / _bin_size;
      bin2 += c.mpos / _bin_size;
    }

    assert( bin1 < _N_bins );
    assert( bin2 < _N_bins );

    // Don't bother marking intra-bin links; these are not informative for clustering.
    if ( bin1 == bin2 ) continue;

    // TEMP: Re-weight links by their distance from the edge of the contig. (this isn't done yet, doesn't seem to have an effect - I need to work on it more!)
    double weight = 1;
    if(0) {
      const int & len1 = _contig_lengths[bin1];
      const int & len2 = _contig_lengths[bin2];

      int dist1 = min( c. pos, abs( len1 - c. pos ) );
      int dist2 = min( c.mpos, abs( len2 - c.mpos ) );
      int dist = dist1 + dist2; // minimum = 0; maximum = (len1+len2)/2
      double dist_norm = 2 * dist / double( len1+len2 ); // minimum = 0; maximum = 1
      weight = 2 * ( 1 - dist_norm ); // minimum = 0; maximum = 2 (max is at dist = 0)
      //PRINT3( dist, dist_norm, weight );
    }

    // Tally the appropriate spots in the 2-d matrix.
    mapped_matrix(bin1,bin2) += weight;
    mapped_matrix(bin2,bin1) += weight;
  }




  if ( verbose ) cout << endl;
  cout << "N aligns read from " << SAM_file << ": " << stepper.N_aligns_read() << endl;


  // Convert all the data from this SAM file to compressed_matrix format, and add it in.
  cout << "Compressing mapped_matrix data..." << endl;
  boost::numeric::ublas::compressed_matrix<int64_t> compressed_matrix = mapped_matrix;
  _matrix = _matrix + compressed_matrix;
}




// Count the number of bins in the human genome.  This only works in non-de novo GLMs.
// Return the total number of bins, plus a vector indicating how many bins are in each chromosome.
vector<int>
GenomeLinkMatrix::BinsPerChromInHumanGenome() const
{
  assert( _species == "human" );
  assert( !DeNovo() );

  vector<int> bins_per_chrom;

  const vector<string> chroms = HumanGenome_chroms();
  const map<string,int> chrom_lens = HumanGenome_chrom_lengths();

  // Loop over all chromosomes in ascending order.
  for ( size_t i = 0; i < chroms.size(); i++ ) {
    int chrom_len = chrom_lens.at( chroms[i] ) + HumanGenome_centro_size;
    int N_bins = ceil( double(chrom_len) / _bin_size );
    bins_per_chrom.push_back( N_bins );
    //cout << chroms[i] << "  has length = " << chrom_len << " and N_bins = " << N_bins << endl;
  }

  assert( bins_per_chrom.size() == HumanGenome_n_chroms );
  return bins_per_chrom;
}





// SetClusters: Assign the informative contigs (_skip_contig=false) into clusters in accordance with bin_to_clusterID.
// If NONINFORMATIVE_RATIO > 0, also assign the skipped contigs to clusters by how well they match the non-skipped contigs.  To be assigned to a cluster, a
// skipped contig needs to link into that cluster with at least NONINFORMATIVE_RATIO times as many links as any other cluster.  So, lower values of
// NONINFORMATIVE_RATIO (above 0) lead to more aggressive clustering of skipped contigs.
void
GenomeLinkMatrix::SetClusters( const vector<int> & bin_to_clusterID, const double NONINFORMATIVE_RATIO )
{
  cout << "SetClusters" << endl;
  assert( (int) bin_to_clusterID.size() == _N_bins );
  assert( NONINFORMATIVE_RATIO == 0 || NONINFORMATIVE_RATIO > 1 ); // negative values and values in the range (0,1] make no sense for this parameter


  // Assign the non-skipped contigs to clusters.
  _clusters = ClusterVec( bin_to_clusterID, !cluster_noninformative_into_singletons ); // TEMP: is this ok?  recent addition
  //PRINT2( _N_bins, _clusters.N_contigs() );
  //for ( int i = 0; i < _N_bins; i++ ) PRINT2( i, bin_to_clusterID[i] );
  //for ( size_t j = 0; j < _clusters.size(); j++ ) cout << "CLUSTER #" << j << " HAS SIZE " << _clusters[j].size() << endl;



  if ( NONINFORMATIVE_RATIO == 0 ) { CanonicalizeClusters(); return; }



  // Now loop over all skipped contigs.  For each skipped contig, determine which cluster it has the largest average linkage to.
  map<int,int> skipped_clusters;
  int N_pass_ratio = 0, N_fail_ratio = 0, N_fail_cluster = 0;

  for ( int i = 0; i < _N_bins; i++ ) {
    if ( bin_to_clusterID[i] != -1 ) continue; // this indicates a skipped contig

    // Determine, for each cluster, the average amount of linkage between this skipped contig and the contigs in that cluster.
    multimap< double, int, greater<double> > cluster_linkages;
    FindClusterLinkages( i, cluster_linkages );

    if ( cluster_linkages.empty() ) { N_fail_cluster++; continue; } // this contig has no links with any established clusters - so don't assign it to a cluster

    // Find the cluster that contains the most average linkage to this contig.
    multimap< double, int, greater<double> >::const_iterator it = cluster_linkages.begin();
    double best_avg_linkage = it->first;
    int best_cluster = it->second;
    it++;

    // Compare this best cluster's linkage to the second-best.
    double second_best_avg_linkage = it->first;
    double ratio = best_avg_linkage / second_best_avg_linkage;
    bool pass_ratio = ratio >= NONINFORMATIVE_RATIO;
    if ( second_best_avg_linkage == 0 ) pass_ratio = ( best_avg_linkage > 0 ); // handle division by 0
    //cout << "CONTIG #" << i << " (len=" << _contig_lengths[i] << ")\tBEST (" << best_cluster << "): " << best_avg_linkage << "\tSECOND BEST: " << second_best_avg_linkage << "\tRATIO: " << ratio << "\t->\t" << ( pass ? "NO" : "YES" ) << endl;

    if ( !pass_ratio ) { N_fail_ratio++; continue; }

    // Record the cluster with the highest average linkage to this skipped contig.
    // Don't yet assign the skipped contig to the cluster (otherwise it will influence the placement of subsequent skipped contigs.)
    skipped_clusters[i] = best_cluster;
    N_pass_ratio++;
  }


  cout << "SUMMARY: NONINFORMATIVE_RATIO = " << NONINFORMATIVE_RATIO << "\t" << N_pass_ratio << " passed ratio, " << N_fail_ratio << " failed ratio, " << N_fail_cluster << " didn't cluster at all" << endl;

  // Now assign the skipped contigs to their clusters.  Skipped contigs with no links to clusters won't be assigned.
  for ( map<int,int>::const_iterator it = skipped_clusters.begin(); it != skipped_clusters.end(); ++it ) {
    //cout << "SetClusters(): ADDING A CONTIG TO CLUSTER " << it->second << endl;
    _clusters[ it->second ].insert( it->first );
  }

  CanonicalizeClusters();
}




// FindClusterLinkages: Make a map of [avg linkage from this contig to a cluster] -> [the cluster ID], sorted in descending order by linkage.
// Uses the variable _clusters in place.
// May return an empty multimap, if this contig has no links to any contigs in any clusters.
void
GenomeLinkMatrix::FindClusterLinkages( const int contig_ID, multimap< double, int, greater<double> > & cluster_linkages ) const
{
  assert( !_clusters.empty() );

  cluster_linkages.clear(); // this is the output object

  for ( size_t j = 0; j < _clusters.size(); j++ ) {

    int64_t total_linkage = 0;
    //int64_t max_linkage = 0;
    //set<int> linkages;

    // Find this cluster's size.  If it's a singleton cluster, skip it (unless cluster_noninformative_into_singletons is set.)
    int cluster_size = _clusters[j].size();
    if ( cluster_size == 1 && !cluster_noninformative_into_singletons ) continue;

    // Calculate the average linkage beween this contig and this cluster. (alternative, commented forms: calculate maximum or median instead)
    for ( set<int>::const_iterator it = _clusters[j].begin(); it != _clusters[j].end(); ++it ) {
      if ( contig_ID == *it ) cluster_size--; // if the contig is in the cluster, normalize the cluster size properly
      else total_linkage += _matrix( contig_ID, *it );
      //PRINT3( contig_ID, *it, _matrix( contig_ID, *it ) );
      //max_linkage = max( max_linkage, _matrix( contig_ID, *it ) );
      //linkages.insert( _matrix( contig_ID, *it ) );
    }
    double avg_linkage = total_linkage / double( cluster_size );

    //int median_linkage = median_of_set( linkages );
    //PRINT3( avg_linkage, median_linkage, max_linkage );

    cluster_linkages.insert( make_pair( avg_linkage, j ) );
    //cout << "Contig #" << contig_ID << " (len=" << _contig_lengths[contig_ID] << ") into cluster #" << j << " (size = " << _clusters[j].size() << "):\tavg_linkage = " << avg_linkage << endl;
  }
}




// cluster_w_len: Helper struct for CanonicalizeClusters, below.
struct cluster_w_len {
  int cluster_ID;
  int64_t total_len;

  // Comparison function: compares total contig length.  This allows a vector<cluster_w_len> to be sorted by total contig length.
  // The directions of the inequalities will cause the longest clusters to end up at the beginning of the sort.
  bool operator<( const cluster_w_len & x ) const
  { return ( total_len > x.total_len ); }
};


// CanonicalizeClusters: Reorder the clusters by total contig length.
void
GenomeLinkMatrix::CanonicalizeClusters()
{
  // Make a vector of cluster_w_len objects.
  vector<cluster_w_len> cluster_lens( _clusters.size() );

  // Fill the cluster_w_len objects with contig length data.
  for ( size_t i = 0; i < _clusters.size(); i++ ) {
    cluster_lens[i].cluster_ID = i;
    cluster_lens[i].total_len = 0;
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      cluster_lens[i].total_len += _contig_lengths[*it];
  }

  // Sort the cluster_w_len objects so that the clusters with the most length end up at the beginning.
  sort( cluster_lens.begin(), cluster_lens.end() );

  // Reorder the clusters in accordance with the sort order of the cluster_w_len objects.
  ClusterVec new_clusters( _clusters.N_contigs() );
  for ( size_t i = 0; i < _clusters.size(); i++ )
    new_clusters.push_back( _clusters[ cluster_lens[i].cluster_ID ] );

  _clusters = new_clusters;
}
