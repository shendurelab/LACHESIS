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


// For documentation, see LinkSizeDistribution.h
#include "LinkSizeDistribution.h"
#include "TextFileParsers.h" // TokenizeFile


#include <assert.h>
#include <cmath> // log2
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip> // setprecision
#include <numeric> // accumulate
#include <algorithm> // max_element

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/SAMStepper.h"


// Boost libraries
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>



/* CONSTANTS */
// These should be LinkSizeDistribution member variables, but the compiler didn't like that idea.

// 2^(1/16) is the geometric "bin size".  In other words, a bin with a lower bound of X will have upper bound X * 2^(1/16).
static const double _SIXTEENTH_ROOT_OF_2 = sqrt( sqrt( sqrt( sqrt( 2 ) ) ) ); // = 1.0442737824274138403219664787399

// Conversion factor for link density.  Setting this to 1e6 = 1000^2 means that the link density is per square Kb.  This has no effect on the results, except
// that it converts the link density numbers into a hopefully more intuitive order of magnitude.
static const double _BIN_NORM = 1e6;



// Calculate the integer log2 of any number (e.g., uint_log2(2) = 1, uint_log2(5) = 2, uint_log2(15) = 3).  Note that uint_log2(0) = 0.
// For maximum speed, this function uses the following method: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
inline unsigned int
uint_log2( const unsigned int x )
{
  static const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
  static const unsigned int S[] = {1, 2, 4, 8, 16};

  unsigned int x2 = x;
  register unsigned int answer = 0;
  if ( x2 & b[4] ) { x2 >>= S[4]; answer |= S[4]; }
  if ( x2 & b[3] ) { x2 >>= S[3]; answer |= S[3]; }
  if ( x2 & b[2] ) { x2 >>= S[2]; answer |= S[2]; }
  if ( x2 & b[1] ) { x2 >>= S[1]; answer |= S[1]; }
  if ( x2 & b[0] ) { x2 >>= S[0]; answer |= S[0]; }
  return answer;
}



// Calculate the fractional part of the log2 of any number in the range [1,2), multiplied by 16 and rounded down (e.g., uint_16log2_frac(1) = 0,
// uint_16log2_frac(1.5) = 9.)  In the name of speed, this function DOES NOT check that its input is in fact in the range [1,2).
// This is used as a recursive helper function for LinkBin(), below.
inline unsigned int
uint_16log2_frac( const double x ) {

  int l = 0;
  double x2 = x * x;
  if ( x2 > 2 ) { x2 /= 2; l += 8; }
  x2 *= x2;
  if ( x2 > 2 ) { x2 /= 2; l += 4; }
  x2 *= x2;
  if ( x2 > 2 ) { x2 /= 2; l += 2; }
  x2 *= x2;
  if ( x2 > 2 ) { x2 /= 2; l += 1; }

  return l;
}






// Constructor for LinkSizeDistribution.  Inputs a set of SAM files and derives the density distribution of link sizes from them.
// Method:
// 1. Make geometrically sized bins of link lengths in the range [MIN_LINK_DIST, MAX_LINK_DIST].
// 2. Determine the range of intra-contig links that we are capable of looking for.
// 3. Find the length of assayable intra-contig sequence in each bin.  This will be used in normalization.
// 4. Read through the SAM file(s) and find the number of links in each bin.
// 5. Normalize to calculate the link density for each intra-contig bin.
// 6. Assume the distribution approximates 1/x for large x, and extrapolate to bins beyond the intra-contig link length.
LinkSizeDistribution::LinkSizeDistribution( const vector<string> & SAM_files )
  : _SAM_files( SAM_files )
{
  cout << "LinkSizeDistribution!" << endl;

  bool verbose = true;

  // Check that the input files all exist.
  assert( !SAM_files.empty() );
  for ( size_t i = 0; i < SAM_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( SAM_files[i] ) );


  // 1. Make geometrically sized bins of link lengths in the range [MIN_LINK_DIST, MAX_LINK_DIST].

  // Make a set of bins that spans the possible range of intra-contig links.  We want the bins to be numerous, but also geometrically sized.
  // We fit 16 bins into each power of 2, so bin 0 goes from [ min, min * 2^1/16 ), bin 1 goes from [ min * 2^1/16, min * 2^1/8 ), etc.
  // Note that, under the the assumption that link density drops as ~1/x, geometrically sized bins will have roughly the same number of links in each of them.
  // Later we'll extrapolate beyond the intra-contig range and make more bins.
  double link_range = log2( double(_MAX_LINK_DIST) / _MIN_LINK_DIST );
  _N_bins = ceil( link_range * 16 );


  // Find the boundaries of the bins.  Each bin has a lower bound equal (before rounding down) to 2^(1/16) =~ 1.044 times the previous bin's lower bound.
  _bin_starts.clear();

  for ( int i = 0; 16*i <= _N_bins; i++ ) {
    double jpower = 1; // jpower = 2^(j/16)
    for ( int j = 0; j < 16 && 16*i+j <= _N_bins; j++ ) {
      _bin_starts.push_back( ( _MIN_LINK_DIST << i ) * jpower );
      jpower *= _SIXTEENTH_ROOT_OF_2;
    }
  }

  assert( _bin_starts.back() == _MAX_LINK_DIST ); // this will fail if _MAX_LINK_DIST is not a power-of-2 multiple of _MIN_LINK_DIST
  assert( (int) _bin_starts.size() == _N_bins + 1 ); // this allows the last bin to have an upper bound


  // 2. Determine the range of intra-contig links that we are capable of looking for.
  // The length of the longest contig in the assembly sets an upper bound on this length.
  vector<int> contig_lens = TargetLengths( SAM_files[0] );
  int N_contigs = contig_lens.size();

  _max_intra_contig_link_dist = *( max_element( contig_lens.begin(), contig_lens.end() ) );
  assert( _max_intra_contig_link_dist < _MAX_LINK_DIST );

  double intra_contig_link_range = log2( double(_max_intra_contig_link_dist) / _MIN_LINK_DIST );
  const int N_intra_contig_bins = ceil( intra_contig_link_range * 16 );

  PRINT7( _MIN_LINK_DIST, _max_intra_contig_link_dist, _MAX_LINK_DIST, link_range, _N_bins, intra_contig_link_range, N_intra_contig_bins );


  // 3. Find the length of assayable intra-contig sequence in each bin.  This will be used in normalization.

  // Specifically, for each bin, calculate the "bin norm", the total length of sequence in the draft assembly that can be assayed for links whose length would
  // place them in that bin.  We'll need these numbers to determine the link densities.
  // The general principle is that links of length L can be found on a contig of length C over a distance of (C-L).  (If C<L, then no links can be found.)
  // Two possible perturbations:
  // -- In reality, there must be room for aligning reads, which have nonzero length.  But we're dealing only with links much longer than a read.
  // -- When aligning to scaffolds rather than contigs, some regions can't be assayed for links because they have gaps.  This is a more serious problem and
  //    I'm not sure how it might perturb the results that we fail to account for gaps.
  vector<int64_t> bin_norms( N_intra_contig_bins, 0 );
  for ( int i = 0; i < N_contigs; i++ ) {

    // For each link bin, calculate the distance on this contig over which links in this bin can be found.
    // Use the minimum link length as indicative of all link lengths in the bin.  This should be ok because the bins are pretty narrow.
    for ( int j = 0; j < N_intra_contig_bins; j++ ) {
      int range = contig_lens[i] - _bin_starts[j];
      if ( range < 0 ) break; // contig is too small for this bin and all remaining bins
      bin_norms[j] += range;
    }
  }


  // 4. Read through the SAM files(s) and find the number of links in each bin.

  // This vector will contain the tally of intra-contig Hi-C link sizes per bin.
  vector<int> N_links( N_intra_contig_bins, 0 );

  // Tally to keep track of how many read pairs pass the filters.
  vector<int> passes( 4, 0 );

  // Loop over the input SAM files.
  for ( size_t i = 0; i < SAM_files.size(); i++ ) {

    cout << "Reading Hi-C data from SAM file " << SAM_files[i] << (verbose ? "\t(dot = 1M alignments)" : "" ) << endl;


    // Set up a SAMStepper object to read in the alignments.
    SAMStepper stepper( SAM_files[i] );
    stepper.FilterAlignedPairs(); // Only look at read pairs where both reads aligned to the assembly.

    // Loop over all pairs of alignments in the SAM file.
    // Note that the next_pair() function assumes that all reads in a SAM file are paired, and the two reads in a pair occur in consecutive order.
    //for ( bam1_t * align = stepper.next_read(); align != NULL; align = stepper.next_read() ) {
    for ( pair< bam1_t *, bam1_t *> aligns = stepper.next_pair(); aligns.first != NULL; aligns = stepper.next_pair() ) {

      if ( verbose && stepper.N_aligns_read() % 1000000 == 0 ) cout << "." << flush;

      const bam1_core_t & c1 = aligns.first->core;
      const bam1_core_t & c2 = aligns.second->core;

      passes[0]++;

      // Ignore reads with mapping quality 0.
      if ( c1.qual == 0 || c2.qual == 0 ) continue;
      passes[1]++;

      // Sanity checks to make sure the read pairs appear as they should in the SAM file.
      assert( c1.tid == c2.mtid );
      assert( c2.tid == c1.mtid );
      assert( c1.pos == c2.mpos );
      assert( c2.pos == c1.mpos );
      assert( c1.tid < N_contigs );
      assert( c2.tid < N_contigs );

      // Only allow intra-contig links.  Note that this is the opposite of the GenomeLinkMatrix and ChromLinkMatrix, in which we only allow intER-contig links.
      if ( c1.tid != c2.tid ) continue;
      passes[2]++;

      // Find the distance implied by this link.
      int dist = abs( c2.pos - c1.pos );
      if ( dist < _MIN_LINK_DIST ) continue; // note that most links will fail this filter unless _MIN_LINK_DIST is set unnecessarily small
      passes[3]++;

      // Convert the distance into a bin.
      int bin = LinkBin( dist );
      //PRINT4( dist, _max_intra_contig_link_dist, bin, N_intra_contig_bins );
      assert( bin != -1 );
      assert( bin < N_intra_contig_bins );
      N_links[bin]++;
    }

    if ( verbose ) cout << endl;
  }

  // Report the results of the pass filtering.
  cout << "Done reading SAM files!  Pass filters:" << endl;
  cout << "\ttotal number of read pairs:\t" << passes[0] << endl;
  cout << "\tboth reads have MQ > 0: \t" << passes[1] << endl;
  cout << "\tlink is intra-contig/scaffold:\t" << passes[2] << endl;
  cout << "\tlink has length >= " << _MIN_LINK_DIST << ":\t" << passes[3] << "\t<- these links are usable in calculating the link size distribution" << endl;



  // 5. Normalize to calculate the link density for each intra-contig bin.
  // This normalization is the product of several calculations:
  // a: Get the number of links found in a bin of range [X,Y).  [Y = X * 2^1/16]
  // b: Divide by the total amount of sequence that can be assayed for this bin.  This gives the expected number of links between a genomic locus l and all
  //    other loci which are at a distance between X and Y from l.
  // c: Divide by the bin size to give the expected number of links between two genomic loci l1, l2 at a specific distance D in the range [X,Y).
  // d: Multiply by 10^6 to get the link density per square Kb - that is, the expected number of links between two 1-Kb regions R1, R2 at a given distance D.
  //    This number is more intuitively human-parseable.
  _link_density = vector<double>( N_intra_contig_bins, 0 );
  for ( int i = 0; i < N_intra_contig_bins; i++ )
    _link_density[i] = double( N_links[i] ) / bin_norms[i] / BinSize(i) * _BIN_NORM;



  // 6. Assume the distribution approximates 1/x for large x, and extrapolate to bins beyond the intra-contig link length.
  // This requires finding the "average" number of links at a large distance.  To do this, we find the set of bins containing the top 1% longest intra-contig
  // links, and we average these bins' link densities.
  int top_bin = N_intra_contig_bins;
  int N_top_links = 0;
  int N_top_links_needed = passes[3] / 100; // 1% of the total number of intra-contig links
  while ( N_top_links < N_top_links_needed )
    N_top_links += N_links[ --top_bin ];

  // Now that we have the bin range corresponding to the 1% longest links, find the average link density of all bins in this range.
  double avg_link_density = 0;
  for ( int i = top_bin; i < N_intra_contig_bins; i++ )
    //PRINT3 (i, _link_density[i], avg_link_density );
    avg_link_density += _link_density[i] * BinSize(i);

  avg_link_density /= ( N_intra_contig_bins - top_bin );
  //PRINT4( N_top_links, N_top_links_needed, top_bin, avg_link_density );


  // Apply this average link density to all bins beyond the last intra-contig bin.
  // Also apply it to the last intra-contig bin - overriding a number generated with what is probably a very low sample size.
  _link_density.resize( _N_bins );
  for ( int i = N_intra_contig_bins - 1; i < _N_bins; i++ )
    _link_density[i] = avg_link_density / BinSize(i);


  /*
  for ( int i = 0; i < N_intra_contig_bins; i++ )
    PRINT5( i, _bin_starts[i], _link_density[i], N_links[i], bin_norms[i] );
  for ( int i = N_intra_contig_bins; i < _N_bins; i++ )
    PRINT3( i, _bin_starts[i], _link_density[i] );
  */
}









// ReadFile: Read a LinkSizeDistribution file that was created with the WriteFile() function, below.
void
LinkSizeDistribution::ReadFile( const string & infile )
{
  // Set dummy values for local variables.  These MUST be overwritten by the input file.
  _SAM_files.clear();
  _N_bins = -1;
  _max_intra_contig_link_dist = -1;

  // Parse the whole input file into tokens.
  cout << "Reading LinkSizeDistribution file " << infile << endl;

  vector< vector<string> > file_as_tokens;
  TokenizeFile( infile, file_as_tokens );


  // Loop over every line of tokens.
  for ( size_t i = 0; i < file_as_tokens.size(); i++ ) {
    const vector<string> & tokens = file_as_tokens[i];


    // Header lines
    if ( tokens[0][0] == '#' ) {

      // Line: "# N_bins = 1"
      if ( tokens.size() == 4 && tokens[1] == "N_bins" )
	_N_bins = boost::lexical_cast<int>( tokens[3] );

      // Line: "# max_intra_contig_link_dist = 1000000"
      if ( tokens.size() == 4 && tokens[1] == "max_intra_contig_link_dist" )
	_max_intra_contig_link_dist = boost::lexical_cast<int>( tokens[3] );

      // Line: "# SAM files used in generating this dataset: [...]"
      else if ( tokens.size() > 8 && tokens[1] == "SAM" )
	for ( size_t j = 8; j < tokens.size(); j++ )
	  _SAM_files.push_back( tokens[j] );

      continue;
    }

    //cout << "Line " << i << "\tN tokens = " << tokens.size() << "\t\t|" << tokens[0] << "|" << tokens[1] << "|" << endl;

    // Non-header lines
    assert( tokens.size() == 2 );
    _bin_starts  .push_back( boost::lexical_cast<int>   ( tokens[0] ) );
    _link_density.push_back( boost::lexical_cast<double>( tokens[1] ) );
  }

  // Make sure all input parameters were read in properly.  If these asserts fail, the input file is missing important header information.
  assert( !_SAM_files.empty() );
  assert( _N_bins != -1 );
  assert( _max_intra_contig_link_dist != -1 );

  // Remove the dummy '0' at the end of the link density vector.
  assert( _link_density.back() == 0 );
  _link_density.pop_back();

  // Make sure the number of bins matches the amount of data we've seen.  Also make sure that MIN_LINK_DIST and MAX_LINK_DIST didn't get changed between runs.
  assert( _N_bins == (int) _link_density.size() );
  assert( _bin_starts[0] == _MIN_LINK_DIST );
  assert( _bin_starts.back() == _MAX_LINK_DIST );
}





void
LinkSizeDistribution::WriteFile( const string & outfile ) const
{
  cout << "Writing a LinkSizeDistribution to file " << outfile << endl;

  ofstream out( outfile.c_str(), ios::out );

  // Write a header explaining this LinkSizeDistribution object.
  out << "# LinkSizeDistribution file - see LinkSizeDistribution.h for documentation of this object type" << endl;
  out << "#\n";
  out << "# The link density for a distance range is defined as the expected number of Hi-C links between two 1-Kb regions at a distance in that range.\n";
  out << "#\n";
  out << "# N_bins = " << _N_bins << endl;
  out << "# max_intra_contig_link_dist = " << _max_intra_contig_link_dist << endl;
  out << "# SAM files used in generating this dataset:";
  for ( size_t i = 0; i < _SAM_files.size(); i++ )
    out << " " << _SAM_files[i];
  out << endl;
  out << "#\n";
  out << "#bin_start\tlink_density" << endl;


  // For each bin, write the bin start and the link density.
  for ( int i = 0; i < _N_bins; i++ )
    out << _bin_starts[i] << '\t' << setprecision(8) << _link_density[i] << endl;

  // Write the last bin stop.  Also write a dummy '0' to maintain this file as a square array that can be read into R.
  out << _MAX_LINK_DIST << "\t0\n";

  out.close();
}




// DrawDotplot: Use QuickDotplot to make a dotplot of this LinkSizeDistribution at out/LinkSizeDistribution.jpg.
// If rescale = true, multiply each bin's link density by the bin size.  This makes a graph that looks roughly flat (y = C) instead of decreasing (y = C/x).
void
LinkSizeDistribution::DrawDotplot( const bool rescale ) const
{
  cout << "Drawing a dotplot of a LinkSizeDistribution" << endl;

  // Make a temp file that will contain this distribution.  This file will be deleted, but QuickDotplot will make its file at out/<temp_file>.jpg.
  string temp_file = "LinkSizeDistribution";
  assert( !boost::filesystem::is_regular_file( temp_file ) );
  ofstream out( temp_file.c_str(), ios::out );

  // Write to the temp file.  Don't draw bins too far beyond the maximum intra-contig link distance.
  for ( int i = 0; i < _N_bins && _bin_starts[i] < 1.5 * _max_intra_contig_link_dist; i++ )
    out << _bin_starts[i] << '\t' << _link_density[i] * ( rescale ? BinSize(i) : 1 ) << endl;
  out.close();

  // Run QuickDotplot to make a JPEG out of this temp file.
  system ( ( "QuickDotplot " + temp_file ).c_str() );

  // Delete the temp file.
  system ( ( "rm " + temp_file ).c_str() );
}



// FindEnrichmentOnContig: Input a contig's length L and its set of intra-contig links.  Determine the local enrichment of links on (and, presumably, in the
// vicinity of) this contig.  This is used in normalization of the gap sizes on either side of this contig.
double
LinkSizeDistribution::FindEnrichmentOnContig( const int L, const vector<int> & links ) const
{
  // We can only examine intra-contig links longer than _MIN_LINK_DIST.  If this contig is shorter than that, return a null enrichment of 1.
  if ( L < _MIN_LINK_DIST ) return 1;


  // Find the number of intra-contig links that are longer than _MIN_LINK_DIST.  Short links are probably spurious ligations, mapping artifacts, etc.
  // And in any event we can't include them in our 'expected number of links' calculations because we haven't recorded their expected density (precisely
  // because most short links are spurious.)
  int total_N_links = 0;
  for ( size_t i = 0; i < links.size(); i++ )
    if ( links[i] >= _MIN_LINK_DIST )
      total_N_links++;

  if ( total_N_links == 0 ) return 0; // might as well save computation time


  // Find the number of intra-contig links (longer than _MIN_LINK_DIST) that we would expect to see on a contig of this size.
  // We calculate the expected number of links in each link-size bin, then add up those numbers to get the total expected number.
  vector<double> N_expected_links;
  FindExpectedIntraContigLinks( L, N_expected_links );
  double total_N_expected_links = accumulate( N_expected_links.begin(), N_expected_links.end(), 0.0 );

  //PRINT4( L, links.size(), total_N_links, total_N_expected_links );
  return total_N_links / total_N_expected_links;
}



const bool allow_longer_links = true; // TEMP: allow links of distance >= L1+L2 (for longer-range calculation).


/* FindDistanceBetweenLinks: Input the lengths of two contigs, a set of Hi-C link distances that describe the position of links between the two contigs,
 * and LDE (link density enrichment), which indicates the relative intra-contig link density on each contig.
 * Determine the distance D to add to each of these link distances to make the total set of distances best match this LinkSizeDistribution.
 *
 * The contigs and links look like this:
 *
 *         read1                                               read2
 *           X|----- dist = d1 -----|         |--- dist = d2 ---|X
 *     ==============================         ====================================
 *             C1 (length L1)       |----D----|         C2 (length L2)
 *
 * The link lengths, as reported to this function, are the values d = d1+d2 for each of the links connecting C1 and C2.
 * Our goal is to estimate the distance D between C1 and C2.  Note that links can only be observed with lengths in the range D <= length <= D+L1+L2.
 *
 * We use a log-likelihood method to determine the best value of D.  (See the function log_likelihood_D() for details of this method.)
 * We expect the log-likelihood as a function of D to have only one major peak, but it will also have a lot of noise, mainly due to binning issues.
 * So we conduct a binary search over the possible values of D, then once we've found the best, we consider all other values in the vicinity.
 * We need to adjust our expectations by the LDE values, due to large-scale variations in link density across the genome.
 *
 *   Thanks to Matthew Snyder and Prof. Joe Felsenstein for help in designing this function. -- Josh
 */
int
LinkSizeDistribution::FindDistanceBetweenLinks( const int L1_0, const int L2_0, const double LDE, const vector<int> & links ) const
{
  // If there are no links, nothing can be done.
  if ( links.empty() ) {
    cout << "WARNING: Called LinkSizeDistribution::FindDistanceBetweenLinks on a pair of contigs with no links.  Why are contigs with no links adjacent to each other?" << endl;
    return INT_MAX;
  }

  // Set L1 <= L2 by convention.
  const int L1 = min( L1_0, L2_0 );
  const int L2 = max( L1_0, L2_0 );

  PRINT3( L1, L2, links.size() );

  // Sanity check on the link lengths.
  for ( size_t i = 0; i < links.size(); i++ ) {
    assert( links[i] > 0 );
    if ( !allow_longer_links ) assert( links[i] <= L1+L2 );
  }


  // As a pre-calculation step, find the logarithms of factorials.
  // This could go in a separate function to save runtime, but I tested it and the savings are pretty minimal.
  vector<double> log_factorial( links.size() + 1, 0 );
  for ( size_t j = 2; j <= links.size(); j++ )
    log_factorial[j] = log_factorial[j-1] + log(j);


  // The maximum value of D is defined by the maximum observed intra-contig link lengths.
  // This guarantees that no link will be considered with a range greater than _MAX_LINK_DIST, which keeps the LinkBin() function safe.
  int max_link = *( max_element( links.begin(), links.end() ) );
  int MAX_D = _MAX_LINK_DIST - max_link;
  //int MAX_D = _max_intra_contig_link_dist - max_link;
  //PRINT( MAX_D );

  int best_D = -1;
  double best_log_likelihood = -INT_MAX;





  // Conduct a binary search on D to find the vicinity of the best log-likelihood LL(D).  Assume that LL(D) is roughly smooth (though with a periodicity
  // associated with bin size.)  Our goal is to find the number D_best that maximizes LL(D).  D_best must be in the range [0, MAX_D).
  // We split this range into 4 equal quartiles: {D_min, D_Q1, D_Q2, D_Q3, D_max}.  We find which one of D_Q1, D_Q2, D_Q3 has the highest log-likelihood;
  // whichever one it is, D_best must be in one of the two quartiles adjacent to it, so our range is cut in half.  We repeat until we've shrunk down the range.
  int D_min = 0;
  int D_max = MAX_D - 1;
  int D_Q2  = ( D_min + D_max ) / 2;
  //double LL_D_min = log_likelihood_D( D_min, L1, L2, LDE, links, log_factorial );
  //double LL_D_max = log_likelihood_D( D_max, L1, L2, LDE, links, log_factorial );
  double LL_D_Q2  = log_likelihood_D( D_Q2,  L1, L2, LDE, links, log_factorial );
  //PRINT6( D_min, LL_D_min, D_max, LL_d_max, D_Q2, LL_D_Q2 );


  while ( D_min + 1 < D_max ) {

    // Find D_Q1, D_Q3.
    int D_Q1 = ( D_min + D_Q2 ) / 2;
    int D_Q3 = ( D_Q2 + D_max ) / 2;

    // Find LL(D_Q1), LL(D_Q3).  Note that the function LL(D) is not actually perfectly smooth; it's periodic with respect to the bin size.  If the range is
    // larger than the bin size, we must make multiple calls to log_likelihood_D and perform some local smoothing in order to get a reliable number.
    int bin_size = D_Q2 * ( _SIXTEENTH_ROOT_OF_2 - 1 );
    double LL_D_Q1 = 0, LL_D_Q3 = 0;
    if ( D_max - D_min > bin_size ) {

      // Take 11 samples, using a step from -5 to +5.
      int step = bin_size / 11;
      for ( int i = -5; i <= 5; i++ ) {
	LL_D_Q1 += log_likelihood_D( D_Q1 + i*step, L1, L2, LDE, links, log_factorial );
	LL_D_Q3 += log_likelihood_D( D_Q3 + i*step, L1, L2, LDE, links, log_factorial );
	//PRINT6( i, step, D_Q1, D_Q1+i*step, D_Q3, D_Q3+i*step );
      }
      LL_D_Q1 /= 11;
      LL_D_Q3 /= 11;

    }
    else {
      LL_D_Q1 = log_likelihood_D( D_Q1, L1, L2, LDE, links, log_factorial );
      LL_D_Q3 = log_likelihood_D( D_Q3, L1, L2, LDE, links, log_factorial );
      //PRINT4( LL_D_Q1, LL_D_Q3 );
    }


    // Determine which of D_Q1, D_Q2, D_Q3 has the highest log-likelihood.  Depending on the result, reset the values of D_min, D_Q2, D_max.
    if ( LL_D_Q1 > LL_D_Q2 && LL_D_Q1 > LL_D_Q3 ) {
      D_max = D_Q2;
      D_Q2  = D_Q1;
      //LL_D_max = LL_D_Q2;
      LL_D_Q2  = LL_D_Q1;
    }
    else if ( LL_D_Q2 >= LL_D_Q1 && LL_D_Q2 > LL_D_Q3 ) {
      D_min = D_Q1;
      D_max = D_Q3;
      //LL_D_min = LL_D_Q1;
      //LL_D_max = LL_D_Q3;
    }
    else {
      D_min = D_Q2;
      D_Q2  = D_Q3;
      //LL_D_min = LL_D_Q2;
      LL_D_Q2  = LL_D_Q3;
    }
  }

  // Record the solution!
  best_D = D_Q2;
  best_log_likelihood = LL_D_Q2;



  cout << "Result of binary search: D = " << best_D << "\twith log-likelihood LL = " << best_log_likelihood << endl;


  if (1) { // TEMP: try an exhaustive search instead
    const int D_step = 1000; // this affects algorithm runtime

    for ( int D = 0; D < MAX_D; D += D_step ) {

      double log_likelihood = log_likelihood_D( D, L1, L2, LDE, links, log_factorial );
      cout << "STUFF:\t" << D << "\t" << log_likelihood << endl;

      // Record the value of D that gives the best log-likelihood.
      if ( log_likelihood > best_log_likelihood ) {
	best_log_likelihood = log_likelihood;
	best_D = D;
      }


    }
  }



  PRINT2( best_D, best_log_likelihood );
  return best_D;
}








// Return the size of bin #bin_ID.
int
LinkSizeDistribution::BinSize( const int bin_ID ) const
{
  assert( bin_ID >= 0 );
  assert( bin_ID < _N_bins );
  return _bin_starts[ bin_ID+1 ] - _bin_starts[ bin_ID ];
}





// LinkBin: convert a link distance to a bin ID.
// Bin #x has a range starting at L = _MIN_LINK_DIST * 2^(x/16).  Therefore, a link of distance L belongs in bin #x = int( 16 * log_2 ( L / _MIN_LINK_DIST ) ).
// This function is called VERY frequently, so it's engineered to be super fast by using bit-shifting and avoiding any direct calls to the log function.
// Also note that if dist >= _MAX_LINK_DIST, the return value will be a bin ID that's higher than _N_bins, and this is not checked for.
int
LinkSizeDistribution::LinkBin( const int dist ) const
{
  if ( dist < _MIN_LINK_DIST ) return -1;

  // First, use bit shifting to find the integral part of log2(dist).  This narrows down the possible bin range to 16 bins.
  unsigned int dist_over_min = dist / _MIN_LINK_DIST; // this division is faster when _MIN_LINK_DIST is a power of 2, because the computer can use bit-shifting
  unsigned int log2_i = uint_log2(dist_over_min);
  //PRINT3( dist, dist_over_min, log2_i );

  // Now, find the fractional part of log2(dist), multiplied by 16.  This determines exactly which bin the link is in.
  int log2_f = uint_16log2_frac( double(dist) / ( _MIN_LINK_DIST << log2_i) );
  return 16 * log2_i + log2_f;

  // OLD VERSION: use a four-step binary search to find exactly which of the 16 possible bins in its range this distance falls into.
  // This sometimes produces slightly different results sometimes due to rounding errors.
  int bin = 16 * log2_i;
  if ( bin+8 < _N_bins && dist >= _bin_starts[bin+8] ) bin += 8;
  if ( bin+4 < _N_bins && dist >= _bin_starts[bin+4] ) bin += 4;
  if ( bin+2 < _N_bins && dist >= _bin_starts[bin+2] ) bin += 2;
  if ( bin+1 < _N_bins && dist >= _bin_starts[bin+1] ) bin += 1;
  //assert( dist >= _bin_starts[bin] );
  //assert( dist <  _bin_starts[bin+1] );
  //PRINT3( dist, bin, _bin_starts[bin] );
  //assert( bin % 16 == log2_f ); // this fails occasionally because of rounding errors

  return bin;
}




// log_likelihood_D: A helper function for FindDistanceBetweenLinks.  Given two contigs of length L1 and L2, calculate the log-likelihood of the observed
// links between them, supposing they're at a distance D from each other.
//
// The algorithmic method is as follows:
// 1. For a given value of D, use the contig lengths and the data from this LinkSizeDistribution to find the expected number of links in each bin.
// 2. Find the actual number of links falling into each bin, by adding D to each of the lengths of the input links.
// 3. Calculate the log-likelihood L(D) of the observations, assuming that the number of links in each bin follows a Poisson distribution.
double
LinkSizeDistribution::log_likelihood_D( const int D, const int L1, const int L2, const double LDE, const vector<int> & links, const vector<double> & log_factorial ) const
{

  const bool verbose = false;
  const int print_winner = -1; // go with 64973 for the first pairing in fly group3

  // 1. For a given value of D, use the contig lengths and the data from this LinkSizeDistribution to find the expected number of links in each bin.
  vector<double> N_expected_links( _N_bins, 0 );
  FindExpectedInterContigLinks( D, L1, L2, LDE, N_expected_links, verbose && 0 );

  double total_N_expected_links = accumulate( N_expected_links.begin(), N_expected_links.end(), 0.0 );
  if ( verbose ) PRINT3( D, total_N_expected_links, links.size() );



  // 2. Find the actual number of links falling into each bin, by adding D to each of the lengths of the input links.
  vector<int> N_observed_links( _N_bins, 0 );
  for ( size_t j = 0; j < links.size(); j++ ) {
    int bin = LinkBin( links[j] + D ); // if we hadn't limited D to MAX_D, this might return a bin ID greater than _N_bins, which would be bad
    if ( bin == -1 ) continue; // link (with D) is too small to fit in any bin
    N_observed_links[bin]++;
  }



  // 3. Calculate the log-likelihood L(D) of the observations, assuming that the number of links at any given size follows a Poisson distribution.
  // The likelihood of the data is the product over all bins of the value:
  // e^(-m) * m^k / k!
  // where m and k are respectively the expected and observed number of links in this bin.  To avoid problems with computation time and underflow, we choose
  // to work with the *log* likelihood, which is the sum (not the product) over all bins of the log of the above expression, i.e.:
  // -m + k ln m - ln(k!)
  double log_likelihood = 0;
  for ( int j = 0; j < _N_bins; j++ ) {
    const double & m = N_expected_links[j];
    const int & k = N_observed_links[j];

    //PRINT4( D, j, m, k );

    // Skip bins with no expectation of links. (Note that all bins within the range of possible link sizes for this D,L1,L2 should have a nonzero expectation
    // of links, unless the bin's link density is 0 due to a total lack of observed intra-contig links.)
    if ( m == 0 ) {
      if ( !allow_longer_links ) assert( k == 0 || _link_density[j] == 0 );
      continue;
    }

    // Calculate the log-likelihood of the observed number of links (k), given the expected number of links (m).
    log_likelihood += -m + k * log(m) - log_factorial[k];
    //PRINT3( m, k, log_likelihood );
  }

  if ( verbose ) cout << "RESULT:\tD=\t" << D << "\tlog_likelihood=\t" << log_likelihood << endl;


  // For the winner, plot the observed and expected distributions.
  if ( D == print_winner && verbose ) {
    for ( int j = 0; j < _N_bins; j++ ) {
      if ( N_expected_links[j] == 0 && N_observed_links[j] == 0 ) continue;
      cout << "WIN AT " << D << "!\t" << j << '\t' << N_expected_links[j] << "\tEXPECTED\n";
      cout << "WIN AT " << D << "!\t" << j << '\t' << N_observed_links[j] << "\tOBSERVED\n";
    }
  }


  assert( log_likelihood <= 0 );
  return log_likelihood;
}





/* FindExpectedIntraContigLinks: A helper function for FindEnrichmentOnContig.  This function determines, for a contig of length L, the expected number of
 * intra-contig Hi-C links in each link-size bin that will be seen.
 *
 * For each bin, the expected number of links is a product of two functions: the link density per bin (i.e., the data stored in this LinkSizeDistribution) and
 * the number of observable links within the contig with sizes that fall in that bin.  The density function of observable link sizes is the tricky part.
 * It's triangular in shape and dependent on _MIN_LINK_DIST, like this:
 *
 *            |
 *            |  .
 *            |  |\           ^
 *    number  |  | \          |
 *        of  |  |  \         |
 *  possible  |  |   \        | height = L - _MIN_LINK_DIST
 *     links  |  |    \       |
 *            |  |     \      |
 *            |  |      \     v
 *            +--+-------+-----
 *             M.L.D.    L
 *
 *                 link size
 *
 ******************************************************/
void
LinkSizeDistribution::FindExpectedIntraContigLinks( const int L, vector<double> & result, const bool verbose ) const
{
  result.resize( _N_bins, 0 );

  if ( verbose ) PRINT2 ( L, _MIN_LINK_DIST );

  for ( int i = 0; i < _N_bins; i++ ) {
    int bin_start = _bin_starts[i], bin_stop = _bin_starts[i+1];
    assert( bin_start < bin_stop );

    // If this contig isn't big enough to include links in this bin, its expectation remains at 0.
    if ( bin_start >= L ) continue;

    // Find the number of observable links of this length (i.e., a vertical slice of the above triangle.)
    int64_t right = min( bin_stop, L );
    int64_t middle_x = ( bin_start + right ) / 2;
    int64_t middle_y = L - middle_x;
    int64_t N_observable_links = ( right - bin_start ) * middle_y;
    assert( N_observable_links >= 0 ); // avoid overflow

    // Calculate the expected number of links by multiplying the number of observable lengths by the density of links in this size range.
    result[i] = N_observable_links * _link_density[i] / _BIN_NORM;
    if ( verbose ) PRINT5( i, bin_start, bin_stop, N_observable_links, result[i] );
  }

}




/* FindExpectedInterContigLinks: A helper function for FindDistanceBetweenLinks.  This function determines, for two contigs of length L1,L2 separated by a
 * putative distance D, the expected number of Hi-C links in each link-size bin that will be seen between the contigs.
 *
 * For each bin, the expected number of links is a product of two functions: the link density per bin (i.e., the data stored in this LinkSizeDistribution) and
 * the number of observable links between the two contigs with sizes that fall in that bin.  The density function of observable link sizes is the tricky part.
 * It's trapezoidal in shape, like this:
 *
 *            |     D+L1    D+L2
 * number of  |      ._______.
 *  possible  |      /       \       ^
 *     links  |     /         \      | height = L1 (because L1 <= L2)
 *            |    /           \     v
 *            +---+-------------+------
 *                D        D+L1+L2
 *
 *                  link size
 *
 ******************************************************/
void
LinkSizeDistribution::FindExpectedInterContigLinks( const int D, const int L1, const int L2, const double LDE, vector<double> & result, const bool verbose ) const
{
  result.resize( _N_bins, 0 );

  if ( verbose ) PRINT4( D, D+L1, D+L2, D+L1+L2 );

  for ( int i = 0; i < _N_bins; i++ ) {
    int bin_start = _bin_starts[i], bin_stop = _bin_starts[i+1];
    assert( bin_start < bin_stop );

    // If this bin doesn't include any lengths of links that could exist between C1 and C2, its expectation remains at 0.
    if ( bin_stop <= D || bin_start >= D+L1+L2 ) continue;

    // Find the number of observable links of this length (i.e., a vertical slice of the above trapezoid.)
    // We do this by determining where the bin falls with respect to the trapezoid's corners, then adding up its area piecemeal.
    int64_t N_observable_links = 0;

    // If the bin falls anywhere along the left slope of the trapezoid...
    if ( bin_start < D+L1 ) {
      int64_t left_slope_min = max( bin_start, D );
      int64_t left_slope_max = min( bin_stop, D+L1 );
      int64_t middle_x = ( left_slope_min + left_slope_max ) / 2;
      int64_t middle_y = middle_x - D;
      if ( verbose ) PRINT5( "LEFT", left_slope_min, left_slope_max, middle_x, middle_y );
      N_observable_links += ( left_slope_max - left_slope_min ) * middle_y;
    }

    // If the bin falls anywhere along the flat middle of the trapezoid...
    if ( bin_stop >= D+L1 && bin_start < D+L2 ) {
      int64_t flat_min = max( bin_start, D+L1 );
      int64_t flat_max = min( bin_stop,  D+L2 );
      if ( verbose ) PRINT3( "FLAT", flat_min, flat_max );
      N_observable_links += ( flat_max - flat_min ) * L1;
    }

    // If the bin falls anywhere along the right slope of the trapezoid...
    if ( bin_stop >= D+L2 ) {
      int64_t right_slope_min = max( bin_start, D+L2 );
      int64_t right_slope_max = min( bin_stop, D+L1+L2 );
      int64_t middle_x = ( right_slope_min + right_slope_max ) / 2;
      int64_t middle_y = D+L1+L2 - middle_x;
      if ( verbose ) PRINT5( "RIGHT", right_slope_min, right_slope_max, middle_x, middle_y );
      N_observable_links += ( right_slope_max - right_slope_min ) * middle_y;
    }

    assert( N_observable_links >= 0 ); // avoid overflow

    // Calculate the expected number of links between C1 and C2 by multiplying the number of observable lengths by the density of links in this size range.
    // Also adjust by the LDE (link density enrichment), to take into account fluctuations in link density across the genome.
    result[i] = N_observable_links * _link_density[i] / _BIN_NORM * LDE;
    if ( verbose ) PRINT5( i, bin_start, bin_stop, N_observable_links, result[i] );

  }

}
