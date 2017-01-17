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
 * LinkSizeDistribution.h
 *
 * A LinkSizeDistribution describes the distribution of Hi-C link sizes in a dataset.
 *
 * The link size distribution is calculated by first examining the intra-contig Hi-C links: that is, links in which both reads map to the same contig in the
 * draft assembly.  The distribution of these link sizes can, with the help of some normalization, be converted into a general distribution of link sizes.
 * Possible sources of error:
 * -- The normalization may be thrown off somewhat by the presence of gaps in scaffolds of the draft assembly.
 * -- The contigs in the draft assembly represent relatively well-behaved, repeat-free regions, so the distribution may not be accurate for trickier regions.
 *
 * Once it's been calculated, the link size distribution can be applied to inter-contig Hi-C links.  Specifically, it can be used to determine a rough gap size
 * between adjacent contigs in an ordering.  See ChromLinkMatrix::SpaceContigs().  The length of intra-contig links is limited by the length of the longest
 * contigs in the draft assembly.  However, it's a well-established fact (cf. Lieberman-Aiden et al., Science 2009) that the density of Hi-C links between two
 * regions A and B on a genome is roughly inversely proportional to the distance between A and B, for large distances.  So we can extrapolate beyond the
 * maximum range of intra-contig links by applying the following equation:
 *
 * Link density(dist > dist_max) = Link density(dist_max) * dist_max / dist
 *
 * With this extrapolation, we can get to arbitrarily high link distances.  However, the distribution still contains little informative value beyond the range
 * of intra-contig links.  The LinkSizeDistribution object should mainly be used for closer-range links (e.g., when ordering and spacing); for more distant
 * links (e.g., when ordering) it's simpler to just use 1/x.
 *
 *
 *
 * There are two main data structures in the LinkSizeDistribution object:
 * 1. _bin_starts: Determines the range of link sizes that go in each bin.
 * 2. _link_density: Describes the link density in each range of link sizes.  Specifically, gives the expected number of links between two 1-Kb regions if the
 *    bin contains the distance between the regions.  Should fall off as roughly 1/distance in accordance with expectations.
 *
 *
 *
 * Josh Burton
 * July 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _LINK_SIZE_DISTRIBUTION__H
#define _LINK_SIZE_DISTRIBUTION__H



#include <string>
#include <vector>
#include <math.h> // sqrt
using namespace std;



class LinkSizeDistribution
{
 public:
  LinkSizeDistribution( const vector<string> & SAM_files );
  LinkSizeDistribution( const string & infile ) { ReadFile( infile ); }


  // ReadFile, WriteFile: Read and write LinkSizeDistribution objects in a simple file format.
  // Once a LinkSizeDistribution is loaded, it's much faster to read and write these files than to parse the SAM files again.
  void ReadFile ( const string & infile );
  void WriteFile( const string & outfile ) const;

  // DrawDotplot: Use QuickDotplot to make a dotplot of this LinkSizeDistribution at out/LinkSizeDistribution.jpg.
  void DrawDotplot( const bool rescale = false ) const;

  // FindEnrichmentOnContig: Input a contig's length L and its set of intra-contig links.  Determine the local enrichment of links on (and, presumably, in the
  // vicinity of) this contig.  This is used in normalization of the gap sizes on either side of this contig.
  double FindEnrichmentOnContig( const int L, const vector<int> & links ) const;

  // FindDistanceBetweenLinks: Input the lengths of two contigs, a set of Hi-C link distances that describe the position of links between the two contigs,
  // and LDE (link density enrichment), which indicates the relative intra-contig link density on these two contigs.
  // Determine the distance D to add to each of these link distances to make the total set of distances best match this LinkSizeDistribution.
  int FindDistanceBetweenLinks( const int L1_0, const int L2_0, const double LDE, const vector<int> & links ) const;

  vector<string> SAM_files() const { return _SAM_files; }


  double log_likelihood_D( const int D, const int L1, const int L2, const double LDE, const vector<int> & links, const vector<double> & log_factorial ) const;

  /* LINK SIZE RANGES */

  // MIN_LINK_DIST: The minimum distance for Hi-C links that we care about here.  This can't be much smaller than the density of restriction sites, or the
  // assumption of constant link density will break down.  Besides, it doesn't have to be very small, since we're going to be using this dataset to look at
  // inter-contig Hi-C links, which by definition are very rarely small.
  // We choose 2^12 = 4,096 because it's (roughly) the density of 6-bp restriction sites.  Also, choosing a power of 2 speeds up the division in dist_bin().
  static const int _MIN_LINK_DIST = 4096;
  // MAX_LINK_DIST: The maximum distance for Hi-C links.  There's no reason why it can't be huge, so it's set to be in the range of entire chromosome sizes.
  // (Link densities above _max_intra_contig_link will be extrapolated anyway.)
  // Must be a power-of-2 multiple of _MIN_LINK_DIST.  We choose 2^16 * _MIN_LINK_DIST = 2^28 = 268,435,456.  This results in _N_bins = 16 * 16 = 256.
  static const int _MAX_LINK_DIST = 1 << 28;


 private:


  /* PRIVATE FUNCTIONS: used locally */
  int BinSize( const int bin_ID ) const; // size in bp of a bin
  int LinkBin( const int dist ) const; // convert a link distance to a bin ID

  void FindExpectedIntraContigLinks( const int L, vector<double> & result, const bool verbose = false ) const;
  void FindExpectedInterContigLinks( const int D, const int L1, const int L2, const double LDE, vector<double> & result, const bool verbose = false ) const;



  /* PRIVATE VARIABLES */

  // max_intra_contig_link_dist: The largest intra-contig link length.  Beyond this length the numbers are just extrapolated as 1/x.
  int _max_intra_contig_link_dist;

  /* THE BINS OF LINK SIZES
   * The bins are in the range [ _MIN_LINK_DIST, MAX_LINK_DIST ).  Their sizes are geometric, so a given bin goes from [X,Y) with Y = X * _SIXTEENTH_ROOT_OF_2.
   */

  int _N_bins;
  vector<int> _bin_starts; // The boundaries of each bin.  This vector has length _N_bins+1 so that the highest bin can have an upper bound.
  vector<double> _link_density; // The bins of link density.  The expected number of links between two regions at a given distance (normaized by _BIN_NORM).


  vector<string> _SAM_files; // the SAM files used to create this distribution

};





#endif
