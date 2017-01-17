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
 * GenomeLinkMatrix.h
 *
 *
 * The GenomeLinkMatrix data structure contains a matrix of data corresponding to a Hi-C heatmap.  It's simply a 2-D matrix of size _N_bins x _N_bins.
 *
 * The main data structure is the 2-D matrix, which is a Boost UBLAS compressed_matrix.  This structure is memory-efficient for large, sparse matrices, which
 * are otherwise O(N^2) in memory and very problematic to use.  It is less memory-efficient than a simple 2-D array for small, dense matrices, but in this case
 * the overall memory usage is quite low.  The only danger is that there may be cases where the matrix is large, yet still dense - perhaps when there is a huge
 * amount of Hi-C data.
 * This class is fundamentally different from ChromLinkMatrix because the main data structure is a matrix of ints instead of a matrix of vectors.
 *
 * There are two different types of GenomeLinkMatrices: "de novo" and "non-de novo".
 * -- In a "de novo" GLM, the Hi-C reads have been aligned to a set of contigs from a de novo assembly in progress (including, possibly, a reference genome
 * broken into equal size chunks, which still counts as "de novo" if the reads are aligned to the broken-up fasta.)  The de novo assembly may or may not be of
 * the human genome.  The constructor calls LoadFromSamDeNovo().  The clustering algorithm puts these contigs together into scaffolds and thus improves the
 * assembly.   The de novo contigs can themselves be aligned onto the reference genome, for the purpose of validation of the GLM clusters (see the TrueMapping
 * class); validation then takes place via ValidateClusters().
 * -- In a "non-de novo" GLM, the Hi-C reads have been aligned to the human reference genome.  The human chromosomes are broken down into simulated contigs of
 * size _bin_size in the function LoadFromSAMNonDeNovo().  The clustering algorithm puts these contigs together into scaffolds without knowing where they came
 * from in the human genome.  Validation is very straightforward because we know the "correct" scaffolding answer.  Non-de novo GLMs always have
 * _species = human.
 *
 *
 * The GenomeLinkMatrix contains algorithms to create chromosomal clusters from Hi-C data.  The method is as follows:
 * -- Load in Hi-C data from a set of SAM files.  In a de novo GLM, each contig becomes a bin.  In a non-de novo GLM, each chromosome gets split up to produce
 *    a maximum bin size, thereby simulating contigs that are (but are not known to be) on the same chromosome.
 * -- Using the SAM files, create a matrix of link density between each pair of bins.
 * -- Cluster these bins based on their link density.  This creates a ClusterVec object.  This is the main algorithm!
 * -- Analyze and validate the clusters.
 *
 * The functions ReadFile() and WriteFile() can be used to read/write GenomeLinkMatrix objects from GenomeLinkMatrix format (*.GLM) files.
 * GLM files can be made into graphics using heatmap.R.
 *
 *
 * Algorithm notes:
 * Short and/or repetitive contigs can screw up the clustering algorithm.  It's best to leave these contigs out of the clustering process, though they can
 * optionally be assigned to clusters after the clusters have been made (see SetClusters).
 * Use SkipShortContigs() to mark contigs for skipping if they are below a given size threshold.
 * Use SkipContigsWithFewREs() to mark contigs for skipping if they don't have enough RE (restriction endonuclease) sites.
 * Use SkipRepeats() to mark contigs for skipping if they are repetitive - i.e., they have a normalized number of Hi-C links that is much greater than average.
 *
 *
 *
 *
 *
 * Josh Burton
 * January 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _GENOME_LINK_MATRIX__H
#define _GENOME_LINK_MATRIX__H


#include "ClusterVec.h"
#include "TrueMapping.h"
#include <string>
#include <vector>
#include <map> // multimap
using namespace std;

// Boost libraries
#include <boost/numeric/ublas/matrix_sparse.hpp>




class GenomeLinkMatrix
{
 public:

  /* CONSTRUCTOR */

  // Make a non-de novo GenomeLinkMatrix with no Hi-C link data.  This isn't good for much other than calling NonDeNovoTrueMapping().
  GenomeLinkMatrix( const string & species, const int bin_size );
  // Load a non-de novo GenomeLinkMatrix (HUMAN ONLY) with a set of SAM files representing human-genome alignments.  Split the GLM into bins of size bin_size.
  GenomeLinkMatrix( const vector<string> & SAM_files, const int bin_size );
  // Load a de novo GenomeLinkMatrix with a set of SAM files representing alignments to contigs.
  GenomeLinkMatrix( const string & species, const vector<string> & SAM_files, const string & RE_sites_file = "" );
  // Load a GenomeLinkMatrix from a file that was previously written with WriteFile().  This may or may not be a de novo GLM.
  GenomeLinkMatrix( const string & LM_file ) { ReadFile( LM_file ); }



  /* QUERY FUNCTIONS */
  int N_bins() const { return _N_bins; }


  /* FILE I/O */

  // ReadFile: Read the data from file GLM_file into this GenomeLinkMatrix.
  void ReadFile( const string & GLM_file );
  // WriteFile: Write the data in this GenomeLinkMatrix to file GLM_file.
  void WriteFile( const string & GLM_file ) const;


  // LoadForSAMDeNovo: A wrapper for LoadFromSAM for de novo GLMs.
  void LoadFromSAMDeNovo( const        string  & SAM_file  );
  void LoadFromSAMDeNovo( const vector<string> & SAM_files );

  // LoadFromSAMNonDeNovo: A wrapper for LoadFromSAM for non-de novo GLMs.
  void LoadFromSAMNonDeNovo( const        string  & SAM_file  );
  void LoadFromSAMNonDeNovo( const vector<string> & SAM_files );


  /* DATA PRE-PROCESSING: NORMALIZATION, RE-ORDERING, ETC. */
  /* Note: if you call one of these, DO NOT call WriteFile() afterward!  GLM files should contain non-normalized, non-reordered data! */

  // NormalizeToDeNovoContigLengths: In a de novo GLM, adjust the matrix data to account for the fact that some bins are smaller.
  // If use_RE_sites, then use the number of RE sites per contig as the denominator for normalization (requires _contig_RE_sites to be loaded.)
  void NormalizeToDeNovoContigLengths( bool use_RE_sites );

  // ReorderContigsByRef: In a de novo GLM, re-order the contigs in _matrix (and _contig_lengths, _contig_RE_sites, _contig_skip) in accordance with this
  // TrueMapping.
  // Also rearrange the TrueMapping itself so that indexing in it will continue to work.
  // This modifies _matrix but does not modify _clusters, so it should be called after all calls to Load...() and Normalize...() but before Cluster().
  void ReorderContigsByRef( TrueMapping & true_mapping );


  /* SKIP FUNCTIONS: These functions mark certain contigs for skipping by setting their _skip_contig bit.
     Contigs marked for skipping will not be used during the clustering algorithm in Cluster(), but may be assigned to clusters afterward by SetClusters(). */

  // SkipShortContigs: Skip contigs of length below min_len.
  void SkipShortContigs( const int & min_len );

  // SkipContigsWithFewREs: Skip contigs with fewer than min_N_REs restriction sites.  This can only be done in a de novo GLM in which RE_sites_file != "".
  void SkipContigsWithFewREs( const int & min_N_REs );

  // SkipRepeats: Skip contigs suspected to be from repetitive regions in the genome.  Contigs are considered repetitive if they have at least (or at MOST, if
  // flip == true) <repeat_multiplicity> times as many links as the average contig.
  // This should be run AFTER normalizing for contig length.
  void SkipRepeats( const double & repeat_multiplicity, const bool flip = false );

  /* MAIN CLUSTERING ALGORITHMS
     Contigs that have been marked as "skipped" by one of the Skip...() functions are not used in clustering.  However, if set_skipped_contigs = true, then
     after clustering, skipped contigs are assigned to clusters by how well they match the non-skipped contigs (see SetClusters()).
     Contigs that have been marked as centromeric (CEN_contigs) will not be merged into the same cluster.
  */

  // Agglomerative Hierarchical clustering
  void AHClustering( const int N_CLUSTERS_MIN, const vector<int> & CEN_contigs, const double MIN_AVG_LINKAGE, const double NONINFORMATIVE_RATIO, const bool DRAW_DOTPLOT, const TrueMapping * true_mapping );

  // Improvements to clustering algorithms.
  void ExcludeLowQualityContigs( const TrueMapping & true_mapping ); // remove from the clusters all contigs whose alignments to reference are sketchy
  void MoveContigsInClusters( const double annealing_factor );



  /* VALIDATION */

  // Create a TrueMapping object for a non-de-novo GenomeLinkMatrix.
  TrueMapping NonDeNovoTrueMapping() const;

  // Some validation exercises for clustering.
  // If true_mapping != NULL (i.e., assembling with a reference) then there will be reference-based validation; otherwise there's just a summary of results.
  void ValidateClusters( const TrueMapping * true_mapping, bool draw_dotplot = false ) const;

  // Report on contigs that appear to be mis-clustered (i.e., they belong to a different chromosome than the plurality of other contigs in their cluster.)
  void ReportMisclusteredContigs( const TrueMapping & true_mapping ) const;


  /* OUTPUT AND REPORTING */

  int64_t NLinks( const int i, const int j ) const { return _matrix(i,j); }

  // SetClusters/GetClusters: Set or return the ClusterVec produced by Cluster().  Note that the clusters are not guaranteed to be in any particular order.
  void SetClusters( const ClusterVec & clusters ) { _clusters = clusters; }
  ClusterVec GetClusters() const;

  // DrawHeatmap: Call WriteFile("heatmap.txt"), then run the R script "heatmap.R", which uses R and ggplot2 to make a heatmap image of this GenomeLinkMatrix.
  void DrawHeatmap( const string & heatmap_file = "" ) const;

  // DrawClusterDotplot: Create a visual dotplot of a set of clusters, using QuickDotplot.  This requires a TrueMapping.
  void DrawClusterDotplot( const TrueMapping & true_mapping ) const;

 private:

  // DeNovo: Return true iff this is a de novo GLM.
  bool DeNovo() const { return _bin_size == 0; }

  // Initialize the matrix of data and allocate memory for it.
  void InitMatrix();

  // LoadRESitesFile: Fill _RE_sites_file, _contig_RE_sites.
  void LoadRESitesFile( const string & RE_sites_file );

  // LoadFromSAM: Fill this GenomeLinkMatrix's matrix with Hi-C data from one or more SAM/BAM files containing human genome-aligned reads.
  // Note that this is NOT the same function as ChromLinkMatrix::LoadFromSAM() because the two objects store Hi-C data differently.
  // DO NOT CALL THIS FUNCTION DIRECTLY - instead call the wrappers LoadFromSAMDeNovo() or LoadFromSAMNonDeNovo(), which fill bins_per_contig.
  void LoadFromSAM( const string & SAM_file, const vector<int> & bins_per_contig );

  // Count the number of bins in the human genome.  This only works in non-de novo GLMs.
  vector<int> BinsPerChromInHumanGenome() const;

  // SetClusters: Assign the informative contigs (_skip_contig=false) into clusters in accordance with bin_to_clusterID.
  // If NONINFORMATIVE_RATIO != 0, also assign the skipped contigs to clusters by how well they match the non-skipped contigs.
  void SetClusters( const vector<int> & bin_to_clusterID, const double NONINFORMATIVE_RATIO );
  // FindClusterLinkages: Make a map of [avg linkage from this contig to a cluster] -> [the cluster ID], sorted in descending order by linkage.
  void FindClusterLinkages( const int contig_ID, multimap< double, int, greater<double> > & cluster_linkages ) const;
  void CanonicalizeClusters(); // reorder the clusters by total contig length



  /* DATA MEMBERS */

  // The species under consideration.  Knowing this helps us avoid confusion.
  string _species;

  // The set of SAM files used to gather this data.
  vector<string> _SAM_files;
  // The RE_sites file name (RE_COUNTS_FILE) used to gather this data.  Left empty for non-de novo GLMs, and for de novo GLMs loaded with RE_sites_file = "".
  string _RE_sites_file;

  // The matrix of Hi-C links.
  int _N_bins; // number of bins
  int _bin_size; // size of bins, in non-de novo GenomeLinkMatrices; set to 0 for de novo GenomeLinkMatrices
  boost::numeric::ublas::compressed_matrix<int64_t> _matrix; // the main data structure!

  // _contig_orig_order: Set by ReorderContigsByRef() so it can remember the original ordering of the contigs and output them properly in GetClusters().
  vector<int> _contig_orig_order;

  // Contig lengths, and numbers of restriction sites per contig.  These are left empty in non-de novo GLMs.
  vector<int> _contig_lengths;
  vector<int> _contig_RE_sites;

  bool _normalized; // has NormalizeToDeNovoContigLengths() been called?

  // contig_skip: Flags indicating which contigs should not be used in clustering (though they may get added in afterward; see SetClusters.)
  // Contigs may be marked for skipping if they are (1) repetitive, as determined by SkipRepeats(); or (2) too short, as determined by SkipShortContigs().
  vector<bool> _contig_skip;

  // _clusters: The output of the SetClusters() function.  Used by GetClusters(), DrawClusterDotplot().
  ClusterVec _clusters;

};




#endif
