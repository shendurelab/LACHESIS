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


/*******************************************************************************
 *
 * ChromLinkMatrix.h
 *
 * The ChromLinkMatrix data structure represents a set of intrachromosomal link data from one cluster.
 * A cluster generally represents a group of contigs which the GenomeLinkMatrix's clustering
 * algorithm has clustered together.  The cluster may also represent a true chromosome, split up
 * into _N_contigs contigs, each of size _contig_size; this called a "non-de novo CLM" and is used
 * for algorithmic testing.
 *
 * A ChromLinkMatrix contains a 2-D array of data, of size (2*_N_contigs) x (2*_N_contigs).  There
 * are two bins for each contig, corresponding to the forward and reverse orientation of each
 * contig.  Each bin contains a vector<int> representing the set of distances of all links between
 * the two contigs, assuming a particular orientation.
 *
 * ChromLinkMatrices take their data from SAM files.  You can load data from SAM files using
 * LoadFromSAM...().
 *
 * The goal of the ChromLinkMatrix class is to find a oriented ordering of contigs (a ContigOrdering
 * object) that is best supported by the Hi-C links. The Make...Order() functions employ a
 * graph-based optimization algorithm to find the best ordering of contigs.  The OrientContigs()
 * function implements the WDAG-based contig orientation algorithm.
 *
 * Josh Burton
 * January 2013
 *
 ******************************************************************************/

#ifndef _CHROM_LINK_MATRIX__H
#define _CHROM_LINK_MATRIX__H

#include <algorithm> // max_element
#include <set>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "ClusterVec.h"
#include "ContigOrdering.h"
#include "LinkSizeDistribution.h"
#include "TrueMapping.h"

using namespace std;

class ChromLinkMatrix {
 public:
  /* CONSTRUCTORS */
  // Null constructor.
  ChromLinkMatrix();
  // Create an empty non-de novo ChromLinkMatrix with a contig size and chromosome length.  This is designed to be used in calls to LoadNonDeNovoCLMsFromSAM().
  ChromLinkMatrix(const string &species,
                  const int contig_size,
                  const int chrom_len);
  // Create an empty ChromLinkMatrix with a specific number of contigs.  This is designed to be used in calls to LoadDeNovoCLMsFromSAM().
  ChromLinkMatrix(const string &species,
                  const int cluster_size);
  // Constructor: Load a de novo ChromLinkMatrix from a set of SAM files, and only include contigs with IDs in cluster #cluster_ID.
  ChromLinkMatrix(const string &species,
                  const vector<string> &SAM_files,
                  const string &RE_sites_file,
                  const ClusterVec &clusters,
                  const int cluster_ID);
  // Load a ChromLinkMatrix from a file that was previously written with WriteFile().
  // ChromLinkMatrix( const string & CLM_file ) : _CP_score_dist(1e7) { ReadFile( CLM_file ); }
 ChromLinkMatrix(const string &CLM_file) : _CP_score_dist(10000000) {
    ReadFile(CLM_file);
  }

  // Destructor: Cleans up.
  ~ChromLinkMatrix();
  /* FILE I/O */

  // ReadFile: Read the data from file CLM_file into this ChromLinkMatrix.
  //void ReadFile(const string &CLM_file);
  void ReadFile(const string &CLM_file);
  // WriteFile: Write the data in this ChromLinkMatrix to file "CLM_file".
  void WriteFile(const string &CLM_file,
                 const bool heatmap = false) const;
  // DrawHeatmap: Call WriteFile("heatmap.txt"), then run the R script "heatmap.R", which uses R and
  // ggplot2 to make a heatmap image of this ChromLinkMatrix.
  void DrawHeatmap(const string &heatmap_file = "") const;

  // LoadFromSAMDeNovo: Fill this de novo ChromLinkMatrix with data from one or more SAM/BAM files,
  // but only for the contigs in cluster #cluster_ID.
  void LoadFromSAMDeNovo(const string &SAM_file,
                         const string &RE_sites_file,
                         const ClusterVec &cluster,
                         const int cluster_ID);
  void LoadFromSAMDeNovo(const vector<string> &SAM_files,
                         const string &RE_sites_file,
                         const ClusterVec &cluster,
                         const int cluster_ID);

  // LoadFromSAMNonDeNovo: Fill this non-de novo ChromLinkMatrix with data from one or more SAM/BAM
  // files describing Hi-C reads aligned to a reference genome.
  void LoadFromSAMNonDeNovo(const string &SAM_file,
                            const int chrom_ID);
  void LoadFromSAMNonDeNovo(const vector<string> &SAM_files,
                            const int chrom_ID);

  /* QUERY FUNCTIONS */

  int N_contigs() const { return _N_contigs; }
  int N_bins() const { return 2*_N_contigs; }
  int contig_size() const { return _contig_size; }
  bool has_links() const;
  int NLinks(const int contig1,
             const int contig2) const { return _matrix[2*contig1][2*contig2].size(); }

  // EmptyRows: Return a vector indicating which contigs in the ChromLinkMatrix have data at all.
  // Contigs in centromeres will end up as false.
  vector<bool> ContigsUsed(const bool flag_adjacent = true) const;

  // ContigOrientLogLikelihood: Return the log-likelihood of observing two contigs in a given
  // orientation, as defined by the links between the contigs.
  double ContigOrientLogLikelihood( const int c1, const bool rc1, const int c2, const bool rc2 ) const;

  // OrderingScore: Find the "score" of this ContigOrdering, indicating how well it matches up with
  // the Hi-C link data.  If oriented = true, takes orientation into account.  oriented=false is
  // faster.  Scores with and without orientation can't be directly compared.  If range_start and
  // range_stop are given, the score is only calculated between contig pairs for which at least one
  // contig falls within the range  [range_start,range_stop).  When testing the effect of adding a
  // single contig or a shred of contigs, use this: it's O(n) instead of O(N^2)!  If range_start is
  // given, but range_stop = -1, the score is only calculated between contig pairs on either side of
  // range_start.
  double OrderingScore(const ContigOrdering &order,
                       const bool oriented = true,
                       const int shred_start = -1,
                       const int shred_stop = -1) const;
  // EnrichmentScore: Find the "enrichment", the degree to which the Hi-C links between contigs are
  // between close contigs.  It's analogous to the concentration  of signal along the main diagonal
  // of the heatmap.
  double EnrichmentScore(const ContigOrdering &order) const;

  /* PRE-PROCESSING FUNCTIONS */
  // FindLongestContig(): Fill _longest_contig.  This function should be called right after loading
  // in _contig_lengths.
  void FindLongestContig() { _longest_contig = *(max_element(_contig_lengths.begin(), _contig_lengths.end())); }
  void SetCPScoreDist(const int CP_score_dist) {
    _CP_score_dist = CP_score_dist;
  }
  // PrefilterLinks: Find contig pairs in which the distribution of Hi-C link positions on the
  // contigs suggest long-range rather than short-range contacts.
  void PrefilterLinks(const set<int> &cluster, const TrueMapping *mapping);

  /* GRAPH ALGORITHM METHODS */
  // MAIN ALGORITHMS
  ContigOrdering MakeTrunkOrder(const int min_N_REs) const;
  ContigOrdering MakeFullOrder(const int min_N_REs, const bool use_CP_score = false) const;
  void SpaceContigs(ContigOrdering &order, const LinkSizeDistribution &link_size_distribution) const;

  vector< vector<int> > FindSpanningTree(const int min_N_REs) const;
  void SmoothThornsInTree(vector< vector<int> > & tree) const;
  ContigOrdering TreeTrunk(const vector< vector<int> > &tree, const bool verbose) const;
  ContigOrdering ReinsertShreds(vector< vector<int> > &tree, const int min_N_REs, const bool use_CP_score = false) const;
  void OrientContigs(ContigOrdering &order) const;

 private:
  // DeNovo: Return true iff this is a de novo CLM.
  bool DeNovo() const { return _contig_size == 0; }
  // Initialize the matrix of data and allocate memory for it.
  void InitMatrix();
  void FreeMatrix();
  // LoadRESitesFile: Fill _contig_RE_sites.
  void LoadRESitesFile(const string & RE_sites_file);
  // AddToMatrix: Add a individual Hi-C link to the matrix.  This function is only used when loading
  // data from SAM files.
  void AddLinkToMatrix(const int contig1,
                       const int contig2,
                       const int read1_dist1,
                       const int read1_dist2,
                       const int read2_dist1,
                       const int read2_dist2);

  // CalculateRepeatFactors: Fill the _repeat_factors vector.  This vector contains the 'factor', or
  // multiplicity, of each contig: the ratio by which the total  density of links in each contig
  // exceeds the density expected by chance.  Used in normalization.  Run this after all the links
  // are loaded in.
  void CalculateRepeatFactors();

  // LinkDensity: Return the number of Hi-C links connecting these two contigs (normalized to contig
  // lengths, if this is a de novo CLM.
  double LinkDensity(const int contig1, const int contig2) const;

  // PlotTree: Use grpahviz to print a spanning tree to a graph image at out/<filename>.  <filename>
  // should end in "png".  Note that graphviz is very slow for large graphs, and the output images
  // themselves are sometimes so large as to cause memory problems.  Hence this is NOT RECOMMENDED
  // for graphs with over 500 vertices.
  void PlotTree(const vector< vector<int> > &tree, const string &filename) const;

  // FindGapSize(): Helper function for SpaceContigs().  Finds the best estimate of the gap size
  // following contig #pos in the ContigOrdering.  Uses not just the links between contigs #pos and
  // #pos+1, but also other links between more distant contigs if their gaps have already been
  // estimated.
  int FindGapSize(const ContigOrdering &order,
                  const int pos,
                  const LinkSizeDistribution &lsd,
                  const vector<double> &enrichments ) const;
  void ReportOrderingSize(const ContigOrdering &order) const;

  /* DATA */
  // The species under consideration.  Knowing this helps us avoid confusion.
  string _species;
  int _N_contigs; // number of contigs; controls size of _matrix
  // Contig lengths
  int _contig_size; // size of contigs (for non-de novo CLMs; 0 otherwise)
  vector<int> _contig_lengths; // sizes of all contigs, normalized to longest contig (for de novo CLMs)
  int _longest_contig; // largest element in _contig_lengths; -1 for non-de novo CLMs
  vector<int> _contig_RE_sites; // number of restriction enzyme (RE) sites per contig
  int _most_contig_REs; // largest element in _contig_RE_sites; -1 for non-de novo CLMs

  /* MAIN DATA STRUCTURE: a matrix with size (2*_N_contigs) x (2*_N_contigs); each element is a vector of read pairs' distances */
  vector<int> ** _matrix;
  bool _matrix_init; // is the matrix initialized? (if not, don't free it!)
  // "Repetitiveness factors" for each contig: the number of total links involving this contig,
  // divided by the average.  Used in normalization.
  vector<double> _repeat_factors;
  // tree: A spanning tree.  An intermediate result, created by MakeTrunkOrder, used by MakeFullOrder.
  mutable vector< vector<int> > _tree;
  // The set of SAM files used to gather this data.
  vector<string> _SAM_files;
  // Maximum distance used in the OrderingScore() function.  Higher values give more precise results but take much more runtime.  Defaults to 10Mb.
  int _CP_score_dist;

  // Friend function declarations (see below for documentation on these functions.)
  friend void LoadNonDeNovoCLMsFromSAM(const string &SAM_file, vector<ChromLinkMatrix *> CLMs);
  friend void LoadNonDeNovoCLMsFromSAM(const vector<string> &SAM_files, vector<ChromLinkMatrix *> CLMs);
  friend void LoadDeNovoCLMsFromSAM(const string &SAM_file,
                                    const string &RE_sites_file,
                                    const ClusterVec &clusters,
                                    vector<ChromLinkMatrix *> CLMs);
  friend void LoadDeNovoCLMsFromSAM(const vector<string> &SAM_files,
                                    const string &RE_sites_file,
                                    const ClusterVec &clusters,
                                    vector<ChromLinkMatrix *> CLMs);
};

// LoadDeNovoCLMsFromSAM: Import one or more SAM/BAM files and create a set of de novo
// ChromLinkMatrices corresponding to a ClusterVec (i.e., a set of contig  clusters that have been
// derived by Lachesis' clustering algorithm; see GenomeLinkMatrix).  As many or as few of the
// ChromLinkMatrix pointers can be non-NULL.  Whichever ones are non-NULL will be assumed to
// represent ChromLinkMatrices for the  cluster ID equal to their index in the vector, and will be
// filled accordingly.  If you are creating a set of ChromLinkMatrices for each chromosome, this is
// much faster than calling LoadFromSAMDeNovo individually for each ChromLinkMatrix  object because
// it only reads through the SAM file(s) once.
void LoadDeNovoCLMsFromSAM(const string &SAM_file,
                           const string &RE_sites_file,
                           const ClusterVec &clusters,
                           vector<ChromLinkMatrix *> CLMs);
void LoadDeNovoCLMsFromSAM(const vector<string> &SAM_files,
                           const string &RE_sites_file,
                           const ClusterVec &clusters,
                           vector<ChromLinkMatrix *> CLMs);

// LoadNonDeNovoCLMsFromSAM: Import one or more SAM/BAM files and create a set of non-de novo
// ChromLinkMatrices corresponding to each chromosome.  As many or as few of the ChromLinkMatrix
// pointers can be non-NULL.  Whichever ones are non-NULL will be assumed to represent
// ChromLinkMatrices with chrID  equal to their index in the vector, and will be filled
// accordingly.  This is "cheating" because we're using the true location of each contig to group
// them into clusters.  Its primary use is for heatmaps (DrawHeatmap).  If you are creating a set of
// ChromLinkMatrices for each chromosome, this is much faster than calling LoadFromSAMNonDeNovo
// individually for each  ChromLinkMatrix object because it only reads through the SAM file(s)
// once.
void LoadNonDeNovoCLMsFromSAM(const string &SAM_file,
                              vector<ChromLinkMatrix *> CLMs);
void LoadNonDeNovoCLMsFromSAM(const vector<string> &SAM_files,
                              vector<ChromLinkMatrix *> CLMs);

#endif
