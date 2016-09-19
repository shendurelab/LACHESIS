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
 * Reporter.h
 *
 * The Reporter class is designed to create a report on the job done by Lachesis.  It outputs metrics of success in human-readable text format and also
 * creates images.
 *
 * In the long term, the Reporter class may be used to write R/ggplot scripts for use in a publication.
 *
 * The Reporter inputs the following data structures:
 * -- The name of any of the SAM files containing alignments to the assembly contigs, which was used in contig clustering and ordering.
 *    The SAM header will be read to get contig names and lengths.
 * -- A TrueMapping, indicating the true locations of assembly contigs on the reference genome (according to a previously computed BLAST alignment)
 * -- A ClusterVec, describing the contig clustering results produced by GenomeLinkMatrix::Cluster()
 * -- Two sets of vector<ContigOrdering>'s, each describing the contig ordering and orienting results produced by ChromLinkMatrix::Optimize() on one or more
 *    clusters.  The first set is trunk orders (i.e., high-quality but incomplete); the second set is full orders (i.e., complete but more risky.)
 *
 *
 *
 * Josh Burton
 * April 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _REPORTER__H
#define _REPORTER__H

#include "TrueMapping.h"
#include "RunParams.h"
#include "ClusterVec.h"
#include "ContigOrdering.h"


// Boost libraries
#include <boost/dynamic_bitset.hpp>

#include <stdlib.h>


#include <iostream>
#include <string>
#include <vector>
using namespace std;





// OrderingFlags: A helper struct for ReporterData.  A set of OrderingFlags describes information about each contig: is it in one of the orderings, is it
// correctly positioned and oriented, etc.
struct OrderingFlags
{
  // DATA FOR CONTIGS IN FULL ORDERINGS: this consists of a bunch of flags for each contig, indicating information about it
  boost::dynamic_bitset<> used;         // this contig is in an ordering
  boost::dynamic_bitset<> high_quality; // this contig is in an ordering and is has orientation quality above a threshold
  boost::dynamic_bitset<> unaligned;    // this contig is in an ordering but doesn't align to any canonical chromosomes
  boost::dynamic_bitset<> chr_mismatch; // this contig is in an ordering and doesn't match the chromosome of the previous contig in the same ordering
  boost::dynamic_bitset<> order_error;  // this contig is in an ordering at a location of an ordering error
  boost::dynamic_bitset<> orient_error; // this contig is in an ordering and is oriented incorrectly
};


// ReporterData
// The Reporter class contains one ReporterData struct (referenced by pointer, so that all the data in it is automatically initialized to 0.)
// The ReporterData class encapsulates all of the raw data points, mostly integer tallies, that are calculated by the Reporter class's Eval...() functions.
// The function Reporter::Report() prints out all the ReporterData in a nice pretty chart.
struct ReporterData
{
  vector<int> cluster_chrom; // for each cluster, the chromosome containing the plurality (by length) of its contigs
  vector<int>     cluster_N_good,   cluster_N_bad,   cluster_N_unaligned;   // for each cluster, the number of contigs (in|not in) the plurality chromosome
  vector<int64_t> cluster_len_good, cluster_len_bad, cluster_len_unaligned; // for each cluster, the length of contigs (in|not in) the plurality chromosome

  boost::dynamic_bitset<> in_cluster; // flags indicating whether each contig is in a cluster

  // OrderingFlags for trunks and full orderings.
  OrderingFlags trunk_flags, order_flags;
};




class Reporter
{
 public:

  /* CONSTRUCTOR: Just load the data structures. */

  Reporter( const RunParams & run_params,
	    const ClusterVec & clusters,
	    const vector<ContigOrdering> & trunks = vector<ContigOrdering>(0, ContigOrdering(0) ),
	    const vector<ContigOrdering> & orders = vector<ContigOrdering>(0, ContigOrdering(0) ) );

  ~Reporter();

  /* EVAL FUNCTIONS: (Mostly) reference-based evaluation.  These put numbers into the ReporterData object.
   * The ReportChart() functions then use these numbers to print pretty output files.
   */

  void Eval() const; // main top-level function; calls other Eval() functions
  void EvalContigUsage() const; // non-reference-based
  void EvalClustering() const;
  void EvalOrderAccuracy( int cluster_ID, bool full_order ) const; // if full_order, then eval _orders[cluster_ID]; otherwise eval _trunks[cluster_ID]
  void EvalGapSizes( int cluster_ID, bool full_order ) const;

  // ReportChart: Print out all the ReporterData info in a pretty chart.
  // There are two versions, depending on whether or not there's a reference available.  The version WithReference should be called after all Eval() functions.
  void ReportChart() const { if ( _true_mapping ) ReportChartWithReference(); else ReportChartNoReference(); }
  void ReportChartWithReference() const;
  void ReportChartNoReference() const;

  /* PLOT FUNCTIONS: These functions produce files at out/<file_head>.jpg. */

  // HistogramOrderingErrors: Make a histogram showing the odds of a contig being mis-ordered, as a function of the contig's length.
  void HistogramOrderingErrors( const bool full_order, const string & file_head ) const;
  // DotplotOrderAccuracy: Make a QuickDotplot of the ordering, highlighting contigs that are mis-ordered.  Uses the output of EvalOrderAccuracy.
  void DotplotOrderAccuracy( const int ordering_ID, const bool full_order, const bool plot_interchrom_contigs, const string & file_head ) const;



 private:

  // RequireReference: Throw a verbose error if _true_mapping == NULL.  This should always be called before anything that uses _true_mapping.
  void RequireReference() const;

  // Helper functions for ReportChart().
  void ReportChartOrderingPercentages( ostream & out ) const;
  void ReportChartOrderingErrors( const bool full_order, ostream & out ) const;

  int N_non_singleton_clusters() const;

  // contig_length_if: Return the length of all contigs with their bit marked as 'true' in the input bitset.
  int64_t contig_length_if( const boost::dynamic_bitset<> & bits ) const;


  /* DATA STRUCTURES.  These are all loaded in by the constructor */

  // Parameters for this run (loaded in from the INI file.)
  const RunParams _run_params;

  // TrueMapping describing contig alignments to referece.  For reference-free assemblies, this is NULL, and many evaluation options/functions are unavailable.
  const TrueMapping * _true_mapping;

  const int _N_contigs;
  const int _N_chroms; // number of chromosomes in reference; will be -1 for a reference-free assembly
  const int _N_clusters, _N_orderings; // the number of orderings may be less than the number of clusters, if orders have not been calculated for all clusters

  const vector<int> _contig_lengths;
  const int64_t _total_contig_length;

  const ClusterVec _clusters;
  const vector<ContigOrdering> _trunks, _orders; // these are empty if !_has_ordering


  // The ReporterData object contains raw data points that are calculated by Eval...() functions.  The numbers are printed out in ReportChart().
  mutable struct ReporterData * _data;
};





// Other reporting functions that are outside the Reporter class.

// MakeWholeAssemblyHeatmap: Make a heatmap of the entire result.
// PLOT_N and USE_RES control how many contigs are considered for plotting and how the link densities are normalized.
void MakeWholeAssemblyHeatmap( const RunParams & run_params, const int PLOT_N = 1000, const bool USE_RES = true );




#endif
