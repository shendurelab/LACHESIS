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


// For documentation, see Reporter.h
#include "Reporter.h"
#include "RunParams.h"
#include "TrueMapping.h"
#include "ClusterVec.h"
#include "ContigOrdering.h"
#include "GenomeLinkMatrix.h" // for MakeWholeAssemblyHeatmap only
#include "TextFileParsers.h" // ParseTabDelimFile()

// C++ includes
#include <assert.h>
#include <math.h> // isfinite
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> // setprecision
#include <algorithm> // sort, max_element
#include <numeric> // accumulate

// Local includes
#include "TimeMem.h"
#include "gtools/N50.h"
#include "gtools/SAMStepper.h" // NTargetsInSAM(), TargetLengths()


// Boost libraries
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>








// This is a local variable for dummy use.  (It makes this module non-thread-safe, but who cares.)
static char _LINE[10000];



// AND(): Calculate the logical AND of two dynamic_bitsets.
boost::dynamic_bitset<>
AND( const boost::dynamic_bitset<> a, const boost::dynamic_bitset<> b )
{
  assert ( a.size() == b.size() );
  boost::dynamic_bitset<> c( a.size(), false );
  for ( size_t i = 0; i < a.size(); i++ )
    if ( a[i] && b[i] )
      c[i] = true;

  return c;
}




// SetToVec: Convert a set<int> into a vector<int> for easy random access.
vector<int>
SetToVec( const set<int> & s )
{
  vector<int> v;
  for ( set<int>::const_iterator it = s.begin(); it != s.end(); ++it )
    v.push_back( *it );
  return v;
}




// Calculate the Pearson correlation coefficient, r, indicating correlation between two data sets.  Ranges from -1 to 1.
double
PearsonCorrelation( const vector<double> & x, const vector<double> & y )
{
  int N = x.size();
  assert( x.size() == y.size() );

  // Calculate the means of x and y.  Also find NaN's/inf's and flag them to be ignored.
  double x_mean = 0, y_mean = 0;
  vector<bool> ignore( N, false );
  int N_ignore = 0;
  for ( int i = 0; i < N; i++ ) {
    if ( !isfinite( x[i] ) || !isfinite( y[i] ) ) {
      ignore[i] = true;
      N_ignore++;
    }
    else {
      x_mean += x[i];
      y_mean += y[i];
    }

  }

  x_mean /= ( N - N_ignore );
  y_mean /= ( N - N_ignore );

  // Calculate the variances of x and y, and the covariance between them.
  double covar = 0;
  double x_var = 0, y_var = 0;

  for ( int i = 0; i < N; i++ ) {
    if ( ignore[i] ) continue;
    double x_diff = x[i] - x_mean;
    double y_diff = y[i] - y_mean;
    covar += x_diff * y_diff;
    x_var += x_diff * x_diff;
    y_var += y_diff * y_diff;
  }

  return covar / sqrt( x_var * y_var );
}





// Constructor!  Load in all the data, and do a bunch of sanity checks.  trunks and orders can both be empty, or they can both be full.
Reporter::Reporter( const RunParams & run_params,
		    const ClusterVec & clusters,
		    const vector<ContigOrdering> & trunks,
		    const vector<ContigOrdering> & orders )
  : _run_params( run_params ),
    _true_mapping( run_params.LoadTrueMapping() ),
    _N_contigs  ( NTargetsInSAM( run_params._SAM_files[0] ) ),
    _N_chroms   ( _true_mapping ? _true_mapping->NTargets() : -1 ),
    _N_clusters ( clusters.size() ),
    _N_orderings( orders.size() ),
    _contig_lengths( TargetLengths( run_params._SAM_files[0] ) ),
    _total_contig_length( accumulate( _contig_lengths.begin(), _contig_lengths.end(), int64_t(0) ) ),
    _clusters( clusters ),
    _trunks( trunks ),
    _orders( orders ),
    _data( new ReporterData() )
{
  cout << "Launching a Reporter for a result with N_contigs = " << _N_contigs << ", N_chroms = " << _N_chroms << ", N_clusters = " << _N_clusters << ", N orderings = " << _N_orderings << ", total assembly length = " << _total_contig_length << endl;

  /* SANITY CHECKS
   *
   * If any of these sanity checks fails, the data structures aren't from the same run.
   *
   */

  // The SAM file should be internally consistent.
  assert( _N_contigs == (int) _contig_lengths.size() );

  // Contig clusters must not contain more than the total number of contigs.
  assert( (int) clusters.SizeSum() <= _N_contigs );

  // The clusters must be non-empty and must contain contig IDs in the appropriate range
  int N_contigs_in_clusters = 0;
  for ( size_t i = 0; i < clusters.size(); i++ ) {
    assert( !clusters[i].empty() );
    assert( *clusters[i].rbegin() < _N_contigs );
    N_contigs_in_clusters += clusters[i].size();
  }


  // If ordering has been run, there must be orderings available.  (If there are no orderings, no ordering evaluation will be done.)
  if ( run_params._do_ordering ) {

    // There may not be as many orderings as there are clusters; all orderings may not have been calculated.  But there can't be more.
    assert( _N_orderings > 0 );
    assert( _N_orderings <= _N_clusters );

    // The trunk orderings and full orderings must be equal in number and must come from the same size clusters.
    assert( (int) _trunks.size() == _N_orderings );
    for ( int i = 0; i < _N_orderings; i++ ) {
      assert( orders[i].N_contigs() == trunks[i].N_contigs() );

      // All of the orderings must match their corresponding clustering by ID.
      if ( (int) clusters[i].size() != orders[i].N_contigs() ) PRINT3( i, clusters[i].size(), orders[i].N_contigs() );
      assert( (int) clusters[i].size() == orders[i].N_contigs() );

      // Trunk orderings must be subsets of full orderings.
      assert( trunks[i].N_contigs_used() <= orders[i].N_contigs_used() );
      for ( int j = 0; j < _trunks[i].N_contigs_used(); j++ )
	assert( orders[i].contig_used( trunks[i].contig_ID(j) ) );

    }
  }


  // Sanity checks for when there's a TrueMapping (i.e., not a reference-free assembly).
  if ( _true_mapping ) {

    // The species names must match.
    assert( _run_params._species == _true_mapping->species() );

    // The SAM file is aligned to the draft assembly, so it must show the same number of targets as the TrueMapping has queries (i.e., assembly contigs).
    assert( _N_contigs == _true_mapping->NQueries() );
  }

}






Reporter::~Reporter()
{
  if ( _true_mapping ) delete _true_mapping;
  delete _data;
}





void
Reporter::Eval() const
{
  // Evaluate the clustering.
  EvalClustering();

  // Evaluate all the orderings, if there are any.
  if ( _N_orderings > 0 ) {
    for ( int i = 0; i < _N_orderings; i++ ) {
      EvalOrderAccuracy( i, false ); // trunks
      EvalOrderAccuracy( i, true ); // full orderings
    }
  }



  // The evaluations past this point are all entirely reference-based.
  if ( _true_mapping == NULL ) return;


  // Evaluate the contig spacing in the full ordering for now.  This can't be done without a reference.
  if ( _orders.size() >= 3 && _orders[3].has_gaps() ) EvalGapSizes( 3, true ); // just one chromosome
  //DotplotOrderAccuracy( 3, true, false, "POA.3.full.no_interchrom" );


  return; // TEMP: skip dotplots for now (save time)

  if ( _N_orderings > 0 ) HistogramOrderingErrors( true, "HOE" );



  // Create colored dotplots to show the accuracy of the orderings.
  for ( int i = 0; i < _N_orderings; i++ ) {
    string i_str = boost::lexical_cast<string>(i);
    DotplotOrderAccuracy( i, false, false, "POA." + i_str + ".trunk.no_interchrom" );
    DotplotOrderAccuracy( i, true,  false, "POA." + i_str + ".full.no_interchrom" );
    DotplotOrderAccuracy( i, false, true,  "POA." + i_str + ".trunk.with_interchrom" );
    DotplotOrderAccuracy( i, true,  true,  "POA." + i_str + ".full.with_interchrom" );
  }

}






// Count which contigs are clustered, and evaluate the clustering by comparing it to the TrueMapping.
// Report on contigs that appear to be mis-clustered (i.e., they belong to a different chromosome than the plurality of other contigs in their cluster.)
// This code is similar to GenomeLinkMatrix::ReportMisclusteredContigs(), but doesn't go into as much detail because only topline numbers are expected, not
// troubleshooting help.
void
Reporter::EvalClustering() const
{
  cout << "EvalClustering" << endl;

  // Mark contigs that are in the cluster.
  _data->in_cluster.resize( _N_contigs );
  for ( int i = 0; i < _N_clusters; i++ )
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      _data->in_cluster[*it] = true;


  // If there's no reference, nothing more can be done.
  if ( _true_mapping == NULL ) return;


  _data->cluster_chrom.clear();

  int N_aligned_total = 0;
  int N_misclustered = 0;
  int64_t len_aligned_total = 0;
  int64_t len_misclustered = 0;


  // For each cluster, calculate the plurality chromosome and figure out the amount of contigs/sequence on that chromosome.
  for ( size_t i = 0; i < _clusters.size(); i++ ) {
    if ( _clusters[i].empty() ) continue;

    // Tallies of which chromosomes contain the contigs in this cluster.
    vector<int> aligns( _N_chroms, 0 );
    vector<int64_t> align_len( _N_chroms, 0 );

    // First, find the chromosome (if any) to which each contig in this cluster aligns.
    bool seen_aligned = false;
    int       N_aligned = 0,   N_unaligned = 0;
    int64_t len_aligned = 0, len_unaligned = 0;
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it ) {
      int chrID = _true_mapping->QTargetID(*it);

      // Kludgey code to collapse/merge chromosomes for the purposes of marking "correct" or "incorrect" chromosomes.
      //if ( chrID == 6 ) chrID = 2; // fly: chr2L -> chr2R
      //if ( chrID == 7 ) chrID = 4; // fly: chr3L -> chr3R
      //if ( chrID == 21 ) chrID = 18; // human: chr22 -> chr19
      //if ( chrID == 20 ) chrID = 19; // human: chr21 -> chr20
      //if ( chrID == 15 || chrID == 21 ) chrID = 18; // postfosmid: chr16, chr22 -> chr19
      //if ( chrID == 16 || chrID == 20 ) chrID = 19; // postfosmid: chr17, chr21 -> chr20

      if ( chrID == -1 ) {
	N_unaligned++;
	len_unaligned += _contig_lengths[*it];
      }
      else {
	N_aligned++;
	len_aligned += _contig_lengths[*it];
	aligns[chrID]++;
	align_len[chrID] += _contig_lengths[*it];
	seen_aligned = true;
      }
    }

    // Now find which chromosome that contains the plurality of these alignments (measured by length, not number of contigs.)
    int N_aligned_correct = 0;
    int64_t len_aligned_correct = 0;
    int cluster_ID = -1;
    for ( size_t j = 0; j < aligns.size(); j++ )
      if ( len_aligned_correct < align_len[j] ) {
        len_aligned_correct = align_len[j];
        N_aligned_correct   = aligns[j];
	cluster_ID = j;
      }

    // If none of the contigs in this cluster align, just go ahead and leave them all marked as incorrectly clustered.
    assert( ( cluster_ID == -1 ) == ( !seen_aligned ) );

    // Add to global tallies.
    N_aligned_total += N_aligned;
    N_misclustered  += N_aligned - N_aligned_correct;
    len_aligned_total += len_aligned;
    len_misclustered  += len_aligned - len_aligned_correct;

    _data->cluster_chrom.push_back( cluster_ID );
    _data->cluster_N_good  .push_back(   N_aligned_correct );
    _data->cluster_len_good.push_back( len_aligned_correct );
    _data->cluster_N_bad   .push_back(   N_aligned -   N_aligned_correct );
    _data->cluster_len_bad .push_back( len_aligned - len_aligned_correct );
    _data->cluster_N_unaligned  .push_back( N_unaligned );
    _data->cluster_len_unaligned.push_back( len_unaligned );
  }


}




// Count which contigs are ordered and oriented in an ordering, and evaluate the accuracy by comparing it to the TrueMapping.
// If full_order, then eval _orders[cluster_ID]; otherwise eval _trunks[cluster_ID].
void
Reporter::EvalOrderAccuracy( int cluster_ID, bool full_order ) const
{
  bool verbose = false;

  // Initialize the variables in the ReporterData struct.
  OrderingFlags & flags = full_order ? _data->order_flags : _data->trunk_flags;
  flags.used        .resize( _N_contigs );
  flags.high_quality.resize( _N_contigs );
  flags.unaligned   .resize( _N_contigs );
  flags.chr_mismatch.resize( _N_contigs );
  flags.order_error .resize( _N_contigs );
  flags.orient_error.resize( _N_contigs );


  // Convert this cluster to a local vector<> for easier random access.
  const vector<int> cluster = SetToVec( _clusters[cluster_ID] );


  const ContigOrdering & scaffold = full_order ? _orders[cluster_ID] : _trunks[cluster_ID];
  //scaffold.Print();
  if ( scaffold.N_contigs() != (int) cluster.size() ) PRINT4( cluster_ID, full_order, scaffold.N_contigs(), cluster.size() );
  assert( scaffold.N_contigs() == (int) cluster.size() );
  int N_contigs_used = scaffold.N_contigs_used();


  // Step through all the contigs in the scaffold.
  for ( int i = 0; i < N_contigs_used; i++ ) {

    // Find the contig position and orientation according to the scaffold.
    int contig0 = cluster[ scaffold.contig_ID(i) ];
    flags.used[contig0] = true;

    // Apply the quality filter (the threshold of which depends on what assembly this is) and determine whether or not this contig passes.
    // If no quality scores are given, the contig automatically passes.
    flags.high_quality[contig0] = !scaffold.has_Q_scores() || scaffold.contig_orient_Q(i) >= _run_params._report_quality_filter;


    // If we have no reference, we can't do any error marking.
    if ( _true_mapping == NULL ) continue;


    // To help call errors, find the previous two contig positions and orientations.
    int contig1 = ( i == 0 ? -1 : cluster[ scaffold.contig_ID(i-1) ] );
    int contig2 = ( i <= 1 ? -1 : cluster[ scaffold.contig_ID(i-2) ] );

    // Find the contig position and orientation of each of these contigs, according to the TrueMapping.
    int chrom0 = _true_mapping->QTargetID   (contig0);
    int start0 = _true_mapping->QTargetStart(contig0);
    int stop0  = _true_mapping->QTargetStop (contig0);
    //if ( _true_mapping->QRC(contig0) ) { int swap = start0; start0 = stop0; stop0 = swap; }

    int chrom1 = -1, start1 = -1, stop1 = -1;
    int chrom2 = -1, start2 = -1, stop2 = -1;
    if ( i >= 1 ) {
      chrom1 = _true_mapping->QTargetID   (contig1);
      start1 = _true_mapping->QTargetStart(contig1);
      stop1  = _true_mapping->QTargetStop (contig1);
      //if ( _true_mapping->QRC(contig1) ) { int swap = start1; start1 = stop1; stop1 = swap; }
    }
    if ( i >= 2 ) {
      chrom2 = _true_mapping->QTargetID   (contig2);
      start2 = _true_mapping->QTargetStart(contig2);
      stop2  = _true_mapping->QTargetStop (contig2);
      //if ( _true_mapping->QRC(contig2) ) { int swap = start2; start2 = stop2; stop2 = swap; }
    }



    // Mark contigs that aren't aligned to (canonical) chromosomes.
    if ( chrom0 == -1 ) flags.unaligned[contig0] = true;

    // Test for inter-chromosome errors.
    // If there is an error - i.e., this contig and its immediate predecessor both align to different canonical chromosomes - mark both contigs as mismatched.
    if ( chrom0 != -1 && chrom1 != -1 && chrom0 != chrom1 ) {
      flags.chr_mismatch[contig0] = true;
      flags.chr_mismatch[contig1] = true;
      if ( verbose ) cout << "ERROR\tCLUSTERING\tCLUSTER:\t" << cluster_ID << endl;
    }


    // Find the relative positions of this contig w.r.t. its predecessor(s).
    // pos =  1 means the subsequent contig is after  the predecessor on the chromosome;
    // pos = -1 means the subsequent contig is before the predecessor on the chromosome;
    // pos =  0 means they overlap.
    // Note that if i <= 1, one or both of these will evaluate to 1, and that's ok.
    int pos01 = ( start0 > start1 && stop0 > stop1 ) ? 1 : ( start0 < start1 && stop0 < stop1 ) ? -1 : 0;
    int pos12 = ( start1 > start2 && stop1 > stop2 ) ? 1 : ( start1 < start2 && stop1 < stop2 ) ? -1 : 0;




    // If this contig and its predecessor both align to the same canonical chromosome...
    if ( i >= 1 && chrom0 != -1 && chrom0 == chrom1 ) {

      // Look for contig orientation errors (i.e., this contig's orientation is not what we would expect based on its position w.r.t. its predecessor.)
      bool expected_rc = _true_mapping->QRC(contig0) ^ scaffold.contig_rc(i);
      if ( ( expected_rc && pos01 == 1 ) || ( !expected_rc && pos01 == -1 ) ) {
	flags.orient_error[contig0] = true;
	if ( verbose ) cout << "ERROR\tORIENTATION\tCLUSTER:\t" << cluster_ID << "\tCONTIG IN CLUSTER:\t" << i << "\tLENGTH:\t" << _contig_lengths[contig0] << "\tCONTIG GLOBAL ID:\t" << contig0 << "\tTRUE POSITION: " << _true_mapping->QTargetName(contig0) << ":" << _true_mapping->QTargetStart(contig0) << "-" << _true_mapping->QTargetStop(contig0) << endl;
      }
      //else
      //cout << "ORIENTING OK!\tQ = " << scaffold.contig_orient_Q(i) << endl;
    }

    // If this contig and its two predecessors both align to the same canonical chromosome...
    if ( i >= 2 && chrom0 != -1 && chrom0 == chrom1 && chrom0 == chrom2 ) {

      // Look for contig ordering errors: three contigs are adjacent in an ordering, but their relative position in the ordering doesn't match the relative
      // positions of their true alignments.
      if ( pos01 != pos12 && pos01 != 0 && pos12 != 0 ) {
	if ( verbose ) cout << "ERROR\tORDERING\tCLUSTER:\t" << cluster_ID << "\tCONTIG IN CLUSTER:\t" << i-1 << "\tLENGTH:\t" << _contig_lengths[contig1] << "\tCONTIG GLOBAL ID:\t" << contig1 << "\tTRUE POSITION: " << _true_mapping->QTargetName(contig1) << ":" << _true_mapping->QTargetStart(contig1) << "-" << _true_mapping->QTargetStop(contig1) << endl;
	flags.order_error[contig1] = true;
      }
    }
  }


  if ( verbose && _true_mapping != NULL )
    for ( int i = 0; i < N_contigs_used; i++ ) {
      int contig = cluster[ scaffold.contig_ID(i) ];
      if ( !_true_mapping->QMaps(contig) ) continue;
      if ( flags.order_error[contig] ) cout << "YESERROR:\t\t" << contig << "\tTRUE POSITION: " << _true_mapping->QTargetName(contig) << ":" << _true_mapping->QTargetStart(contig) << "-" << _true_mapping->QTargetStop(contig) << endl;
      if ( !flags.order_error[contig] ) cout << "NOERROR:\t\t" << contig << "\tTRUE POSITION: " << _true_mapping->QTargetName(contig) << ":" << _true_mapping->QTargetStart(contig) << "-" << _true_mapping->QTargetStop(contig) << endl;
    }
}





// Evaluate the contig spacing.
// This function should be called after calling EvalOrderAccuracy on the same ordering, so that the _data flags are filled.
// TODO: still working on getting contig spacing to work (ChromLinkMatrix::SpaceContigs) thus still working on this file
void
Reporter::EvalGapSizes( int cluster_ID, bool full_order ) const
{
  RequireReference(); // can't eval gap sizes without a reference

  // Convert this cluster to a local vector<> for easier random access.
  const vector<int> cluster = SetToVec( _clusters[cluster_ID] );

  // Get the scaffold.
  const ContigOrdering & scaffold = full_order ? _orders[cluster_ID] : _trunks[cluster_ID];
  int N_contigs_used = scaffold.N_contigs_used();


  // Get the flags indicating errors.  If EvalOrderAccuracy() hasn't been called on this ordering, these asserts will fail.
  const OrderingFlags & flags = full_order ? _data->order_flags : _data->trunk_flags;
  assert( !flags.used.empty() );
  assert( !flags.unaligned   .empty() );
  assert( !flags.chr_mismatch.empty() );
  assert( !flags.order_error .empty() );
  assert( !flags.orient_error.empty() );


  // Find the GC content and mappability for each contig.
  // We'll correlate these measures with the gap size estimation error of each contig, and attempt to determine the cause of errors.
  vector<double> contig_GC          = ParseTabDelimFile<double>( "../human/prefosmid/assembly.GCpct", 1 );
  vector<double> contig_mappability = ParseTabDelimFile<double>( "../human/prefosmid/assembly.MQ>0pct", 1 );
  assert( (int) contig_GC         .size() == _N_contigs );
  assert( (int) contig_mappability.size() == _N_contigs );
  //vector<double> contig_RF = ParseTabDelimFile<double>( "RFs.txt", 0 ); // TEMP
  vector<double> contig_LDE = ParseTabDelimFile<double>( "enrichments.txt", 0 ); // LDE = link density enrichment (on contig, compared to expectations)


  vector<double> dists, dists_est, errors, GCs1, GCs2, GCs_comb, maps1, maps2, maps_comb, LDEs1, LDEs2, LDEs_comb, RFs_comb;

  // Step through all the pairs of adjacent contigs in the scaffold.
  for ( int i = 0; i+1 < N_contigs_used; i++ ) {
    int contig1 = cluster[ scaffold.contig_ID(i) ];
    int contig2 = cluster[ scaffold.contig_ID(i+1) ];
    int dist_est = scaffold.gap_size(i);
    if ( dist_est == INT_MAX ) continue; // this can happen if the gaps were estimated by SpaceContigs() using a subset of the SAM files

    // Ignore pairs of adjacent contigs unless both of them map to the reference, on the same chromosome, with no apparent ordering or orienting errors.
    if ( flags.unaligned   [contig1] || flags.unaligned   [contig2] ) continue;
    if ( flags.chr_mismatch[contig1] || flags.chr_mismatch[contig2] ) continue;
    if ( flags.order_error [contig1] || flags.order_error [contig2] ) continue;
    if ( flags.orient_error[contig1] || flags.orient_error[contig2] ) continue;

    // TEMP: Only look at big contigs.
    // threshold types: OR < additive < multiplicative < AND
    // OR = 0.08; additive = 0.17; multiplicative = 0.20; AND = 0.28
    //if ( _contig_lengths[contig1] < 500000 && _contig_lengths[contig2] < 500000 ) continue;
    //if ( _contig_lengths[contig1] + _contig_lengths[contig2] < 700000 ) continue;
    //if ( _contig_lengths[contig1] / 1000 * _contig_lengths[contig2] / 1000 < 80000 ) continue;
    //if ( _contig_lengths[contig1] < 200000 || _contig_lengths[contig2] < 200000 ) continue;


    // Find the true distance between these contigs, according to their reference alignments.
    assert( _true_mapping->QTargetID( contig1 ) == _true_mapping->QTargetID( contig2 ) );
    int start1 = _true_mapping->QTargetStart( contig1 );
    int stop1  = _true_mapping->QTargetStop ( contig1 );
    int start2 = _true_mapping->QTargetStart( contig2 );
    int stop2  = _true_mapping->QTargetStop ( contig2 );

    int dist = max( start2 - stop1, start1 - stop2 );
    if ( dist < 100 ) dist = 100; // if the contigs overlap, define the distance as 100


    //PRINT7( contig1, contig2, start1, stop1, start2, stop2, dist );
    cout << "EvalGapSizes: Contigs #" << i << ", " << i+1 << " in the scaffold\t(global IDs: " << contig1 << ", " << contig2 << ", lengths = " << _contig_lengths[contig1] << ", " << _contig_lengths[contig2] << ")\tare at a true distance of\t" << dist << "\tvs. predicted distance of\t" << dist_est << endl;

    // TODO: compare these true distances to the distances found in ChromLinkMatrix::SpaceContigs, and investigate differences
    // grep COMPARISON b | cut -f2,3 > c ; QuickDotplot c
    cout << "COMP1\t" << log(dist+1) << "\t" << log(dist_est+1) << endl;

    double error = double( dist_est+1 ) / (dist+1); // "error" = the multiplicative factor of wrongness by which the estimate exceeds the true distance


    double GC1  = contig_GC         [contig1], GC2  = contig_GC[contig2];
    double map1 = contig_mappability[contig1], map2 = contig_mappability[contig2];
    double LDE1 = contig_LDE        [contig1], LDE2 = contig_LDE[contig2];
    //double RF1  = contig_RF         [contig1], RF2  = contig_RF[contig2];
    double GC_combined  = ( GC1  * _contig_lengths[contig1] + GC2  * _contig_lengths[contig2] ) / ( _contig_lengths[contig1] + _contig_lengths[contig2] );
    double map_combined = ( map1 * _contig_lengths[contig1] + map2 * _contig_lengths[contig2] ) / ( _contig_lengths[contig1] + _contig_lengths[contig2] );
    double LDE_combined = ( LDE1 * _contig_lengths[contig1] + LDE2 * _contig_lengths[contig2] ) / ( _contig_lengths[contig1] + _contig_lengths[contig2] );
    //double RF_combined  = ( RF1  * _contig_lengths[contig1] + RF2  * _contig_lengths[contig2] ) / ( _contig_lengths[contig1] + _contig_lengths[contig2] );
    LDE_combined = log( LDE_combined );



    cout << "COMP2\t" <<  GC_combined << '\t' << log(error) << endl;
    cout << "COMP3\t" << map_combined << '\t' << log(error) << endl;
    cout << "COMP4\t" << LDE_combined << '\t' << log(error) << endl;
    //cout << "COMP4\t" <<  RF_combined << '\t' << log(error+1) << endl;
    dists_est.push_back( log(dist_est+1) );
    dists    .push_back( log(dist+1) );
    errors   .push_back( log(error) );
    GCs1     .push_back( GC1 );
    GCs2     .push_back( GC2 );
    GCs_comb .push_back( GC_combined );
    maps1    .push_back( map1 );
    maps2    .push_back( map2 );
    maps_comb.push_back( map_combined );
    LDEs1    .push_back( LDE1 );
    LDEs2    .push_back( LDE2 );
    LDEs_comb.push_back( LDE_combined );
    //RFs_comb .push_back( RF_combined );
  }

  int N_data = dists.size();
  PRINT( N_data );
  PRINT( PearsonCorrelation( dists, dists_est ) );
  PRINT( PearsonCorrelation( errors, GCs_comb ) );
  PRINT( PearsonCorrelation( errors, maps_comb ) );
  PRINT( PearsonCorrelation( errors, LDEs_comb ) );
  //PRINT( PearsonCorrelation( errors, RFs_comb ) );

}




// ReportChartWithReference: Print out all the ReporterData info in a pretty chart.  This is the version for an assembly with a reference available.
void
Reporter::ReportChartWithReference() const
{
  RequireReference();

  assert( _N_clusters == (int) _data->cluster_chrom.size() );
  assert( _N_clusters == (int) _data->cluster_N_good.size() );
  assert( _N_clusters == (int) _data->cluster_N_bad.size() );
  assert( _N_clusters == (int) _data->cluster_N_unaligned.size() );
  assert( _N_clusters == (int) _data->cluster_len_good.size() );
  assert( _N_clusters == (int) _data->cluster_len_bad.size() );
  assert( _N_clusters == (int) _data->cluster_len_unaligned.size() );


  const double epsilon = 1.0e-12; // small number to add to denominators to avoid division by 0

  // Open the ReportChart output file.
  string report_chart_file = _run_params._out_dir + "/REPORT.txt";
  cout << "Writing a ReportChart to " << report_chart_file << endl;
  ofstream out( report_chart_file.c_str() );


  // Write the run parameters to the file.
  _run_params.PrintParams( out );
  out << "\n\n\n\n\n";


  // Print basic info about this assembly.
  out << "ReportChart!" << endl << endl;
  out << "Info about input assembly:" << endl;
  if ( _run_params._sim_bin_size > 0 ) out << "REFERENCE GENOME CHOPPED UP INTO " << _run_params._sim_bin_size << "-BP BINS" << endl;
  else                                 out << "DE NOVO ASSEMBLY, reference genome available for validation" << endl;
  out << "Species: " << _run_params._species << endl;
  out << "N chromosomes in reference, including non-canonical:\t" << _N_chroms << endl;
  out << "N contigs:\t" << _N_contigs << "\t\tTotal length:\t" << _total_contig_length << "\t\tN50:\t" << N50( _contig_lengths ) << endl;
  out << "N clusters (derived):\t\t" << _N_clusters << endl;
  out << "N non-singleton clusters:\t" << N_non_singleton_clusters() << endl;
  out << "N orderings found:\t\t" << _N_orderings << endl;
  out << endl << endl;


  out << "############################\n";
  out << "#                          #\n";
  out << "#    CLUSTERING METRICS    #\n";
  out << "#                          #\n";
  out << "############################\n\n\n";

  out << setprecision(4);

  // Derive and print basic numbers about clustering.
  int       N_in_clusters = _data->in_cluster.count();
  int64_t len_in_clusters = contig_length_if( _data->in_cluster );
  double   pct_N_in_clusters = 100.0 *   N_in_clusters / _N_contigs;
  double pct_len_in_clusters = 100.0 * len_in_clusters / _total_contig_length;
  out << "Number of contigs in clusters:\t" <<   N_in_clusters << "\t\t(" << pct_N_in_clusters << "% of all contigs)" << endl;
  out << "Length of contigs in clusters:\t" << len_in_clusters << "\t(" << pct_len_in_clusters << "% of all sequence length)" << endl;
  out << endl;

  int N_singletons = 0;
  int64_t len_singletons = 0;

  // Print a chart with info about each cluster, including what portion of the cluster (length, contigs) aligns to the plurality chromosome or fails to align.
  string horiz_line = "+---------+----------------+---------+---------------+------------------+-------------+-------------------+----------------------+\n";
  out << horiz_line;
  out << "| CLUSTER |   PLURALITY    | Number of contigs in cluster...            | Length of contigs in cluster...                        |\n";
  out << "| NUMBER  |   CHROMOSOME   |  TOTAL  |   UNALIGNED   | WRONG CHROMOSOME |    TOTAL    |     UNALIGNED     |   WRONG CHROMOSOME   |\n";
  out << horiz_line;
  for ( int i = 0; i < _N_clusters; i++ ) {

    // Skip singleton clusters.
    if ( _clusters[i].size() == 1 ) {
      N_singletons++;
      len_singletons += _contig_lengths[ *(_clusters[i].begin()) ];
      continue;
    }


    int chrom_ID = _data->cluster_chrom[i];
    // Re-format the chromosome name to appear nicely in the output chart: <= 14 characters and centered in whitespace.
    string chrom_name = ( chrom_ID == -1 ? "none" : _true_mapping->TargetName(chrom_ID) );
    chrom_name = chrom_name.substr(0,14);
    for ( size_t j = 0; j + chrom_name.size() < 14; j++ ) chrom_name += ' '; // note that chrom_name.size() increases with each iteration


    // Calculate the total cluster length.
    boost::dynamic_bitset<> in_cluster( _N_contigs );
    int cluster_N = 0;
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it ) {
      cluster_N++;
      in_cluster[*it] = true;
    }
    int64_t cluster_len = contig_length_if( in_cluster );

    // For each chromosome, calculate percentages, and print a line to the chart.
    double pct_N_bad         = 100.0 * _data->cluster_N_bad  [i] / ( _data->cluster_N_good  [i] + _data->cluster_N_bad  [i] + epsilon );
    double pct_len_bad       = 100.0 * _data->cluster_len_bad[i] / ( _data->cluster_len_good[i] + _data->cluster_len_bad[i] + epsilon );
    double pct_N_unaligned   = 100.0 * _data->cluster_N_unaligned[i] / cluster_N;
    double pct_len_unaligned = 100.0 * _data->cluster_len_unaligned[i] / cluster_len;

    sprintf( _LINE, "|%6d   | %14s |%7ld  |%5d (%6.2f%%)| %5d  (%6.2f%%) | %11ld |%9ld (%6.3f%%)| %9ld  (%6.3f%%) |\n", i, chrom_name.c_str(), _clusters[i].size(), _data->cluster_N_unaligned[i], pct_N_unaligned, _data->cluster_N_bad[i], pct_N_bad, cluster_len, _data->cluster_len_unaligned[i], pct_len_unaligned, _data->cluster_len_bad[i], pct_len_bad );
    out << _LINE;
  }
  out << horiz_line;

  // Print a final line with totals.
  int     total_N_good   = accumulate( _data->cluster_N_good  .begin(), _data->cluster_N_good  .end(), 0 );
  int64_t total_len_good = accumulate( _data->cluster_len_good.begin(), _data->cluster_len_good.end(), int64_t(0) );
  int     total_N_bad    = accumulate( _data->cluster_N_bad   .begin(), _data->cluster_N_bad   .end(), 0 );
  int64_t total_len_bad  = accumulate( _data->cluster_len_bad .begin(), _data->cluster_len_bad .end(), int64_t(0) );
  int     total_N_unaligned   = accumulate( _data->cluster_N_unaligned  .begin(), _data->cluster_N_unaligned.end(), 0 );
  int64_t total_len_unaligned = accumulate( _data->cluster_len_unaligned.begin(), _data->cluster_len_unaligned.end(), int64_t(0) );
  double pct_N_bad   = 100.0 * total_N_bad   / ( total_N_good   + total_N_bad + epsilon );
  double pct_len_bad = 100.0 * total_len_bad / ( total_len_good + total_len_bad + epsilon );
  double pct_N_unaligned   = 100.0 * total_N_unaligned   / N_in_clusters;
  double pct_len_unaligned = 100.0 * total_len_unaligned / len_in_clusters;

  sprintf( _LINE, "|       TOTAL              |%7d  |%5d (%6.2f%%)| %5d  (%6.2f%%) | %11ld |%9ld (%6.3f%%)| %9ld  (%6.3f%%) |\n",
	   N_in_clusters, total_N_unaligned, pct_N_unaligned, total_N_bad, pct_N_bad, len_in_clusters, total_len_unaligned, pct_len_unaligned, total_len_bad, pct_len_bad );
  out << _LINE;
  out << horiz_line;

  if ( N_singletons > 0 )
    out << endl << "NOTE: The above chart omits " << N_singletons << " `singleton' clusters, with only one contig each (total length = " << len_singletons << ")" << endl;

  out << endl << endl;



  if ( _N_orderings == 0 ) return;


  out << "############################\n";
  out << "#                          #\n";
  out << "#     ORDERING METRICS     #\n";
  out << "#                          #\n";
  out << "############################\n\n\n";


  // Derive and print basic numbers about ordering.
  ReportChartOrderingPercentages(out);

  // This horizontal line matches the chart lines in ReportChartOrderingErrors().
  horiz_line = "+-------------+---------------------------+--------------+--------------+\n";

  // Print the chart for ordering errors.
  out << horiz_line;
  out << "| TRUNK/FULL? |  ERROR TYPE               |  % CONTIGS   |   % LENGTH   |\n";
  out << horiz_line;
  ReportChartOrderingErrors( false, out ); // trunk only
  out << horiz_line;
  ReportChartOrderingErrors( true, out ); // full ordering
  out << horiz_line;
  out << endl << endl;


  // Print a footer with useful definitions.
  out << "###################################\n\n";
  out << "Definitions:" << endl;
  out << "'Plurality chromosome':\t\tThe plurality chromosome for a cluster is the chromosome that contains the most contigs (by length!) in that cluster.\n\t\t\t\tContigs are considered mis-clustered if they align to a canonical chromosome other than the plurality chromosome.\n";
  out << "'Inter-chromosome error':\tTwo contigs are adjacent in an ordering but align to different canonical chromosomes.\n";
  out << "'Ordering error':\t\tThree contigs are adjacent in an ordering, but their relative position in the ordering doesn't match the relative positions of their true alignments.\n";
  out << "'Orientation error':\t\tA contig's assigned orientation in its ordering doesn't match the orientation it should have, imputed from its and its neighbors' true alignment position.\n";
  out << endl << endl;


  out.close();
}

// ReportChartNoReference: Print out all the ReporterData info in a pretty chart.  This is the version for a reference-free assembly.
void
Reporter::ReportChartNoReference() const
{
  assert( !_run_params._use_ref );

  // Open the ReportChart output file.
  string report_chart_file = _run_params._out_dir + "/REPORT.txt";
  cout << "Writing a ReportChart to " << report_chart_file << endl;
  ofstream out( report_chart_file.c_str() );


  // Write the run parameters to the file.
  _run_params.PrintParams( out );
  out << "\n\n\n\n\n";


  // Print basic info about this assembly.
  out << "ReportChart!" << endl << endl;
  out << "Info about input assembly:" << endl;
  out << "DE NOVO ASSEMBLY, with no reference genome (less validation available)" << endl;
  out << "Species: " << _run_params._species << endl;
  out << "N contigs:\t" << _N_contigs << "\t\tTotal length:\t" << _total_contig_length << "\t\tN50:\t" << N50( _contig_lengths ) << endl;
  out << "N clusters (derived):\t" << _N_clusters << endl;
  out << "N non-singleton clusters:\t" << N_non_singleton_clusters() << endl;
  out << "N orderings found:\t" << _N_orderings << endl;
  out << endl << endl;


  out << "############################\n";
  out << "#                          #\n";
  out << "#    CLUSTERING METRICS    #\n";
  out << "#                          #\n";
  out << "############################\n\n\n";

  out << setprecision(4);

  // Derive and print basic numbers about clustering.
  int       N_in_clusters = _data->in_cluster.count();
  int64_t len_in_clusters = contig_length_if( _data->in_cluster );
  double   pct_N_in_clusters = 100.0 *   N_in_clusters / _N_contigs;
  double pct_len_in_clusters = 100.0 * len_in_clusters / _total_contig_length;
  out << "Number of contigs in clusters:\t" <<   N_in_clusters << "\t\t(" << pct_N_in_clusters << "% of all contigs)" << endl;
  out << "Length of contigs in clusters:\t" << len_in_clusters << "\t(" << pct_len_in_clusters << "% of all sequence length)" << endl;
  out << endl;

  // Print a chart with info about each cluster.  It's not a ton of info because there's no reference.
  string horiz_line = "+----------+-----------+-------------+\n";
  out << horiz_line;
  out << "|  CLUSTER | NUMBER OF |  LENGTH OF  |\n";
  out << "|  NUMBER  |  CONTIGS  |   CONTIGS   | \n";
  out << horiz_line;
  for ( int i = 0; i < _N_clusters; i++ ) {

    // Calculate the total cluster length.
    boost::dynamic_bitset<> in_cluster( _N_contigs );
    for ( set<int>::const_iterator it = _clusters[i].begin(); it != _clusters[i].end(); ++it )
      in_cluster[*it] = true;
    int64_t cluster_len = contig_length_if( in_cluster );

    sprintf( _LINE, "| %6d   | %7ld   | %11ld |\n", i, _clusters[i].size(), cluster_len );
    out << _LINE;
  }
  out << horiz_line;

  // Print a final line with totals.
  sprintf( _LINE, "|   TOTAL  | %7d   | %11ld |\n", N_in_clusters, len_in_clusters );
  out << _LINE;
  out << horiz_line;

  out << endl << endl;



  if ( _N_orderings == 0 ) return;


  out << "############################\n";
  out << "#                          #\n";
  out << "#     ORDERING METRICS     #\n";
  out << "#                          #\n";
  out << "############################\n\n\n";

  // Derive and print basic numbers about ordering.
  ReportChartOrderingPercentages(out);

  out.close();
}







// RequireReference: Throw a verbose error if _true_mapping == NULL.  This should always be called before anything that uses _true_mapping.
void
Reporter::RequireReference() const
{
  if ( _true_mapping != NULL ) return;

  cerr << "ERROR: Reporter module tried to call an evaluation function requiring a TrueMapping object, but there's no TrueMapping because this is a reference-free assembly (USE_REFERENCE=0)." << endl;
  exit(1);
}




// HistogramOrderingErrors: Make a histogram showing the odds of a contig being mis-ordered, as a function of the contig's length.
void
Reporter::HistogramOrderingErrors(const bool full_order,
                                  const string &file_head) const {
  RequireReference();
  // For contig 'lengths', either use actual lengths or the number of RE sites.
  bool use_RE_lens = true;
  vector<int> lens = use_RE_lens ? ParseTabDelimFile<int>(_run_params.DraftContigRESitesFilename(), 1) : _contig_lengths;
  // If this assert fails, the input contig_lengths file is inconsistent with the original dataset (wrong number of contigs).
  assert((int)lens.size() == _N_contigs);
  cout << ": HistogramOrderingErrors with file_head = " << file_head << endl;
  const OrderingFlags &flags = full_order ? _data->order_flags : _data->trunk_flags;

  // Make logarithmic-sized bins for contig lengths.  Bin #N will contain data for contigs with sizes in the range 2^N <= size < 2^(N+1).
  int longest_contig = *(max_element( lens.begin(), lens.end()));
  int N_bins = 1;
  while (longest_contig >>= 1) {
    N_bins++;
  }
  // PRINT2(longest_contig, N_bins);

  vector<int64_t> N_total(N_bins, 0), len_total(N_bins, 0), N_order_error(N_bins, 0), N_orient_error(N_bins, 0), N_o_and_o_error(N_bins, 0);

  // Now loop over all contigs and gather data.
  for (int i = 0; i < _N_contigs; i++) {
    // Skip contigs that aren't used or that don't align; these don't have the option of being ordering errors.
    if (!flags.used[i] || flags.unaligned[i]) {
      continue;
    }

    // For each contig, look up its length, and determine the bin it belongs in.
    int len = lens[i];
    int len_bin = 0;
    while (len >>= 1) {
      len_bin++;
    }
    assert(len_bin < N_bins);

    // Look up the contig's error(s) and fill the data structures.
    N_total[len_bin]++;
    len_total[len_bin] += lens[i];
    if (flags.order_error[i])  N_order_error[len_bin]++;
    if (flags.orient_error[i]) N_orient_error[len_bin]++;
    if (flags.order_error[i] && flags.orient_error[i]) {
      N_o_and_o_error[len_bin]++;
    }
  }

  // Make the histogram.
  ofstream out( file_head.c_str() );
  out << "#log2(contig_N_RE_sites)\tN_contigs\ttotal_contig_N_RE_sites\tpct_ok\tpct_order_error\tpct_orient_error\tpct_order_and_orient_errors" << endl;

  for (int i = 0; i < N_bins; i++) {
    double pct_order_error = 100.0 * N_order_error[i] / N_total[i];
    double pct_orient_error = 100.0 * N_orient_error[i] / N_total[i];
    double pct_o_and_o_error = 100.0 * N_o_and_o_error[i] / N_total[i];
    // PRINT5( i, N_total[i], pct_order_error, pct_orient_error, pct_o_and_o_error );
    out << i
	<< '\t' << N_total[i]
	<< '\t' << len_total[i]
	<< '\t' << (100 - pct_order_error - pct_orient_error + pct_o_and_o_error)
	<< '\t' << (pct_order_error - pct_o_and_o_error)
	<< '\t' << (pct_orient_error - pct_o_and_o_error)
	<< '\t' << pct_o_and_o_error
	<< endl;
  }
  out.close();
}

// DotplotOrderAccuracy: Make a QuickDotplot of the ordering, highlighting contigs that are mis-ordered.  Uses the output of EvalOrderAccuracy.
void
Reporter::DotplotOrderAccuracy( const int ordering_ID, const bool full_order, const bool plot_interchrom_contigs, const string & file_head ) const
{
  RequireReference();

  assert( ordering_ID >= 0 && ordering_ID < _N_orderings );

  cout << "DotplotOrderAccuracy with ordering_ID = " << ordering_ID << ", full_order = " << int(full_order) << "; file_head = " << file_head << endl;

  // Pick the ordering to analyze.
  const ContigOrdering & scaffold = full_order ? _orders[ordering_ID] : _trunks[ordering_ID];
  const OrderingFlags & flags = full_order ? _data->order_flags : _data->trunk_flags;
  //scaffold.Print();


  // Get the full ordering of the contigs on reference.
  vector<int> contig_order = _true_mapping->QueriesToGenomeOrder();


  // Convert this ordering to a local vector<> for easier random access.
  vector<int> cluster = SetToVec( _clusters[ordering_ID] );
  assert( scaffold.N_contigs() == (int) cluster.size() );
  int N_contigs_used = scaffold.N_contigs_used();

  // Open the dotplot file for writing.  This filename is the hard-wired input to the script QuickDotplot.POA.R.
  ofstream out( "QuickDotplot.POA.txt", ios::out );

  // Step through all the contigs in the scaffold.  For each contig, determine its coordinates and write them to the file.
  for ( int i = 0; i < N_contigs_used; i++ ) {

    // Find the contig position and orientation according to the scaffold.
    int contig = cluster[ scaffold.contig_ID(i) ];

    // Find this contig's true position on the genome (i.e., its index among all sorted contigs, including those not in this cluster.
    int true_chrom = _true_mapping->QTargetID( contig );
    int true_pos = contig_order[contig];

    //PRINT5( contig, flags.used[contig], flags.chr_mismatch[contig], flags.order_error[contig], flags.orient_error[contig] );
    //PRINT3( i, contig, contig_order[contig] );

    // If this contig isn't on the plurality chromosome, we may want to skip it, to prevent the output dotplot from gettin' all scrunched.
    if ( !plot_interchrom_contigs )
      if ( true_chrom != _data->cluster_chrom[ordering_ID] ) continue;


    // Find which of the 3 types of errors (if any) are occurring at this contig.
    // There are 3 types of errors, although there are fewer than 2^3 = 8 possible error states, because chrom mismatch can't coexist with the other two.
    int e1 = flags.chr_mismatch[contig];
    int e2 = flags.order_error [contig];
    int e3 = flags.orient_error[contig];
    int error_state = (e1 << 2) + (e2 << 1) + e3;

    // Convert this set of error states into a descriptive string, used for labeling in the dotplot.
    string error_str;
    switch( error_state ) {
    case 0: error_str = "awesome"; break;
    case 1: error_str = "err_orient"; break;
    case 2: error_str = "err_order"; break;
    case 3: error_str = "err_order_orient"; break;
    case 4: error_str = "err_interchrom"; break;
    case 5: error_str = "err_interchrom_orient"; break;
    case 6: error_str = "err_interchrom_order"; break;
    case 7: error_str = "err_interchrom_order_orient"; break;
    }


    // Plot a dotplot point.  The x-axis is the true position; the y-axis is the position in this scaffold.
    out << true_pos << '\t' << i << '\t' << error_str << endl;

    //PRINT5( i, contig, true_chrom, true_pos, error_str );
  }

  out.close();


  // Run the QuickDotplot script to generate a dot plot image, which gets placed at out/<file_head>.jpg.
  system( "QuickDotplot.POA.R" ); // hard-wired to input the file 'QuickDotplot.POA.txt' and create 'out/POA.jpg'
  if ( file_head != "POA" ) {
    string cmd = "cp QuickDotplot.POA.txt " + file_head + ".txt";
    //system( cmd.c_str() );
    cmd = "cp out/POA.jpg out/" + file_head + ".jpg";
    system( cmd.c_str() );
  }
}






// Helper function for the ReportChart() functions.
void
Reporter::ReportChartOrderingPercentages( ostream & out ) const
{
  int       N_in_clusters = _data->in_cluster.count();
  int64_t len_in_clusters = contig_length_if( _data->in_cluster );

  int N_in_orders    = _data->order_flags.used.count();
  int N_in_trunks    = _data->trunk_flags.used.count();
  int N_HQ_in_orders = _data->order_flags.high_quality.count();
  int N_HQ_in_trunks = _data->trunk_flags.high_quality.count();
  int64_t len_in_orders    = contig_length_if( _data->order_flags.used );
  int64_t len_in_trunks    = contig_length_if( _data->trunk_flags.used );
  int64_t len_HQ_in_orders = contig_length_if( _data->order_flags.high_quality );
  int64_t len_HQ_in_trunks = contig_length_if( _data->trunk_flags.high_quality );
  double      pct_N_in_orders = 100.0 * N_in_orders / N_in_clusters;
  double      pct_N_in_trunks = 100.0 * N_in_trunks / N_in_orders;
  double   pct_N_HQ_in_orders = 100.0 * N_HQ_in_orders / N_in_orders;
  double   pct_N_HQ_in_trunks = 100.0 * N_HQ_in_trunks / N_in_trunks;
  double    pct_len_in_orders = 100.0 * len_in_orders / len_in_clusters;
  double    pct_len_in_trunks = 100.0 * len_in_trunks / len_in_orders;
  double pct_len_HQ_in_orders = 100.0 * len_HQ_in_orders / len_in_orders;
  double pct_len_HQ_in_trunks = 100.0 * len_HQ_in_trunks / len_in_trunks;
  out << "Number of contigs in orderings:\t" <<   N_in_orders << "\t\t(" << pct_N_in_orders << "% of all contigs in clusters, "
      << ( 100.0 * N_in_orders/ _N_contigs )  << "% of all contigs)" << endl;
  out << "Length of contigs in orderings:\t" << len_in_orders << "\t(" << pct_len_in_orders << "% of all length in clusters, "
      << ( 100.0 * len_in_orders/ _total_contig_length )  << "% of all sequence length)" << endl;
  out << "Number of contigs in trunks:\t"    <<   N_in_trunks << "\t\t(" << pct_N_in_trunks << "% of contigs in orderings)" << endl;
  out << "Length of contigs in trunks:\t"    << len_in_trunks << "\t(" << pct_len_in_trunks << "% of length in orderings)" << endl;
  out << endl;
  out << "Fraction of contigs in orderings with high orientation quality:\t" << N_HQ_in_orders << " (" << pct_N_HQ_in_orders << "%), with length " << len_HQ_in_orders << " (" << pct_len_HQ_in_orders << "%)" << endl;
  out << "Fraction of contigs in trunks    with high orientation quality:\t" << N_HQ_in_trunks << " (" << pct_N_HQ_in_trunks << "%), with length " << len_HQ_in_trunks << " (" << pct_len_HQ_in_trunks << "%)" << endl;
  out << endl;
}




// Helper function for ReportChartWithReference().
void
Reporter::ReportChartOrderingErrors( const bool full_order, ostream & out ) const
{
  // Get the dataset: is this the full order or not?
  OrderingFlags flags = full_order ? _data->order_flags : _data->trunk_flags;
  string label = full_order ? "FULL ORDER" : "TRUNK";

  // Derive numbers.
  int N_contigs_in_chroms = flags.used.count() - flags.unaligned.count();
  int64_t len_contigs_in_chroms = contig_length_if( flags.used ) - contig_length_if( flags.unaligned );
  double   pct_N_chr_mismatches = 100.0 * flags.chr_mismatch.count() / N_contigs_in_chroms;
  double pct_len_chr_mismatches = 100.0 * contig_length_if( flags.chr_mismatch ) / len_contigs_in_chroms;
  double   pct_N_order_errors = 100.0 * flags.order_error.count() / N_contigs_in_chroms;
  double pct_len_order_errors = 100.0 * contig_length_if( flags.order_error ) / len_contigs_in_chroms;
  double   pct_N_orient_errors = 100.0 * flags.orient_error.count() / N_contigs_in_chroms;
  double pct_len_orient_errors = 100.0 * contig_length_if( flags.orient_error ) / len_contigs_in_chroms;
  double   pct_N_HQ_order_errors = 100.0 * AND( flags.high_quality, flags.order_error ).count() / flags.high_quality.count();
  double pct_len_HQ_order_errors = 100.0 * contig_length_if( AND( flags.high_quality, flags.order_error ) ) / contig_length_if( flags.high_quality );
  double   pct_N_HQ_orient_errors = 100.0 * AND( flags.high_quality, flags.orient_error ).count() / flags.high_quality.count();
  double pct_len_HQ_orient_errors = 100.0 * contig_length_if( AND( flags.high_quality, flags.orient_error ) ) / contig_length_if( flags.high_quality );


  // Print the chart, not including horizontal lines.
  sprintf( _LINE, "|  %-10s | Inter-chromosome errors   |    %6.2f%%   |    %6.2f%%   |\n", label.c_str(), pct_N_chr_mismatches, pct_len_chr_mismatches );
  out << _LINE;
  sprintf( _LINE, "|  %-10s | Ordering errors           |    %6.2f%%   |    %6.2f%%   |\n", label.c_str(), pct_N_order_errors, pct_len_order_errors );
  out << _LINE;
  sprintf( _LINE, "|  %-10s | Ordering errors    (HQ)   |    %6.2f%%   |    %6.2f%%   |\n", label.c_str(), pct_N_HQ_order_errors, pct_len_HQ_order_errors );
  out << _LINE;
  sprintf( _LINE, "|  %-10s | Orientation errors        |    %6.2f%%   |    %6.2f%%   |\n", label.c_str(), pct_N_orient_errors, pct_len_orient_errors );
  out << _LINE;
  sprintf( _LINE, "|  %-10s | Orientation errors (HQ)   |    %6.2f%%   |    %6.2f%%   |\n", label.c_str(), pct_N_HQ_orient_errors, pct_len_HQ_orient_errors );
  out << _LINE;
}



int
Reporter::N_non_singleton_clusters() const
{
  assert( _N_clusters == (int) _clusters.size() );

  int N = 0;
  for ( int i = 0; i < _N_clusters; i++ )
    if ( _clusters[i].size() > 1 )
      N++;

  return N;
}



// contig_length_if: Return the length of all contigs with their bit marked as 'true' in the input bitset.
int64_t
Reporter::contig_length_if( const boost::dynamic_bitset<> & bits ) const
{
  assert( (int) bits.size() == _N_contigs );

  int64_t length = 0;
  for ( int i = 0; i < _N_contigs; i++ )
    if ( bits[i] )
      length += _contig_lengths[i];

  return length;
}






// Make a heatmap of the entire result.  That is, take contigs that have been both clustered and ordered by Lachesis, and plot an NxN heatmap of the Hi-C link
// densities between these contigs.
// Only the PLOT_N longest contigs will be considered for plotting (less, if total number of contigs < PLOT_N) and these will only be plotted if they are
// clustered and ordered.  If USE_RES = true, length (and Hi-C link normalization) is measured by RE sites; otherwise it's measured by contig length.
// This is a useful reference-free visual evaluation.
// Method:
// 1. Load the ClusterVec and the ContigOrderings to determine how these contigs have been clustered and ordered by Lachesis.
// 2. Convert these data structures into a single vector that maps the original order of contigs onto an order that describes how Lachesis has ordered them.
// 3. Load a GenomeLinkMatrix to get the quantity of Hi-C links between all pairs of contigs.
// 4. Make a heatmap of the GLM's data, with the contigs reordered as in Step 2.
// Uses the R script heatmap.MWAH.R.  Makes an output file heatmap.txt and an auxiliary output file heatmap.chrom_breaks.txt.
void
MakeWholeAssemblyHeatmap( const RunParams & run_params, const int PLOT_N, const bool USE_RES )
{
  cout << "MakeWholeAssemblyHeatmap!" << endl;


  // 1. Load the ClusterVec and the ContigOrderings to determine how these contigs have been clustered and ordered by Lachesis.
  ClusterVec clusters( run_params._out_dir + "/main_results/clusters.txt" );

  vector<ContigOrdering> orders;
  for ( size_t i = 0; i < clusters.size(); i++ )
    orders.push_back( ContigOrdering( run_params._out_dir + "/main_results/group" + boost::lexical_cast<string>(i) + ".ordering" ) );

  // Also get contig lengths.
  vector<int> contig_lengths = TargetLengths( run_params._SAM_files[0] );
  int N_contigs = contig_lengths.size();
  vector<int> contig_N_REs = ParseTabDelimFile<int>( run_params.DraftContigRESitesFilename(), 1 );
  assert( N_contigs == (int) contig_lengths.size() );

  // Find which contigs are the PLOT_N longest.  Note that even these contigs will not be plotted if they are not also clustered and ordered.
  vector<int> contig_lengths_sorted = USE_RES ? contig_N_REs : contig_lengths;
  sort( contig_lengths_sorted.begin(), contig_lengths_sorted.end(), greater<int>() );
  int length_cutoff = contig_lengths_sorted[ PLOT_N < N_contigs ? PLOT_N : N_contigs - 1 ];

  vector<bool> is_long( N_contigs, false );
  for ( int i = 0; i < N_contigs; i++ ) {
    int len = USE_RES ? contig_N_REs[i] : contig_lengths[i];
    if ( len >= length_cutoff )
      is_long[i] = true;
  }



  // 2. Convert these data structures into a single vector that maps the original order of contigs onto an order that describes how Lachesis has ordered them.
  // This necessarily entails some arbitrary symmetry-breaking because there's no inherent ordering to the clusters or orientation for each ordering.
  // We resolve these by taking the default: the ClusterVec is already sorted by decreasing total contig length, and the ContigOrderings are canonicalized.
  // Note that contigs that are clustered but not ordered, or aren't clustered at all, are left out of the new ordering.  This will make a cleaner heatmap.
  vector<int> old_to_new( N_contigs, -1 ); // lookup table
  vector<int> new_to_old( N_contigs, -1 ); // reverse lookup table
  int new_ID = 0;
  int64_t total_len = 0;

  // Make an output file that delineates chromosomes.
  ofstream out2( "heatmap.chrom_breaks.txt" );
  out2 << "0" << endl;

  for ( size_t i = 0; i < clusters.size(); i++ ) {

    //cout << "New cluster starts at #" << new_ID << endl;

    // Convert this cluster from a set<int> into a vector<int> for easy random access.
    vector<int> cluster;
    for ( set<int>::const_iterator it = clusters[i].begin(); it != clusters[i].end(); ++it )
      cluster.push_back( *it );

    // Loop over all the contigs in this cluster ordering, and add them to the lookup table.
    for ( int j = 0; j < orders[i].N_contigs_used(); j++ ) {

      // The ContigOrdering stores the ID of each contig within the local ordering; we must convert this to a global ID.
      int local_ID = orders[i].contig_ID(j);
      int global_ID = cluster[ local_ID ];

      if ( !is_long[global_ID] ) continue;

      total_len += contig_lengths[global_ID];

      // Add to the lookup tables.
      old_to_new[ global_ID ] = new_ID;
      new_to_old[ new_ID ] = global_ID;
      new_ID++;

      //cout << "Cluster #" << i << ", contig #" << j << " = " << local_ID << "\t" << global_ID << endl;

    }

    out2 << ( double(new_ID) - 0.5 ) << endl;

    /*
    // Find all the contigs in this cluster that haven't been added to the ordering, and add them in at the end.
    for ( size_t j = 0; j < cluster.size(); j++ )
      if ( old_to_new[ cluster[j] ] == -1 ) {
	old_to_new[ cluster[j] ] = new_ID;
	new_to_old[ new_ID ] = cluster[j];
	new_ID++;
      }
    */

  }
  out2.close();

  int N_contigs_used = new_ID;
  cout << "Total number of contigs that are sufficiently long (among the top " << PLOT_N << ") and are clustered and ordered, and will therefore go into the histogram = " << N_contigs_used << " (length = " << total_len << ")" << endl;


  // 3. Load a GenomeLinkMatrix to get the quantity of Hi-C links between all pairs of contigs.
  GenomeLinkMatrix glm( run_params._out_dir + "/cached_data/all.GLM" );

  // Normalize this GenomeLinkMatrix.
  glm.NormalizeToDeNovoContigLengths( USE_RES );



  // 4. Make a heatmap of the GLM's data, with the contigs reordered as in Step 2.
  cout << "Writing a heatmap to heatmap.txt" << endl;
  ofstream out( "heatmap.txt" );
  out << "X\tY\tZ\n";

  for ( int i = 0; i < N_contigs_used; i++ )
    for ( int j = 0; j < N_contigs_used; j++ ) {
      if ( i == j ) continue; // leave the main diagonal blank - there's no data there anyway

      int i_old = new_to_old[i];
      int j_old = new_to_old[j];

      double data = glm.NLinks( i_old, j_old );
      //if ( data == 0 ) continue; // don't print a lack of data, including along the main diagonal

      // Rescale the data to put it into a more intuitive range for visual clarity.
      data /= 1000;

      //PRINT5( i, j, i_old, j_old, data );
      out << i << '\t' << j << '\t' << data << endl;
    }

  out.close();

  cout << "Making out/heatmap.MWAH.jpg... this may take a while" << endl;
  system( "heatmap.MWAH.R" );

  cout << "Done with that heatmap!" << endl;
}
