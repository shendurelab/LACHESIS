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
 * Lachesis.cc
 *
 * This is the top-level module for the Lachesis program.  This class reads in an INI file and executes the Lachesis algorithms at a high level.
 * From here, the program accesses the other modules:
 *
 * RunParams: The parameters for a Lachesis run, as read from the INI file.
 * GenomeLinkMatrix: A matrix of Hi-C links between all contigs.  Contains contig clustering algorithms.
 * ChromLinkMatrix: A matrix of Hi-C links between contigs in a group.  Contains contig ordering, orienting, and spacing algorithms.
 * LinkSizeDistribution: The density of Hi-C links as a function of distance.  Used for contig spacing.
 * ClusterVec: A vector< set<int> > describing a clustering result.
 * ContigOrdering: An ordering of contigs in a group, eventually including orientation and spacing.
 * TrueMapping: The true location of each contig on the reference assembly, if there is one.  Used for reference-based validation.
 * Reporter: Tools to evaluate the Lachesis result and produce the REPORT.txt file.
 * TextFileParsers: A set of useful functions to parse text files.
 *
 *
 *
 *
 * Josh Burton
 * December 2012
 *
 *************************************************************************************************************************************************************/



// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// STL declarations
#include <ctime>
#include <cerrno>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Local includes
#include "RunParams.h"
#include "GenomeLinkMatrix.h"
#include "ChromLinkMatrix.h"
#include "LinkSizeDistribution.h"
#include "ClusterVec.h"
#include "ContigOrdering.h"
#include "TrueMapping.h"
#include "Reporter.h"










// Run the Lachesis clustering algorithm.
void
LachesisClustering( const RunParams & run_params )
{
  cout << "\n\t|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|\n\t|                                 |\n\t|       LACHESIS CLUSTERING       |\n\t|                                 |\n\t|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|\n\n";

  // Make the directories needed for this run, if necessary.
  system ( ( "mkdir -p " + run_params._out_dir + "/cached_data" ).c_str() );
  system ( ( "mkdir -p " + run_params._out_dir + "/main_results" ).c_str() );

  // Set up a TrueMapping object.
  TrueMapping * true_mapping = run_params.LoadTrueMapping();
  const bool postfosmid = false; // placeholder for parameters optimized for the post-fosmid case

  GenomeLinkMatrix * glm;

  // Look for the *.GLM file, which describes the data in a GenomeLinkMatrix.
  // If the OVERWRITE_GLM flag is not set, and if the file exists (because of a previous run), read the data from it to make a GenomeLinkMatrix object.
  // Otherwise, create the data by reading the SAM files, which takes longer.
  string GLM_file = run_params._out_dir + "/cached_data/all.GLM";
  if ( !boost::filesystem::is_regular_file( GLM_file ) || run_params._overwrite_GLM ) {
    glm = new GenomeLinkMatrix( run_params._species, run_params._SAM_files, run_params.DraftContigRESitesFilename() );
    glm->WriteFile( GLM_file );
  }
  else
    glm = new GenomeLinkMatrix( GLM_file );

  // Pre-processing.
  glm->NormalizeToDeNovoContigLengths( true );
  if ( true_mapping && run_params._sim_bin_size == 0 ) glm->ReorderContigsByRef( *true_mapping );

  glm->SkipContigsWithFewREs( run_params._cluster_min_RE_sites );
  glm->SkipRepeats( postfosmid ? 1.2 : run_params._cluster_max_link_density );

  if ( run_params._cluster_draw_heatmap ) glm->DrawHeatmap( "heatmap.jpg" );


  glm->AHClustering( run_params._cluster_N, run_params._cluster_CEN_contig_IDs, 0, run_params._cluster_noninformative_ratio, run_params._cluster_draw_dotplot, true_mapping );

  // Improve the clustering results, in the postfosmid case.
  if ( postfosmid ) glm->MoveContigsInClusters( 1.2 );
  //glm->UndoMisjoins();

  // If only using high-quality (i.e., well-aligning to reference) contigs, throw out the low-quality contigs at the last minute.
  //glm->ExcludeLowQualityContigs( true_mapping );

  // Report on the clustering.  If there is a TrueMapping, perform reference-based validation.
  glm->ValidateClusters( true_mapping, run_params._cluster_draw_dotplot );

  ClusterVec clusters = glm->GetClusters();
  clusters.WriteFile( run_params._out_dir + "/main_results/clusters.txt" );
  clusters.WriteFile( run_params._out_dir + "/main_results/clusters.by_name.txt", run_params.LoadDraftContigNames() );


  if ( true_mapping ) delete true_mapping; // cleanup
}




// Run the Lachesis ordering and orienting algorithms.
void
LachesisOrdering( const RunParams & run_params )
{
  cout << "\n\t|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|\n\t|                                 |\n\t|        LACHESIS ORDERING        |\n\t|                                 |\n\t|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|\n\n";

  // Load the clusters of contigs.
  string clusters_file = run_params._out_dir + "/main_results/clusters.by_name.txt";
  //string clusters_file = run_params._out_dir + "/main_results/clusters.merged_with_RAD.by_name.txt.split"; // TEMP
  if ( !boost::filesystem::is_regular_file( clusters_file ) ) {
    cerr << "ERROR: Can't find file '" << clusters_file << "' with clusters.  Maybe you're trying to run ordering without having run clustering first? (DO_CLUSTERING = 0, DO_ORDERING = 1)" << endl;
    exit(1);
  }
  const ClusterVec clusters( clusters_file, run_params.LoadDraftContigNames() );
  assert( (int) clusters.size() >= run_params._cluster_N );


  system ( ( "mkdir -p " + run_params._out_dir + "/cached_data" ).c_str() ); // make the directory, if necessary

  // Look for the complete set of ChromLinkMatrix files (*.CLM).  If the CLM files don't all already exist (or if the OVERWRITE_CLMS flag is set), create them.
  // This requires loading the SAM files, which is time-consuming, so we only do it if we have to.
  // But creating the set of CLMs all at once only requires reading through the SAM files once, so it's much faster than creating them all individually.
  for ( size_t i = 0; i < clusters.size(); i++ ) {
    string CLM_file = run_params._out_dir + "/cached_data/group" + boost::lexical_cast<string>( i ) + ".CLM";

    cout << "01 HERE.\n";
    if ( !boost::filesystem::is_regular_file( CLM_file ) || run_params._overwrite_CLMs ) {
      cout << "Need to read SAM files and create ChromLinkMatrix files at " << run_params._out_dir << "/cached_data/group*.CLM.  This will take a while." << endl;

      // Initialize a set of ChromLinkMatrix objects, each with the proper number of contigs.
      vector<ChromLinkMatrix *> CLMs( clusters.size() );
      for ( size_t j = 0; j < clusters.size(); j++ )
	CLMs[j] = new ChromLinkMatrix( run_params._species, clusters[j].size() );

      // Read all of the SAM files and fill all of the ChromLinkMatrices.
      LoadDeNovoCLMsFromSAM( run_params._SAM_files, run_params.DraftContigRESitesFilename(), clusters, CLMs );

      // Write the ChromLinkMatrices to files.
      for ( size_t j = 0; j < clusters.size(); j++ ) {
	string new_file_head = run_params._out_dir + "/cached_data/group" + boost::lexical_cast<string>( j );
	string new_CLM_file = new_file_head + ".CLM";
	CLMs[j]->WriteFile( new_CLM_file );
	delete CLMs[j];
      }

      break;
    }
  }



  // Loop over all clusters.  For each cluster, load a ChromLinkMatrix object and use it to order and orient the contigs.
  for ( size_t i = 0; i < clusters.size(); i++ ) {
    cout << ": Ordering on cluster #" << i << endl;

    // Read in the ChromLinkMatrix from the *.CLM file.  This file should have
    // been created above if it didn't already exist.
    string i_str = boost::lexical_cast<string>(i);
    string clm_input = run_params._out_dir + "/cached_data/group" + i_str + ".CLM";
    cout << "TESTME: " + clm_input + "\n";
    ChromLinkMatrix clm(clm_input);

    //clm.PrefilterLinks( clusters[i], run_params.LoadTrueMapping() );

    // Main algorithms to find the orderings in this chromosome: first the
    // 'trunk' ordering, then the full ordering.  Each of the orderings is also oriented.
    //clm.DrawHeatmap( "heatmap." + boost::lexical_cast<string>(i) + ".jpg" );
    ContigOrdering trunk = clm.MakeTrunkOrder(run_params._order_min_N_REs_in_trunk);
    ContigOrdering order = clm.MakeFullOrder (run_params._order_min_N_REs_in_shreds);
    string trunk_file = run_params._out_dir + "/cached_data/group"  + i_str + ".trunk.ordering";
    trunk.WriteFile( trunk_file, clusters[i], run_params.LoadDraftContigNames());
    string ordering_file = run_params._out_dir + "/main_results/group" + i_str + ".ordering";
    order.WriteFile(ordering_file, clusters[i], run_params.LoadDraftContigNames());
    if (run_params._use_ref && run_params._order_draw_dotplots) {
      string dotplot_file = "clm." + i_str + ".dotplot.txt";
      order.DrawDotplotVsTruth(clusters[i], *(run_params.LoadTrueMapping()), dotplot_file);
    }
  }


}






// Run the Lachesis reporting functions, centered around the Reporter class.
void LachesisReporting(const RunParams &run_params) {
  // Load the clusters of contigs.
  string clusters_file = run_params._out_dir + "/main_results/clusters.by_name.txt";
  //string clusters_file = run_params._out_dir + "/main_results/clusters.merged_with_RAD.by_name.txt.split"; // TEMP
  if ( !boost::filesystem::is_regular_file( clusters_file ) ) {
    cerr << "ERROR: Can't find file '" << clusters_file << "' with clusters.  Maybe you're trying to run reporting without having run clustering first? (DO_CLUSTERING = 0, DO_REPORTING = 1)" << endl;
    exit(1);
  }
  ClusterVec clusters( clusters_file, run_params.LoadDraftContigNames() );


  vector<ContigOrdering> trunks, orders;

  // Load from file the ContigOrdering objects describing the trunks and full orderings (including orientations.)
  for ( size_t i = 0; i < clusters.size(); i++ ) {
    string trunk_file = run_params._out_dir + "/cached_data/group"  + boost::lexical_cast<string>(i) + ".trunk.ordering";
    string order_file = run_params._out_dir + "/main_results/group" + boost::lexical_cast<string>(i) + ".ordering";
    if ( !boost::filesystem::is_regular_file( order_file ) ) {
      cerr << "WARNING: Can't find file '" << order_file << "' with contig orderings.  Maybe you're trying to run reporting without having run ordering first? (DO_ORDERING = 0, DO_REPORTING = 1)  Reporting can still proceed but will report on the clustering only." << endl;
      trunks.clear();
      orders.clear();
      break;
    }

    trunks.push_back( ContigOrdering( trunk_file ) );
    orders.push_back( ContigOrdering( order_file ) );
  }

  // If any groups are marked for exclusion, blank out their orderings.
  if ( !trunks.empty() && !orders.empty() )
    for ( size_t i = 0; i < run_params._report_excluded_groups.size(); i++ ) {
      int ID = run_params._report_excluded_groups[i];
      assert( ID >= 0 );
      if ( ID >= (int) clusters.size() ) cerr << "ERROR: REPORT_EXCLUDED_GROUPS contains group number " << ID << ", which is higher than the total number of groups here (" << clusters.size() << ")" << endl;
      assert( ID < (int) clusters.size() );
      trunks[ID].Clear();
      orders[ID].Clear();
    }



  // Set up the Reporter object with all of these data structures.
  Reporter reporter( run_params, clusters, trunks, orders );
  reporter.Eval();

  // Write a report file.
  reporter.ReportChart();

  // Optional: Make a heatmap of the whole assembly.  This is a usual reference-free visual evaluation.
  if ( run_params._report_draw_heatmap ) MakeWholeAssemblyHeatmap( run_params );
}









int main(int argc, char * argv[]) {
  /* This is stupid
     system ( "cat splash_screen.txt" );
  */
  const char* ini_file;
  cout << endl << endl;

  // If an INI file was not specified, print syntax and exit.
  if (argc != 2) {
    cout << "Syntax: Lachesis <ini_file>" << endl;
    cout << "For a sample ini_file, see Lachesis.ini." << endl << endl;
    cout << "Defaulting to test_case.ini\n";
    ini_file = "INIs/test_case.ini";
  } else {
    ini_file = argv[1];
  }

  // Input the Lachesis.ini file and find run parameters.
  const RunParams run_params(ini_file);

  // Run the steps of the Lachesis ordering!

  if ( run_params._do_clustering ) LachesisClustering( run_params );
  if ( run_params._do_ordering )   LachesisOrdering  ( run_params );
  if ( run_params._do_reporting )  LachesisReporting ( run_params );

  cout << ": Done!" << endl;
  return 0;
}
