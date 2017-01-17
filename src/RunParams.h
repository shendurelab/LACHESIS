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
 * RunParams.h
 *
 * A RunParams object contains the run parameters for a Lachesis run.  That is, it tells Lachesis where to look for input files (draft assembly, aligned Hi-C
 * reads), what algorithmic steps to do, and what heuristic parameters to use.
 *
 *
 *
 * Josh Burton
 * June 2013
 *
 *************************************************************************************************************************************************************/


#ifndef _RUN_PARAMS__H
#define _RUN_PARAMS__H

#include "TrueMapping.h"

#include <vector>
#include <string>
using namespace std;



class RunParams
{
 public:

  RunParams( const string & ini_file ) { ParseIniFile( ini_file ); }

  void ParseIniFile( const string & ini_file );

  // Get the set of contig/chromosome names in the reference assembly fasta.  This fails if USE_REFERENCE = 0.
  // After the first call, this vector is cached for faster retrieval in the future.
  vector<string> * LoadRefGenomeContigNames() const;

  // Get the set of contig names in the draft assembly fasta.  After the first call, this vector is cached for faster retrieval in the future.
  vector<string> * LoadDraftContigNames() const;

  // Get the filename that contains the number of restriction enzyme sites for each contig in the draft assembly.  This might require generating the file,
  // by calling the script CountMotifsInFasta.pl.
  string DraftContigRESitesFilename() const;

  // Load a TrueMapping using the files in this RunParams object.
  // If _use_ref == false, returns a NULL pointer; otherwise returns a new object that must be 'delete'd later to save memory.
  TrueMapping * LoadTrueMapping() const;

  // Report the values of each parameter in this RunParams object (as they appeared in the ini file.)
  void PrintParams( ostream & out = cout ) const;


 public:
  /* These variables are all loaded in by ParseIniFile(). */

  string _species; // effect: certain Lachesis functions are only available if _species == "human"; if _species == "fly", the reference genome is modified
  string _out_dir; // output directory


  // Draft assembly files.
  string _draft_assembly_fasta; // FASTA file of draft assembly
  string _SAM_dir; // directory containing SAM/BAM files
  vector<string> _SAM_files; // the set of SAM/BAM files
  string _RE_site_seq; // used in LoadDraftContigRESites(), which passes it to CountMotifsInFasta.pl

  // Reference assembly files (optional).
  bool _use_ref;
  int _sim_bin_size;
  string _ref_assembly_fasta, _BLAST_file_head;

  // Options for what steps of Lachesis to run.
  bool _do_clustering, _do_ordering, _do_reporting;
  bool _overwrite_GLM, _overwrite_CLMs;

  // Heuristic parameters for clustering.
  int _cluster_N, _cluster_min_RE_sites;
  vector<int> _cluster_CEN_contig_IDs;
  double _cluster_max_link_density, _cluster_noninformative_ratio;
  bool _cluster_draw_heatmap, _cluster_draw_dotplot;

  // Heuristic parameters for ordering.
  int _order_min_N_REs_in_trunk, _order_min_N_REs_in_shreds;
  bool _order_draw_dotplots;

  // Heuristic parameters for reporting.
  vector<int> _report_excluded_groups; // groups chosen not to be included in reporting numbers (e.g., small, chimeric groups)
  int _report_quality_filter;
  bool _report_draw_heatmap;

 private:
  // A listing of all of the lines from the ini file that were used in the creation of this RunParams object.
  vector<string> _params;


  /* These functions and variables are used during ParseIniFile(). */


  // ReportParseFailure: A helper function that gives useful verbose output whenever ParseIniFile fails to parse a line.
  void ReportParseFailure( const string & err_description ) const;

  // ConvertOrFail: Recast this string as object of type T (T=int, double, bool).  If it can't be converted, throw a ReportParseFailure.
  template<class T> T ConvertOrFail( const string & token ) const;

  // VerifySAMFileHeaders: Some sanity checks on the _SAM_files.
  void VerifySAMFileHeaders() const;

  // Variables used in ParseIniFile and its helper functions.
  string _ini_file;
  static const size_t _LINE_LEN = 500000;
  char _line[_LINE_LEN];
  int _line_N;


  // Cached stuff.  This stuff starts out empty.
  mutable vector<string> _ref_contig_names;   // filled by GetRefGenomeContigNames()
  mutable vector<string> _draft_contig_names; // filled by GetDraftContigNames()

};



#endif
