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


// For documentation, see RunParams.h
#include "RunParams.h"
#include "TextFileParsers.h" // ParseTabDelimFile, GetFastaNames

#include <assert.h>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iostream>


// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "gtools/SAMStepper.h" // NTargetsInSAM()

// Local includes
#include "TimeMem.h"











// ReportParseFailure: A helper function that gives useful verbose output whenever ParseIniFile fails to parse a line.
void
RunParams::ReportParseFailure( const string & err_description ) const
{
  cerr << endl
       << "ERROR: Parsing failure in Lachesis INI file!" << endl
       << "FILE: " << _ini_file << endl
       << "LINE " << _line_N << ": '" << _line << "'" << endl
       << err_description << endl << endl;
  exit(1);
}



// ConvertOrFail: Recast this string as object of type T (T=int, double, bool).  If it can't be converted, throw a ReportParseFailure.
// Note that for T=bool, the only strings that are guaranteed to recast correctly are "0" and "1".
template<class T> T
RunParams::ConvertOrFail( const string & value ) const
{
  try {
    return boost::lexical_cast<T> (value);
  }

  catch ( boost::bad_lexical_cast& ) {

    // This function should only be called with T = int, double, or bool.
    string type_name;
    if      ( typeid(T) == typeid(int)    ) type_name = "int";
    else if ( typeid(T) == typeid(double) ) type_name = "double";
    else if ( typeid(T) == typeid(bool)   ) type_name = "bool";
    else ReportParseFailure ( "RunParams::ConvertOrFail called with unrecognized typeid?!?  This is a programming error, not a runtime error." );

    ReportParseFailure( "Value '" + value + "' should be of type '" + type_name + "' but isn't." );
    exit(1); // control never reaches here
  }
}



// VerifySAMFileHeaders: Some sanity checks on the _SAM_files.
void
RunParams::VerifySAMFileHeaders() const
{
  int N_targets;

  // Require the SAM file headers to agree on the number of targets they have.
  for ( size_t i = 0; i < _SAM_files.size(); i++ ) {
    if ( i == 0 ) N_targets = NTargetsInSAM( _SAM_files[i] );
    else assert( N_targets == NTargetsInSAM( _SAM_files[i] ) );
  }
}



// ParseIniFile: Read a Lachesis.ini-type file and fill the variables in this RunParams object.
// Note that the input file must be carefully formatted, with each (non-commented) line of the form "key = value" or "key = value1 value2 [...]".
// If any line in the input file is not of this form, ParseIniFile will choke and complain.
void
RunParams::ParseIniFile( const string & ini_file )
{
  cout << "RunParams is parsing INI file at " << ini_file << endl;
  if ( !boost::filesystem::is_regular_file( ini_file ) ) {
    cerr << "ERROR: Can't find file " << ini_file << endl;
    exit(1);
  }

  const bool verbose = false;


  /* QUALITY CONTROL VARIABLES.  The following set of variables are designed to make sure the INI file is properly formatted as we read it in. */

  // This is the order in which RunParams expects to see the parameters.  It's also the order in which the parameters appear in the default INI files.
  // It's important to enforce the order because some parameters depend on earlier ones (e.g., SAM_DIR must be loaded so we know where to look for SAM_FILES.)
  // If you want to permanently add or remove any parameters, make sure to update N_keys as well as keys_order_array.
  const int N_keys = 28;
  const char * keys_order_array[] = { "SPECIES", "OUTPUT_DIR",
				      "DRAFT_ASSEMBLY_FASTA", "SAM_DIR", "SAM_FILES", "RE_SITE_SEQ",
				      "USE_REFERENCE", "SIM_BIN_SIZE", "REF_ASSEMBLY_FASTA", "BLAST_FILE_HEAD",
				      "DO_CLUSTERING", "DO_ORDERING", "DO_REPORTING", "OVERWRITE_GLM", "OVERWRITE_CLMS",
				      "CLUSTER_N", "CLUSTER_CONTIGS_WITH_CENS", "CLUSTER_MIN_RE_SITES", "CLUSTER_MAX_LINK_DENSITY",
				      "CLUSTER_NONINFORMATIVE_RATIO", "CLUSTER_DRAW_HEATMAP", "CLUSTER_DRAW_DOTPLOT",
				      "ORDER_MIN_N_RES_IN_TRUNK", "ORDER_MIN_N_RES_IN_SHREDS", "ORDER_DRAW_DOTPLOTS",
				      "REPORT_EXCLUDED_GROUPS", "REPORT_QUALITY_FILTER", "REPORT_DRAW_HEATMAP" };
  const vector<string> keys_order( keys_order_array, keys_order_array + N_keys );

  // For certain keys, we can have any (nonzero) number of values appear after the key.  Mark these keys.  For all other keys, exactly one value is required.
  map<string, bool> variable_N_values;
  for ( int i = 0; i < N_keys; i++ ) {
    string key = keys_order[i];
    if ( key == "SAM_FILES" || key == "CLUSTER_CONTIGS_WITH_CENS" || key == "REPORT_EXCLUDED_GROUPS" ) variable_N_values[key] = true;
    else variable_N_values[key] = false;
  }




  vector<string> tokens;
  int current_key_ID = 0;

  _ini_file = ini_file;
  ifstream in( _ini_file.c_str(), ios::in );
  _line_N = 0;

  // Read the file line-by-line.
  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;
    _line_N++;


    // Remove trailing (but not leading) whitespace from the line.
    size_t stop = string(_line).find_last_not_of( " \t" );
    if ( stop == string::npos ) continue; // no non-whitespace in this line, so skip it
    _line[stop+1] = '\0';

    // Skip commented lines.  Note that a line is only considered commented if it starts with a '#' symbol.  Mid-line comments in the ini file aren't allowed.
    if ( _line[0] == '#' ) continue;

    // Split each line into tokens, using whitespace (i.e., each set of one or more consecutive spaces/tabs) as delimiters.
    boost::split( tokens, _line, boost::is_any_of(" \t"), boost::token_compress_on );

    // Require each line to be of the form "key = value" (with room for more than one value.)
    if ( tokens.size() < 3 ) ReportParseFailure( "There must be at least three whitespace-delimited tokens (`e.g., key = value`)." );
    if ( tokens[1] != "=" )  ReportParseFailure( "The second whitespace-delimited token must be an equals sign (=)." );

    string key = tokens[0];
    string value = tokens[2]; // this may be the first of multiple values

    // Require the key to be one of the keys that RunParams expects to see.
    if ( variable_N_values.find( key ) == variable_N_values.end() )
      ReportParseFailure( "Key '" + key + "' isn't recognized as an informative piece of information for Lachesis." );

    // If this key expects only one value (most do), check that only one value is given.
    if ( !variable_N_values.at(key) && tokens.size() != 3 )
      ReportParseFailure( "Key '" + key + "' should be followed by only one value, but multiple values are given." );

    // Require the key to be the next key in the ordered list.
    if ( key != keys_order[current_key_ID] )
      ReportParseFailure( "Key '" + key + "' appears too early in the INI file; we expect to see the key '" + keys_order[current_key_ID] + "' instead.\nMake sure not to change the order in which the keys appear in the INI file." );
    current_key_ID++;



    // Finally, set the local variable that corresponds to this key.
    // For values that are numbers, check to make sure the value is actually a number.
    // For values that describe input files or directories, check to make sure the input file/directory actually exists.
    if ( key == "SPECIES" ) {
      // During software development, Lachesis used to accept only a small set of species names as input.  There's no reason to enforce this anymore.
      // However, the name "human" is still required in some places for non-de-novo CLMs, and "fly" still causes assumptions about which chromosomes to add or
      // drop (see LoadTrueMapping() below.)
      /*
      vector<string> allowable_species;
      allowable_species.push_back("human");
      allowable_species.push_back("mouse");
      allowable_species.push_back("fly");
      if ( find( allowable_species.begin(), allowable_species.end(), value ) == allowable_species.end() )
	ReportParseFailure( "'" + value + "' is not recognized by Lachesis.  Accepted species are: " + boost::algorithm::join( allowable_species, ", " ) );
      */
      _species = value;
    }
    else if ( key == "OUTPUT_DIR" ) {
      _out_dir = value;
      if ( boost::filesystem::is_regular_file( value ) )
	ReportParseFailure( "File '" + value + "' already exists as a non-directory file." );
    }
    else if ( key == "DRAFT_ASSEMBLY_FASTA" ) {
      _draft_assembly_fasta = value;
      if ( !boost::filesystem::is_regular_file( value ) && !boost::filesystem::is_regular_file( value + ".names" ) ) // technically this should also check for the existence of <assembly-fasta>.counts_<RE_SITE_SEQ>.txt, but RE_SITE_SEQ hasn't been loaded in yet
	ReportParseFailure( "Can't find file '" + value + "'." );
    }
    else if ( key == "SAM_DIR" ) {
      _SAM_dir = value;
      if ( !boost::filesystem::is_directory( value ) )
	ReportParseFailure( "Can't find directory '" + value + "'." );
    }
    else if ( key == "SAM_FILES" ) {
      for ( size_t i = 2; i < tokens.size(); i++ ) {
	_SAM_files.push_back( _SAM_dir + "/" + tokens[i] );
	if ( !boost::filesystem::is_regular_file( _SAM_dir + "/" + tokens[i] ) )
	  ReportParseFailure( "Can't find file '" + tokens[i] + "' in directory SAM_DIR = '" + _SAM_dir + "'." );
      }
      VerifySAMFileHeaders();
    }
    else if ( key == "RE_SITE_SEQ" ) _RE_site_seq = value;
    else if ( key == "USE_REFERENCE" ) _use_ref = ConvertOrFail<bool>( value );
    else if ( key == "SIM_BIN_SIZE" ) {
      _sim_bin_size = ConvertOrFail<int>( value );
      if ( _use_ref )
	if ( _sim_bin_size != 0 && _species != "human" )
	  ReportParseFailure( "SIM_BIN_SIZE can only be nonzero if species = 'human'." );
    }
    else if ( key == "REF_ASSEMBLY_FASTA" ) {
      _ref_assembly_fasta = value;
      if ( _use_ref )
	if ( !boost::filesystem::is_regular_file( value ) && !boost::filesystem::is_regular_file( value + ".names" ) )
	  ReportParseFailure( "Can't find file '" + value + "'." );
    }
    else if ( key == "BLAST_FILE_HEAD" ) {
      _BLAST_file_head = value;
      if ( _use_ref && _sim_bin_size == 0 && !boost::filesystem::is_regular_file( value + ".1.blast.out" ) && !boost::filesystem::is_regular_file( _out_dir + "/cached_data/TrueMapping.assembly.txt" ) )
	ReportParseFailure( "BLAST_FILE_HEAD = " + value + " doesn't seem to point to any usable BLAST alignment files, and no cached data is available in OUTPUT_DIR." );
    }
    else if ( key == "DO_CLUSTERING" )  _do_clustering  = ConvertOrFail<bool>( value );
    else if ( key == "DO_ORDERING" )    _do_ordering    = ConvertOrFail<bool>( value );
    else if ( key == "DO_REPORTING" )   _do_reporting   = ConvertOrFail<bool>( value );
    else if ( key == "OVERWRITE_GLM" )  _overwrite_GLM  = ConvertOrFail<bool>( value );
    else if ( key == "OVERWRITE_CLMS" ) _overwrite_CLMs = ConvertOrFail<bool> ( value );
    else if ( key == "CLUSTER_N" )                    _cluster_N                    = ConvertOrFail<int>   ( value );
    else if ( key == "CLUSTER_CONTIGS_WITH_CENS" ) {
      _cluster_CEN_contig_IDs.clear();
      if ( value != "-1" ) // if the first listed value is -1, don't do anything - leave the _cluster_CEN_contig_IDs vector empty
      for ( size_t i = 2; i < tokens.size(); i++ ) {
	int cID = ConvertOrFail<int>( tokens[i] );
	if ( cID < 0 ) ReportParseFailure( "CLUSTER_CONTIGS_WITH_CENS can't contain any contig IDs less than 0 (unless the list consists entirely of `-1', indicating an empty list.)" );
	assert( cID >= 0 );
	_cluster_CEN_contig_IDs.push_back( cID );
      }
      if ( (int) _cluster_CEN_contig_IDs.size() > _cluster_N ) ReportParseFailure( "CLUSTER_CONTIGS_WITH_CENS = contains more contig IDs than the number of clusters CLUSTER_N." );
    }
    else if ( key == "CLUSTER_MIN_RE_SITES" )         _cluster_min_RE_sites         = ConvertOrFail<int>   ( value );
    else if ( key == "CLUSTER_MAX_LINK_DENSITY" )     _cluster_max_link_density     = ConvertOrFail<double>( value );
    else if ( key == "CLUSTER_NONINFORMATIVE_RATIO" ) {
      _cluster_noninformative_ratio = ConvertOrFail<double>( value );
      if ( _cluster_noninformative_ratio != 0 && _cluster_noninformative_ratio <= 1 )
	ReportParseFailure( "CLUSTER_NONINFORMATIVE_RATIO must either be 0 or >1." );
    }
    else if ( key == "CLUSTER_DRAW_HEATMAP" )         _cluster_draw_heatmap         = ConvertOrFail<bool>  ( value );
    else if ( key == "CLUSTER_DRAW_DOTPLOT" )         _cluster_draw_dotplot         = ConvertOrFail<bool>  ( value );
    else if ( key == "ORDER_MIN_N_RES_IN_TRUNK" )     _order_min_N_REs_in_trunk     = ConvertOrFail<int>   ( value );
    else if ( key == "ORDER_MIN_N_RES_IN_SHREDS" )    _order_min_N_REs_in_shreds    = ConvertOrFail<int>   ( value );
    else if ( key == "ORDER_DRAW_DOTPLOTS" )          _order_draw_dotplots          = ConvertOrFail<bool>  ( value );
    else if ( key == "REPORT_EXCLUDED_GROUPS" ) {
      _report_excluded_groups.clear();
      if ( value != "-1" ) // if the first listed value is -1, don't do anything - leave the _report_excluded_groups vector empty
	for ( size_t i = 2; i < tokens.size(); i++ ) {
	  int group_ID = ConvertOrFail<int>( tokens[i] );
	  //if ( group_ID < 0 || group_ID >= _cluster_N ) ReportParseFailure( "Group ID '" + tokens[i] + "' is not in the range [0,N_clusters)." );
	  _report_excluded_groups.push_back( group_ID );
	}
    }
    else if ( key == "REPORT_QUALITY_FILTER" ) _report_quality_filter = ConvertOrFail<int>( value );
    else if ( key == "REPORT_DRAW_HEATMAP" )   _report_draw_heatmap   = ConvertOrFail<bool>( value );


    // Record this line.
    _params.push_back( _line );


    // Verbose output (optional).
    if ( !verbose ) continue;
    cout << "LINE #" << _line_N << endl;
    cout << "KEY: " << key << endl;
    cout << "VALUES:";
    for ( size_t i = 2; i < tokens.size(); i++ )
      cout << '\t' << tokens[i];
    cout << endl << endl;
  }


  assert( current_key_ID == N_keys );
}




// Get the set of contig/chromosome names in the reference assembly fasta.  This fails if USE_REFERENCE = 0.
// After the first call, this vector is cached for faster retrieval in the future.
vector<string> *
RunParams::LoadRefGenomeContigNames() const
{
  // If this function hasn't already been run, parse the fasta file and get the reference assembly's contig/chromosome names.
  // Note that this reads the file <_ref_assembly_fasta>.names, and will create the file if necessary.
  if ( _ref_contig_names.empty() ) _ref_contig_names = GetFastaNames( _ref_assembly_fasta );
  return &_ref_contig_names;
}




// Get the set of contig names in the draft assembly.  After the first call, this vector is cached for faster retrieval in the future.
vector<string> *
RunParams::LoadDraftContigNames() const
{
  // If this function hasn't already been run, parse the fasta file and get the draft contig names.
  // Note that this reads the file <_draft_assembly_fasta>.names, and will create the file if necessary.
  if ( _draft_contig_names.empty() ) _draft_contig_names = GetFastaNames( _draft_assembly_fasta );
  return &_draft_contig_names;
}



// Get the filename that contains the number of restriction enzyme sites for each contig in the draft assembly.  This might require generating the file,
// by calling the script CountMotifsInFasta.pl.
string
RunParams::DraftContigRESitesFilename() const
{
  string RE_sites_file = _draft_assembly_fasta + ".counts_" + _RE_site_seq + ".txt";

  // If this function hasn't already been run, the file may not exist, in which case the script CountMotifsInFasta.pl needs to be run.
  if ( !boost::filesystem::is_regular_file( RE_sites_file ) ) {
    string cmd = "CountMotifsInFasta.pl " + _draft_assembly_fasta + " " + _RE_site_seq;
    system( cmd.c_str() );
    assert( boost::filesystem::is_regular_file( RE_sites_file ) ); // if this fails, the CountMotifsInFasta.pl script didn't run correctly
  }

  return RE_sites_file;
}



// Load a TrueMapping using the files in this RunParams object.
// If _use_ref == false, returns a NULL pointer; otherwise returns a new object that must be 'delete'd later to save memory.
TrueMapping *
RunParams::LoadTrueMapping() const
{
  if ( !_use_ref ) return NULL;

  assert( !_SAM_files.empty() );

  // Call these functions to fill the cache variables _ref_contig_names and _draft_contig_names, respectively.
  LoadRefGenomeContigNames();
  LoadDraftContigNames();

  // Create a TrueMapping object, which records where the contigs are truly located on the reference.
  // If the draft assembly consists of simulated bins from the reference assembly, use a special constructor that deduces the true location of each bin.
  // Otherwise, pass the BLAST file head into the constructor so it can use those alignments.

  TrueMapping * mapping =
    _sim_bin_size != 0 ?
    new TrueMapping( _species, _sim_bin_size, _draft_contig_names, _ref_contig_names )
    :
    new TrueMapping( _species, _draft_contig_names, _ref_contig_names, _BLAST_file_head, _out_dir, _SAM_files[0] );


  // Remove heterochromatic regions from the fly reference.
  if ( _species == "fly" ) {
    mapping->RemoveTarget( "2LHet" );
    mapping->RemoveTarget( "2RHet" );
    mapping->RemoveTarget( "3LHet" );
    mapping->RemoveTarget( "3RHet" );
    mapping->RemoveTarget( "XHet" );
    mapping->RemoveTarget( "YHet" );
    mapping->RemoveTarget( "U" );
    mapping->RemoveTarget( "Uextra" );

    // Also merge chromosomes 2 and 3.
    mapping->MergeTargets( "2L", "2R", "2" );
    mapping->MergeTargets( "3L", "3R", "3" );
  }

  return mapping;
}






// Report the values of each parameter in this RunParams object (as they appeared in the ini file.)
void
RunParams::PrintParams( ostream & out ) const
{
  for ( size_t i = 0; i < _params.size(); i++ )
    out << _params[i] << endl;
}
