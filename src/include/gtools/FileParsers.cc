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


// For documentation, see FileParsers.h
#include "FileParsers.h"

#include <assert.h>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

// Local modules in the gtools library.
#include "ChromInterval.h"
#include "VCF_variant_info.h"
#include "HumanGenome.h"

// Boost libraries.
#include <boost/algorithm/string.hpp> // split, to_upper
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/regex.hpp>

// Modules in ~/include.
#include "../TimeMem.h"


static const unsigned _LINE_LEN = 50000;
char _line[_LINE_LEN];


static const bool EXCLUDE_HELA_LOH_ARMS = false; // TEMP




// GetEnv: Return the shell environment value for a variable.  If the variable doesn't exist, return an empty string.
// This function is a wrapper to the function getenv(const char *).
string
GetEnv( const char * var )
{
  const char * value = getenv(var);
  if ( value == NULL ) return "";
  return string( value );
}







// ParseBED: Parse a BED file and return a vector of chrom_intervals
//     representing the first three columns in the file (chrom, start, stop).
// If chrom != "" (the default value), only return intervals on that chrom.
vector<chrom_interval>
ParseBED( const string & BED_file, const string & chrom )
{
  assert( boost::filesystem::is_regular_file( BED_file ) );
  vector<chrom_interval> intervals;

  // Read the file line-by-line.
  ifstream in( BED_file.c_str(), ios::in );
  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;

    if ( _line[0] == '#' ) continue; // ignore commented lines

    // Create the chrom_interval object by parsing the BED line.
    chrom_interval ci = chrom_interval_from_BED_line( _line );
    if ( chrom != "" && chrom != ci.chrom() ) continue;
    intervals.push_back( ci );
  }

  in.close();


  return intervals;
}



// ParseBEDgraph: Like ParseBED, but return a vector<pair<ci,double> > that
//     indicates the value in the fourth column of the BED/BEDgraph file.
vector< pair<chrom_interval,double> >
ParseBEDgraph( const string & BED_file, const string & chrom )
{
  assert( boost::filesystem::is_regular_file( BED_file ) );
  vector< pair<chrom_interval,double> > intervals;

  // Read the file line-by-line.
  ifstream in( BED_file.c_str(), ios::in );
  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;

    if ( _line[0] == '#' ) continue; // ignore commented lines

    // Create the chrom_interval object by parsing the BED line.
    chrom_interval ci = chrom_interval_from_BED_line( _line );
    if ( chrom != "" && chrom != ci.chrom() ) continue;

    // Break up the line into tab-delimited tokens.
    // The fourth token is the value we want.
    vector<string> tokens;
    boost::split( tokens, _line, boost::is_any_of("\t") );
    assert( tokens.size() >= 4 );
    double value = boost::lexical_cast<double>( tokens[3] );

    intervals.push_back( make_pair( ci, value ) );
  }

  in.close();


  return intervals;
}





// ParseAndMergeBED: Parse a BED/BEDgraph file, merge the input chrom_intervals
//     into equal-sized "windows", and return those windows sorted by chrom.
// The output is a map of chromosome name to the set of chrom_intervals on that
// chromosome.
// The typical usage of this is for SUNK windows.  These windows are typically
// ~2k each, which is too small to be used individually for variant or LOH
// analysis.  So as a pre-processing step, we merge sets of
// <n_intervals_per_window> consecutive "SUNK intervals" into a "SUNK window".
map< string, vector<chrom_interval> >
ParseAndMergeBED( const string & BED_file, const int n_intervals_per_window )
{
  // Output: a set of "windows", sorted by chromosome.
  map< string, vector<chrom_interval> > windows;

  vector<chrom_interval> empty_vec;
  chrom_interval window;
  int n_in_window = 0;


  // Read the bedgraph file line-by-line.
  ifstream in( BED_file.c_str(), ios::in );
  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;

    // Parse this line into a chrom_interval.
    chrom_interval interval = chrom_interval_from_BED_line( _line );

    // At the very first interval only, open a new window.
    if ( n_in_window == 0 ) {
      windows.insert( make_pair( interval.chrom(), empty_vec ) );

      window = interval;
      n_in_window = 1;
    }

    // At every interval past the first, the window may need to be closed,
    // either because it's "full" or because we've changed chromosomes.
    // Check for this.
    else if ( n_in_window == n_intervals_per_window || window.chrom() != interval.chrom() ) {

      //cout << "Adding window: " << window << "\twith " << n_in_window << " intervals" << endl;
      windows[window.chrom()].push_back( window );

      if ( window.chrom() != interval.chrom() ) windows.insert( make_pair( interval.chrom(), empty_vec ) );

      // Open a new window, starting with this interval.
      window = interval;
      n_in_window = 1;
    }

    // If the window is still open, add the current interval to it.
    else {
      assert( window.stop <= interval.start );
      window.stop = interval.stop;
      n_in_window++;
    }

  }

  in.close();

  return windows;
}





// ParseCNFile
// Parse a BEDgraph file describing a genome-wide copy number (CN) profile.
// Return a map of chrom_intervals to CN calls (as opposed to ParseBEDgraph,
// which returns a vector<pair>.)  It has to be a multimap because for reasons
// I don't entirely understand, the STL map<> container doesn't properly
// differentiate between different chrom_intervals.
multimap<chrom_interval, int>
ParseCNFile( const string & CN_calls_file, const string & chrom )
{
  multimap<chrom_interval, int> CN_calls;


  // Open the file and read it line-by-line.
  ifstream in( CN_calls_file.c_str(), ios::in );

  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;

    if ( _line[0] == '#' ) continue; // ignore comments.

    // Create the chrom_interval object by parsing the BED line.
    chrom_interval ci = chrom_interval_from_BED_line( _line );
    if ( chrom != "" && chrom != ci.chrom() ) continue;

    // Break up the line into tab-delimited tokens.
    // The fourth token is the value we want.
    vector<string> tokens;
    boost::split( tokens, _line, boost::is_any_of("\t") );
    assert( tokens.size() >= 4 );
    int value = boost::lexical_cast<int>( tokens[3] );
    CN_calls.insert( make_pair( ci, value ) );
  }

  in.close();

  return CN_calls;
}




// ParseSimHapMatrix: Parse a HapMatrix of simulated clone-call data, as
//     generated by the Python script make_matrix_file.py.
void
ParseSimHapMatrix
( const string & matrix_file, // input file
  int & n_frags, // number of clones in output
  int & n_loci, // number of variant loci
  int & frag_size, // number of variant calls in each fragment
  vector<string> & frag_data, // raw data from each fragment
  vector<int> & frag_offsets, // offset location of each fragment
  vector<boost::tribool> & frag_truth, // tribool::indeterminate = clone error
  string & loci_truth // the true haplotype of all loci
  )
{
  frag_size = -1;
  int frag_ID = 0;

  // Read the file line-by-line.
  ifstream matrix( matrix_file.c_str(), ios::in );
  while ( 1 ) {
    matrix.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( matrix.fail() ) break;

    // Search for the line with "N FRAGMENTS" on it.  There will be one vertex
    // in the HapGraph for each fragment.
    boost::cmatch match;
    if ( boost::regex_search( _line, match, boost::regex("N FRAGMENTS: (\\d+)") ) ) {
      n_frags = boost::lexical_cast<int>( match[1] );
      frag_data   .resize( n_frags );
      frag_offsets.resize( n_frags );
      frag_truth  .resize( n_frags );
    }

    // Search for a line with "N LOCI" on it.  This is used as a sanity check.
    else if ( boost::regex_search( _line, match, boost::regex("N LOCI: (\\d+)") ) )
      n_loci = boost::lexical_cast<int>( match[1] );

    // Search for the line with "FRAGMENT SIZE" on it.  This will be used as a
    // sanity check.
    else if ( boost::regex_search( _line, match, boost::regex("FRAGMENT SIZE: (\\d+)") ) )
      frag_size = boost::lexical_cast<int>( match[1] );

    // Search for the line with "TRUTH HAPLOTYPE" on it.  This is the true
    // haplotype of all loci - used for validation only!
    else if ( boost::regex_search( _line, match, boost::regex("TRUTH HAPLOTYPE: (\\d+)") ) )
      loci_truth = match[1];

    // Search for all lines with data about a fragment.
    else if ( ! boost::regex_search( _line, boost::regex( "\\#" ) ) ) {
      assert( boost::regex_match( _line, match, boost::regex( "(\\d+)\\s+(\\d+)\\s+(\\w+)" ) ) );

      // In the simulated matrix file format, each fragment line contains an
      // integer offset (to keep the matrix sparse), followed by a set of 1's
      // and 0's indicating the value of the variant at each location in the
      // fragment, followed by a character from [0,1,E] indicating the truth
      // identity of this fragment.
      int offset = boost::lexical_cast<int>( match[1] );
      string values = match[2];

      assert( values.find_first_not_of( "01" ) == string::npos );
      assert( (int) values.size() == frag_size );
      assert( offset + frag_size <= n_loci );

      // Make a list of fragment offsets and contents.
      frag_data[frag_ID] = values;
      frag_offsets[frag_ID] = offset;

      // Grab the truth bit (the haplotype this fragment is really from.)
      // An "E" here represents a clone error.  Use a tribool to store this.
      assert( match[3] == "0" || match[3] == "1" || match[3] == "E" );
      if      ( match[3] == "1" ) frag_truth[frag_ID] = true;
      else if ( match[3] == "0" ) frag_truth[frag_ID] = false;
      else                        frag_truth[frag_ID] = boost::indeterminate;

      assert( frag_ID++ < n_frags );
    }
  }



  // Verify that all necessary pieces of data have been found.
  assert( n_frags > 0 );
  assert( n_loci > 0 );
  assert( frag_size > 0 );
  assert( frag_ID == n_frags );
  assert( (int) loci_truth.size() == n_loci );
}





// ParseRealHapMatrix: Parse a HapMatrix of real clone-call data, as generated
//     by the Python script VCFtoHaploMatrix.py or a related script.
// The main output data is the sets var_calls and clone_calls.
// Each element in var_calls is a variant, mapped to a vector of the clones in
// which it appears and the way in which it is called (true=alt, false=ref) in
// each clone.
// clone_calls is a reverse mapping (clone->call instead of variant->call).
void
ParseRealHapMatrix
( const string & matrix_file, // input file
  int & n_clones, // number of clones in output
  int & n_loci, // number of variant loci
  map< string, vector< pair<int, bool> > > & var_calls, // variant -> call
  vector< map<int,string> > & clone_calls, // clone -> set of variant calls
  vector<chrom_interval> & clone_intervals, // chromosomal location of clones
  vector<double> & clone_qscores // quality score of each clone, from the file
  )
{
  assert( boost::filesystem::is_regular_file( matrix_file ) );


  // These variables are used locally, for internal consistency checks.
  set<string> variants;
  int total_n_calls = 0;
  int n_calls = 0;
  int frag_ID = 0;

  vector<chrom_interval> HeLa_LOH_arms = HumanGenome_GetHeLaLOHChromosomeArms();



  // Read the file line-by-line.
  ifstream matrix( matrix_file.c_str(), ios::in );
  while ( 1 ) {
    matrix.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( matrix.fail() ) break;

    // Search for the line with "N CLONES" on it.  For the purposes of this
    // algorithm, a "clone" is a type of fragment, and so there will be one
    // vertex in the HapGraph for each clone.
    boost::cmatch match;
    if ( boost::regex_search( _line, match, boost::regex("N CLONES: (\\d+)") ) ) {
      n_clones = boost::lexical_cast<int>( match[1] );
      clone_calls.resize( n_clones );
      clone_intervals.resize( n_clones );
      clone_qscores  .resize( n_clones );
    }

    // Search for a line with "N VARIANTS" on it.  For the purposes of this
    // algorithm, a "variant" is a locus.
    else if ( boost::regex_search( _line, match, boost::regex("N VARIANTS: (\\d+)") ) )
      n_loci = boost::lexical_cast<int>( match[1] );

    // Search for a line with "TOTAL N CLONE CALLS" on it.  This number will be
    // used as a sanity check.
    else if ( boost::regex_search( _line, match, boost::regex("TOTAL N CLONE CALLS OF VARIANTS: (\\d+)") ) )
      total_n_calls = boost::lexical_cast<int>( match[1] );

    // Search for all lines with data about a fragment.
    // In the real matrix-file format created by VCFtoHaploMatrix.py, each line
    // is a series of tab-delimited tokens.  The first token is a clone's tag.
    // and the subsequent token[s] are all of the form "<var>,<call>", where
    // <var> is a variant tag, and <call> is either 0 or 1, meaning that this
    // clone calls this variant as either ref or non-ref.
    else if ( ! boost::regex_search( _line, boost::regex( "\\#" ) ) ) {

      // Break up the line into tab-delimited tokens.
      vector<string> tokens;
      boost::split( tokens, _line, boost::is_any_of("\t") );
      assert( tokens.size() >= 3 );

      // The first token gives the clone "tag", in the form: chr_start_end
      // Parse this tag to find the clone location on the chromosome.
      vector<string> clone_tokens;
      boost::split( clone_tokens, tokens[0], boost::is_any_of("_") );
      assert( clone_tokens.size() == 3 );
      clone_intervals[frag_ID] =
	chrom_interval( clone_tokens[0],
			boost::lexical_cast<int>( clone_tokens[1] ),
			boost::lexical_cast<int>( clone_tokens[2] ) );
      clone_intervals[frag_ID].ID = frag_ID;

      if ( EXCLUDE_HELA_LOH_ARMS ) {
	bool overlap = false;
	for ( size_t j = 0; j < HeLa_LOH_arms.size(); j++ )
	  if ( HeLa_LOH_arms[j].overlaps( clone_intervals[frag_ID] ) ) {
	    cout << "\tSKIPPING... LOH" << endl;
	    overlap = true;
	    break;
	  }
	if ( overlap ) continue;
      }


      // The second token gives the clone weight.
      try {
	clone_qscores[frag_ID] = boost::lexical_cast<double>( tokens[1] );
      }
      catch( boost::bad_lexical_cast & ) {
	cerr << "WARNING: Using an obsolete HapMatrix file!" << endl;
	cerr << "PROBLEM LINE: " << _line << endl;
      }



      // Each token after the second describes a variant called by this clone.
      for ( size_t i = 2; i < tokens.size(); i++ ) {
	string t = tokens[i];
	assert( t[t.size()-2] == ',' );
	assert( t[t.size()-1] == '0' || t[t.size()-1] == '1' );
	bool call = t[t.size()-1] == '1';
	string var = t.substr( 0, t.size()-2 );
	//cout << "CLONE: " << tokens[0] << ", VARIANT:" << var << ", CALL = " << int(call) << endl;

	// Add this variant and call to the data structures.
	n_calls++;
	var_calls[var].push_back( make_pair( frag_ID, call ) );
	variants.insert( var );

	// Split this variant name to get the identity of the ref/alt alleles.
	// Record, for each clone, the set of variants and their descriptions,
	// as well as which allele appears on this fragment.
	vector<string> var_tokens;
	boost::split( var_tokens, var, boost::is_any_of("_") );
	assert( var_tokens.size() == 4 );
	int var_pos = boost::lexical_cast<int>( var_tokens[1] );
	assert( var_pos >= clone_intervals[frag_ID].start );
	assert( var_pos <  clone_intervals[frag_ID].stop  );
	string ref_alt_call = var_tokens[2] + "_" + var_tokens[3] + "_" + t[t.size()-1];
	clone_calls[frag_ID].insert( make_pair( var_pos, ref_alt_call ) );


      }

      assert( frag_ID++ < n_clones );

    }

  }


  //cout << Time() << ": ParseRealHapMatrix:\tN FRAGS (vertices): " << n_clones << ", N LOCI: " << n_loci << ", TOTAL N CALLS: " << total_n_calls << endl;

  if ( EXCLUDE_HELA_LOH_ARMS ) {
    n_clones = frag_ID;
    n_loci = variants.size();
  }
  else {
    // Sanity checks.
    // These will fail if the header info in the matrix file is incorrect.
    assert( frag_ID == n_clones );
    assert( total_n_calls > 0 );
    assert( n_calls == total_n_calls );
    assert( n_loci == (int)variants .size() );
    assert( n_loci == (int)var_calls.size() );
  }
}





// ParseVCFLine: Helper function for ParseVCF.
// Parse a line in a VCF file.
// If this line represents a good variant, fill the var object and return true.
// Return false for commented lines, variants with undefined reference alleles,
// triallelics, and variants that fail the VCF_input_filter.
bool
ParseVCFLine( const char * line, const VCF_input_filter & filter, VCF_variant_info & var )
{
  if ( line[0] == '#' ) return false; // ignore commented lines

  // TODO: Try a boost::split instead of a regex to speed this up

  // Match this line to a regex of a VCF line.
  // There are two distinct VCF formats, created by GATK and samtools,
  // respectively, but both formats produce the same set of 11 tokens:
  // 0: The complete line
  // 1: Chromosome
  // 2: Position
  // 3: rsID (or '.' if not in dbSNP)
  // 4. reference base
  // 5. alternate base
  // 6. QUAL column (ignored in favor of PL in the format field)
  // 7. Pass filter ('.', 'PASS', or 'LowQual')
  // 8. The INFO field: a set of key-value pairs in the form FOO=bar;BAZ=...
  // 9. The FORMAT field: a set of 2-letter keys in the form "AD:DP:..:etc."
  // 10. The values corresponding to the keys to the FORMAT field.
  // We use one regex for both GATK and samtools formats, and we REQUIRE that
  // every non-commented line of a VCF file matches our regex completely.
  // The only thing different between the two formats is the way in which the
  // read depths are calculated.
  boost::cmatch match;

  if ( !boost::regex_match( line, match, boost::regex("([\\w\\.]+)\\t(\\w+)\\t([rs\\d\\;\\.]+)\\t([ACGTN]+)\\t([ACGTN\\.\\,]+)\\t([\\d\\.]+)\\t(\\.|PASS|LowQual)\\t(\\S+)\\t([\\w\\:]+)\\t(\\S+)") ) ) {

    // No match.  Bail out.
    cerr << "FileParsers::ParseVCFLine: LINE DOES NOT MATCH ANY KNOWN VCF FORMAT:\n" << line << endl;
    exit(1);
  }


  assert( match.size() == 11 );
  assert( match[0] == line );


  // Apply the VCF_input_filter and filter out variants.
  string chrom = match[1];
  if ( filter.chrom != "" && filter.chrom != chrom ) return false;
  bool dbSNP = ( match[3] != "." );
  if ( filter.dbSNP != -1 && filter.dbSNP != dbSNP ) return false;


  // Parse the FORMAT and values fields - i.e., the "AD:GT:DP:GQ:PL" and
  // their values.  Make a map<> for these key-value pairs.
  map<string,string> format;
  {
    string keys0 = match[9], values0 = match[10];
    vector<string> keys, values;
    boost::split( keys,   keys0,   boost::is_any_of(":") );
    boost::split( values, values0, boost::is_any_of(":") );
    assert( keys.size() == values.size() );
    for ( size_t i = 0; i < keys.size(); i++ ) {
      assert( keys[i].size() == 2 );
      format[keys[i]] = values[i];
    }
  }



  // Determine the genotype call.  Filter out variants by genotype if chosen.
  int call = ( format["GT"][0] == '1' ) + ( format["GT"][2] == '1' );
  if ( filter.genotype != -1 && filter.genotype != call ) return false;


  // Find the ref/alt base(s).  They will be single characters unless this is
  // an indel.
  // Skip vars with undefined reference alleles, and also triallelics.
  string ref_base = match[4];
  string alt_base = match[5];
  boost::to_upper( ref_base );
  boost::to_upper( alt_base );
  if ( ref_base.find_first_not_of( "ACGT" ) != string::npos ) return false;
  if ( alt_base.find_first_not_of( "ACGT" ) != string::npos ) return false;


  // Parse out the ref/alt depths.  The exact method for doing this depends
  // on the VCF format, which we distinguish by which FORMAT fields appear.
  int ref_depth, alt_depth;

  // In GATK-produced VCF files, we use the FORMAT field called "AD".
  if ( format.find( "AD" ) != format.end() ) {
    vector<string> depths;
    boost::split( depths, format["AD"], boost::is_any_of(",") );

    // The AD value contains two elements: ref depth and alt depth.
    assert( depths.size() == 2 );
    ref_depth = boost::lexical_cast<int>( depths[0] );
    alt_depth = boost::lexical_cast<int>( depths[1] );
  }
  else {

    // In samtools-produced VCF files, we use the INFO field called "DP4".
    // To save runtime, parse out the elements following "DP4=" using string
    // functions instead of regex functions.
    string INFO = match[8];
    int str_start = INFO.find( "DP4=" ) + 4;
    int str_stop  = INFO.find( ';', str_start );
    string DP4 = INFO.substr( str_start, str_stop - str_start );

    // The DP4 value has 4 elements: ref depth fwd/rc, and alt depth fwd/rc.
    vector<string> dp4;
    boost::split( dp4, DP4, boost::is_any_of(",/") );
    assert( dp4.size() == 4 );
    ref_depth = boost::lexical_cast<int>( dp4[0] ) + boost::lexical_cast<int>( dp4[1] );
    alt_depth = boost::lexical_cast<int>( dp4[2] ) + boost::lexical_cast<int>( dp4[3] );

  }

  // Parse out the log-probabilities of each genotype call.  Higher numbers
  // indicate lower probabilities.
  // For a quality score, we use the differential between the log-probability
  // of the genotype call that was made and the next best log-probability.
  vector<string> logprobs;
  boost::split( logprobs, format["PL"], boost::is_any_of(",") );
  assert( logprobs.size() == 3 );
  int call_logprob = boost::lexical_cast<int>( logprobs[call] );
  int call_bestother = INT_MAX;
  for ( int i = 0; i < 3; i++ ) {
    if ( i == call ) continue;
    int logprob = boost::lexical_cast<int>( logprobs[i] );
    call_bestother = min( call_bestother, logprob );
  }
  int qual_PL = max( call_bestother - call_logprob, 0 );


  // Create a VCF_variant_info object and fill it with the data we've found.
  var = VCF_variant_info();
  var.chrom = chrom;
  if ( var.chrom[0] != 'c' ) var.chrom = "chr" + var.chrom;
  var.pos = boost::lexical_cast<int>( match[2] );
  var.dbSNP = dbSNP;
  var.in1KG = boost::indeterminate; // this can be set later by calling Set1KGFlags()
  var.qual = qual_PL;
  var.ref_base = ref_base;
  var.alt_base = alt_base;
  var.ref_depth = ref_depth;
  var.alt_depth = alt_depth;
  var.call = call;

  // Some variables in the VCF_variant_info object aren't known from the VCF
  // but may be figured out later.  For now, set these variables' values to
  // "unknown".
  var.CN = -1; // copy number
  var.in_repeat = boost::indeterminate; // is this var in a repeat?

  return true;
}





/* ParseVCF: Parse one or more VCF files and return a struct containing various
 *     bits of data.  The regex here is currently designed to handle VCFs
 *     created with GATK (including all-positions) and samtools.
 * The VCF_input_filter object allows the user to filter out variants that
 * don't match a set of criteria (e.g., genotype call, which chromosome, etc.)
 *
 * This function has been tested on the following VCF files:
 * /net/shendure/vol7/HELA/110831.HELA_ABCDEF.srt.rmdup.GATK.recal.realigned.vars.raw.vcf (GATK)
 * /net/shendure/vol7/HELA/120110.HELA_ALL.UG.1kG_CEU_YRI.vcf (GATK)
 * /net/shendure/vol10/HELA/VARIANTS/120525.HELA_ALL.SNP.SAMTOOLS_CURRENT.vcf (samtools)
 *
 * The Python equivalent of this code is at VCFParsers.py::parse_VCF_line.
 *
 *****************************************************************************/
vector<VCF_variant_info>
ParseVCF( const vector<string> & VCF_files, const struct VCF_input_filter & filter )
{
  bool verbose = false;

  for ( size_t i = 0; i < VCF_files.size(); i++ )
    assert( boost::filesystem::is_regular_file( VCF_files[i] ) );

  vector<VCF_variant_info> output;

  // Open each VCF file in turn, and read it line-by-line.
  for ( size_t i = 0; i < VCF_files.size(); i++ ) {

    ifstream in( VCF_files[i].c_str() );

    while ( 1 ) {
      in.getline( _line, _LINE_LEN );
      assert( strlen(_line)+1 < _LINE_LEN );
      if ( in.fail() ) break;

      // ParseVCFLine is where all the hardcore parsing happens.
      struct VCF_variant_info var;
      if ( !ParseVCFLine( _line, filter, var ) ) continue;

      if ( verbose ) cout << "ParseVCF: File " << VCF_files[i] << " reports variant #" << output.size() << ": " << var.all_info() << endl;

      output.push_back(var);
    }

    in.close();
  }


  return output;
}


// ParseVCF on only one VCF file.
vector<VCF_variant_info>
ParseVCF( const string & VCF_file, const struct VCF_input_filter & filter )
{
  return ParseVCF( vector<string>( 1, VCF_file ), filter );
}



// ParseVCF, but only return variants on a given chromosome.
// If chrom == "", this simplifies to the default call to ParseVCF.
vector<struct VCF_variant_info>
ParseVCF( const vector<string> & VCF_files, const string & chrom )
{
  struct VCF_input_filter filter;
  filter.chrom = chrom;
  return ParseVCF( VCF_files, filter );
}


// ParseVCF on only one VCF file, and with only one chromosome.
vector<VCF_variant_info>
ParseVCF( const string & VCF_file, const string & chrom )
{
  struct VCF_input_filter filter;
  filter.chrom = chrom;
  return ParseVCF( vector<string>( 1, VCF_file ), filter );
}





// Set1KGFlags: Take a vector<VCF_variant_info> created by ParseVCF, and set
//     the in1KG flags of the variants in accordance with which variants appear
//     in a second file (presumably of variants only in 1KG).
int
Set1KGFlags( vector<VCF_variant_info> & variants, const string & VCF_1KG_file, const string & chrom )
{
  // Make a lookup table of variant position to VCF_variant_info.  This allows
  // us to take a variant position and lookup or modify that variant.
  map<int, VCF_variant_info *> var_lookup;
  for ( unsigned i = 0; i < variants.size(); i++ ) {
    assert( variants[i].chrom == chrom );
    variants[i].in1KG = false; // i.e., make them not indeterminate anymore
    var_lookup[ variants[i].pos ] = &variants[i];
  }



  // Read the VCF file of variants that are only in 1KG.
  // We don't need to do a full parsing; just look at the first two tokens on
  // each line.
  assert( boost::filesystem::is_regular_file( VCF_1KG_file ) );
  ifstream in( VCF_1KG_file.c_str() );

  int n_in1KG = 0;

  while ( 1 ) {
    in.getline( _line, _LINE_LEN );
    assert( strlen(_line)+1 < _LINE_LEN );
    if ( in.fail() ) break;

    // Split this line into tokens.
    // We only care about the first two tokens, which indicate chrom and pos.
    vector<string> tokens;
    boost::split( tokens, _line, boost::is_any_of("\t") );

    if ( tokens[0] != chrom ) continue; // not this chromosome

    // Find the VCF_variant_info at this position (if there is one there.)
    int pos = boost::lexical_cast<int>( tokens[1] );
    map<int,VCF_variant_info *>::const_iterator it = var_lookup.find(pos);
    if ( it == var_lookup.end() ) continue;

    // Mark this VCF_variant_info as in 1KG.
    it->second->in1KG = true;
    n_in1KG++;
  }

  return n_in1KG;
}





// Parse1KGFreqs: Parse one or more VCF files from 1000 Genomes, and report
//     the frequency of each variant that appears.  The output is a map of
//     variant tag (format: <chrom>_<pos>_<ref-base>_<alt-base>) to frequency.
map<string,double>
Parse1KGFreqs( const vector<string> & VCFs_1KG )
{
  set<string> vars; // variants seen
  map<string,int> n_alt; // number of samples in which a variant appears
  map<string,int> n_total; // number of samples in which a variant could appear

  // Open each VCF file and read it line-by-line.
  for ( unsigned i = 0; i < VCFs_1KG.size(); i++ ) {
    //cout << "Reading file: " << VCFs_1KG[i] << endl;
    ifstream in( VCFs_1KG[i].c_str() );

    while ( 1 ) {
      in.getline( _line, _LINE_LEN );
      assert( strlen(_line)+1 < _LINE_LEN );
      if ( in.fail() ) break;

      if ( _line[0] == '#' ) continue; // ignore commented lines

      // Break up the VCF line into tab-delimited tokens.
      vector<string> tokens;
      boost::split( tokens, _line, boost::is_any_of("\t") );
      assert( tokens.size() == 8 );

      // Skip lines with multiple ref and/or alt alleles.
      if ( tokens[3].size() > 1 || tokens[4].size() > 1 ) continue;

      // Make the variant tag.
      string tag = tokens[0] + "_" + tokens[1] + "_" + tokens[3] + "_" + tokens[4];

      // Match the last token (the one that contains a bunch of pieces of info)
      // to a regex, to get the variant and total coverage.
      // Match this line to a regex of a VCF line.
      boost::cmatch match;
      if ( !boost::regex_search( tokens[7].c_str(), match, boost::regex("AC=(\\d+);\\.*AN=(\\d+)") ) ) {
	cerr << "Parse1KGFreqs: FAIL LINE:\n" << _line << endl;
	exit(1);
      }

      int alt   = boost::lexical_cast<int>( match[1] );
      int total = boost::lexical_cast<int>( match[2] );

      // Add this variant to the datasets.
      if ( vars.find(tag) == vars.end() ) {
	vars.insert(tag);
	n_alt  [tag] = 0;
	n_total[tag] = 0;
      }

      n_alt  [tag] += alt;
      n_total[tag] += total;
      //cout << "tag:\t" << tag << "\tn_alt = " << alt << "\tn_total = " << total << endl;
    }

    in.close();
  }


  // For each variant, divide its number of appearances by its number of
  // possible appearances to get its frequency.
  map<string,double> freqs;

  for ( set<string>::const_iterator it = vars.begin(); it != vars.end(); ++it )
    freqs[*it] = double( n_alt[*it] ) / n_total[*it];

  return freqs;
}
