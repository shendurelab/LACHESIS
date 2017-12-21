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


// For documentation, see TextFileParsers.h
#include "TextFileParsers.h"

// C libraries
#include <assert.h>

// STL declarations
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iostream>

// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>



static const unsigned LINE_LEN = 1000000000;
char LINE[LINE_LEN];



// Helper function and sanity check.  Make sure the variable LINE doesn't overrun the available LINE_LEN.
void check_line_len() {
  if ( strlen(LINE)+1 >= LINE_LEN ) {
    cerr << "Line too long: " << LINE << endl;
    assert( strlen(LINE)+1 < LINE_LEN );
  }
}



// TokenizeFile: Split up a file into lines, and split each line into tokens using whitespace (' ' or '\t') as delimiters.
// Return all tokens as strings, in the output variable tokens.  There are no guarantees about the number of lines or the number of tokens per line.
// Aside from the line delimiter ('\n') and the token delimiters (' ', '\t') this function makes no assumptions whatsoever about the file contents.
// For large files, this function is somewhat slower than parsing files locally because it requires the creation of a large data structure.
// If compress = true, use the token_compress_on flag to compress multiple consecutive whitespace delimiters into one.
void
TokenizeFile( const string & infile, vector< vector<string> > & tokens, const bool & compress, const string & delimiters )
{
  assert( boost::filesystem::is_regular_file( infile ) );
  tokens.clear();

  vector<string> tokens_in_line;

  // Read the file line-by-line.
  ifstream in( infile.c_str(), ios::in );
  while ( 1 ) {
    in.getline( LINE, LINE_LEN );
    check_line_len();
    if ( in.fail() ) break;

    // Convert each line into a set of tokens by splitting on whitespace.
    boost::split( tokens_in_line, LINE, boost::is_any_of(delimiters), compress ? boost::token_compress_on : boost::token_compress_off );
    tokens.push_back( tokens_in_line );
  }

}



// TokenizeCSV: Like TokenizeFile, but recognize as delimiters the regex /\,\s+/ (i.e., a comma followed by any amount of whitespace).
// This is computationally inefficient as written.
void
TokenizeCSV( const string & infile, vector< vector<string> > & tokens )
{
  assert( boost::filesystem::is_regular_file( infile ) );
  tokens.clear();

  vector<string> tokens_in_line;

  // Read the file line-by-line.
  ifstream in( infile.c_str(), ios::in );
  while ( 1 ) {
    in.getline( LINE, LINE_LEN );
    check_line_len();
    if ( in.fail() ) break;

    // Convert each line into a set of comma-delimited tokens by splitting on commas, then removing whitespace from the start of all tokens after the first.
    boost::split( tokens_in_line, LINE, boost::is_any_of(",") );
    for ( size_t i = 0; i < tokens_in_line.size(); i++ ) {
      if ( tokens_in_line[i].empty() ) continue;
      while ( tokens_in_line[i][0] == ' ' || tokens_in_line[i][0] == '\t' ) tokens_in_line[i] = tokens_in_line[i].substr(1);
    }
    tokens.push_back( tokens_in_line );
  }
}




// ParseTabDelimFile: Parse a tab-delimited file.  Return a vector of the <column_ID>'th token (zero-indexed) on each line, recast as objects of class T.
// T can be any type accepted by boost::lexical_cast<>, but for each type there must be a template instantiaion for this function (see below).
// Syntax: "vector<int> v = ParseTabDelimFile<int>( infile )"
// Throws an error if any row in this file has no more than <column_ID> tokens, or if one of the tokens can't be recast as a T.
template<class T> vector<T>
ParseTabDelimFile( const string & infile, const size_t column_ID )
{
  assert( boost::filesystem::is_regular_file( infile ) );

  vector<T> output;
  vector<string> tokens;

  // Read the file line-by-line.
  ifstream in( infile.c_str(), ios::in );
  while ( 1 ) {
    in.getline( LINE, LINE_LEN );
    check_line_len();
    if ( in.fail() ) break;

    // Get the tokens in this line.
    boost::split( tokens, LINE, boost::is_any_of(" \t") );
    assert( tokens.size() > column_ID );

    // Find the token that we want and recast it as class T.
    output.push_back( boost::lexical_cast<T>( tokens[column_ID] ) );
  }

  return output;
}


// Template instantiations for ParseTabDelimFile.
template vector<int>    ParseTabDelimFile( const string & infile, const size_t column_ID );
template vector<double> ParseTabDelimFile( const string & infile, const size_t column_ID );
template vector<string> ParseTabDelimFile( const string & infile, const size_t column_ID );





// GetFastaNames: Input a FASTA filename and return the set of contig names in that FASTA.
// This function uses ParseTabDelimFile() on <fasta-file>.names, and if necessary it calls MakeFastaNamesFile() first to create <fasta-file>.names.
vector<string>
GetFastaNames( const string & fasta_file )
{
  string names_file = fasta_file + ".names";
  if( !boost::filesystem::is_regular_file( names_file ) ) {
    cout << "Calling MakeFastaNamesFile on " << names_file << endl;
    assert( boost::filesystem::is_regular_file( fasta_file ) );
    MakeFastaNamesFile( fasta_file );
  }
  vector<string> names = ParseTabDelimFile<string>( names_file, 0 );

  // Sanity check.
  for ( size_t i = 0; i < names.size(); i++ )
    if ( names[i][0] == '>' )
      cout << "WARNING: In fasta names file " << names_file << ", name `" << names[i] << "' starts with a `>' symbol; this is weird and possibly obsolete." << endl;

  return names;
}




// GetFastaSizes: Input a FASTA filename and return the set of contig lengths in that FASTA.
// This function uses TokenizeFile() on <fasta-file>.FastaSize, and if necessary it runs the command FastaSize to create the file <fasta-file>.FastaSize.
vector<int>
GetFastaSizes( const string & fasta_file )
{
  string FastaSize_file = fasta_file  + ".FastaSize";

  // Make the FastaSize file, if necessary.
  if ( !boost::filesystem::is_regular_file( FastaSize_file ) ) {
    string cmd = "FastaSize " + fasta_file + " > " + FastaSize_file;
    system( cmd.c_str() );
    assert( boost::filesystem::is_regular_file( FastaSize_file ) );
  }

  vector<int> contig_sizes;

  vector< vector<string> > FastaSize_tokens;
  TokenizeFile( fasta_file + ".FastaSize", FastaSize_tokens, true );
  int N_contigs = FastaSize_tokens.size() - 1; // the FastaSize file has one line for each contig, plus a summary line at the end

  // Get the contig lengths from the fasta.FastaSize file.  For each line in the FastaSize file (before the last), there should be 3 tokens: a blank one
  // (whitespace), a number indicating length, and a chromosome name.
  for ( int i = 0; i < N_contigs; i++ ) {
    vector<string> & t = FastaSize_tokens[i];
    assert( t.size() == 3 );
    assert( t[0] == "" );
    contig_sizes.push_back( boost::lexical_cast<int>( t[1] ) );
  }

  return contig_sizes;
}




// MakeFastaNamesFile: Input a FASTA filename.  Create a file at <fasta-file>.names, containing all of the contig names in the FASTA, without the leading '>'.
// This is a wrapper for the Unix command: grep "\>" fasta-file | cut -c2- > fasta-file.names
// After running this, the contig names can be read in via: ParseTabDelimFile<string>( fasta-file.names, 0 )
void
MakeFastaNamesFile( const string & fasta_file )
{
  string cmd = "grep '>' " + fasta_file + " | cut -c2- > " + fasta_file + ".names";
  int ret = system( cmd.c_str() );
  if ( ret > 0 ) cout << "WARNING: MakeFastaNamesFile: system() return nonzero value " << ret << endl;
}






// BlastAlignmentVec: Helper struct for the function ParseBlastAlignmentFiles.
// Describes a set of alignments from a query to a set of target sequences.  Used to find where exactly this query is on its target.
struct BlastAlignmentVec
{

  // Minimum length of alignments to be considered for the 'main part' of this alignment - unless there are none this long, in which case take the longest
  static const int MIN_ALIGN_LEN = 5000;


  // Constructor:
  BlastAlignmentVec( int N_targets )
  {
    align_len_to_start.resize( N_targets );
    align_rc.resize( N_targets, false );
    align_starts.resize( N_targets );
    align_stops .resize( N_targets );
  }

  // Add an alignment on one target sequence.
  void Add( int target, int start, int stop ) {

    align_starts[target].push_back(start);
    align_stops [target].push_back(stop);

    // Iff stop < start, the alignment is RC.  We note the orientation if this is the first alignment added.
    if ( align_len_to_start[target].empty() ) align_rc[target] = ( stop < start );

    // Find the length and the actual start site (the lower of start and stop).  Record these data in the map, which automatically sorts by length.
    int start0 = min( start, stop );
    int len = max( start, stop ) - start0;

    align_len_to_start[target].insert( make_pair( len, start0 ) );
  }


  // Determine where this query is on the given target sequence, assuming it is on that sequence (the target is chosen by TabulateAlignsToTarget()).
  // We employ a "growing from a seed" algorithm: start with the region corresponding to one alignment, and keep expanding it until it's as big as the query.
  void FindAlignRegionOnTarget( const int target_ID, const int query_len, int & region_start, int & region_stop ) const {

    // No target -> no start or stop
    if ( target_ID == -1 ) { region_start = -1; region_stop = -1; return; }
    assert( target_ID >= 0 && target_ID < (int) align_len_to_start.size() );

    // Get the set of alignments on this target.  Note that in the map, they are indicated by length and by start position, and sorted by decreasing length.
    const map<int,int,greater<int> > & aligns = align_len_to_start[target_ID];

    // Find the genomic region corresponding to the first (longest) alignment.
    region_start = aligns.begin()->second; // start
    region_stop  = aligns.begin()->second + aligns.begin()->first; // stop = start + len
    assert( region_stop - region_start > 0 );

    // Starting with this region as a 'seed', extend it to form a single contiguous interval that includes as many of the alignments on the target sequence as
    // possible, while not exceeding the overall length of the query sequence.
    // Note that, because the alignment may be gapped and may be a slightly different length on the query and the target, the strict requirement that the
    // region length not exceed the query size is overly restrictive and may cause a slightly lower length.  I don't really care.
    for ( map<int,int>::const_iterator it = aligns.begin(); it != aligns.end(); ++it ) {
      int start = min( region_start, it->second );
      int stop  = max( region_stop,  it->second + it->first ); // stop = start + len

      // Find the total region size that would ensue if this alignment were included.  If it's not too big, include this region.
      int region_len_if_included = stop - start;
      if ( region_len_if_included > query_len ) {
	//cout << "Excluding:\t" << it->second << "\t" << it->second+it->first << endl;
	continue;
      }
      //cout << "Including:\t" << it->second << "\t" << it->second+it->first << "\t\tNew region:\t" << start << "-" << stop << endl;
      region_stop  = stop;
      region_start = start;
    }


    if (0) { // verbosity
      int region_size = region_stop - region_start;
      double size_frac = double( region_size ) / query_len;
      cout << "FINAL REGION:\t" << region_start << "\t" << region_stop << "\t\tSIZE:\t" << region_size << "\tquery_len = " << query_len << "\tsize_frac = " << size_frac << endl;
    }


    // If the first (longest) alignment was rc, flip start and stop to indicate this.
    if ( align_rc[target_ID] ) { int swap = region_start; region_start = region_stop; region_stop = swap; }
  }


  // DATA
  vector< vector<int> > align_starts, align_stops;
  vector< map<int,int,greater<int> > > align_len_to_start; // vector: for each target, a list of alignment lengths keyed to alignment start locations
  vector<bool> align_rc; // vector: for each target, the orientation of the first (longest) alignment on that target

};



// TabulateAlignsToTarget: Helper function for ParseBlastAlignmentFiles.
// Counts up all of the bp in a query, examines which target(s) it's aligned to, and assigns it a plurality target.
// Also calculates 2 quality scores for "unique alignability" and "target specificity".  Prints all of this data to the output stream.
void
TabulateAlignsToTarget( const int query_ID, const vector<int> & align_on_query, const BlastAlignmentVec & BLAST_aligns, int N_targets, ostream & out )
{
  assert( query_ID >= 0 );

  int query_len = align_on_query.size();

  // Tally up data from the align_on_query vector.
  // The code for the align_on_query vector is: -1: no alignment seen (yet); X>=0: aligns to exactly one location on target #X; -2: multiple alignment
  int N_unaligned = 0, N_nonunique = 0;
  vector<int> N_unique_on_target( N_targets, 0 );
  for ( int i = 0; i < query_len; i++ )
    switch ( align_on_query[i] ) {
    case -1: N_unaligned++; break;
    case -2: N_nonunique++; break;
    default: N_unique_on_target[ align_on_query[i] ]++; break;
    }

  // Find the target with the most alignment from this query.
  vector<int>::iterator max = max_element( N_unique_on_target.begin(), N_unique_on_target.end() );
  int N_on_best_target = *max;
  int best_target = max - N_unique_on_target.begin();
  if ( N_on_best_target == 0 ) best_target = -1; // if this query doesn't uniquely align anywhere, mark it as entirely unmapped

  // Find the boundaries of this alignment on the target sequence.
  int start;
  int  stop;
  BLAST_aligns.FindAlignRegionOnTarget( best_target, query_len, start, stop );

  // Calculate the quality scores.
  int N_aligned = query_len - N_unaligned - N_nonunique;
  double unique_alignability = N_aligned / double ( query_len > 0 ? query_len : 1 );
  double target_specificity = N_on_best_target / double ( N_aligned > 0 ? N_aligned : 1 );
  //PRINT6( query_ID, query_len, N_unaligned, N_nonunique, best_target, N_on_best_target );
  //PRINT3( query_ID, unique_alignability, target_specificity );




  // Write to the output file.
  out << query_ID << '\t' << best_target << '\t' << start << '\t' << stop << '\t' << unique_alignability << '\t' << target_specificity << endl;

  //for ( int i = 0; i * 10000 < query_len; i++ ) {
  //  out << "ALIGNABILITY\t" << unique_alignability << endl;
  //  out << "SPECIFICITY\t"  << target_specificity << endl;
  //}

}




// ParseBlastAlignmentFiles: Input a set of BLAST files describing a set of queries aligning to targets.  Also input query lengths and target names.
// Determine the full extent of each query sequence's position on target, and write the results to outfile.
//
// For each query sequence, find the best target (chromosome) - that is, the one containing the plurality of aligned sequence.
// Then, once we've chosen that target, we employ a "growing from a seed" algorithm: start with the region corresponding to the one best alignment, then keep
// expanding it to include other alignments, until it's as big as the query.
//
// Also derive two "mapping quality scores" for each query:
// Unique alignability: What fraction of the bases in this contig appear in exactly one alignment to reference?
// Target specificity: Of the bases that align uniquely, what fraction aligns to the plurality target?
void
ParseBlastAlignmentFiles( const vector<string> & BLAST_files, const vector<int> & query_lengths, const vector<string> & target_names, const string & outfile )
{
  ifstream in;
  char line[LINE_LEN];
  vector<string> tokens;


  ofstream out( outfile.c_str(), ios::out );

  int N_queries = query_lengths.size();
  int N_targets = target_names.size();

  // Make a lookup table of target name -> ID.
  map<string,int> target_names_to_IDs;
  for ( size_t i = 0; i < target_names.size(); i++ )
    target_names_to_IDs[ target_names[i] ] = i;


  // Make a header for the cache file.  This doesn't get parsed anywhere but it's useful for human readability.
  out << "# This file was created by the function ParseBlastAlignmentFiles in TextFileParsers.cc" << endl;
  out << "#" << endl;
  out << "# N assembly contigs = " << N_queries << endl;
  out << "# N reference contigs = " << N_targets << endl;
  out << "#" << endl;
  out << "# There is one row for each query, containing six numbers:" << endl;
  out << "# query_ID\tbest_target\tstart_on_target\tstop_on_target\tunique_alignability\ttarget_specificity" << endl;


  int file_ID = 0; // index in BLAST_files
  int query_ID = -1; // will get incremented to 0 for first query
  int N_hits = 0;
  vector<int> align_on_query; // will record which parts of the query sequence are aligned, and to which target
  BlastAlignmentVec BLAST_aligns( N_targets ); // will record the locations of all BLAST alignments of this query onto the target sequences


  // Loop through the lines of the BLAST output file, and if necessary through multiple BLAST output files, until we've seen enough queries.
  while ( query_ID+1 < N_queries ) {

    // Open a new file if necessary.
    if ( !in.is_open() ) {
      assert( N_hits == 0 ); // If this fails, the BLAST file has incorrectly reported its number of hits.

      string BLAST_chunk_file = BLAST_files[file_ID++];
      //cout << "Reading from input file: " << BLAST_chunk_file << " (query_ID = " << query_ID+1 << ")" << endl;
      assert( boost::filesystem::is_regular_file( BLAST_chunk_file ) );
      in.open( BLAST_chunk_file.c_str(), ios::in );
    }

    // Read the BLAST output file line-by-line.  When we're done with the file, close it and prepare to open the next one.
    in.getline( line, LINE_LEN );
    check_line_len();
    if ( in.fail() ) { in.close(); continue; }


    boost::split( tokens, line, boost::is_any_of(" \t") );
    assert( tokens.size() > 0 );

    // Parse commented lines to get metadata
    if ( tokens[0] == "#" ) {

      assert( tokens.size() >= 3 ); // if this assert fails, the BLAST file is formatted in an unexpected way, maybe due to a blastn version mismatch

      // "BLASTN" line indicates the start of a new query
      if ( tokens[1] == "BLASTN" ) {

	// Tabulate data from previous query.
	if ( query_ID != -1 ) TabulateAlignsToTarget( query_ID, align_on_query, BLAST_aligns, N_targets, out );

	// Start a new query!
	query_ID++;
	align_on_query = vector<int>( query_lengths[query_ID], -1 );
	BLAST_aligns = BlastAlignmentVec( N_targets );
      }

      // "Query:" line gives query name; make sure it matches the name in the query_names (the queries must appear in the same order!)
      // Actually, don't make sure of this; this is too stringent because different BLAST files have slightly different formats
      // We want to get the token in the line that represents the first word after "Query:", but taking into account the possibility of multiple spaces.
      // This is too hard so let's just comment it out.
      else if ( tokens[1] == "Query:" ) {} // assert( tokens[2] == _query_names[query_ID] );

      // "XXX hits found" line indicates the number of lines that will follow that describe alignments
      else if ( tokens[2] == "hits" ) N_hits = boost::lexical_cast<int>( tokens[1] );

      continue;
    }


    // If control reaches here, this line is non-commented and thus describes a hit on query_ID.
    N_hits--;

    assert( tokens.size() >= 10 );
    if ( target_names_to_IDs.find( tokens[1] ) == target_names_to_IDs.end() ) {
      cout << "ERROR: Target (chromosome) name mismatch.  Name `" << tokens[1] << "' in BLAST file " << BLAST_files[--file_ID] << " not found in reference file." << endl;
      system( ( "rm " + outfile ).c_str() );
      exit(1);
    }

    // Find where this query aligns.
    int target_ID = target_names_to_IDs[ tokens[1] ];
    int start_on_Q = boost::lexical_cast<int>( tokens[6] );
    int  stop_on_Q = boost::lexical_cast<int>( tokens[7] );
    int start_on_T = boost::lexical_cast<int>( tokens[8] );
    int  stop_on_T = boost::lexical_cast<int>( tokens[9] );
    assert( start_on_Q < stop_on_Q );
    assert(  stop_on_Q <= query_lengths[query_ID] );  // if this fails, the query sequences in the BLAST file don't match the target sequences in the SAM file
    assert( start_on_T >= 0 );
    assert(  stop_on_T >= 0 ); // note: the sign of ( start_on_T - stop_on_T ) may be positive or negative, and determines the orientation of the alignment

    // Mark each base on this query as aligning to this target sequence.
    // The code for the align_on_query vector is: -1: no alignment seen (yet); X>=0: aligns to exactly one location on target #X; -2: multiple alignment
    for ( int i = start_on_Q; i < stop_on_Q; i++ ) {
      if ( align_on_query[i] == -1 ) align_on_query[i] = target_ID;
      else align_on_query[i] = -2;
    }

    // Record this alignment in the BLAST_aligns struct.
    // If this target turns out to be the 'best' one, we'll use the alignments onto it to determine where the query sequence is on the target.
    BLAST_aligns.Add( target_ID, start_on_T, stop_on_T );


  }


  // Don't forget the last listed query!
  assert( query_ID + 1 == N_queries ); // if this assert fails, the SAM file you used to find the query_lengths may have an internal inconsistency
  TabulateAlignsToTarget( query_ID, align_on_query, BLAST_aligns, N_targets, out );

  in.close();


}
