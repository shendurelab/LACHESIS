// For documentation, see TextFileParsers.h
#include "TextFileParsers.h"


// C libraries
#include <assert.h>

// STL declarations
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>



static const unsigned LINE_LEN = 1000000;
char LINE[LINE_LEN];




// TokenizeFile: Split up a file into lines, and split each line into tokens using whitespace (' ' or '\t') as delimiters.
// Return all tokens as strings, in the output variable tokens.
// Aside from the line delimiter ('\n') and the token delimiters (' ', '\t') this function makes no assumptions whatsoever about the file contents.
// For large files, this function is somewhat slower than parsing files locally because it requires the creation of a large data structure.
// If compress = true, use the token_compress_on flag to compress multiple consecutive whitespace delimiters into one.
void
TokenizeFile( const string & infile, vector< vector<string> > & tokens, const bool & compress )
{
  assert( boost::filesystem::is_regular_file( infile ) );
  tokens.clear();
  
  vector<string> tokens_in_line;
  
  // Read the file line-by-line.
  ifstream in( infile.c_str(), ios::in );
  while ( 1 ) {
    in.getline( LINE, LINE_LEN );
    if ( strlen(LINE)+1 >= LINE_LEN ) {
      cerr << "Line too long: " << LINE << endl;
      assert( strlen(LINE)+1 < LINE_LEN );
    }
    if ( in.fail() ) break;
    
    // Convert each line into a set of tokens by splitting on whitespace.
    boost::split( tokens_in_line, LINE, boost::is_any_of(" \t"), compress ? boost::token_compress_on : boost::token_compress_off );
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
    assert( strlen(LINE)+1 < LINE_LEN );
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
  assert( boost::filesystem::is_regular_file( fasta_file ) );
  
  string names_file = fasta_file + ".names";
  if( !boost::filesystem::is_regular_file( names_file ) ) {
    cout << "Calling MakeFastaNamesFile on " << names_file << endl;
    MakeFastaNamesFile( fasta_file );
  }
  return ParseTabDelimFile<string>( names_file, 0 );
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
