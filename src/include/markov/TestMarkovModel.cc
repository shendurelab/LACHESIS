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


/******************************************************************************
 *
 * TestMarkovModel.cc
 *
 *
 * Josh Burton
 *
 *****************************************************************************/


// C libraries
#include <assert.h>
#include <stdio.h>

// STL declarations
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;


// Local declarations
#include "MarkovChain.h"
#include "HMM.h"
#include "SymbolSet.h"

#include <boost/filesystem.hpp> // is_regular_file



// Define the allowable states of the model.
// This enum creates convenient numbered labels for each of the states.
enum { NEUTRAL, CONSERVED };





void TestMarkovChain()
{
  const int N_STATES = 2;
  const double SWITCH_PROB = 0.06;

  cout << "Creating a " << N_STATES << "-state Markov chain.\nEach state has a " << 100*SWITCH_PROB << "% chance of switching to the other at any given time." << endl << endl;

  MarkovChain GC_content_MC( N_STATES );

  // Set basic initiation and transition probabilities.
  GC_content_MC.SetInitProbsUniform();
  GC_content_MC.SetTransProbsUniformSwitchProb( SWITCH_PROB );

  // Simulate a 100-base sequence.
  const int SEQ_LENGTH = 100;
  vector<int> seq = GC_content_MC.GenerateChain( SEQ_LENGTH, false );

  cout << "Random (non-seeded) sequence: ";
  for ( size_t i = 0; i < seq.size(); i++ )
    cout << seq[i];



  cout << endl << endl;

}







// Read a file where each line is of the form "STRING ... INT".
// The values are normalized so they sum to 1.
// Used on anc_rep_counts.txt and codon1_2_counts.txt.
void ParseCountColumnFile( const string & filename, vector<string> & keys,
			   vector<double> & values )
{
  assert( boost::filesystem::is_regular_file( filename ) );
  keys.clear();
  values.clear();


  vector<int> values_raw; // not normalized
  string key;
  int value;
  int value_sum = 0;

  // Read the file line-by-line.
  ifstream file( filename.c_str() );
  char line[1000];
  while ( file.getline( line, 1000 ) ) {

    // Get the key-value pair on each line, and store it.
    istringstream iss(line);
    iss >> key;
    iss >> value;
    value_sum += value;
    keys.push_back( key );
    values_raw.push_back( value );
  }


  // Normalize the values.
  values.resize( values_raw.size(), 0 );
  for ( size_t i = 0; i < values.size(); i++ )
    values[i] = values_raw[i] / double( value_sum );



  file.close();
}





// Parse a multiple alignment file and return the sequence of observed symbols.
// This function is hard-wired to work with ENm008.aln - so it expects three
// aligned sequences, called "hg18", "canFam2", and "mm9".
vector<int>
ParseMultipleAlignFile( const string & filename, const SymbolSet & symbol_set )
{
  vector<int> observations;

  // Hard-wired sequence number and nametags.
  int n_seqs = 3;
  vector<string> seqs( n_seqs, "" );
  vector<string> seq_names;
  seq_names.push_back( "hg18" );
  seq_names.push_back( "canFam2" );
  seq_names.push_back( "mm9" );


  // Read through the file and look for the sequence nametags.
  ifstream file( filename.c_str() );
  char line_c[5000];
  string line;
  while ( file.getline( line_c, 5000 ) ) {
    line = line_c;

    // See if this line matches any of our nametags.
    int match = -1;
    for ( int i = 0; i < n_seqs; i++ )
      if ( line.substr( 0, seq_names[i].size() ) == seq_names[i] )
	match = i;
    if ( match == -1 ) continue;

    // Get the portion of the line after the nametag: i.e. the sequence itself.
    line = line.substr( seq_names[match].size() );
    string::size_type not_whitespace = line.find_first_not_of( " \t" );
    line.erase( 0, not_whitespace );

    // Remember this sequence.
    seqs[match] = line;

    // Assume that the sequences always come in order, so if we have a sequence
    // of the last type, we've got a complete set.
    if ( match+1 == n_seqs ) {
      int match_len = line.size();
      for ( int i = 0; i < n_seqs; i++ )
	assert( match_len == int( seqs[i].size() ) );

      // At each step of the match, combine characters from each sequence in
      // order to determine the observed symbol.
      string observation( n_seqs, '\0' );
      for ( int i = 0; i < match_len; i++ ) {
	for ( int j = 0; j < n_seqs; j++ )
	  observation[j] = seqs[j][i];

	observations.push_back( symbol_set.as_int( observation ) );
      }
    }

  }


  return observations;
}







int main( int argc, char * argv[] )
{
  //cout << Time() << ": TestMarkovModel!" << endl;

  TestMarkovChain();


  // Load the alignment counts files.  Use these to create a set of observable
  // symbols, and also to determine the states' emission probabilities.
  // CS = conserved state; NS = neutral/non-conserved state
  vector<string> NS_syms;
  vector<double> NS_probs;
  ParseCountColumnFile( "anc_rep_counts.txt", NS_syms, NS_probs );
  vector<string> CS_syms;
  vector<double> CS_probs;
  ParseCountColumnFile( "codon1_2_counts.txt", CS_syms, CS_probs );

  // Create a SymbolSet, encapsulating all the observed symbols seen in
  // the alignment count files.
  assert( NS_syms == CS_syms );
  SymbolSet symbol_set( NS_syms );


  // Initialize an HMM with 2 states and all the symbols of the SymbolSet.
  HMM hmm( 2, symbol_set.size() );


  // Set initiation probabilities, as defined in the Assignment instructions.
  vector<double> init_probs( 2, 0 );
  init_probs[NEUTRAL  ] = 0.95;
  init_probs[CONSERVED] = 0.05;
  hmm.SetInitProbs( init_probs );


  // Set transition probabilities, as defined in the Assignment instructions.
  // These state names come from an enum in HMM.h.
  vector< vector<double> > trans_probs( 2, vector<double>( 2, 0 ) );
  trans_probs[NEUTRAL  ][NEUTRAL  ] = 0.95;
  trans_probs[NEUTRAL  ][CONSERVED] = 0.05;
  trans_probs[CONSERVED][NEUTRAL  ] = 0.1;
  trans_probs[CONSERVED][CONSERVED] = 0.9;
  hmm.SetTransProbs( trans_probs );


  // Set emission probabilities, in accordance with the data in the alignment
  // counts files.
  vector< vector<double> > emiss_probs;
  emiss_probs.push_back( NS_probs );
  emiss_probs.push_back( CS_probs );
  hmm.SetSymbolEmissProbs( emiss_probs );



  // Parse the multiple alignment file and find the sequence of observed
  // symbols.  Load this sequence into the HMM.
  vector<int> obs = ParseMultipleAlignFile( "ENm008.aln", symbol_set );
  hmm.SetObservations( obs );






  cout << "<gs540_hw assignment=\"8\" name=\"Josh Burton\" email=\"jnburton@uw.edu\">" << endl;
  cout << "\t<results>" << endl;



  // Run Viterbi training.
  int N_viterbi = 1; // Answer converges after 6 iterations
  for ( int i = 1; i <= N_viterbi; i++ ) {
    //cout << Time() << ": Begin Viterbi training: iteration " << i << " of " << N_viterbi << endl;
    vector<int> states;
    hmm.ViterbiTraining(states);
    hmm.Print();
  }



  cout << "\t</results>" << endl;
  cout << "\t</program>" << endl;
  cout << "</gs540_hw>" << endl;



  // TEMP
  cout << "Drawing PNG at test.png..." << endl;
  hmm.DrawPNGAtState( "test1", 0 );
  hmm.DrawPNGAtState( "test2", 2 );
  hmm.DrawPNGAtState( "test3", 4 );
  hmm.DrawPNGAtState( "test4", 498138 );
  hmm.DrawPNGAtState( "test5", 498140 );



  //cout << Time() << ": Done!" << endl;
  return 0;
}
