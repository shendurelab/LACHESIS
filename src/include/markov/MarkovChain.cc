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


// For general documentation, see MarkovChain.h
#include "MarkovChain.h"

// C and C++ includes
#include <assert.h>
#include <math.h> // fabs exp, log
#include <stdio.h>
#include <stdlib.h> // srand48, drand48

#include <iostream>
#include <vector>
#include <numeric> // accumulate



// log_prob_index
// Given a log-probability vector and a double between 0 and 1, determine which
// element in the vector the double "points" to.
size_t
log_prob_index( const vector<double> & logprobs, double threshold )
{
  assert( threshold >= 0 && threshold < 1 );

  double total = 0;
  for ( size_t i = 0; i < logprobs.size(); i++ ) {
    total += exp( logprobs[i] );
    if ( total >= threshold ) return i;
  }

  // control should never reach here
  cerr << "FAIL: MarkovChain::log_prob_index" << endl;
  cerr << "Debug info:" << endl;
  cerr << "total = " << total << ", threshold = " << threshold << endl;
  for ( size_t i = 0; i < logprobs.size(); i++ ) {
    cerr << "logprobs[" << i << "] = " << logprobs[i] << endl;
  }
  assert(0);
  return -1;
}





MarkovChain::~MarkovChain() {}



// Generate a random sequence of states based on the initation and transition
// probabilities.  This function is the main aspect of the MarkovChain class
// that sets it apart from a generic MarkovModel.
// If set_seed = true (default), (re)sets a random seed.
// Returns a vector of states with length chain_length.
vector<int>
MarkovChain::GenerateChain( const int chain_length, const bool set_seed ) const
{
  vector<int> chain;

  if ( set_seed ) srand48( time(0) ); // set random seed

  // Choose an initial state at random.
  int state = log_prob_index( _init_probs, drand48() );

  // Each state transitions randomly into the next state.
  // Record each state, until the requested number of states has been observed.
  for ( int i = 0; i < chain_length; i++ ) {
    chain.push_back(state);
    state = log_prob_index( _trans_probs[state], drand48() );
  }

  return chain;
}






// Print out useful information about the current state of this HMM.
// NOTE: When printing probabilities, convert from logarithmic scale back to
// regular scale.
void
MarkovChain::Print( ostream & out ) const
{
  out << "MarkovChain::Print is not yet implemented." << endl; // TEMP
}
