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


// For general documentation, see MarkovModel.h
#include "MarkovModel.h"

// C and C++ includes
#include <assert.h>
#include <math.h> // fabs, exp, log
#include <stdio.h>
#include <iostream>
#include <vector>





// Assert that this vector represents a set of probabilities.
// I.e.: the numbers sum to 1, and none of them is less than 0.
// Local-use function.
void assert_prob_vector( const vector<double> & v, const int size ) {

  if ( size != 0 )
    assert( int( v.size() ) == size );

  double prob_sum = 0;
  for ( int i = 0; i < size; i++ ) {
    if ( v[i] < 0 ) {
      cerr << "MarkovModel::assert_prob_vector: This vector has an element less than 0: v[" << i << "] = " << v[i] << endl;
      assert( v[i] >= 0 );
    }
    prob_sum += v[i];
  }

  if ( fabs( prob_sum - 1 ) > 1e-6 ) { // allow floating-point rounding error
    cerr << "MarkovModel::assert_prob_vector:  The elements in this probability vector do not add to 1.  Instead they add to: " << prob_sum << endl;
    cerr << "Probability vector:" << endl;
    for ( int i = 0; i < size; i++ )
      cerr << "v[" << i << "] = " << v[i] << endl;
    assert( prob_sum == 1 );
  }
}




// Constructor.
MarkovModel::MarkovModel( const int N_states )
  : _N_states( N_states )
{
  // Clear the flags for loaded data.
  _has_init_probs   = false;
  _has_trans_probs  = false;

  // Initialize probability matrices.
  _init_probs .resize( _N_states );
  _trans_probs.resize( _N_states );

  for ( int i = 0; i < N_states; i++ )
    _trans_probs[i].resize( _N_states );
}




MarkovModel::~MarkovModel() {}





// Load the initial state probabilities.
void
MarkovModel::SetInitProbs( const vector<double> & probs )
{
  // Before using this probability set, verify that it makes sense.
  assert_prob_vector( probs, _N_states );

  // Convert the probabilities to log scale.
  for ( int i = 0; i < _N_states; i++ )
    _init_probs[i] = log( probs[i] );

  _has_init_probs = true;
}


// Load the initial state probabilities.
// Assume an equal likelihood of starting in any state.
// This functionis a wrapper to SetInitProbs().
void
MarkovModel::SetInitProbsUniform()
{
  SetInitProbs( vector<double>( _N_states, 1.0 / _N_states ) );
}




// Load the state-to-state transition probabilities.
// The input is an N_states x N_states matrix; probs[i][j] describes the
// probability of state i going to state j.
void
MarkovModel::SetTransProbs( const vector< vector<double> > & probs )
{
  // Before using this probability set, verify that it makes sense.
  assert( int( probs.size() ) == _N_states );
  for ( int i = 0; i < _N_states; i++ )
    assert_prob_vector( probs[i], _N_states );

  // Convert the probabilities to log scale.
  for ( int i = 0; i < _N_states; i++ )
    for ( int j = 0; j < _N_states; j++ )
      _trans_probs[i][j] = log( probs[i][j] );

  _has_trans_probs = true;
}


// Load the state-to-state transition probabilities.
// Assume each state has an equal chance of transitioning to each other state.
// (This is commonly used as an initial assumption before training.)
// This function is a wrapper to SetTransProbs().
void
MarkovModel::SetTransProbsUniformSwitchProb( const double switch_prob )
{
  assert( switch_prob >= 0 && switch_prob <= 1 );
  double diag_prob = 1 - (_N_states-1) * switch_prob;
  assert( diag_prob > 0 ); // if this fails, switch_prob is too high

  // Initialize the NxN matrix.
  vector< vector<double> > trans_probs( _N_states, vector<double>( _N_states, switch_prob ) );
  // Modify the diagonal.
  for ( int i = 0; i < _N_states; i++ )
    trans_probs[i][i] = diag_prob;

  SetTransProbs( trans_probs );
}








// Get the frequency of one of the states, as determined by the MarkovModel.
// This frequency is a DERIVED value and can only be found if one of the
// AdjustProbsTo...() functions has already been run.
double
MarkovModel::GetStateFreq( const int state_ID ) const
{
  assert( state_ID < _N_states );
  assert( !_state_freqs.empty() ); // vector is filled by AdjustProbsTo...()
  return _state_freqs[state_ID];
}


// Get the transition frequency from one state (state A) to another (state B):
// that is, the probability that the model, if in state A, will transition to
// state B.
// This probability may be prior (if LoadTransitionProbs has been run, but
// AdjustProbsTo...() has not) or posterior (if AdjustProbsTo...() has been
// run).
double
MarkovModel::GetTransitionFreq( const int stateA, const int stateB ) const
{
  assert( _has_trans_probs );
  assert( stateA < _N_states );
  assert( stateB < _N_states );
  return exp( _trans_probs[stateA][stateB] );
}
