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
 * MarkovModel
 *
 * This class describes a Markov model.  A Markov model, in general form, is
 * a system that steps through a series of timepoints and is in one of a set
 * of possible states at each timepoint.  Between timepoints, it changes states
 * according to the Markov property: the state at each timepoint depends only
 * on the state at the previous timepoint.
 *
 *
 * The MarkovModel class is a superclass and cannot be instantiated directly.
 * You probably want to use one of its subclasses:
 *
 * -- MarkovChain: In a Markov chain, the states are directly observable.
 *    Markov chains are typically used to produce randomized output based on
 *    transition and emission probabilities that are known (possibly through
 *    training.)
 * -- HMM: In a Hidden Markov Model (HMM), the states can't be directly
 *    observed; all you see is a series of observations, each of which suggests
 *    but does not guarantee what the underlying state is.  HMMs are typically
 *    used to find the underlying states of a system with noisy observed data.
 *
 *****************************************************************************/


#ifndef _MARKOV_MODEL_H
#define _MARKOV_MODEL_H

#include <vector>
#include <iostream>
using namespace std;



// Assert that this vector represents a set of probabilities.
void assert_prob_vector( const vector<double> & v, const int size = 0 );



class MarkovModel
{

 public:

  // Constructor function.
  // This function does not fill all of the model parameters.  You must
  // subsequently call the data-loading functions, below.
  MarkovModel( const int N_states );

  // Destructor.
  ~MarkovModel();

  /* DATA-LOADING FUNCTIONS */
  /* Each of these must be called before running the model. */

  // Initial state probabilities
  void SetInitProbs( const vector<double> & probs );
  void SetInitProbsUniform();
  // State-to-state transition probabilities
  void SetTransProbs( const vector< vector<double> > & probs );
  void SetTransProbsUniformSwitchProb( const double switch_prob );


  /* FUNCTIONS TO GET DATA FROM THE MARKOV MODEL */

  // These functions return frequencies in regular (non-log) form.
  double GetStateFreq( const int state_ID ) const;
  double GetTransitionFreq( const int stateA, const int stateB ) const;

  int NStates() const { return _N_states; }
  virtual void Print( ostream & out ) const = 0;


 protected:


  /* MARKOV MODEL PARAMETERS */

  const int _N_states; // number of states in the model

  // Flags for whether or not types of data have been loaded in.
  bool _has_trans_probs, _has_init_probs;

  // Model probabilities.
  // NOTE: THESE ARE ALL STORED AS LOGARITHMS.
  vector<double> _init_probs; // initial state probabilities
  vector< vector<double> > _trans_probs; // state-to-state transition probs


  /* DERIVED DATA */

  // Posterior probabilities of each state (NOT stored in log form).
  // These are calculated during a call to AdjustProbsTo...().
  vector<double> _state_freqs;


};


#endif
