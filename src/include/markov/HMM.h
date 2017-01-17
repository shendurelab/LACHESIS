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
 * HMM: Hidden Markov Model
 *
 * The HMM is a special case of the MarkovModel (see MarkovModel.h for further
 * documentation.)  The HMM is applied to situations in which a sequence of
 * data is observed, and each timepoint of data is assumed to be emitted by a
 * state whose identity is hidden.
 *
 *
 * There are two different fundamental types of HMMs: those in which a finite
 * set of possible observations (or "symbols") may be observed at each time
 * point, and those in which the set of possible observations is continuous and
 * thus effectively infinite.  We call these discrete and continuous HMMs,
 * respectively.
 *
 *
 * An HMM contains several types of probabilities, and the types differ between
 * discrete and continuous HMMs.  All HMMs include initation and transition
 * probabilities, as inherited from the MarkovModel class, which still govern
 * the model's random path through state space.  Discrete HMMs contain an
 * "observations" vector giving the IDs of the sequence of observed symbols, as
 * well as a "symbol emission probabilities" matrix that gives the probability
 * of each state emitting each symbol.  When discrete HMMs are trained, their
 * symbol emission probabilities may change.  Continuous HMMs, on the other
 * hand, contain an unchanging matrix of "time emission probabilities" that
 * describe the likelihood of each state emitting the observed data at each
 * timepoint.
 *
 *
 *
 * Josh Burton
 * June 2012
 *
 *****************************************************************************/


#ifndef _HIDDEN_MARKOV_MODEL_H
#define _HIDDEN_MARKOV_MODEL_H

#include "MarkovModel.h" // superclass
#include "WDAG.h"
#include <vector>
#include <string>
using namespace std;





class HMM : public MarkovModel
{

 public:

  // Constructor function.
  // IMPORTANT: If you set N_symbols = 0, this will be a continuous HMM;
  // otherwise, it will be a discrete HMM with N symbols.
  // This function does not fill all of the model parameters.  You must
  // subsequently call the data-loading functions, below.
  HMM( const int N_states, const int N_symbols );

  // Destructor.
  ~HMM();

  /* DATA-LOADING FUNCTIONS */

  // You must call some of these functions before running the HMM.
  // For discrete HMMs, call SetSymbolEmissProbs and SetObservations.
  // For continuous HMMs, call SetTimeEmissProbs.
  // For both, you must also call the other MarkovModel data-loading functions:
  // SetInitProbs and SetTransProbs.

  // In a discrete HMM: Probabilities of each state emitting each symbol.
  // Input these as non-logarithms.  The vector is by states, then symbols.
  void SetSymbolEmissProbs( const vector< vector<double> > & probs );
  // In a discrete HMM: The set of observed symbols.
  void SetObservations( const vector<int> & symbols );
  // In a continuous HMM: Log-likelihoods of each state emitting the observed
  // data at each timepoint.  Input these as logarithms.
  // Note that the vector is by timepoints, then states.
  void SetTimeEmissProbs( const vector< vector<double> > & probs );

  bool HasAllData() const;

  /* HMM ALGORITHMS */

  // Viterbi or Baum-Welch training to improve the transition probabilities.
  // For iterative training, run these functions repeatedly.
  // Return true if any of the probabilities change.
  bool ViterbiTraining( vector<int> & predicted_states );
  bool BaumWelchTraining( double & log_like );


  /* OUTPUT FUNCTIONS UNIQUE TO THE HMM CLASS */

  size_t NTimepoints() const;
  void WriteWDAGToFile( const string & file ) const { to_WDAG().WriteToFile(file); }
  void DrawPNGAtState( const string & PNG_file_head, const size_t T, const size_t depth = 2 ) const;
  void Print( ostream & out = cout ) const;


  //private: // TEMP

  /* HELPER FUNCTIONS */

  // Is this a discrete or a continuous HMM?
  bool is_discrete_HMM() const { return _N_symbols != 0; }

  // Create a WDAG representing this HMM with its current parameters.
  WDAG to_WDAG() const;

  // Change the transition and/or emission probabilities in accordance with
  // calculated data.  This is the final step of both Viterbi and Baum-Welch.
  // Return true if any of the probabilities change.
  bool AdjustProbsToViterbi( const vector<string> & best_path, vector<int> & states );
  bool AdjustProbsToBaumWelch( const WDAG & wdag );




  /* HIDDEN MARKOV MODEL PARAMETERS */

  // Number of obervable symbols.  If 0, this is a continuous HMM.
  const int _N_symbols;

  // Flags for whether or not types of data have been loaded in.
  // (There are more of these in the MarkovModel class.)
  bool _has_symbol_emiss_probs, _has_observations; // used in discrete HMMs
  bool _has_time_emiss_probs; // used in continuous HMMs

  // Flags for whether or not algorithms have been run.
  bool _ran_viterbi, _ran_baum_welch;

  // NOTE: THESE ARE STORED AS LOGARITHMS
  vector< vector<double> > _symbol_emiss_probs; // used in discrete HMMs
  vector< vector<double> >   _time_emiss_probs; // used in continuous HMMs

  // Set of observed symbols.  Used in discrete HMMs.
  vector<int> _observations;


};


#endif
