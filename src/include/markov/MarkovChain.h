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
 * MarkovChain
 *
 * See MarkovModel.h for documentation.
 *
 *****************************************************************************/


#ifndef _MARKOV_CHAIN_H
#define _MARKOV_CHAIN_H

#include "MarkovModel.h" // superclass
#include <iostream>
#include <vector>
#include <string>
using namespace std;





class MarkovChain : public MarkovModel
{

 public:

  // Constructor function.
  // This function does not fill all of the model parameters.
  // You must call the functions SetInitProbs and SetTransProbs,
  // which are inherited from the MarkovModel class.
  MarkovChain( const int N_states ) :
    MarkovModel( N_states ) {}

  // Destructor.
  ~MarkovChain();

  // Generate a random sequence of symbols based on the initation, transition,
  // and emission probabilities.  This function is the main aspect of the
  // MarkovChain class that sets it apart from a generic MarkovModel.
  // If set_seed = true (default), (re)sets a random seed.
  // Returns a vector of symbols with length chain_length.
  vector<int> GenerateChain( const int chain_length, const bool set_seed = true ) const;


  void Print( ostream & out = cout ) const;


 private:

};


#endif
