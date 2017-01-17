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
 * WDAG: Weighted Directed Acyclic Graph
 *
 *****************************************************************************/


#ifndef _WDAG_H
#define _WDAG_H

#include <assert.h>
#include <cmath> // exp, log
#include <iostream>
#include <map>
#include <vector>
#include <string>
using namespace std;



// A node in a WDAG.
class WDAGNode {

  // WDAG and HMM can see all of the WDAGNode member variables.
  friend class WDAG;
  friend class HMM;

 public:

  // Add an edge that enters this WDAGNode.
  // Note that each edge has a list of the nodes that enter it only - so
  // children know about their parents, but not vice versa.
  void AddEdge( const WDAGNode * parent, const string & E_name, const double & weight ) {
    assert( parent != NULL );
    assert( this != NULL );
    _n_parents++;
    _parents.push_back( parent );
    _in_e_names.push_back( E_name );
    _in_e_weights.push_back( weight );
  }

  int ID() const { return _ID; }

 private:

  int _ID; // unique identifier


  // Information about this node's parents and its entering edges.
  int _n_parents;
  vector<const WDAGNode *> _parents; // This node's immediate parents
  vector<string> _in_e_names; // Names of in-edges (same order as 'parents')
  vector<double> _in_e_weights; // Weights of in-edges (same order as 'parents')


  // Information about the best path entering this node.
  double _best_weight; // Weight of highest-weight path entering this node
  int _best_parent_id; // Index of parent that gives the highest-weight path (-1 if start of path)

  // Forward and backward probabilities (in log space).  Used in Baum-Welch
  // training only, during the function FindPosteriorProbs().
  double _fw_prob, _bw_prob;


};




class WDAG {

  friend class HMM;

 public:

  // Constructor.
  // This function does not fill all of the model parameters.  You must
  // subsequently call the data-loading functions, below.
  WDAG();

  // Destructor.
  ~WDAG();


  // Load an entire WDAG from a file in the format of GS 540 Assignment 2.
  void ReadFromFile( const string & filename );
  void WriteToFile( const string & filename ) const;

  // Reserve memory for N nodes.
  void Reserve( const int N );

  // Add a node to the WDAG.  At first, this node has no parents.
  // Returns the newly added node, in case you need to use it for something.
  WDAGNode * AddNode();
  int N() const { return _N; }
  WDAGNode * GetNode( const int i ) const { return _nodes[i]; }
  WDAGNode * GetLastNode() const { return _nodes.back(); }

  // Set required start/end nodes.  DO NOT set these to WDAGNodes outside of this WDAG.
  void SetReqStart( const WDAGNode * node ) { _req_start = node; }
  void SetReqEnd  ( const WDAGNode * node ) { _req_end   = node; }

  // Run dynamic programming on this WDAG to find out the weights of all paths.
  void FindBestPath(); // Viterbi
  void FindPosteriorProbs(); // Baum-Welch

  // Reporting functions.  Call a Find function before running these; otherwise they will give no output (or fail asserts).

  // Requires FindBestPath().
  double BestWeight() const { return _best_path_weight; }
  void ReportBestPath( ostream & out = cout ) const;
  vector<int> BestNodeIDs() const;

  // Requires FindPosteriorProbs().
  double Alpha() const { return _req_end  ->_fw_prob; }
  double Beta () const { return _req_start->_bw_prob; }




 private:


  int _N; // number of nodes
  vector<WDAGNode *> _nodes; // nodes are keyed by ID and allocated

  // Required start/end nodes, if any.
  // These are not allocated - they're just pointers.
  const WDAGNode * _req_start, * _req_end;

  /* RESULTS
     These are created by FindBestPath() (Viterbi) but not by FindPosteriorProbs() (Baum-Welch) */
  double _best_path_weight;
  vector<string> _best_edges;
  vector<const WDAGNode *> _best_nodes;

};






// Log of zero.
#include <limits.h>
#define LOG_ZERO (-INFINITY)


// Formula to add two numbers in natural log space.
// If z = x + y, then ln(z) = lnsum( ln(x) + ln(y) ).
// Taken from:  http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
double lnsum( const double & ln_x, const double & ln_y );




#endif
