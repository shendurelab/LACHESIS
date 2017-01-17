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


// For general documentation, see WDAG.h
#include "WDAG.h"


// C++ modules.
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>


// Boost libraries.
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp> // is_regular_file
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/log1p.hpp>





WDAG::WDAG()
{
  _N = 0;
  _req_start = _req_end = NULL;
  _best_path_weight = 0;
}



WDAG::~WDAG()
{
  // Delete all WDAGNodes.
  for ( int i = 0; i < _N; i++ )
    delete _nodes[i];
}



// Load an entire WDAG from a file in the format of GS 540 Assignment 2.
void
WDAG::ReadFromFile( const string & filename )
{
  assert( boost::filesystem::is_regular_file( filename ) );

  ifstream file( filename.c_str() );
  char line[100];
  while ( file.getline( line, 100 ) ) {
    vector<string> tokens;
    boost::split( tokens, line, boost::is_any_of("\t") );

    // Identify and parse a line representing a vertex.
    if ( tokens[0] == "V" ) {

      WDAGNode * node = AddNode(); // TODO: add a name based on tokens[1]
      if ( tokens.back() == "START" )
	_req_start = node;
      if ( tokens.back() == "END" )
	_req_end = node;

    }


    // Identify and parse a line representing an edge.
    else if ( tokens[0] == "E" ) {

      string E_name = tokens[1];
      string V_in = tokens[2];
      string V_out = tokens[3];
      double weight = boost::lexical_cast<double>( tokens[4] );

      // Find the parent and child node.  They should both exist already.
      WDAGNode * parent= _nodes[ boost::lexical_cast<int>( V_in .substr(1) ) ];
      WDAGNode * child = _nodes[ boost::lexical_cast<int>( V_out.substr(1) ) ];

      // Tell the child node about its parent.
      // Because children deserve to know the truth about their parents.
      child->AddEdge( parent, E_name, weight );
    }
  }



  file.close();

}


// Write a WDAG to a file in the format of GS 540 Assignment 2.
void
WDAG::WriteToFile( const string & filename ) const
{
  ofstream file( filename.c_str() );

  // Output a line for each node/vertex.
  // Format: "V node_ID"
  for ( int i = 0; i < _N; i++ ) {
    file << "V " << _nodes[i]->_ID;
    if ( _nodes[i] == _req_start )
      file << " START";
    if ( _nodes[i] == _req_end )
      file << " END";
    file << endl;
  }

  // Output a line for each edge.
  // Format: "E edge_name parent_ID child_ID node_weight"
  for ( int i = 0; i < _N; i++ ) {
    const WDAGNode * node = _nodes[i];
    for ( int j = 0; j < node->_n_parents; j++ ) {
      file << "E " << node->_in_e_names[j] << " " << node->_parents[j]->_ID
	   << " " << node->_ID << " " << node->_in_e_weights[j] << endl;
    }
  }

  file.close();
}




// Reserve memory for N nodes.
void
WDAG::Reserve( const int N )
{
  _nodes.reserve( N );
}




// Add a node to the WDAG.  At first, this node has no parents.
// Returns the newly added node, in case you need to use it for something.
WDAGNode *
WDAG::AddNode()
{
  // Initialize a new WDAGNode, with a new unique ID, and add it to the
  // data structures.
  WDAGNode * node = new WDAGNode();
  _nodes.push_back( node );
  node->_ID = _N;
  _N++;

  // Leave placeholders for the node's parents.
  node->_n_parents = 0;
  node->_parents.clear();

  // Leave empty placeholders for the WDAG dynamic programming algorithm.
  node->_best_weight = 0;
  node->_best_parent_id = -1;
  node->_in_e_names.clear();
  node->_in_e_weights.clear();
  node->_fw_prob = 0;
  node->_bw_prob = 0;
  return node;
}






// Run dynamic programming on this WDAG to find out the weights of all paths.
// The WDAG should already contain nodes, loaded with vertex and edge info.
// Find the highest-weight path in the WDAG and report it.
// Note that this algorithm will not work properly unless parents are
// guaranteed to appear before children in _nodes.
void
WDAG::FindBestPath()
{
  _best_path_weight = 0;
  const WDAGNode * best_path_end_node = _nodes[0]; // default solution

  // Process the nodes in the same order in which they appeared in the file.
  for ( int i = 0; i < _N; i++ ) {
    WDAGNode * node = _nodes[i];

    // If constraining the path start, disallow the null path in all nodes
    // except the specified start.
    if ( _req_start != NULL && _req_start != node )
      node->_best_weight = LOG_ZERO;


    // Look at all the parents of this node, and find which one gives the
    // highest-weight path.
    for ( int j = 0; j < node->_n_parents; j++ ) {
      const WDAGNode * parent = node->_parents[j];
      double weight = parent->_best_weight + node->_in_e_weights[j];

      // When we've found the parent that gives the highest-weight path, update
      // the info in this node.
      if ( node->_best_weight < weight ) {
	node->_best_weight = weight;
	node->_best_parent_id = j;
      }

      // Keep track of the highest weight we've seen anywhere.
      if ( node->_best_weight > _best_path_weight ) {
	_best_path_weight = node->_best_weight;
	best_path_end_node = node;
      }
    }

  }


  // If we are constraining the path end, override the best path we've seen
  // and replace it with the best path ending at the specified end.
  if ( _req_end != NULL ) {
    best_path_end_node = _req_end;
    _best_path_weight = _req_end->_best_weight;
  }


  // We've stepped through the entire graph and found the final node of the
  // best path.  Trace backward from this node and find the complete path.
  _best_edges.clear();
  _best_nodes.clear();
  _best_nodes.push_back( best_path_end_node );

  const WDAGNode * node = best_path_end_node;
  assert( node != NULL );
  int parent = node->_best_parent_id;
  while ( parent != -1 ) {
    _best_edges.push_back( node->_in_e_names[parent] );
    node = node->_parents[parent];
    _best_nodes.push_back( node );
    parent = node->_best_parent_id;
  }


  // Flip the path into the proper direction.
  reverse( _best_nodes.begin(), _best_nodes.end() );
  reverse( _best_edges.begin(), _best_edges.end() );
}





// Find the forward and backward probabilities of each node.  Used in the
// Baum-Welch algorithm.
void
WDAG::FindPosteriorProbs()
{
  // We can only do this in a WDAG with a definitive start and end point.
  assert( _req_start != NULL );
  assert( _req_end != NULL );

  // Set the initial (log) probabilities of all nodes.
  // log prob = 0 implies prob = 1; log prob = LOG_ZERO implies prob = 0
  for ( int i = 0; i < _N; i++ ) {
    WDAGNode * node = _nodes[i];

    if ( node == _req_start ) node->_fw_prob = 0;
    else                      node->_fw_prob = LOG_ZERO;
    if ( node == _req_end   ) node->_bw_prob = 0;
    else                      node->_bw_prob = LOG_ZERO;
  }



  // Process the nodes in forward order, and calculate forward probabilities.
  // This assumes, as always, that parents appear before children.
  for ( int i = 0; i < _N; i++ ) {
    WDAGNode * node = _nodes[i];

    // Look at all the edges entering this node, and use this to find the
    // posterior probability of this node.
    for ( int j = 0; j < node->_n_parents; j++ ) {
      const WDAGNode * parent = node->_parents[j];

      // Get the weight from this parent and add it to the posterior prob.
      double weight = parent->_fw_prob + node->_in_e_weights[j];
      node->_fw_prob = lnsum( node->_fw_prob, weight );
    }
  }


  // Now process the nodes in reverse order, and calculate backward probs.
  // Now parents appear AFTER children.
  for ( int i = _N-1; i >= 0; i-- ) {
    WDAGNode * node = _nodes[i];

    // Look at all the edges entering this node, and modify the posterior
    // probabilities of this node's *parents* accordingly.
    // In this way we will eventually step through all edges, always examining
    // all the edges leaving each parent after examining its children.
    for ( int j = 0; j < node->_n_parents; j++ ) {
      WDAGNode * parent = (WDAGNode *) node->_parents[j];

      // Get the weight from this node and add it to the PARENT's prob.
      double weight = node->_bw_prob + node->_in_e_weights[j];
      parent->_bw_prob = lnsum( parent->_bw_prob, weight );
    }
  }



  // Verify that the total forward and backward probabilities are identical
  // (or at least within the range of floating-point error.
  double error = ( Alpha() - Beta() ) / Alpha();
  assert( error < 1e-6 );

}




// Run this after running FindBestPath().
void
WDAG::ReportBestPath( ostream & out ) const
{
  assert( !_best_nodes.empty() );

  out << "WDAG::ReportBestPath: EDGES:" << endl;
  for ( size_t i = 0; i < _best_edges.size(); i++ )
    out << "\t edge[" << i << "]\t" << _best_edges[i] << endl;
  out << endl;
  out << "WDAG::ReportBestPath: NODES:" << endl;
  for ( size_t i = 0; i < _best_nodes.size(); i++ )
    out << "\t node[" << i << "]\t" << _best_nodes[i]->_ID << endl;
  out << endl;
}




// Run this after running FindBestPath().
// Returns the IDs of the best WDAGNodes.
vector<int>
WDAG::BestNodeIDs() const
{
  vector<int> best_node_IDs;
  for ( size_t i = 0; i < _best_nodes.size(); i++ )
    best_node_IDs.push_back( _best_nodes[i]->_ID );
  return best_node_IDs;
}



// Formula to add two numbers in natural log space.
// If z = x + y, then ln(z) = lnsum( ln(x) + ln(y) ).
double
lnsum( const double & ln_x, const double & ln_y )
{
  // If either number is 0 (in normal space), return the other number.
  if ( ln_x == LOG_ZERO ) return ln_y;
  if ( ln_y == LOG_ZERO ) return ln_x;

  if ( ln_x > ln_y )
    return ln_x + boost::math::log1p( exp( ln_y - ln_x ) );
  else
    return ln_y + boost::math::log1p( exp( ln_x - ln_y ) );
}
