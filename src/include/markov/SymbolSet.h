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
 * SymbolSet: A collection of all of the observable symbols that exist in an
 * HMM sequence.  Each symbol is stored as a string.
 *
 *
 *****************************************************************************/


#ifndef _SYMBOL_SET_H
#define _SYMBOL_SET_H

#include <assert.h>
#include <map>
#include <string>
#include <vector>
using namespace std;


class SymbolSet {

 public:

  // Constructors.
  SymbolSet() {}
  SymbolSet( const int n_symbols );
  SymbolSet( const vector<string> & syms );

  void push_back( const string & sym );

  // Convert a symbol to a printable int/char object, or back.
  // Add 28 to any character values, so as to avoid the wackiness of the
  // special ASCII symbols that are less than 28.
  int as_int( const string & sym ) const { assert( _sym_to_ID.find(sym) != _sym_to_ID.end() ); return _sym_to_ID.at(sym); }
  string symbol( int i ) const { return _ID_to_sym.at(i); }
  char as_char( const string & sym ) const { return as_int(sym) + 28; }
  string as_symbol( const unsigned char & c ) const { return symbol( c-28 ); }

  size_t size() const { return _ID_to_sym.size(); }

  // Print a set of symbols, with no spaces in between.
  string SymbolSequence( const vector<int> & symbols ) const;

  // Return a string with characters representing all observable symbols.
  string all_symbols() const;

 private:

  // Lookup table for ID -> symbol.
  vector<string> _ID_to_sym;

  // Lookup table for symbol -> ID.
  map<string,int> _sym_to_ID;

};


#endif
