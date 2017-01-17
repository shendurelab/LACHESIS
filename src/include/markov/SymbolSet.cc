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


// For general documentation, see SymbolSet.h
#include "SymbolSet.h"

//#include <assert.h>
#include <string>
#include <stdio.h>


SymbolSet::SymbolSet( const int n_symbols )
{
  for ( int i = 0; i < n_symbols; i++ ) {
    char symbol[100];
    sprintf( symbol, "%d", i );
    _ID_to_sym.push_back( symbol );
    _sym_to_ID[ symbol ] = i;
  }

}



// Constructor: Construct the whole SymbolSet at once.
SymbolSet::SymbolSet( const vector<string> & syms )
{
  _ID_to_sym = syms;

  for ( size_t i = 0; i < syms.size(); i++ )
    _sym_to_ID[ syms[i] ] = i;
}








void
SymbolSet::push_back( const string & sym )
{
  _sym_to_ID[sym] = _ID_to_sym.size();
  _ID_to_sym.push_back( sym );
}




// Print a set of symbols, with no spaces in between.
string
SymbolSet::SymbolSequence( const vector<int> & symbols ) const
{
  string seq;
  for ( size_t i = 0; i < symbols.size(); i++ )
    seq += symbol( symbols[i] );
  return seq;
}


// Return a string with characters representing all observable symbols.
string
SymbolSet::all_symbols() const
{
  string symbols( size(), '\0' );
  for ( size_t i = 0; i < size(); i++ )
    symbols[i] = i + 28;

  return symbols;
}
