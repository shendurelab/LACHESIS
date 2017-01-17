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
 * ParsedArgs.h
 *
 * A tool for parsing command-line arguments of the form:
 * ExecutableName ARG1=value ARG2=value ...
 *
 * The main function is ParseArgs(), which inputs the command-line arguments
 * (argc, argv) and returns a ParsedArgs object containing the arguments as
 * key-value pairs.
 *
 *
 * To create a ParsedArgs, put this at the beginning of your main() function:
 * ParsedArgs args( argc, argv );
 * args.Require( "ARG1" ); // the user MUST submit a value for ARG1
 * args.RequireOrDefault( "ARG2", "value" ); // the user can submit a value
 *                             // for ARG2; but if not, a default value is used
 *
 * To get values out of the ParsedArgs object, write the following:
 * string value = args[ARG1];
 * You can also convert the value into non-string types:
 * int value = args.ValueAsInt( ARG1 );
 * double value = args.ValueAsDouble( ARG1 );
 * bool value = args.ValueAsBool( ARG1 );
 *
 *
 * Syntax notes:
 * -- The argument-value pairs *must* be of the form "ARG=value", with an
 *    equals sign in the middle (and ONLY there) indicating the break between
 *    argument and value.  There should be no spaces on either side of the
 *    equals sign.
 * -- The argument and value may be any alphanumeric strings.  By convention,
 *    the argument is in all caps, but this is not enforced.
 * -- There must be at least one whitespace character in between each
 *    argument-value pair.
 * -- The return value is a ParsedArgs object: documentation below.
 *
 *
 *
 * Josh Burton
 * January 2011
 *
 *****************************************************************************/



#ifndef _PARSED_ARGS
#define _PARSED_ARGS

#include <map>
#include <string>
#include <sstream>


/* ParsedArgs
 *
 * A ParsedArgs object is essentially a map<string,string> that also contains
 * a few functions to convert map values to int or double.
 * To get arguments out of a ParsedArgs, use the following syntax:
 *
 * ParsedArgs args( argc, argv );
 * args[KEY]; // returns the value associated with KEY, as a string.
 * args.ValueAsInt( KEY ); // recasts the value of KEY as an integer.
 * args.ValueAsDouble( KEY ); // recasts the value of KEY as a double.
 * args.ValueAsBool( KEY ); // recasts the value of KEY as a Boolean.
 *
 * Calling these functions on a non-existent KEY will cause a run-time abort.
 *
 * You can force your user to submit KEYs with the function Require( KEY ).
 * Alternatively, you can supply a default value if necessary, with
 * RequireOrDefault( KEY, value ).
 *
 */
class ParsedArgs {

 public:
  ParsedArgs() {}


  void Require( const string & key ) const {
    map<string,string>::const_iterator value_iter = args_map.find( key );
    if ( value_iter == args_map.end() ) {
      cerr << "ERROR: ParsedArgs: User is expected to input a value for argument key `" << key << "'" << endl;
      exit(1);
    }
  }

  void RequireOrDefault( const string & key, const string & value ) {
    map<string,string>::iterator value_iter = args_map.find( key );
    if ( value_iter == args_map.end() )
      args_map[key] = value;
  }


  string operator[]( const string & key ) const {
    map<string,string>::const_iterator value_iter = args_map.find( key );
    if ( value_iter == args_map.end() ) KeyFail( key );
    return value_iter->second;
  }

  int ValueAsInt( const string & key ) const {
    string value = (*this)[key];
    int i;
    istringstream buffer( value );
    buffer >> i;
    if ( !buffer.eof() )
      ConversionFail( value, "ValueAsInt" );
    return i;
  }

  double ValueAsDouble( const string & key ) const {
    string value = (*this)[key];
    double d;
    istringstream buffer( value );
    buffer >> d;
    if ( !buffer.eof() )
      ConversionFail( value, "ValueAsDouble" );
    return d;
  }

  bool ValueAsBool( const string & key ) const {
    string value = (*this)[key];
    // Convert the string to lowercase.
    transform( value.begin(), value.end(), value.begin(), ::tolower );
    // Try as hard as possible to interpret this string as a Boolean.
    if ( value == "false" ) return false;
    if ( value == "true"  ) return true;
    if ( value == "f" ) return false;
    if ( value == "t"  ) return true;
    if ( value == "0" ) return false;
    if ( value == "1" ) return true;
    cerr << "ERROR: ParsedArgs: Can't process argument value `" << (*this)[key] << "' as a Boolean" << endl;
    exit(1); // failed
    return false;
  }


  friend ParsedArgs ParseArgs( int argc, char * argv[] ); // the below function

 private:

  // Error message and abort, if the user asks for a non-existent key.
  void KeyFail( const string & key ) const {
      cerr << "ERROR: ParsedArgs: Can't find requested argument key `" << key << "'" << endl;
      exit(1);
  }

  // Error message and abort, if the user submits a non-numeric "number".
  void ConversionFail( const string & value, const string & function ) const {
      cerr << "ERROR: ParsedArgs: Can't call function `" << function
	<< "' on non-numeric value `" << value << "'" << endl;
      exit(1);
  }

  map<string,string> args_map;
};



ParsedArgs
ParseArgs( int argc, char * argv[] )
{
  ParsedArgs args;

  // Parse each argument individually.  The arguments have already been
  // separated by whitespace (by the shell.)
  // Note that the first argument - argv[0] - is just the command name.
  for ( int i = 1; i < argc; i++ ) {

    // Find the location of the = sign.
    char * equals_loc = strstr( argv[i], "=" );
    if ( equals_loc == NULL ) {
      cerr << "ERROR: ParsedArgs: Couldn't find a `=' token in argument `" << argv[i] << "'" << endl;
      exit(1);
    }

    // Check that there is no more than one = sign.
    char * equals_loc_2 = strstr( equals_loc+1, "=" );
    if ( equals_loc_2 != NULL ) {
      cerr << "ERROR: ParsedArgs: Found multiple `=' tokens in argument `" << argv[i] << "'" << endl;
      exit(1);
    }

    // Grab the strings on either side of the = sign.
    string key = argv[i];
    key.resize( equals_loc - argv[i] );
    string value = equals_loc + 1;

    args.args_map[key] = value;
  }

  return args;
}

#endif
