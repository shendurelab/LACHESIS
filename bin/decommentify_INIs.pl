#!/usr/bin/perl -w
use strict;

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


# decommentify_INI.pl
#
# A simple script to remove all comments (and blank lines) from an INI file.  This action can be undone with the script commentify_INI.pl.
#
# Josh Burton
# August 2013


# Check for input arguments.
if ( @ARGV != 1 ) {
    print STDERR "\nSyntax: $0 <INI-file>   (writes to stdout)\n\n";
    exit;
}


# Open the filename provided.
die "Can't find input file `$ARGV[0]`: $!\n" unless -e $ARGV[0];
open IN, '<', $ARGV[0] or die;

# Filter out empty lines and lines starting with a comment ('#').  Print all other lines.
map { print unless /^[\n\#]/ } <IN>;

close IN;
