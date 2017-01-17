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
 * TimeMem.h
 *
 * This module contains two functions: Time() and MemUsage().
 *
 *
 * Usage example:
 *
 * cout << Time() << ": This program is using " << MemUsage() << " KB of memory." << endl;
 *
 *****************************************************************************/



#ifndef __TIME_MEM_H
#define __TIME_MEM_H

#include <string>
using namespace std;

// Return the current timestamp, in the format:
// Thu Jan  6 16:34:55 2011
string Time();


// Return memory being used by this process in KB.
// If VM = true, returns virtual memory usage; otherwise (default) return resident set usage.
// Returns 0 if unable to read from /proc/self/stat.
double MemUsage( const bool VM = false );


#endif
