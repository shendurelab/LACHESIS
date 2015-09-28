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


// For general documentation, see TimeMem.h

#include "TimeMem.h"
#include <string>
#include <fstream>
#include <unistd.h> // sysconf, _SC_PAGE_SIZE


// Return the current timestamp, in the format:
// Thu Jan  6 16:34:55 2011
string Time()
{
  time_t rawtime;
  time ( &rawtime );
  string time_str = ctime( &rawtime );
  time_str.resize( time_str.size()-1 ); // chomp the newline
  return time_str;
}






// Return memory being used by this process in KB.
// If VM = true, returns virtual memory usage; otherwise (default) return resident set usage.
// Returns 0 if unable to read from /proc/self/stat.
// This function has been taken and modified from: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
double MemUsage( const bool VM )
{
   using std::ios_base;
   using std::ifstream;
   using std::string;
   
   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);
   
   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;
   
   // the two fields we want
   //
   unsigned long vsize;
   long rss;
   
   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
   
   stat_stream.close();
   
   if ( VM ) {
     double vm_usage = vsize / 1024.0; // this gives virtual memory usage
     return vm_usage;
   }
   
   else {
     long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
     double resident_set = rss * page_size_kb;
     return resident_set;
   }
}

