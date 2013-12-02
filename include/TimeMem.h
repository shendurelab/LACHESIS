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
