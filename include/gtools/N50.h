// N50
// Function to calculate the N50 (length-weighted median) of a vector.
// Return 0 if the vector is empty.
// This code was adapted from math/Functions.cc in the ALLPATHS source code.
// The function is templatized so it works on int, int64_t, float, double
// (template instantiations are in N50.cc.)


#ifndef __N50_H
#define __N50_H

#include <vector>
using namespace std;

template<class T> T N50( const vector<T>& v );

#endif
