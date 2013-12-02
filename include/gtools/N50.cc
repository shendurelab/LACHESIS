// For general documentation, see N50.h
#include "N50.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <set> // needed for greater<T>()
#include <algorithm> // sort
#include <numeric> // accumulate


// Template instantiations.
template int     N50( const vector<int>     & v );
template int64_t N50( const vector<int64_t> & v );
template float   N50( const vector<float>   & v );
template double  N50( const vector<double>  & v );




template<class T>
T
N50( const vector<T>& v )
{
  if ( v.empty() ) return 0;
  
  int64_t sum = accumulate( v.begin(), v.end(), (int64_t) 0 );
  int64_t half = 0;
  
  // Copy v onto a local vector and sort it.
  vector<T> sorted_v = v;
  sort( sorted_v.begin(), sorted_v.end(), greater<T>() );
  
  // Step through sorted v until reaching the halfway point.
  for ( size_t i = 0; i < sorted_v.size(); i++ ) {
    half += sorted_v[i];
    if ( 2 * half == sum && i < sorted_v.size() - 1 ) 
      return (sorted_v[i] + sorted_v[i+1])/2;
    if ( 2 * half >= sum ) return sorted_v[i];
  }
  
  assert(0); // never executed
  return 0;
}



