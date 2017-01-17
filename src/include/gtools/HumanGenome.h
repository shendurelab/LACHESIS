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
 * HumanGenome
 *
 * This module serves as a way of loading data about the human genome.
 * It contains four functions, each of which returns a const data structure:
 *
 * vector<string> HumanGenome_chroms()
 *     Chromosome names in standard order.
 * map<string,int> HumanGenome_chrom_IDs()
 *     Lookup table of chrom names to their indices in HumanGenome_chroms().
 * map<string,int> HumanGenome_chrom_lengths()
 *     Map of chromosome names to chromosome sizes (not including centromeres).
 * map<string,int> HumanGenome_centromere_locs()
 *     Map of chromosome names to centromere locations.
 *
 * We also assume centromeres are exactly 3MB, which they are in hg19.
 *
 *
 *
 *
 * Josh Burton
 * February 2012
 *
 *****************************************************************************/


#ifndef _HUMAN_GENOME__H
#define _HUMAN_GENOME__H

#include <map>
#include <string>
#include <vector>
#include <fstream>
using namespace std;


// Local includes
#include "ChromInterval.h"



static const size_t HumanGenome_n_chroms = 24; // including X,Y
const int HumanGenome_centro_size = 3000000;




// The set of chromosome names, in standard order.
const vector<string> HumanGenome_chroms();


// A lookup table of chrom names to their indices in HumanGenome_chroms().
const map<string,int> HumanGenome_chrom_IDs();


// Chromosome lengths taken from files in human_ann.  Lengths do not include
// 3MB of centromeres.
const map<string,int> HumanGenome_chrom_lengths();


// Centromere locations.  Centromeres must be exactly 3MB in length, as in hg19.
const map<string,chrom_interval> HumanGenome_centromere_intervals();
// Integer positions indicate the exact middle of centromeres.
const map<string,int> HumanGenome_centromere_locs();




// Return a hard-wired set of ChromIntervals corresponding to the chromosome
// arms that are observed to have partial (distal) or complete LOH in HeLa.
// Complete LOH: 5p, 6q, Xp, Xq
// Near-complete distal LOH: 3q, 6p, 13q, 22q
// Partial distal LOH: 2q, 11q, 19p
// ("Near-complete" may actually be complete, if the centromere is mis-placed.)
// I've called LOH regions manually, by loading HET_freq tracks into IGV.
vector<chrom_interval> HumanGenome_GetHeLaLOHChromosomeArms();




#endif
