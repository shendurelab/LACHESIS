// For documentation, see HumanGenome.h
#include "HumanGenome.h"



// Local includes
#include "ChromInterval.h"


// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/lexical_cast.hpp>


// The set of chromosome names, in standard order.
const vector<string>
HumanGenome_chroms()
{
  vector<string> chroms;
  chroms.push_back("chr1");
  chroms.push_back("chr2");
  chroms.push_back("chr3");
  chroms.push_back("chr4");
  chroms.push_back("chr5");
  chroms.push_back("chr6");
  chroms.push_back("chr7");
  chroms.push_back("chr8");
  chroms.push_back("chr9");
  chroms.push_back("chr10");
  chroms.push_back("chr11");
  chroms.push_back("chr12");
  chroms.push_back("chr13");
  chroms.push_back("chr14");
  chroms.push_back("chr15");
  chroms.push_back("chr16");
  chroms.push_back("chr17");
  chroms.push_back("chr18");
  chroms.push_back("chr19");
  chroms.push_back("chr20");
  chroms.push_back("chr21");
  chroms.push_back("chr22");
  chroms.push_back("chrX");
  chroms.push_back("chrY");
  return chroms;
}


// A lookup table of chrom names to their indices in HumanGenome_chroms().
const map<string,int>
  HumanGenome_chrom_IDs()
{
  vector<string> chroms = HumanGenome_chroms();
  map<string,int> chrom_IDs;
  for ( unsigned i = 0; i < chroms.size(); i++ )
    chrom_IDs[ chroms[i] ] = i;
  return chrom_IDs;
}



// Chromosome lengths taken from files in human_ann.  Lengths do not include
// 3MB of centromeres.
const map<string,int>
HumanGenome_chrom_lengths()
{
  map<string,int> chrom_lens;
  chrom_lens["chr1"] = 246250621;
  chrom_lens["chr2"] = 240199373;
  chrom_lens["chr3"] = 195022430;
  chrom_lens["chr4"] = 188154276;
  chrom_lens["chr5"] = 177915260;
  chrom_lens["chr6"] = 168115067;
  chrom_lens["chr7"] = 156138663;
  chrom_lens["chr8"] = 143364022;
  chrom_lens["chr9"] = 138213431;
  chrom_lens["chr10"] = 132534747;
  chrom_lens["chr11"] = 132006516;
  chrom_lens["chr12"] = 130851895;
  chrom_lens["chr13"] = 112169878;
  chrom_lens["chr14"] = 104349540;
  chrom_lens["chr15"] = 99531392;
  chrom_lens["chr16"] = 87354753;
  chrom_lens["chr17"] = 78195210;
  chrom_lens["chr18"] = 75077248;
  chrom_lens["chr19"] = 56128983;
  chrom_lens["chr20"] = 60025520;
  chrom_lens["chr21"] = 45129895;
  chrom_lens["chr22"] = 48304566;
  chrom_lens["chrX"] = 152270560;
  chrom_lens["chrY"] = 56373566;
  return chrom_lens;
}



static const unsigned LINE_LEN = 1000;



// Centromere locations.  Centromeres must be exactly 3MB in length, as in hg19.
const map<string,chrom_interval>
HumanGenome_centromere_intervals()
{
  map<string,chrom_interval> intervals;
  
  // This file gives centromere locations and is available from the UCSC
  // Genome Browser.
  const string centro_file = "../human_ann/centromeres.txt";
  
  // Read the file line-by-line.
  char line[LINE_LEN];
  ifstream in( centro_file.c_str(), ios::in );
  while ( 1 ) {
    in.getline( line, LINE_LEN );
    assert( strlen(line)+1 < LINE_LEN );
    if ( in.fail() ) break;
    
    if ( line[0] == '#' ) continue; // ignore commented lines
    
    // Split up the line into tokens.
    vector<string> tokens;
    boost::split( tokens, line, boost::is_any_of("\t") );
    int centro_start = boost::lexical_cast<int>( tokens[2] );
    int centro_stop  = boost::lexical_cast<int>( tokens[3] );
    assert( centro_stop - centro_start == HumanGenome_centro_size );
    intervals[ tokens[1] ] = chrom_interval( tokens[1], centro_start, centro_stop );
    
  }
  
  
  return intervals;
}


// Centromere locations.  Centromeres must be exactly 3MB in length, as in hg19.
// Integer positions indicate the exact middle of centromeres.
const map<string,int>
HumanGenome_centromere_locs()
{
  const map<string,chrom_interval> intervals = HumanGenome_centromere_intervals();
  
  map<string,int> interval_centers;
  for ( map<string,chrom_interval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it )
    interval_centers[it->first] = ( it->second.start + it->second.stop ) / 2;
  
  return interval_centers;
}




// Return a hard-wired set of ChromIntervals corresponding to the chromosome
// arms that are observed to have partial (distal) or complete LOH in HeLa.
// Complete LOH: 5p, 6q, Xp, Xq
// Near-complete distal LOH: 3q, 6p, 13q, 22q
// Partial distal LOH: 2q, 11q, 19p
// ("Near-complete" may actually be complete, if the centromere is mis-placed.)
// I've called LOH regions manually, by loading HET_freq tracks into IGV.
vector<chrom_interval>
HumanGenome_GetHeLaLOHChromosomeArms()
{
  vector<chrom_interval> arms;
  // Centromere locations come from centromeres.txt.
  // Telomere locations come from the human fasta.
  arms.push_back( chrom_interval( "chr2", 106690500, 243199373 ) ); // 2q
  arms.push_back( chrom_interval( "chr3", 94505000, 198022430 ) ); // 3q
  arms.push_back( chrom_interval( "chr5", 0, 46405641 ) ); // 5p
  arms.push_back( chrom_interval( "chr6", 0, 57182000 ) ); // 6p
  arms.push_back( chrom_interval( "chr6", 58830166, 171115067 ) ); // 6q
  arms.push_back( chrom_interval( "chr11", 102210000, 135006516 ) ); // 11q
  arms.push_back( chrom_interval( "chr13", 19200000, 115169878 ) ); // 13q
  arms.push_back( chrom_interval( "chr19", 0, 12885000 ) ); // 19p
  arms.push_back( chrom_interval( "chr22", 17515000, 51304566 ) ); // 22q
  arms.push_back( chrom_interval( "chrX", 0, 58632012 ) ); // Xp
  arms.push_back( chrom_interval( "chrX", 61632012, 155270560 ) ); // Xq
  return arms;
}
