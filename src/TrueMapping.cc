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


// For documentation, see TrueMapping.h
#include "TrueMapping.h"
#include "TextFileParsers.h" // ParseBlastAlignmentFiles


#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <numeric> // accumulate
#include <algorithm> // count, max_element, min_element

// Boost libraries
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/HumanGenome.h"
#include "gtools/SAMStepper.h" // TargetLengths, TargetNames





static const unsigned LINE_LEN = 10000;












// Load in a set of BLAST alignments for a de novo GLM.  Also load in query and target names.  The dummy_SAM_file is for getting contig lengths.
// Files should exist: <BLAST_file_head>.blast.out, <BLAST_file_head>.*.blast.out (for * = 1,2,3,...)
TrueMapping::TrueMapping( const string & species, const vector<string> & query_names, const vector<string> & target_names, const string & BLAST_file_head, const string & out_dir, const string & dummy_SAM_file )
  : _species( species ),
    _query_names( query_names ),
    _target_names( target_names )
{
  cout << "Creating a TrueMapping" << endl;

  // These parameters are taken from the lengths of the names vectors.
  assert( NQueries() > 0 );
  assert( NTargets() > 0 );

  // Now that we know how many query and target contigs there are, reserve memory in local data structures.
  // All query contigs are initially marked as unaligned.  This will be changed if an alignment is seen for it in the BLAST file.
  _target.resize( NQueries(), -1 );
  _start .resize( NQueries(), -1 );
  _stop  .resize( NQueries(), -1 );
  _qual_alignability.resize( NQueries(), 0 );
  _qual_specificity .resize( NQueries(), 0 );

  // Determine the cache file name and the set of available BLAST files for ReadBlastAlignsFromFileSet, below.
  // The BLAST files take the form <BLAST_file_head>.*.blast.out, with * = 1,2,...
  string cache_file = out_dir + "/cached_data/TrueMapping.assembly.txt";
  vector<string> BLAST_files;
  for ( int i = 1;; i++ ) {
    string file = BLAST_file_head + "." + boost::lexical_cast<string>(i) + ".blast.out";
    if ( boost::filesystem::is_regular_file( file ) ) BLAST_files.push_back( file );
    else break;
  }

  //PRINT( BLAST_files.size() );

  // ReadBlastAlignsFromFileSet: Read the alignments of assembly contigs onto the reference genome.
  // This is initially done by parsing a set of BLAST files (which takes a lot of runtime due to file I/O) and then carefully expanding the BLAST alignments
  // into whole-contig alignments.  The results of this method are written to a cache file.  If the cache filename already exists, we can save time by reading
  // the alignments directly from it.  Either way, the function will fill local variables.
  ReadBlastAlignsFromFileSet( species, dummy_SAM_file, BLAST_files, cache_file );

  // Count the number of unaligned sequences.
  int N_unaligned = count( _target.begin(), _target.end(), -1 );

  cout << "TrueMapping contains " << NQueries() << " query sequences aligned to " << NTargets() << " target sequences (incl. " << N_unaligned << " unaligned)" << endl;
}






// Create a TrueMapping for a de novo GLM made by chopping up a reference genome into bins of size BIN_SIZE.  The species must be human (for now).
TrueMapping::TrueMapping( const string & species, const int BIN_SIZE, const vector<string> & query_names, const vector<string> & target_names )
  : _species( species ),
    _query_names( query_names ),
    _target_names( target_names )
{
  // Initially, clear the set of alignments.
  _target.clear();
  _start.clear();
  _stop.clear();
  _qual_alignability.clear();
  _qual_specificity.clear();

  // These parameters are taken from the lengths of the names vectors.
  assert( NQueries() > 0 );
  assert( NTargets() > 0 );



  // Use the human chromosome lengths to determine the true position of each bin.
  assert( _species == "human" );
  vector<string> human_chroms = HumanGenome_chroms(); // chrom names start with "chr"; doesn't include non-canonical chromosomes

  for ( size_t i = 0; i < human_chroms.size(); i++ ) {

    // Find the length of this chromosome, and the number of bins on it.
    assert( human_chroms[i] == "chr" + _target_names[i] );
    int chrlen = HumanGenome_chrom_lengths().at( human_chroms[i] ) + 3e6; // the 3MB accounts for missing centromere lengths in HumanGenome_chrom_lengths()
    int N_bins = ( chrlen + BIN_SIZE - 1 ) / BIN_SIZE;
    int centro_loc = HumanGenome_centromere_locs().at( human_chroms[i] );

    // For each of these bins, create an alignment to reference that simply describes where the bin is from.
    for ( int j = 0; j < N_bins; j++ ) {

      // If this is on top of a centromere, don't mark it as aligned.
      if ( j * BIN_SIZE - 1.5e6 < centro_loc && (j+1) * BIN_SIZE + 1.5e6 > centro_loc && 0 ) {
	_target.push_back(-1);
	_start .push_back(-1);
	_stop  .push_back(-1);
	//PRINT5( i, j, BIN_SIZE, human_chroms[i], centro_loc );
	continue;
      }

      _target.push_back( i );
      _start .push_back( j * BIN_SIZE );
      _stop  .push_back( min( (j+1) * BIN_SIZE, chrlen ) ); // don't overflow beyond the edge of the chromosome!
    }

  }



  // Sanity check.  If this fails, the bins don't add up right.
  if ( _target.size() != _query_names.size() ) {
    PRINT2( _target.size(), _query_names.size() );
    assert( _target.size() == _query_names.size() );
  }


  // Mark empty quality scores.
  _qual_alignability.resize( NQueries(), 0 );
  _qual_specificity .resize( NQueries(), 0 );
}




// Create a TrueMapping for a non-de novo GLM.
TrueMapping::TrueMapping( const string & species, const int & bin_size, const vector<string> & chrom_names, const map<string,int> & chrom_lengths )
  : _species( species ),
    _target_names( chrom_names )
{
  _query_names.clear();

  // Determine the number of bins in each chromosome of this genome.
  for ( int i = 0; i < NTargets(); i++ ) {
    int chrom_len = chrom_lengths.at( _target_names[i] );
    int N_bins = ceil( double(chrom_len) / bin_size );

    // Each of the bins in a chromosome has a "true mapping" onto a specific location on that chromosome.
    for ( int j = 0; j < N_bins; j++ ) {
      _target.push_back(i);
      _start.push_back( j * bin_size );
      _stop .push_back( min( (j+1) * bin_size, chrom_len ) );
      _qual_alignability.push_back(1);
      _qual_specificity .push_back(1);

      // Make up a descriptive string "name" for this bin.
      _query_names.push_back( _target_names[i] + "_BIN" + boost::lexical_cast<string>(j) );
    }
  }

}






// Find the target ID corresponding to a target name.  If the target name is not found, throw an exception.
int
TrueMapping::TargetID( const string & target_name ) const
{
  // Use an O(n) linear search.  This is slow, but it doesn't matter unless many repeated calls to TargetID are made, which shouldn't be an issue.
  for ( int i = 0; i < NTargets(); i++ )
    if ( _target_names[i] == target_name )
      return i;

  // If control reaches here, none of the target IDs matched the input name.
  cerr << "TrueMapping::TargetID: Couldn't find a target with the name `" << target_name << "'" << endl;
  assert(0);
}



// Print: Print human-readable output of all alignments.  If genome_order == true, sort output lines in reference order; otherwise, sort by query ID.
void
TrueMapping::Print( const bool genome_order ) const
{
  cout << "TrueMapping::Print" << endl;

  // Determine the order in which to print the query sequences.
  vector<int> order( NQueries(), -1 );
  if ( genome_order ) {
    vector<int> query_to_genome_order = QueriesToGenomeOrder();
    for ( int i = 0; i < NQueries(); i++ )
      order[ query_to_genome_order[i] ] = i;
  }
  else
    for ( int i = 0; i < NQueries(); i++ )
      order[i] = i;


  // Print!
  for ( int i = 0; i < NQueries(); i++ ) {
    int ii = order[i];
    cout << ii << "\t" << _query_names[ii] << "\t" << QTargetName(ii) << "\t" << _start[ii] << "\t" << _stop[ii] << endl;
  }
}




// PrintSeqLengthOnTargets: Print a chart report about the total length of query sequence aligned to each target sequence.
// Note that "the total length of query sequence" may include portions of the query sequences that aren't actually aligned to the target in question.
void
TrueMapping::PrintSeqLengthOnTargets( const string & dummy_SAM_file, ostream & out ) const
{
  cout << "TrueMapping::PrintSeqLengthOnTargets" << endl;

  // Get query sequence lengths.
  vector<int> query_lengths = TargetLengths( dummy_SAM_file );
  assert( (int) query_lengths.size() == NQueries() );

  //for ( int i = 0; i < NQueries(); i++ )
  //PRINT4( i, query_lengths[i], QQAlignability(i), QQSpecificity(i) );

  vector<int>       N_queries   ( NTargets(), 0 ); // number of queries on this target
  vector<int64_t> len_queries   ( NTargets(), 0 ); // length of queries on this target
  vector<int>       N_queries_HQ( NTargets(), 0 ); // number of queries on this target (high-quality only, as determined by the thresholds in QHighQuality)
  vector<int64_t> len_queries_HQ( NTargets(), 0 ); // length of queries on this target (high-quality only, as determined by the thresholds in QHighQuality)
  int       N_unaligned = 0;
  int64_t len_unaligned = 0;

  // Loop over all queries, and find what target they align to.
  for ( int i = 0; i < NQueries(); i++ ) {
    int target_ID = QTargetID(i);
    int align_len = query_lengths[i];

    // If the query doesn't map to any target, record it separately.
    if ( target_ID == -1 ) {
      N_unaligned++;
      len_unaligned += align_len;
      continue;
    }

    //int align_len = abs( _start[i] - _stop[i] ); // use this definition to exclude portions of the query sequence not aligned to the target

    // Fill data structures.
    N_queries  [target_ID]++;
    len_queries[target_ID] += align_len;
    if ( QHighQuality(i) ) {
      N_queries_HQ  [target_ID]++;
      len_queries_HQ[target_ID] += align_len;
    }
  }

  for ( int i = 0; i < NTargets(); i++ )
    out << "Target " << i << ": " << TargetName(i) << "\tN = " << N_queries[i] << "\tlen = " << len_queries[i] << "\tN_HQ = " << N_queries_HQ[i] << "\tlen_HQ = " << len_queries_HQ[i] << endl;

  out << "UNALIGNED CONTIGS: N = " << N_unaligned << "\tlen = " << len_unaligned << endl;
}





// Return a mapping of query contig ID to chromosome ID in the genome (as implied in the fasta order.)
// Query contigs that don't map to any target are given the number -1 in the output vector.
vector<int>
TrueMapping::QueriesToChromIDs() const
{
  // If this is human, use the mapping defined in HumanGenome.h.
  if ( _species == "human" ) return QueriesToHumanChromIDs();

  // Convert the local mapping of query contigs onto targets so that it maps onto human chromosome IDs (including -1's for non-canonical chromosomes.)
  vector<int> query_to_fasta_order;
  for ( int i = 0; i < NQueries(); i++ ) {
    if ( _target[i] == -1 ) query_to_fasta_order.push_back( -1 ); // this query didn't map
    else                    query_to_fasta_order.push_back( _target[i] );
  }

  return query_to_fasta_order;
}




// Mapping of contig ID to genome ordering (e.g., the first contig on the first chromosome in the fasta maps to 0, etc.)
// If you call this, consider calling ReorderQueries() afterward.
// Query contigs that don't map to any target are put at the end of the ordering, in increasing order of contig ID.
vector<int>
TrueMapping::QueriesToGenomeOrder() const
{
  if ( _species == "human" ) return QueriesToHumanGenome();


  vector<int> genome_order( NQueries(), -1 );
  int genome_pos = 0;

  // Loop over each chromosome in the genome, in the order in which they appear in the fasta.
  for ( int chrID = 0; chrID < NTargets(); chrID++ ) {
    if ( chrID == -1 ) continue; // no queries align to this target

    // For each chromosome, find all of the query contigs that align to this chromosome.
    // Put these queries' IDs into a multimap, indexed by their alignment position (defined rather broadly as the middle position of the alignment.)
    multimap<int, int> align_pos_to_query_ID;

    for ( int i = 0; i < NQueries(); i++ )
      if ( _target[i] == chrID )
	align_pos_to_query_ID.insert( make_pair( (_start[i]+_stop[i]) / 2, i ) );

    // Iterate through the multimap; this returns all of the query indices, sorted by their alignment position.
    // Assign each successive query index to an increasing genome ID.
    for ( multimap<int,int>::const_iterator it = align_pos_to_query_ID.begin(); it != align_pos_to_query_ID.end(); ++it )
      genome_order[it->second] = genome_pos++;
  }


  // The query contigs that do not map to any targets have not yet been assigned an order.  Do that now.
  for ( int i = 0; i < NQueries(); i++ )
    if ( genome_order[i] == -1 ) {
      assert( _target[i] == -1 );
      genome_order[i] = genome_pos++;
    }

  assert( genome_pos == NQueries() );

  return genome_order;
}



// Find the first contig on a chromosome.  Useful for non-de novo assemblies.
int
TrueMapping::FirstContigOnChrom( const int chrom_ID ) const
{
  assert( chrom_ID < NTargets() );
  vector<int>::const_iterator it = find( _target.begin(), _target.end(), chrom_ID );
  return int( it - _target.begin() );
}



// Return a mapping of query contig ID to chromosome ID in the human genome.
// This assumes the chromosomes are named in accordance with the standard in HumanGenome.h.
// Query contigs that don't map, or that map to non-canonical chromosomes, are given the number -1 in the output vector.
vector<int>
TrueMapping::QueriesToHumanChromIDs() const
{
  assert( _species == "human" );

  // Make a mapping of target IDs (in this TrueMapping's index) to human chromosome IDs.
  vector<int> target_to_chrom = TargetToHumanChromIDs();

  // Now, convert the local mapping of query contigs onto targets so that it maps onto human chromosome IDs (including -1's for non-canonical chromosomes.)
  vector<int> query_to_chrom;
  for ( int i = 0; i < NQueries(); i++ ) {
    if ( _target[i] == -1 ) query_to_chrom.push_back( -1 ); // this query didn't map
    else                    query_to_chrom.push_back( target_to_chrom[ _target[i] ] ); // this query mapped, but maybe onto a non-canonical chromosome
  }

  return query_to_chrom;
}



// Return a mapping of query contig ID to overall ordering in the human genome (e.g., the first query contig on chr1 maps to 0, and so on.)
// This assumes the chromosomes are named in accordance with the standard in HumanGenome.h.
// Query contigs that don't map, or that map to non-canonical chromosomes, are put at the end of the ordering, in increasing order of contig ID.
vector<int>
TrueMapping::QueriesToHumanGenome() const
{
  assert( _species == "human" );

  // Make a mapping of target IDs (in this TrueMapping's index) to human chromosome IDs.
  vector<int> target_to_chrom = TargetToHumanChromIDs();

  vector<int> genome_order( NQueries(), -1 );
  int genome_pos = 0;

  // Loop over each chromosome in the human genome.  Index the chroms by chrID (e.g., 0 = "chr1".)
  for ( size_t chrID = 0; chrID < HumanGenome_n_chroms; chrID++ ) {

    // For each chromosome, find all of the query contigs that align to this chromosome.
    // Put these queries' IDs into a multimap, indexed by their alignment position (defined rather broadly as the middle position of the alignment.)
    multimap<int, int> align_pos_to_query_ID;

    for ( int i = 0; i < NQueries(); i++ )
      if ( _target[i] != -1 ) // skip unmapped contigs for now
	if ( target_to_chrom[ _target[i] ] == (int) chrID )
	  align_pos_to_query_ID.insert( make_pair( (_start[i]+_stop[i]) / 2, i ) );

    // Iterate through the multimap; this returns all of the query indices, sorted by their alignment position.
    // Assign each successive query index to an increasing genome ID.
    for ( multimap<int,int>::const_iterator it = align_pos_to_query_ID.begin(); it != align_pos_to_query_ID.end(); ++it )
      genome_order[it->second] = genome_pos++;
  }


  // The remaining query contigs that do not map to any targets (or map to non-canonical targets) have not yet been assigned an order.  Do that now.
  for ( int i = 0; i < NQueries(); i++ )
    if ( genome_order[i] == -1 )
      genome_order[i] = genome_pos++;

  assert( genome_pos == NQueries() );

  return genome_order;
}





// Remove a target sequence from consideration by marking all query sequences that map to it as unmapped.  (Don't remove its name from the lists though.)
// This function is useful if a reference fasta contains unannotated or unplaced scaffolds (e.g., the "het" and "U" reference contigs in Drosophila.)
void
TrueMapping::RemoveTarget( const int target_ID )
{
  assert( target_ID >= 0 );
  assert( target_ID < NTargets() );

  bool verbose = false;

  if ( verbose ) cout << "TrueMapping::RemoveTarget on target " << target_ID << "\t(name: " << _target_names[target_ID] << ")" << flush;

  int N_removed = 0;

  for ( int i = 0; i < NQueries(); i++ )
    if ( _target[i] == target_ID ) {
      _target[i] = -1;
      _start [i] = -1;
      _stop  [i] = -1;
      _qual_alignability[i] = 0;
      _qual_specificity [i] = 0;
      N_removed++;
    }

  if ( verbose ) cout << "\t" << N_removed << " queries now marked as unmapped." << endl;
}






// Merge two target sequences - that is, mark them as both the same sequence in actuality.  Assign a new name to the combined sequence.
// This function is useful if a reference fasta contains two arms of a chromosome as separate target sequences (e.g., 2L/2R and 3L/3R in Drosophila.)
void
TrueMapping::MergeTargets( const int target_ID_1, const int target_ID_2, const string & merged_name )
{
  assert( target_ID_1 >= 0 );
  assert( target_ID_2 >= 0 );
  assert( target_ID_1 < NTargets() );
  assert( target_ID_2 < NTargets() );
  if ( _target_names[target_ID_1] != merged_name )
    assert( find( _target_names.begin(), _target_names.end(), merged_name ) == _target_names.end() ); // new name can't already be the name of a target

  bool verbose = false;

  if ( verbose ) cout << "TrueMapping::MergeTargets is merging two targets (#" << target_ID_1 << ", name: " << _target_names[target_ID_1] << "; #" << target_ID_2 << ", name: " << _target_names[target_ID_2] << ") into a single target (new #" << target_ID_1 << ", name: " << merged_name << ")" << endl;

  int N_merged_1 = 0, N_merged_2 = 0;

  for ( int i = 0; i < NQueries(); i++ ) {

    // The new target inherits the same ID as the old target 1.  This matters only for internal bookkeeping.
    if ( _target[i] == target_ID_1 )
      N_merged_1++;

    // Queries that had previously aligned to the old target 2 must be moved over.
    else if ( _target[i] == target_ID_2 ) {
      _target[i] = target_ID_1;
      N_merged_2++;
    }
  }

  if ( verbose ) cout << "\t\t" << N_merged_1 << " and " << N_merged_2 << " queries were merged from the old targets onto the new one." << endl;

  // Rename the new target.  (Target 2 can keep its name in the vector; it doesn't matter because now nothing aligns to it.)
  _target_names[target_ID_1] = merged_name;
}










// Reorder the query sequences.  This accepts output from QueriesToGenomeOrder() and QueriesToHumanGenome().
void
TrueMapping::ReorderQueries( const vector<int> & new_ordering )
{
  assert( (int) new_ordering.size() == NQueries() );

  // Verify that the new_ordering is one-to-one.
  {
    vector<int> ordering_rev( NQueries(), -1 );
    for ( int i = 0; i < NQueries(); i++ ) {
      assert( ordering_rev[ new_ordering[i] ] == -1 );
      ordering_rev[ new_ordering[i] ] = i;
    }
  }


  vector<string> new_query_names( NQueries(), "" );
  vector<int> new_target( NQueries(), -1 );
  vector<int> new_start ( NQueries(), -1 );
  vector<int> new_stop  ( NQueries(), -1 );
  vector<double> new_QA ( NQueries(), 0 );
  vector<double> new_QS ( NQueries(), 0 );

  for ( int i = 0; i < NQueries(); i++ ) {
    int j = new_ordering[i];
    assert( j >= 0 && j < NQueries() );
    assert( new_target[j] == -1 ); // this will fail if the new_ordering is not one-to-one

    new_query_names[j] = _query_names[i];
    new_target[j] = _target[i];
    new_start [j] = _start [i];
    new_stop  [j] = _stop  [i];
    new_QA[j] = _qual_alignability[i];
    new_QS[j] = _qual_specificity [i];

  }

  _query_names = new_query_names;
  _target = new_target;
  _start  = new_start;
  _stop   = new_stop;
  _qual_alignability = new_QA;
  _qual_specificity  = new_QS;
}





// Helper function for the TrueMapping constructor.
// Read alignment info from the BLAST files and write them in a simple format to TrueMapping_file.  If TrueMapping_file already exists, just read from it
// directly, to save runtime.  Either way, load the alignment data into this TrueMapping object.
void
TrueMapping::ReadBlastAlignsFromFileSet( const string & species, const string & dummy_SAM_file, const vector<string> & BLAST_files, const string & TrueMapping_file )
{

  // If the cache file doesn't already exist, we must create it.
  // Parse the alignments, calculate the best target and the quality metrics, and create a cache file.  Runtime on human: ~1 min.
  if( !boost::filesystem::is_regular_file( TrueMapping_file ) ) {


    // Get query sequence lengths.
    assert( boost::filesystem::is_regular_file( dummy_SAM_file ) );
    vector<int> query_lengths = TargetLengths( dummy_SAM_file );

    // Sanity checks.
    assert( _query_names == TargetNames( dummy_SAM_file ) ); // if this fails, the SAM file isn't aligned to the draft assembly fasta file

    if ( query_lengths.empty() ) {
      cout << "ERROR: SAM file '" << dummy_SAM_file << "' seems to have no SAM/BAM header." << endl;
      assert( !query_lengths.empty() );
    }

    if ( (int) query_lengths.size() != NQueries() ) {
      cout << "ERROR: List of query names (from <DRAFT_ASSEMBLY_FASTA>.names) and SAM file '" << dummy_SAM_file << "' seem to be from different datasets." << endl;
      PRINT2( query_lengths.size(), NQueries() );
      assert( (int) query_lengths.size() == NQueries() );
    }

    // Do the parsing!
    cout << "Parsing BLAST files to find contig alignments to reference; will cache results at " << TrueMapping_file << endl;
    ParseBlastAlignmentFiles( BLAST_files, query_lengths, _target_names, TrueMapping_file );
  }



  // The cache file now exists.  Read it line-by-line.  Each (non-commented) line in the file describes a contig.
  assert( boost::filesystem::is_regular_file( TrueMapping_file ) );

  vector< vector<string> > tokens;
  TokenizeFile( TrueMapping_file, tokens );

  int query_ID = 0;

  for ( size_t i = 0; i < tokens.size(); i++ ) {

    const vector<string> & line = tokens[i];
    if ( line[0][0] == '#' ) continue; // skip commented lines
    assert( line.size() == 6 );

    //PRINT6( line[0], line[1], line[2], line[3], line[4], line[5] );
    assert( query_ID == boost::lexical_cast<int>( line[0] ) );
    _target[query_ID] = boost::lexical_cast<int>( line[1] );
    _start [query_ID] = boost::lexical_cast<int>( line[2] );
    _stop  [query_ID] = boost::lexical_cast<int>( line[3] );
    _qual_alignability[query_ID] = boost::lexical_cast<double>( line[4] );
    _qual_specificity [query_ID] = boost::lexical_cast<double>( line[5] );

    query_ID++;
  }


  // Sanity check.  If this fails, the save file may not match the dataset.  Specifically, if query_ID == 0, the save file may be empty; just delete it.
  if ( query_ID != NQueries() ) {
    PRINT2( query_ID, NQueries() );
    cout << "ERROR: Assembly size in draft assembly fasta doesn't seem to match the cached file at " << TrueMapping_file << ".  Maybe the cached file is empty or corrupted due to an earlier aborted run?  If so, just delete it and let Lachesis re-create it." << endl;
  }
  assert( query_ID == NQueries() );

}








// Make a mapping of *target* ID to chromosome ID in the human genome.
// This assumes the chromosomes are named in accordance with the standard in HumanGenome.h.
vector<int>
TrueMapping::TargetToHumanChromIDs() const
{
  assert( _species == "human" );

  const vector<string> human_chrom_names = HumanGenome_chroms();
  const vector<string>::const_iterator not_found = human_chrom_names.end();


  // Initialize the vector to have all -1's.  For chromosome names that aren't canonical human-genome chroms, the value will stay at -1.
  vector<int> target_to_chrom( NTargets(), -1 );


  // For each chromosome name seen, look for a human-genome chromosome of the same name.
  for ( int i = 0; i < NTargets(); i++ ) {

    // Allow the prefix "chr" to appear or not - e.g., "X" and "chrX" are both acceptable.
    string chrom = _target_names[i];
    vector<string>::const_iterator  human_chrom = find( human_chrom_names.begin(), not_found, chrom );
    if ( human_chrom == not_found ) human_chrom = find( human_chrom_names.begin(), not_found, ( "chr" + chrom ) );
    if ( human_chrom == not_found ) continue;

    int human_chrID = human_chrom - human_chrom_names.begin(); // e.g., "0" for chr1, "22" for chrX
    target_to_chrom[i] = human_chrID;
  }

  return target_to_chrom;
}
