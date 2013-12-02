// For documentation, see TrueMapping.h
#include "TrueMapping.h"
#include "TextFileParsers.h" // TODO: make the stuff here text-file-parseable


#include <assert.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <numeric> // accumulate
#include <algorithm> // count, max_element, min_element

// Boost libraries
#include <boost/algorithm/string.hpp> // split
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

// Modules in ~/include (must add -L~/include and -lJ<module> to link)
#include "TimeMem.h"
#include "gtools/HumanGenome.h"
#include "gtools/SAMStepper.h" // TargetLengths





static const unsigned LINE_LEN = 10000;





// Helper function for the TrueMapping constructor.
// Read a *.fasta.names file and find out all the contig names.  Output a vector<string> of contig names and also a lookup table of contig name to contig ID.
void
read_fasta_names( const string & fasta_names_file,
		  vector<string> & contig_names, // output
		  map<string,int> & contig_IDs ) // output
{
  char line[LINE_LEN];
  vector<string> tokens;
  
  contig_names.clear();
  contig_IDs.clear();
  
  // Read the file line-by-line.
  ifstream in;
  in.open( fasta_names_file.c_str(), ios::in );
  
  while ( 1 ) {
    in.getline( line, LINE_LEN );
    assert( strlen(line)+1 < LINE_LEN );
    if ( in.fail() ) break;
    
    // Each line in the file is of the format ">name ..."; that is, the first token starts with a '>' and there may or may not be other tokens.
    boost::split( tokens, line, boost::is_any_of(" \t") );
    if ( tokens[0][0] != '>' ) {
      cout << "ERROR: FASTA names file " << fasta_names_file << " should only have lines beginning with '>'.  Saw the line: " << line << endl;
      assert( tokens[0][0] == '>' );
    }
    string name = tokens[0].substr(1);
    
    // Fill data structures.
    contig_IDs[name] = contig_names.size();
    contig_names.push_back(name);
  }
  
  in.close();
}




// BlastAlignmentVec: Helper struct for the TrueMapping constructor.
// Describes a set of alignments from a query to a set of target sequences, will be used to find where exactly this query is on its target.
struct BlastAlignmentVec
{
  
  // Minimum length of alignments to be considered for the 'main part' of this alignment - unless there are none this long, in which case take the longest
  static const int MIN_ALIGN_LEN = 5000;
  
  
  // Constructor: 
  BlastAlignmentVec( int N_targets )
  {
    align_len_to_start.resize( N_targets );
    align_rc.resize( N_targets, false );
    align_starts.resize( N_targets );
    align_stops .resize( N_targets );
  }
  
  // Add an alignment on one target sequence.
  void Add( int target, int start, int stop ) {
    
    align_starts[target].push_back(start);
    align_stops [target].push_back(stop);
    
    // Iff stop < start, the alignment is RC.  We note the orientation if this is the first alignment added.
    if ( align_len_to_start[target].empty() ) align_rc[target] = ( stop < start );
    
    // Find the length and the actual start site (the lower of start and stop).  Record these data in the map, which automatically sorts by length.
    int start0 = min( start, stop );
    int len = max( start, stop ) - start0;
    
    align_len_to_start[target].insert( make_pair( len, start0 ) );
  }
  
  
  // Determine where this query is on the given target sequence, assuming it is on that sequence (the target is chosen by TabulateAlignsToTarget()).
  // We employ a "growing from a seed" algorithm: start with the region corresponding to one alignment, and keep expanding it until it's as big as the query.
  void FindAlignRegionOnTarget( const int target_ID, const int query_len, int & region_start, int & region_stop ) const {
    
    // No target -> no start or stop
    if ( target_ID == -1 ) { region_start = -1; region_stop = -1; return; }
    assert( target_ID >= 0 && target_ID < (int) align_len_to_start.size() );
    
    // Get the set of alignments on this target.  Note that in the map, they are indicated by length and by start position, and sorted by decreasing length.
    const map<int,int,greater<int> > & aligns = align_len_to_start[target_ID];
    
    // Find the genomic region corresponding to the first (longest) alignment.
    region_start = aligns.begin()->second; // start
    region_stop  = aligns.begin()->second + aligns.begin()->first; // stop = start + len
    assert( region_stop - region_start > 0 );
    
    // Starting with this region as a 'seed', extend it to form a single contiguous interval that includes as many of the alignments on the target sequence as
    // possible, while not exceeding the overall length of the query sequence.
    // Note that, because the alignment may be gapped and may be a slightly different length on the query and the target, the strict requirement that the
    // region length not exceed the query size is overly restrictive and may cause a slightly lower length.  I don't really care.
    for ( map<int,int>::const_iterator it = aligns.begin(); it != aligns.end(); ++it ) {
      int start = min( region_start, it->second );
      int stop  = max( region_stop,  it->second + it->first ); // stop = start + len
      
      // Find the total region size that would ensue if this alignment were included.  If it's not too big, include this region.
      int region_len_if_included = stop - start;
      if ( region_len_if_included > query_len ) {
	//cout << "Excluding:\t" << it->second << "\t" << it->second+it->first << endl;
	continue;
      }
      //cout << "Including:\t" << it->second << "\t" << it->second+it->first << "\t\tNew region:\t" << start << "-" << stop << endl;
      region_stop  = stop;
      region_start = start;
    }
    
    
    if (0) { // verbosity
      int region_size = region_stop - region_start;
      double size_frac = double( region_size ) / query_len;
      cout << "FINAL REGION:\t" << region_start << "\t" << region_stop << "\t\tSIZE:\t" << region_size << "\tquery_len = " << query_len << "\tsize_frac = " << size_frac << endl;
    }
    
    
    // If the first (longest) alignment was rc, flip start and stop to indicate this.
    if ( align_rc[target_ID] ) { int swap = region_start; region_start = region_stop; region_stop = swap; }
  }
  
  
  // DATA
  vector< vector<int> > align_starts, align_stops;
  vector< map<int,int,greater<int> > > align_len_to_start; // vector: for each target, a list of alignment lengths keyed to alignment start locations
  vector<bool> align_rc; // vector: for each target, the orientation of the first (longest) alignment on that target
  
};






// Load in a BLAST output file.
TrueMapping::TrueMapping( const string & species, const string & query_names_file, const string & target_names_file, const string & BLAST_file_head, const string & out_dir, const string & dummy_SAM_file )
  : _species( species )
{
  // Initially, clear the list of query and target names, and the set of alignments.
  _query_names.clear();
  _target_names.clear();
  _target.clear();
  _start.clear();
  _stop.clear();
  _qual_alignability.clear();
  _qual_specificity.clear();
  
  
  
  cout << Time() << ": TrueMapping\t<-\t" << BLAST_file_head << endl;
  assert( boost::filesystem::is_regular_file( query_names_file ) );
  assert( boost::filesystem::is_regular_file( target_names_file ) );
  assert( boost::filesystem::is_regular_file( dummy_SAM_file ) );
  
  // Read the query_names and target_names files, respectively.  Fill _query_names and _target_names with the names of queries in the respective fastas.
  map<string,int> query_names_to_IDs, target_names_to_IDs;
  read_fasta_names(  query_names_file,  _query_names,  query_names_to_IDs );
  read_fasta_names( target_names_file, _target_names, target_names_to_IDs );
  
  // Get query sequence lengths.
  vector<int> query_lengths = TargetLengths( dummy_SAM_file );
  
  // Sanity checks.
  if ( query_lengths.empty() ) {
    cout << "ERROR: SAM file '" << dummy_SAM_file << "' seems to have no SAM/BAM header." << endl;
    assert( !query_lengths.empty() );
  }
  
  if ( (int) query_lengths.size() != NQueries() ) {
    cout << "ERROR: DRAFT_ASSEMBLY_NAMES_FILE (" << query_names_file << ") and SAM file '" << dummy_SAM_file << "' seem to be from different datasets." << endl;
    PRINT2( query_lengths.size(), NQueries() );
    assert( (int) query_lengths.size() == NQueries() );
  }
  
  
  // Now that we know how many query and target contigs there are, reserve memory in local data structures.
  // All query contigs are initially marked as unaligned.  This will be changed if an alignment is seen for it in the BLAST file.
  _target.resize( NQueries(), -1 );
  _start .resize( NQueries(), -1 );
  _stop  .resize( NQueries(), -1 );
  _qual_alignability.resize( NQueries(), 0 );
  _qual_specificity .resize( NQueries(), 0 );
  
  
  
  ifstream in;
  char line[LINE_LEN];
  vector<string> tokens;
  
  
  /* For each query sequence, find the best target (chromosome) - that is, the one containing the plurality of aligned sequence.
   * Also derive two "mapping quality scores" for each query:
   * Unique alignability: What fraction of the bases in this contig appear in exactly one alignment to reference?
   * Target specificity: Of the bases that align uniquely, what fraction aligns to the plurality target.
   * To calculate these numbers, we must read in the *complete* set of BLAST alignments, which is in multiple files
   * These files take the form <BLAST_file_head>.*.blast.out, with * = 1,2,... until all query contigs have been seen.
   */
  
  
  // Filename for 'saved' file.  The first time this function is run on a given assembly, it creates this file to store queries' alignments and qualities.
  // Then in future runs, this file can be loaded, to save runtime.  There are asserts to prevent file mismatches.
  //string save_file = BLAST_file_head + ".TrueMapping_align_data.txt";
  string save_file = out_dir + "/cached_data/TrueMapping.assembly.txt";
  
  if( boost::filesystem::is_regular_file( save_file ) ) {
    
    // Read the save file line-by-line.
    // Each (non-commented) line in the file describes a contig.
    in.open( save_file.c_str(), ios::in );
    
    int query_ID = 0;
    
    while ( 1 ) {
      in.getline( line, LINE_LEN );
      assert( strlen(line)+1 < LINE_LEN );
      if ( in.fail() ) break;
      
      // Skip commented lines.
      if ( line[0] == '#' ) continue;
      
      // Parse the line into its six tokens, as written by TabulateAlignsToTarget().
      boost::split( tokens, line, boost::is_any_of("\t") );
      assert( tokens.size() == 6 );
      assert( boost::lexical_cast<int>( tokens[0] ) == query_ID );
      
      _target[query_ID] = boost::lexical_cast<int>( tokens[1] );
      _start [query_ID] = boost::lexical_cast<int>( tokens[2] );
      _stop  [query_ID] = boost::lexical_cast<int>( tokens[3] );
      _qual_alignability[query_ID] = boost::lexical_cast<double>( tokens[4] );
      _qual_specificity [query_ID] = boost::lexical_cast<double>( tokens[5] );
      
      query_ID++;
    }
    
    // Sanity check.  If this fails, the save file may not match the dataset.  Specifically, if query_ID == 0, the save file may be empty; just delete it.
    if ( query_ID != NQueries() ) PRINT2( query_ID, NQueries() );
    assert( query_ID == NQueries() );
    
  }
  
  
  
  else { // No save file exists.  Parse the alignments, calculate the best target and the quality metrics, and create a save file.  Runtime on human: ~1 min.
    
    cout << Time() << ": Parsing BLAST output to find contig alignments to reference; will cache results at " << save_file << endl;
    ofstream out( save_file.c_str(), ios::out );
    // Make a header for this save file.
    out << "# This file was created by the TrueMapping constructor in TrueMapping.cc" << endl;
    out << "#" << endl;
    out << "# species = " << species << endl;
    out << "# NQueries() = " << NQueries() << endl;
    out << "# NTargets() = " << NTargets() << endl;
    out << "# BLAST_file_head = " << BLAST_file_head << endl;
    out << "#" << endl;
    out << "# There is one row for each query, containing six numbers:" << endl;
    out << "# query_ID\tbest_target\tstart_on_target\tstop_on_target\tunique_alignability\ttarget_specificity" << endl;
    
    
    
    int BLAST_file_ID = 1; // first file is <BLAST_file_head>.1.blast.out
    int query_ID = -1; // will get incremented to 0 for first query
    int N_hits = 0;
    vector<int> align_on_query; // will record which parts of the query sequence are aligned, and to which target
    BlastAlignmentVec BLAST_aligns( NTargets() ); // will record the locations of all BLAST alignments of this query onto the target sequences
    
    
    // Loop through the lines of the BLAST output file, and if necessary through multiple BLAST output files, until we've seen enough queries.
    while ( query_ID+1 < NQueries() ) {
      
      // Open a new file if necessary.
      if ( !in.is_open() ) {
	assert( N_hits == 0 ); // If this fails, the BLAST file has incorrectly reported its number of hits.
	
	string BLAST_chunk_file = BLAST_file_head + "." + boost::lexical_cast<string>(BLAST_file_ID) + ".blast.out";
	cout << Time() << ": Reading from input file: " << BLAST_chunk_file << " (query_ID = " << query_ID+1 << ")" << endl;
	assert( boost::filesystem::is_regular_file( BLAST_chunk_file ) );
	in.open( BLAST_chunk_file.c_str(), ios::in );
	
	BLAST_file_ID++; // increment for next time.
      }
      
      // Read the BLAST output file line-by-line.  When we're done with the file, close it and prepare to open the next one.
      in.getline( line, LINE_LEN );
      assert( strlen(line)+1 < LINE_LEN );
      if ( in.fail() ) { in.close(); continue; }
      
      
      boost::split( tokens, line, boost::is_any_of(" \t") );
      assert( tokens.size() > 0 );
      
      // Parse commented lines to get metadata
      if ( tokens[0] == "#" ) {
	
	assert( tokens.size() >= 3 ); // if this assert fails, the BLAST file is formatted in an unexpected way, maybe due to a blastn version mismatch
	
	// "BLASTN" line indicates the start of a new query
	if ( tokens[1] == "BLASTN" ) {
	  
	  // Tabulate data from previous query.
	  if ( query_ID != -1 ) TabulateAlignsToTarget( query_ID, align_on_query, BLAST_aligns, out );
	  
	  // Start a new query!
	  query_ID++;
	  align_on_query = vector<int>( query_lengths[query_ID], -1 );
	  BLAST_aligns = BlastAlignmentVec( NTargets() );
	}
	
	// "Query:" line gives query name; make sure it matches the name in the query_names_file (the queries must appear in the same order!)
	else if ( tokens[1] == "Query:" ) assert( tokens.back() == _query_names[query_ID] );
	
	// "XXX hits found" line indicates the number of lines that will follow that describe alignments
	else if ( tokens[2] == "hits" ) N_hits = boost::lexical_cast<int>( tokens[1] );
	
	continue;
      }
      
      
      // If control reaches here, this line is non-commented and thus describes a hit on query_ID.
      N_hits--;
      
      assert( tokens.size() >= 10 );
      assert( tokens[0] == _query_names[query_ID] );
      
      // Find where this query aligns.
      int target_ID = target_names_to_IDs[ tokens[1] ];
      int start_on_Q = boost::lexical_cast<int>( tokens[6] );
      int  stop_on_Q = boost::lexical_cast<int>( tokens[7] );
      int start_on_T = boost::lexical_cast<int>( tokens[8] );
      int  stop_on_T = boost::lexical_cast<int>( tokens[9] );
      assert( start_on_Q < stop_on_Q );
      assert(  stop_on_Q <= query_lengths[query_ID] );
      assert( start_on_T >= 0 );
      assert(  stop_on_T >= 0 ); // note: the sign of ( start_on_T - stop_on_T ) may be positive or negative, and determines the orientation of the alignment
      
      // Mark each base on this query as aligning to this target sequence.
      // The code for the align_on_query vector is: -1: no alignment seen (yet); X>=0: aligns to exactly one location on target #X; -2: multiple alignment
      for ( int i = start_on_Q; i < stop_on_Q; i++ ) {
	if ( align_on_query[i] == -1 ) align_on_query[i] = target_ID;
	else align_on_query[i] = -2;
      }
      
      // Record this alignment in the BLAST_aligns struct.
      // If this target turns out to be the 'best' one, we'll use the alignments onto it to determine where the query sequence is on the target.
      BLAST_aligns.Add( target_ID, start_on_T, stop_on_T );
      
      
    }
    
    
    // Don't forget the last listed query!
    assert( query_ID + 1 == NQueries() );
    TabulateAlignsToTarget( query_ID, align_on_query, BLAST_aligns, out );
    
    in.close();
    
  }
  
  
  // Count the number of unaligned sequences.
  int N_unaligned = count( _target.begin(), _target.end(), -1 );
  
  cout << Time() << ": TrueMapping contains " << NQueries() << " query sequences aligned to " << NTargets() << " target sequences (incl. " << N_unaligned << " unaligned)" << endl;
  
  
}






// Create a TrueMapping for a de novo GLM made by chopping up a reference genome into bins of size BIN_SIZE.  The species must be human (for now).
TrueMapping::TrueMapping( const string & species, const int BIN_SIZE, const string & query_names_file, const string & target_names_file )
  : _species( species )
{
  // Initially, clear the list of query and target names, and the set of alignments.
  _query_names.clear();
  _target_names.clear();
  _target.clear();
  _start.clear();
  _stop.clear();
  _qual_alignability.clear();
  _qual_specificity.clear();
  
  
  // Read the query_names and target_names files, respectively.  Fill _query_names and _target_names with the names of contigs in the respective fastas.
  map<string,int> query_names_to_IDs, target_names_to_IDs; // not used
  read_fasta_names(  query_names_file,  _query_names,  query_names_to_IDs );
  read_fasta_names( target_names_file, _target_names, target_names_to_IDs );
  
  
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
  cout << Time() << ": TrueMapping::Print" << endl;
  
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
  cout << Time() << ": TrueMapping::PrintSeqLengthOnTargets" << endl;
  
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
// Counts up all of the bp in a query, examines which target(s) it's aligned to, and assigns it a plurality target (modifies _target).
// Also calculates 2 quality scores for "unique alignability" and "target specificity", and prints them to output file.
void
TrueMapping::TabulateAlignsToTarget( const int query_ID, const vector<int> & align_on_query, const BlastAlignmentVec & BLAST_aligns, ostream & out )
{
  assert( query_ID >= 0 && query_ID < NQueries() );
  
  int query_len = align_on_query.size();
  
  // Tally up data from the align_on_query vector.
  // The code for the align_on_query vector is: -1: no alignment seen (yet); X>=0: aligns to exactly one location on target #X; -2: multiple alignment
  int N_unaligned = 0, N_nonunique = 0;
  vector<int> N_unique_on_target( NTargets(), 0 );
  for ( int i = 0; i < query_len; i++ )
    switch ( align_on_query[i] ) {
    case -1: N_unaligned++; break;
    case -2: N_nonunique++; break;
    default: N_unique_on_target[ align_on_query[i] ]++; break;
    }
  
  // Find the target with the most alignment from this query.
  vector<int>::iterator max = max_element( N_unique_on_target.begin(), N_unique_on_target.end() );
  int N_on_best_target = *max;
  int best_target = max - N_unique_on_target.begin();
  if ( N_on_best_target == 0 ) best_target = -1; // if this query doesn't uniquely align anywhere, mark it as entirely unmapped
  
  // Find the boundaries of this alignment on the target sequence.
  int start;
  int  stop;
  BLAST_aligns.FindAlignRegionOnTarget( best_target, query_len, start, stop );
  
  // Calculate the quality scores.
  int N_aligned = query_len - N_unaligned - N_nonunique;
  double unique_alignability = N_aligned / double ( query_len > 0 ? query_len : 1 );
  double target_specificity = N_on_best_target / double ( N_aligned > 0 ? N_aligned : 1 );
  //PRINT6( query_ID, query_len, N_unaligned, N_nonunique, best_target, N_on_best_target );
  //PRINT3( query_ID, unique_alignability, target_specificity );
  
  
  
  // Set info about this query.
  _target[query_ID] = best_target;
  _start [query_ID] = start;
  _stop  [query_ID] = stop;
  _qual_alignability[query_ID] = unique_alignability;
  _qual_specificity [query_ID] = target_specificity;
  
  
  // Write to the output file.
  out << query_ID << '\t' << best_target << '\t' << start << '\t' << stop << '\t' << unique_alignability << '\t' << target_specificity << endl;
  
  //for ( int i = 0; i * 10000 < query_len; i++ ) {
  //  out << "ALIGNABILITY\t" << unique_alignability << endl;
  //  out << "SPECIFICITY\t"  << target_specificity << endl;
  //}
  
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
  



