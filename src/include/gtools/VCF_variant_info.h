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
 * VCF_variant_info
 *
 * This module contains a single, simple struct: VCF_variant_info.
 *
 * This struct describes a variant and contains all the information in a VCF
 * line about the variant.
 * VCF_variant_info objects are created by ParseVCF(), in FileParsers.h.
 *
 *
 * Josh Burton
 * April 2012
 *
 *****************************************************************************/


#ifndef _VCF_VARIANT_INFO__H
#define _VCF_VARIANT_INFO__H

#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>

// TODO: consider splitting into 2 classes: one that's just "variant", with
// chrom, pos, ref/alt base; and one that's "variant info"

struct VCF_variant_info {

  /* DATA IN VCF */
  string chrom; // e.g., "chr1"
  int pos; // on chrom, 0-indexed
  bool dbSNP; // is it in dbSNP?
  boost::tribool in1KG; // is it in 1KG? (this can be set by calling FileParsers::Set1KGFlags)
  double qual; // "PL" variant quality score (more useful than "QUAL" column)
  string ref_base, alt_base; // identity of ref/alt base(s); will consist of ACGT and will be single-letter except for indels
  int ref_depth, alt_depth; // ref/alt allele depth
  int call; // genotype call in the VCF file: 0 = hmz ref; 1 = het; 2 = hmz alt - note that this is made by GATK and may be fallible!

  /* AUXILIARY DATA, NOT IN VCF - but may be loaded later */
  int CN; // copy number of region around variant; -1 if unknown
  boost::tribool in_repeat; // boost::indeterminate if unknown


  /* FUNCTIONS */

  // Simple processing functions.
  string name() const { return chrom + ":" + boost::lexical_cast<string>(pos); }
  string tag()  const { return chrom + "_" + boost::lexical_cast<string>(pos) + "_" + ref_base + "_" + alt_base; }
  bool is_SNV() const { return ref_base.size() == 1 && alt_base.size() == 1; }
  bool is_indel() const { return !is_SNV(); }

  int total_depth() const { return ref_depth + alt_depth; }
  string call_word() const { return call == 0 ? "HMZ_REF" : call == 1 ? "HET" : "HMZ_ALT"; }
  bool call_het() const { return call == 1; }
  // RAF/AAF/MAF: Reference/alternate/minor allele frequency.
  // We calculate these from allele depths, not from the VCF's "AB=" field
  // (which is sometimes RAF and sometimes AAF, and is always imprecise.)
  // MAF is the smaller of RAF and AAF.
  double RAF() const { return double(ref_depth) / ( ref_depth + alt_depth ); }
  double AAF() const { return double(alt_depth) / ( ref_depth + alt_depth ); }
  double MAF() const { return double( min( ref_depth, alt_depth ) ) / (ref_depth + alt_depth ); }

  // Info dump.
  string all_info() const { return name() + "\tdbSNP? " + boost::lexical_cast<string>(dbSNP) + "\tin1KG? " + boost::lexical_cast<string>(in1KG) + "\tQUAL: " + boost::lexical_cast<string>(qual) + "\t" + ref_base + "->" + alt_base + " (depths: " + boost::lexical_cast<string>(ref_depth) + "," + boost::lexical_cast<string>(alt_depth) + ")\tCALL: " + call_word();
  }
};




#endif
