#include "ssw_cpp.h"

#include "global_variables.h"
#include "bitseq_utils.h"
#include "process_cigar.h"
#include "output_utils.h"

using namespace std;
using namespace boost;

void processSeedMotifWise(tuple<int, int> seed_position, int chunk_start, int &motif_length, int &seed_type, string &sequence_id, string &sequence,
                          int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                          boost::dynamic_bitset<> &N_bset, ofstream &out,
                          StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment,
                          vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);

int longestContinuousMatches(boost::dynamic_bitset<> &bset);
