#include "ssw_cpp.h"

#include "global_variables.h"
#include "bitseq_utils.h"
#include "process_cigar.h"
#include "output_utils.h"
#include "parse_smallmotif_seed.h"

using namespace std;
using namespace boost;

void processSeed(tuple<int, int> seed_position, int seq_start, int &motif_length, int &seed_type, string &sequence_id, string &sequence,
                 int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                 boost::dynamic_bitset<> &N_bset, ostream &out, vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                 vector<boost::dynamic_bitset<>*> &MATRIX, StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter,
                 StripedSmithWaterman::Alignment &alignment, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);
