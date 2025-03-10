#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "ssw_cpp.h"

using namespace std;
using namespace boost;


void processSeed(tuple<int, int> seed_position, int seq_start, int &motif_length, int &seed_type, string &sequence_id, string &sequence,
                 int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                 boost::dynamic_bitset<> &N_bset, int &continuous_threshold, ostream &out, vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                 vector<boost::dynamic_bitset<>*> &MATRIX, StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter,
                 StripedSmithWaterman::Alignment &alignment, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);

int longestContinuousMatches(boost::dynamic_bitset<> &bset);

void longRepeatSSWAlignment(string &alignment_cigar, int seed_sequence_length, string &seed_sequence, int motif_length, string &motif,
                            StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter,
                            StripedSmithWaterman::Alignment &alignment, int slice_length);
