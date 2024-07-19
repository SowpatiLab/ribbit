#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "ssw_cpp.h"

using namespace std;
using namespace boost;


void processSeed(tuple<int, int> seed_position, int &motif_length, int &seed_type, string &sequence_id, string &sequence, int &sequence_length, 
                 boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                 boost::dynamic_bitset<> &N_bset, int &continuous_threshold, ostream &out, vector<boost::dynamic_bitset<>*> &MATRIX,
                 StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment);

int longestContinuousMatches(boost::dynamic_bitset<> &bset);