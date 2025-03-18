#include <boost/dynamic_bitset.hpp>

#include "global_variables.h"

using namespace std;

vector<tuple<int, int, int, int>> processShiftXORsPerfect(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                          int &window_length);

bool retainNestedSeed(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                      int nested_midx, int parent_midx, int bset_size);

bool retainIdenticalSeeds(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                           int nested_midx, int parent_midx, int bset_size);
