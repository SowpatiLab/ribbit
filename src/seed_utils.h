#include "global_variables.h"
#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace boost;

void filterPerfectSeeds(vector<tuple<int, int, int, int>> &seed_positions_perfect,
                        vector<tuple<int, int, int, int>> &seed_positions_substut);

void filterShortSeeds(vector<tuple<int, int, int, int>> &seeds);

bool retainNestedSeed(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                      int nested_midx, int parent_midx, int bset_size);

bool retainIdenticalSeeds(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                           int nested_midx, int parent_midx, int bset_size);
