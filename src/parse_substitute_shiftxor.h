#include "global_variables.h"
#include "parse_perfect_shiftxor.h"
#include "seed_utils.h"

using namespace std;

vector<tuple<int, int, int, int, int, int, int>> processShiftXORswithSubstitutions(vector<boost::dynamic_bitset<>> &motif_bsets, vector<boost::dynamic_bitset<>> &perfect_bsets,
                                                                    vector<boost::dynamic_bitset<>> &anchored_bsets, boost::dynamic_bitset<> &N_bset,
                                                                    vector<tuple<int, int, int, int,int,int,int>> &seed_positions_perfect);

void filterPerfectSeeds(vector<tuple<int, int, int, int, int, int, int>> &seed_positions_perfect,
                        vector<tuple<int, int, int, int, int, int, int>> &seed_positions_substut);
