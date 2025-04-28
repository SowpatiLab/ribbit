#include "global_variables.h"
#include "parse_perfect_shiftxor.h"

using namespace std;

vector<tuple<int, int, int, int>> processShiftXORswithSubstitutions(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                                    int &window_length, int &window_bitcount_threshold,
                                                                    vector<tuple<int, int, int, int>> &seed_positions_perfect);
                                                                    
void filterPerfectSeeds(vector<tuple<int, int, int, int>> &seed_positions_perfect,
                        vector<tuple<int, int, int, int>> &seed_positions_substut);
