#include <boost/dynamic_bitset.hpp>

#include "global_variables.h"
#include "merge_types.h"
#include "binomial_thresholds.h"
#include "seed_utils.h"

using namespace std;

// converting the shift XOR bitsets to anchored shift XOR bitsets
void generateAnchorShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                               vector<boost::dynamic_bitset<>> &lsxor_anchor_bsets, int anchor_size);

// converting the shift XOR bitsets to anchored shift XOR bitsets
void generatePerfectShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                             vector<boost::dynamic_bitset<>> &lsxor_perfect_bsets, int anchor_size);

vector<tuple<int,int,int,int,int,int,int>> processShiftXORsAnchored(vector<boost::dynamic_bitset<>> &lshift_anchored_bsets,
                                                        vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                                                        vector<boost::dynamic_bitset<>> &lshift_perfect_bsets,
                                                        boost::dynamic_bitset<> &N_bset,
                                                        vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                                                        vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut,
                                                        int chunk_start);
