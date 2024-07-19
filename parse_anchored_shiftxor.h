#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

using namespace std;

// converting the shift XOR bitsets to anchored shift XOR bitsets
void generateAnchoredShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                               vector<boost::dynamic_bitset<>> &lsxor_anchor_bsets, int min_shift, int nshifts, int anchor_size);


vector<tuple<int,int,int,int>> processShiftXORsAnchored(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                        int &window_length, int &window_bitcount_threshold, int nshifts, int &min_motif_size,
                                                        int &min_shift, vector<tuple<int, int, int, int>> &seed_positions_perfect,
                                                        vector<tuple<int, int, int, int>> &seed_positions_substut);