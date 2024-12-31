#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

using namespace std;

vector<tuple<int, int, int, int>> processShiftXORswithSubstitutions(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                                    int &window_length, int &window_bitcount_threshold,
                                                                    vector<tuple<int, int, int, int>> &seed_positions_perfect);
