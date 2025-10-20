#include "global_variables.h"
#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace boost;

void filterPerfectSeeds(vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                        vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut);

void filterShortSeeds(vector<tuple<int,int,int,int,int,int,int>> &seeds);

bool retainNestedSeed(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                      int nested_midx, int parent_midx, int bset_size);

bool retainIdenticalSeeds(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                           int nested_midx, int parent_midx, int bset_size);

void getBitCount(boost::dynamic_bitset<> &bset, int start_pos, int end_pos, int &motif_bitcount);

int longestContinuousMatches(boost::dynamic_bitset<> &bset);

int longestContinuousMatches(boost::dynamic_bitset<> &bset, int start_pos, int end_pos);

void previouslyIdentifiedMotif(int seed_start, int seed_end, int motif_length, int chunk_start,
                               vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                               string &motif);

void checkAtomicity(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                    vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets);

void filterLowerMatchOverlapSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds,
                                  vector<boost::dynamic_bitset<>> &motif_bsets, vector<boost::dynamic_bitset<>> &perfect_bsets,
                                  vector<boost::dynamic_bitset<>> &anchored_bsets);

void filterNearAtomicSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                           vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets);

void processOverlappingSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                             vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets, int bset_size,
                             vector<set<int>> &skip_atomicity);

void MergeIdenticalMotifSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                              vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets);
