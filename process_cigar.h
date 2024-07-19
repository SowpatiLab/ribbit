#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

tuple<vector<int>, string, float> processCIGARWithPruning(int seed_start, int seed_sequence_length, string &cigar,
                                                          string &seed_sequence, int motif_length);

tuple<vector<int>, string, float> processCIGARMotifWise(int seed_start, int seed_sequence_length, string &cigar,
                                                        string &seed_sequence, int motif_length);