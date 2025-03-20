#include "global_variables.h"

using namespace std;

void processCIGARWithPruning(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length, 
                             int &repeat_start, int &repeat_end, int &alignment_length, int &match_units, string &new_cigar,
                             double &purity, double &avg_motifpurity, int &avg_motifindels);

void processCIGARMotifWise(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length,
                           int &repeat_start, int&repeat_end, int &alignment_length, string &cigar_string, double &purity,
                           int &interruptions, double &avg_motifpurity, int &avg_motifindels, int &avg_matchlen);

void motifwiseParameters(string &cigar, int motif_length, double &avg_motifpurity, int &avg_motifindels);

tuple<vector<int>, vector<char>> cigarSplit(string cigar);
