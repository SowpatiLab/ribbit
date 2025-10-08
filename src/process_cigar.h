#include "global_variables.h"
#include "cigar_utils.h"

using namespace std;

void processCIGARWithPruning(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length, 
                             int &repeat_start, int &repeat_end, int &alignment_length, int &match_units, string &new_cigar,
                             double &purity, double &avg_motifpurity, int &avg_motifindels);

void processCIGARMotifWise(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length,
                           int &repeat_start, int&repeat_end, int &alignment_length, string &cigar_string, double &purity,
                           int &substitutions, int &indels, double &avg_motifpurity, int &avg_motifindels, int &avg_matchlen);

void trimCigarMotifPurity(string &cigar, int &motif_length, int & repeat_start, int &repeat_end, int &alignment_length,
                          double &purity, double &motifwise_purity);


void trimCigarMotifPurity(string &cigar, int &motif_length, int &repeat_start, int &repeat_end, int &alignment_length, double &purity);
