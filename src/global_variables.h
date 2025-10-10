#include <vector>
#include <set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/multiprecision/cpp_int.hpp>
#include <time.h>
#include <mutex>
#include <cmath>
#include <zlib.h>
#include <numeric>
#include <iomanip>

using namespace boost::multiprecision;

#ifndef GLOBAL_VARIABLES_H
#define GLOBAL_VARIABLES_H

using namespace std;

extern uint32_t **REPEAT_CLASSES;
extern int NUM_MOTIFS;

extern int *WINDOW_LENGTHS;
extern int *WINDOW_THRESHOLDS;
extern int *SEEDLEN_CUTOFF;

extern int *MOTIF_START;
extern int *MOTIF_END;
extern int *MOTIF_UNITS;
extern int *MOTIF_GAPS;
extern int *MOTIF_GAPSIZE;
extern uint32_t *MOTIF_NEXT;

extern string SEQUENCE_ID;
extern string SEQUENCE;

extern int SPLIT_LENGTH;
extern int SPLIT_OVERLAP;

extern bool CIGAROUTPUT;

extern int RANK_P;
extern int RANK_Q;
extern int RANK_S;
extern int RANK_F;
extern int RANK_C;
extern int RANK_A;
extern int RANK_N;

extern int MINIMUM_MLEN;
extern int MAXIMUM_MLEN;
extern int SMALL_MLEN_LIMIT;
extern int NMLENS;
extern int MINIMUM_SHIFT;
extern int MAXIMUM_SHIFT;
extern int NSHIFTS;
extern int THREADS;
extern std::mutex MTX;

extern unordered_map<int, int> MINIMUM_LENGTH;
extern unordered_map<int, int> MINIMUM_UNITS;
extern unordered_map<int, int> PERFECT_UNITS;

extern unordered_map<int, unordered_map<uint256_t, int>> ATOMICITY_MAP; // store atomicity of motifs
extern unordered_map<int, unordered_map<uint256_t, string>> MOTIFS_MAP; // Store motif to int32

extern bool LENGTH_CUTOFF_MODE;

// cutoffs for different measures of purity
extern double PURITY_THRESHOLD;
extern double MOTIFPURITY_THRESHOLD;
extern time_t START_TIME;

#endif // MATRIX_H
