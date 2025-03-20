#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_int.hpp>
#include <mutex>
#include <time.h>

#ifndef GLOABL_VARIABLES_H
#define GLOABL_VARIABLES_H

using namespace std;

extern uint32_t **REPEAT_CLASSES;
extern int NUM_MOTIFS;

extern int *MOTIF_UNITS;
extern int *MOTIF_START;
extern int *MOTIF_END;
extern int *MOTIF_GAPS;
extern int *MOTIF_GAPSIZE;
extern uint32_t *MOTIF_NEXT;

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

extern bool LENGTH_CUTOFF_MODE;

// cutoffs for different measures of purity
extern double PURITY_THRESHOLD;
extern double MOTIFPURITY_THRESHOLD;
extern time_t START_TIME;

#endif // MATRIX_H
