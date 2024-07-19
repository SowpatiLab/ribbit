#include <cstdint>
#include <unordered_map>
#include <time.h>

#include "global_variables.h"

using namespace std;

// Define rclasses matrix
uint32_t **REPEAT_CLASSES = nullptr;
int NUM_MOTIFS;

int *MOTIF_FREQUENCY = nullptr;
int *MOTIF_UNITS = nullptr;
int *MOTIF_START = nullptr;
int *MOTIF_END = nullptr;
int *MOTIF_GAPS = nullptr;
int *MOTIF_GAPSIZE = nullptr;
uint32_t *MOTIF_NEXT = nullptr;

// Define rclasses matrix
int MINIMUM_MLEN = 2;
int MAXIMUM_MLEN = 100;

int RANK_P = 5;
int RANK_Q = 4;
int RANK_S = 3;
//int RANK_F = 2;
int RANK_C = 1;
int RANK_A = 0;
int RANK_N = -1;

unordered_map<int, int> MINIMUM_LENGTH;
unordered_map<int, int> MINIMUM_UNITS;
unordered_map<int, int> PERFECT_UNITS;

bool LENGTH_CUTOFF_MODE = true;
bool PURITY_CUTOFF_MODE = true;

// cutoffs for different measures of purity
float PURITY_THRESHOLD = 0.85; 
int   MISMATCHES_THRESHOLD = 0;
int   INTERRUPTIONS_THRESHOLD = 0;
time_t START_TIME = time(0);