#include "global_variables.h"

using namespace std;
using namespace boost::multiprecision;

// Define rclasses matrix
uint32_t **REPEAT_CLASSES = nullptr;
int NUM_MOTIFS;

int *WINDOW_LENGTHS = nullptr;
int *WINDOW_THRESHOLDS = nullptr;
int *SEEDLEN_CUTOFF = nullptr;

int *MOTIF_START = nullptr;
int *MOTIF_END = nullptr;
int *MOTIF_UNITS = nullptr;
int *MOTIF_GAPS = nullptr;
int *MOTIF_GAPSIZE = nullptr;
uint32_t *MOTIF_NEXT = nullptr;

// Define rclasses matrix
int MINIMUM_MLEN = 2;
int MAXIMUM_MLEN = 100;
int SMALL_MLEN_LIMIT = 6;
int NMLENS;
int MINIMUM_SHIFT = 2;
int MAXIMUM_SHIFT = 100;
int NSHIFTS = 100;
int THREADS = 1;
std::mutex MTX;

string SEQUENCE_ID = "";
string SEQUENCE = "";

extern int SPLIT_LENGTH = 10000000;
extern int SPLIT_OVERLAP = 2000;
extern int CHUNK_START = 0;

bool CIGAROUTPUT = false;

int RANK_P = 5;
int RANK_Q = 4;
int RANK_S = 3;
int RANK_F = 2;
int RANK_C = 1;
int RANK_A = 0;
int RANK_N = -1;

unordered_map<int, int> MINIMUM_LENGTH;
unordered_map<int, int> MINIMUM_UNITS;
unordered_map<int, int> PERFECT_UNITS;

unordered_map<int, unordered_map<uint256_t, int>> ATOMICITY_MAP; // store atomicity of motifs
unordered_map<int, unordered_map<uint256_t, string>> MOTIFS_MAP; // Store motif to int32

bool LENGTH_CUTOFF_MODE = true;

// cutoffs for different measures of purity
double PURITY_THRESHOLD = 0.8;
double MOTIFPURITY_THRESHOLD = 0.8;
int    INTERRUPTIONS_THRESHOLD = 0;
time_t START_TIME = time(0);
