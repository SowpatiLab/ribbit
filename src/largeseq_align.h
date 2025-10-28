#include "global_variables.h"
#include "ssw_cpp.h"
#include "cigar_utils.h"

string alignLargeSequence(string &sequence, string &motif, int motif_length, StripedSmithWaterman::Aligner &aligner,
                          StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment);
