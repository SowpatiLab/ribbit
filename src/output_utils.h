#include "global_variables.h"
#include "process_cigar.h"
#include "cigar_utils.h"
#include "ssw_cpp.h"

#include <cmath>

using namespace std;
using namespace boost;

void addLocusToOutput(string sequence_id, int repeat_start, int repeat_end, string motif, double purity, string cigar_string,
                      int motif_length, int repeat_length, int repeat_units, ostream* out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                      int &recursion_level, vector<tuple<string, int, int, string, double, string, int, int, int>> &new_repeat_loci);

void printRepeatsToOutput(ostream* out, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                          int end_index);
