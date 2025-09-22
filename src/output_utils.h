#include "global_variables.h"
#include "process_cigar.h"
#include "cigar_utils.h"
#include "ssw_cpp.h"

#include <cmath>

using namespace std;
using namespace boost;

void addLocusToOutput(string &sequence_id, int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                      int atomicity, int repeat_length, int repeat_units, ofstream &out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);
