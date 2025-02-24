#include <iostream>
#include <fstream>
#include <unordered_map>
#include <mutex>
#include <boost/dynamic_bitset.hpp>

#include <cstdint>
#include <numeric>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <boost/multiprecision/cpp_int.hpp>

#include "global_variables.h"

using namespace std;
using namespace boost;

void addLocusToOutput(string &sequence_id, int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                      int atomicity, int repeat_length, int repeat_units, ostream &out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);