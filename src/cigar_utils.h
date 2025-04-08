#include "global_variables.h"
#include "process_cigar.h"

using namespace std;
using namespace boost;

string buildCigar(vector<int> &clens, vector<char> &ctypes);
string buildCigar(tuple<vector<int>, vector<char>> &cigar_values);

tuple<vector<int>, vector<char>> extractDownCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end);
tuple<vector<int>, vector<char>> extractUpCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end);
tuple<vector<int>, vector<char>> extractRegionCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end);

void extendTillMatch(string &non_olcigar, string &cigar, int &position, bool end);
void extendTillMatch(tuple<vector<int>, vector<char>> &non_olcigar, string &cigar, int &position, bool end);

void separateCigarsOverlappingLoci(int upstart, int upend, tuple<vector<int>, vector<char>> &up_splitcigar_values,
                                   int dnstart, int dnend, tuple<vector<int>, vector<char>> &dn_splitcigar_values,
                                   tuple<vector<int>, vector<char>> &up_olcigar, tuple<vector<int>, vector<char>> &dn_olcigar,
                                   tuple<vector<int>, vector<char>> &up_non_olcigar, tuple<vector<int>, vector<char>> &dn_non_olcigar);

int getAlignmentLength(string &cigar);
int getAlignmentLength(vector<int> &clens, vector<char> &ctypes);
int getAlignmentLength(tuple<vector<int>, vector<char>> &cigar_values);

int getAlignmentLengthWithDeletions(string &cigar);

void getMatches(string &cigar, int &matches, int &longest_match);
void getMatches(tuple<vector<int>, vector<char>> &cigar, int &matches, int &longest_match);
void getMatches(vector<int> &clens, vector<char> &ctypes, int &matches, int &longest_match);
int getMatches(string &cigar);
int getMatches(tuple<vector<int>, vector<char>> &cigar);
int getMatches(vector<int> &clens, vector<char> &ctypes);

