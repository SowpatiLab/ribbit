#include "global_variables.h"

using namespace std;

unordered_map<int, pair<int, int>> getWindowThresholds(int minimum_mlen, int maximum_mlen);
long double probWithRunApprox(int n, int x, int r, long double p);
int minimumNumberOfSuccesses(int n, int r, long double p, unordered_map<int, int> &threshold_bits);
void calculateWindowThresholds();