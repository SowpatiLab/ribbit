#include "global_variables.h"

using namespace std;

void mergeAllLists(vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                   vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut,
                   int from_index_perfect, int from_index_substut, vector<int> &last_types,
                   vector<int> &last_indices, int seed_start);
