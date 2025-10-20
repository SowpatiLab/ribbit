#include "merge_types.h"

using namespace std;


void mergeAllLists(vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                   vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut,
                   int from_index_perfect, int from_index_substut, vector<int> &last_types,
                   vector<int> &last_indices, int seed_start) {
    /*
     *  merges the list of seeds of across all types
     *  @param seed_positions_perfect vector of perfect seed positions
     *  @param seed_positions_substut vector of repeat positions with mismatches allowed
     *  @param seed_positions_anchored vector of repeat positions with indels allowed
     *  @param from_index_perfect index of perfect seed to resume merging from
     *  @param from_index_substut index of substitute seeds to resume merging from
     *  @param last_types index of substitute seeds to resume merging from
     *  @param last_indices index of substitute seeds to resume merging from
     *  @param seed_start start of the seed
     *  @returns void
     */

    int perfect_start_bool = false, substut_start_bool = false;
    int perfect_index = from_index_perfect, substut_index = from_index_substut;
    int perfect_end, substut_end;
    int perfect_type, substut_type;

    if (seed_positions_perfect.size() == 0) {
        perfect_start_bool = true;
    }

    if (seed_positions_substut.size() == 0) {
        substut_start_bool = true;
    }

    while (!(perfect_start_bool && substut_start_bool)) {
        if (substut_start_bool) {
            while (perfect_index >= 0 || !perfect_start_bool) {
                perfect_end = get<1> (seed_positions_perfect[perfect_index]);
                perfect_type = get<3> (seed_positions_perfect[perfect_index]);
                if (perfect_end >= seed_start) {
                    if (perfect_type != RANK_N) {
                        last_types.push_back(RANK_P);
                        last_indices.push_back(perfect_index);
                    }
                    perfect_index -= 1;
                }
                if (perfect_index < 0 || perfect_end < seed_start) {
                    perfect_start_bool = true; break;
                }
            }
        }

        else if (perfect_start_bool) {
            while (substut_end >= 0 || !substut_start_bool) {
                substut_end = get<1> (seed_positions_substut[substut_index]);
                substut_type = get<3> (seed_positions_substut[substut_index]);
                if (substut_end >= seed_start) {
                    if (substut_type != RANK_N) {
                        last_types.push_back(RANK_S);
                        last_indices.push_back(substut_index);
                    }
                    substut_index -= 1;
                }
                if (substut_index < 0 || substut_end < seed_start) {
                    substut_start_bool = true; break;
                }
            }
        }

        else {
            perfect_end = get<1> (seed_positions_perfect[perfect_index]);
            substut_end = get<1> (seed_positions_substut[substut_index]);
            perfect_type = get<3> (seed_positions_perfect[perfect_index]);
            substut_type = get<3> (seed_positions_substut[substut_index]);

            if (substut_end > perfect_end) {
                if (substut_type != RANK_N) {
                    last_types.push_back(RANK_S);
                    last_indices.push_back(substut_index);
                }
                substut_index -= 1;
            }

            else if (substut_end <= perfect_end) {
                if (perfect_type != RANK_N) {
                    last_types.push_back(RANK_P);
                    last_indices.push_back(perfect_index);
                }
                perfect_index -= 1;
            }

            if (perfect_index < 0 || perfect_end < seed_start) {
                perfect_start_bool = true;
            }

            if (substut_index < 0 || substut_end < seed_start) {
                substut_start_bool = true;
            }
        }
    }
}
