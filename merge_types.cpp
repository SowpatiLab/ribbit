#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

#include "merge_types.h"
#include "global_variables.h"

using namespace std;
void mergeAllLists(vector<tuple<int,int,int,int>> &seed_positions_perfect,
                   vector<tuple<int,int,int,int>> &seed_positions_substut,
                   vector<tuple<int,int,int,int>> &seed_positions_anchored,
                   int from_index_perfect, int from_index_substut, vector<int> &last_types,
                   vector<int> &last_indices, int seed_start) {

    vector<int> last_subperf_types;
    vector<int> last_subperf_indices;
    int perfect_start_bool = false, substut_start_bool = false;
    int perfect_index = from_index_perfect, substut_index = from_index_substut;
    int perfect_end, substut_end;
    int perfect_type, substut_type;

    if (seed_positions_perfect.size() == 0) {
        perfect_start_bool = true;
    }

    while (!(perfect_start_bool && substut_start_bool)) {
        // cout << "Loop-1: " << perfect_index << "\t" << substut_index << "\n";
        if (substut_start_bool) {
            while (perfect_index >= 0 || !perfect_start_bool) {
                perfect_end = get<1> (seed_positions_perfect[perfect_index]);
                perfect_type = get<3> (seed_positions_perfect[perfect_index]);
                if (perfect_end >= seed_start) {
                    if (perfect_type != RANK_N) {
                        last_subperf_types.push_back(RANK_P);
                        last_subperf_indices.push_back(perfect_index);
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
                        last_subperf_types.push_back(RANK_S);
                        last_subperf_indices.push_back(substut_index);
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
                    last_subperf_types.push_back(RANK_S);
                    last_subperf_indices.push_back(substut_index);
                }
                substut_index -= 1;
            }

            else if (substut_end <= perfect_end) {
                if (perfect_type != RANK_N) {
                    last_subperf_types.push_back(RANK_P);
                    last_subperf_indices.push_back(perfect_index);
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

    // cout << "Subperf lengths: " << last_subperf_indices.size() << "\t" << last_subperf_types.size() << "\n";

    int subperf_start_bool = false, anchored_start_bool = false;
    int subperf_index = last_subperf_indices.size()-1, anchored_index = seed_positions_anchored.size()-1;
    int subperf_end, anchored_end;
    int subperf_type, anchored_type; int idx;

    if (seed_positions_anchored.size() == 0) {
        for (int _=0; _<last_subperf_indices.size(); _++) last_indices.push_back(last_subperf_indices[_]);
        for (int _=0; _<last_subperf_types.size(); _++) last_types.push_back(last_subperf_types[_]);
    }
    else if (last_subperf_indices.size() == 0) {
        while (anchored_end >= 0 || !anchored_start_bool) {
            anchored_end = get<1> (seed_positions_anchored[anchored_index]);
            anchored_type = get<3> (seed_positions_anchored[anchored_index]);
            if (anchored_end >= seed_start) {
                if (anchored_type != RANK_N) {
                    last_types.push_back(RANK_A);
                    last_indices.push_back(anchored_index);
                }
                anchored_index -= 1;
            }
            if (anchored_index < 0 || anchored_end < seed_start) {
                anchored_start_bool = true; break;
            }
        }
    }
    else {
        while (!(subperf_start_bool && anchored_start_bool)) {
            // cout << "Loop-2: " << subperf_index << "\t" << anchored_index << "\n";
            if (anchored_start_bool) {
                while (subperf_index >= 0 || !subperf_start_bool) {
                    subperf_type = last_subperf_types[subperf_index];
                    idx = last_subperf_indices[subperf_index];
                    if (subperf_type == RANK_P) { subperf_end = get<1> (seed_positions_perfect[idx]); }
                    else if (subperf_type == RANK_S) { subperf_end = get<1> (seed_positions_substut[idx]); }
                    if (subperf_end >= seed_start) {
                        last_types.push_back(subperf_type);
                        last_indices.push_back(idx);
                        subperf_index -= 1;
                    }
                    if (subperf_index < 0 || subperf_end < seed_start) {
                        subperf_start_bool = true; break;
                    }
                }
            }

            else if (subperf_start_bool) {
                while (anchored_end >= 0 || !anchored_start_bool) {
                    anchored_end = get<1> (seed_positions_anchored[anchored_index]);
                    anchored_type = get<3> (seed_positions_anchored[anchored_index]);
                    if (anchored_end >= seed_start) {
                        if (anchored_type != RANK_N) {
                            last_types.push_back(RANK_A);
                            last_indices.push_back(anchored_index);
                        }
                        anchored_index -= 1;
                    }
                    if (anchored_index < 0 || anchored_end < seed_start) {
                        anchored_start_bool = true; break;
                    }
                }
            }

            else {
                subperf_type = last_subperf_types[subperf_index];
                idx = last_subperf_indices[subperf_index];

                if (subperf_type == RANK_P) { subperf_end = get<1> (seed_positions_perfect[idx]); }
                else if (subperf_type == RANK_S) { subperf_end = get<1> (seed_positions_substut[idx]); }                    
                anchored_end = get<1> (seed_positions_anchored[anchored_index]);

                if (anchored_end > subperf_end) {
                    last_types.push_back(RANK_A); last_indices.push_back(anchored_index);
                    anchored_index -= 1;
                }

                else if (anchored_end <= subperf_end) {
                    last_types.push_back(subperf_type);
                    last_indices.push_back(idx);
                    subperf_index -= 1;
                }

                if (subperf_index < 0 || subperf_end < seed_start) {
                    subperf_start_bool = true;
                }

                if (anchored_index < 0 || anchored_end < seed_start) {
                    anchored_start_bool = true;
                }
            }
        }
    }
}



