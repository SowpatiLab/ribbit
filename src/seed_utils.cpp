#include "seed_utils.h"

using namespace std;
using namespace boost;

void filterPerfectSeeds(vector<tuple<int, int, int, int>> &seed_positions_perfect,
                        vector<tuple<int, int, int, int>> &seed_positions_substut) {
    /*
     *  filter the perfect seeds from the substitute seeds
     *  @param seed_positions_perfect vector of identified perfect repeat seeds
     *  @param seed_positions_substut vector of identified repeat seeds with allowed substitutions
     *  @return none
     */

    int start_index = 0;
    int perfect_start, perfect_end, perfect_mlen, perfect_type;
    int substut_start, substut_end, substut_mlen, substut_type;
    for (int i = 0; i < seed_positions_perfect.size(); i++) {
        // for each perfect seed check if it is contained within any of the substitute seeds
        perfect_start = get<0>(seed_positions_perfect[i]);
        perfect_end  = get<1>(seed_positions_perfect[i]);
        perfect_mlen = get<2>(seed_positions_perfect[i]);
        perfect_type = get<3>(seed_positions_perfect[i]);

        for (int j = start_index; j < seed_positions_substut.size(); j++) {
            substut_start = get<0>(seed_positions_substut[j]);
            substut_end  = get<1>(seed_positions_substut[j]);
            substut_mlen = get<2>(seed_positions_substut[j]);
            substut_type = get<3>(seed_positions_substut[j]);

            // saving the index for the next perfect seed search.
            // the seed position should be before the start of the current perfect seed
            if (substut_end < perfect_start - MAXIMUM_MLEN) start_index = j;

            // ignore if the perfect seed or the substitute seed is invalid
            if (perfect_type == RANK_N || substut_type == RANK_N) continue;

            if (perfect_mlen == substut_mlen) {
                if (substut_start <= perfect_start && substut_end >= perfect_end &&
                    perfect_start - substut_start < perfect_mlen && substut_end - perfect_end < perfect_mlen) {
                    // if the perfect seed is contained within the substitute seed and
                    // the substitute seed is not longer than the perfect seed by more than 1 motif length on either side
                    // then we discard the perfect seed
                    seed_positions_perfect[i] = tuple<int, int, int, int>{perfect_start, perfect_end, perfect_mlen, RANK_N};
                }
            }

            if (substut_end > perfect_end + perfect_mlen) break;
        }
    }
}


void filterShortSeeds(vector<tuple<int, int, int, int>> &seeds) {
    /*
     *  marks seeds which are shorter than the cutoff length as invalid
     *  @param seeds vector of seed positions
     *  @return none
     */

    int seed_start, seed_end, seed_mlen, seed_type;
    tuple<int, int, int, int> seed;
    int seedlen_cutoff = SEEDLEN_CUTOFF[0];
    for (int i = 0; i < seeds.size(); i++) {
        seed = seeds[i];

        seed_type = get<3>(seed);
        if (seed_type == RANK_N) continue;

        seed_start = get<0>(seed);
        seed_end = get<1>(seed);
        seed_mlen = get<2>(seed);

        if (THREADS > 1) MTX.lock();
        seedlen_cutoff = SEEDLEN_CUTOFF[seed_mlen - MINIMUM_MLEN];
        if (THREADS > 1) MTX.unlock();

        if (seed_end - seed_start < seedlen_cutoff) {
            seeds[i] = tuple<int, int, int, int>{seed_start, seed_end, seed_mlen, RANK_N};
        }
    }
}


bool retainNestedSeed(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                      int nested_midx, int parent_midx, int bset_size) {
    /*
     *  compares the number of matches in the nested and the parent bitsets and decides to retain the nested repeat
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param start start of the nested locus
     *  @param end end of the nested locus
     *  @param nested_midx index for the shift XOR bitset of the nested repeat
     *  @param parent_midx index for the shift XOR bitset of the parent repeat
     *  @param bset_size size of the shift XOR bitset
     *  @return bool if the nested repeat should be retained or not
    */
    int nested_count = 0, parent_count = 0;
    for(int i=start; i<end; i++) {
        if (motif_bsets[nested_midx][bset_size - 1 - i] == 1) nested_count += 1;
        if (motif_bsets[parent_midx][bset_size - 1 - i] == 1) parent_count += 1;
    }

    if (nested_count < parent_count) { return false; }
    else { return true; }
}


bool retainIdenticalSeeds(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                           int nested_midx, int parent_midx, int bset_size) {
    /*
     *  compares the number of matches in both bitsets and decides which one to retain
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param start start of the nested locus
     *  @param end end of the nested locus
     *  @param nested_midx index for the shift XOR bitset of the nested repeat
     *  @param parent_midx index for the shift XOR bitset of the parent repeat
     *  @param bset_size size of the shift XOR bitset
     *  @return bool if the nested repeat should be retained or not
    */
    int nested_count = 0, parent_count = 0;
    for(int i=start; i<end; i++) {
        if (motif_bsets[nested_midx][bset_size - 1 - i] == 1) nested_count += 1;
        if (motif_bsets[parent_midx][bset_size - 1 - i] == 1) parent_count += 1;
    }

    if      (nested_count < parent_count)  { return false; }
    else if (nested_count == parent_count) { return nested_midx < parent_midx; }
    else    { return true; }
}
