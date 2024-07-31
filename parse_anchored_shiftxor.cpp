/*
 * Different methods for parsing shift XOR and identification of tandem repeats
*/


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

#include "global_variables.h"
#include "parse_anchored_shiftxor.h"
#include "parse_seed.h"
#include "merge_types.h"

using namespace std;

void generateAnchoredShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                               vector<boost::dynamic_bitset<>> &lsxor_anchor_bsets, int anchor_size) {

    /*
     *  generates shift XOR bitsets only retaining the anchors
     *  @param lshift_xor_bsets vector of left shift XOR bitsets of all shifts
     *  @param N_bset bitset with information of N positions
     *  @param lsxor_anchor_bsets vector of the left shift anchor bitsets
     *  @param anchor_size the length of the anchor size
    */

    int bset_size = N_bset.size();
    int anchor_start = -1;
    int motif_length;
    for (int lsxor_idx=0; lsxor_idx < NSHIFTS; lsxor_idx++) {
        motif_length = MINIMUM_SHIFT + lsxor_idx;
        boost::dynamic_bitset<> anchor_bset(bset_size, 0ull);
        for (int xor_idx = bset_size-1; xor_idx >= lsxor_idx + MINIMUM_SHIFT; xor_idx--) {

            if (lshift_xor_bsets[lsxor_idx][xor_idx] == 1) {
                if (anchor_start == -1) anchor_start = xor_idx;
            }

            else {
                if (anchor_start - xor_idx >= anchor_size && anchor_start - xor_idx < 2*motif_length) {
                    // the anchors to be retained should be at least of the minimum anchor size mentioned
                    // and not more than twice of the motif length it being tagged in because this will retain all the
                    // perfect repeats of that motif length
                    anchor_bset.set(xor_idx+1, anchor_start - xor_idx, 1);
                } anchor_start = -1;
            }
        }

        lsxor_anchor_bsets.push_back(anchor_bset);
        anchor_start = -1;
    }
}


bool retainNestedSeedAnchored(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                              int nested_midx, int parent_midx, int bset_size) {

    int nested_count = 0, parent_count = 0;
    for(int i=start; i<end; i++) {
        if (motif_bsets[nested_midx][bset_size - 1 - i] == 1) nested_count += 1;
        if (motif_bsets[parent_midx][bset_size - 1 - i] == 1) parent_count += 1;
    }

    if (nested_count < parent_count) { return false; }
    else { return true; }
}

bool retainIdeniticalSeedAnchored(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                              int nested_midx, int parent_midx, int bset_size) {

    int nested_count = 0, parent_count = 0;
    for(int i=start; i<end; i++) {
        if (motif_bsets[nested_midx][bset_size - 1 - i] == 1) nested_count += 1;
        if (motif_bsets[parent_midx][bset_size - 1 - i] == 1) parent_count += 1;
    }

    if (nested_count < parent_count) { return false; }
    else if (nested_count == parent_count) { return nested_midx < parent_midx; }
    else { return true; }
}


int clearSeedsAnchored(int from_index, vector<tuple<int, int, int, int>> &seed_positions, int seed_start) {

    vector<int> remove_seeds;
    int last_end, last_type;

    for (int i=from_index; i>=0; i--) {
        // starting from the last seed and decrementing in indices
        last_end   = get<1> (seed_positions[i]);
        last_type  = get<3> (seed_positions[i]);

        if (last_type == RANK_N) { remove_seeds.push_back(i); }
        if (last_end < seed_start) { break; }
    }

    from_index -= remove_seeds.size();
    for (int i=0; i<remove_seeds.size(); i++) {
        // removing the redundant seeds
        // because the seeds are being removed in reverse order the index of the next
        // seed to be removed is not changed
        seed_positions.erase(seed_positions.begin() + remove_seeds[i]);
    }

    return from_index;
}


tuple<int,int> addSeedToSeedPositionsAnchored(int seed_start, int seed_end, int motif_length, vector<tuple<int, int, int, int>> &seed_positions_perfect,
                                              vector<tuple<int, int, int, int>> &seed_positions_substut, vector<tuple<int, int, int, int>> &seed_positions_anchored,
                                              int* seedlen_cutoffs, vector<boost::dynamic_bitset<>> &motif_bsets, int bset_size, tuple<int,int> from_indices, int seed_type) {
    /*
     *  add seed to the existing seed positions list
     *  @param seed_start start position of the seed
     *  @param seed_end end position of the seed
     *  @param motif_length motif length of the TR seed
     *  @param min_shift the minimu size of the shift
     *  @param seed_positions vector of existing seed positions
     *  @param motif_bsets shift XOR bitsets of all the motif sizes
     *  @param bset_size the total size of a shift XOR bitset
     *  @return none add the seed to seed_position
    */

    int last_start, last_end, last_rend, last_mlen;
    int last_length, last_rlen, last_type;       // coordinate variables for existing seeds

    int from_index_perfect = get<0> (from_indices);
    int from_index_substut = get<1> (from_indices);
    for (int i=from_index_perfect; i<seed_positions_perfect.size(); i++) {
        last_start   = get<0> (seed_positions_perfect[i]);

        // go to the point where the start of the last seed is beyond the current seed
        // this logic will leave us with the last seed that is alteast overlapping at least by 1 base at the end
        if (last_start > seed_end) { break; }
        else if (from_index_perfect == seed_positions_perfect.size() - 1) { break; }
        else { from_index_perfect += 1;  }
    }

    for (int i=from_index_substut; i<seed_positions_substut.size(); i++) {
        last_start   = get<0> (seed_positions_substut[i]);

        // go to the point where the start of the last seed is beyond the current seed
        // this logic will leave us with the last seed that is alteast overlapping at least by 1 base at the end
        if (last_start > seed_end) { break; }
        else if (from_index_substut == seed_positions_substut.size() - 1) { break; }
        else { from_index_substut += 1;  }
    }

    if (seed_end-seed_start < seedlen_cutoffs[motif_length-MINIMUM_MLEN]) { return tuple<int,int>{from_index_perfect, from_index_substut}; }

    vector<int> last_types, last_indices;
    mergeAllLists(seed_positions_perfect, seed_positions_substut, seed_positions_anchored,
                  from_index_perfect, from_index_substut, last_types, last_indices, seed_start);   

    int seed_rend   = seed_end + motif_length;
    int seed_length = seed_end - seed_start;
    int seed_rlen   = seed_length + motif_length;

    // indices for different shifts in motif_bsets
    int  seed_midx = motif_length - MINIMUM_SHIFT;
    int  last_midx = 0;
    int merge_start = 0, merge_end = 0, overlap_length = 0;

    vector<int> identical, nestedin, overlap;
    vector<int> parentof_subperf_factor, parentof_subperf_multiple, parentof_subperf_nonfactor;
    vector<int> parentof_subperf_factorsizes, parentof_subperf_nonfactorsizes;
    vector<int> parentof_subperf_factortypes, parentof_subperf_multipletypes, parentof_subperf_nonfactortypes;
    vector<int> parentof_anchored_factor, parentof_anchored_nonfactor;
    tuple<int,int> from_indices_new = {0, 0};

    for (int _=0; _<last_indices.size(); _++) {

        // starting from the last seed and decrementing in indices
        int i = last_indices[_];
        if (last_types[_] == RANK_P) {
            last_start = get<0> (seed_positions_perfect[i]);
            last_mlen  = get<2> (seed_positions_perfect[i]);
            last_end   = get<1> (seed_positions_perfect[i]);
            last_rend  = get<1> (seed_positions_perfect[i]) + last_mlen;
            last_type  = get<3> (seed_positions_perfect[i]);
        }
        else if (last_types[_] == RANK_S) {
            last_start = get<0> (seed_positions_substut[i]);
            last_mlen  = get<2> (seed_positions_substut[i]);
            last_end   = get<1> (seed_positions_substut[i]);
            last_rend  = get<1> (seed_positions_substut[i]) + last_mlen;
            last_type  = get<3> (seed_positions_substut[i]);
        }
        else if (last_types[_] == RANK_A) {
            last_start = get<0> (seed_positions_anchored[i]);
            last_mlen  = get<2> (seed_positions_anchored[i]);
            last_end   = get<1> (seed_positions_anchored[i]);
            last_rend  = get<1> (seed_positions_anchored[i]) + last_mlen;
            last_type  = get<3> (seed_positions_anchored[i]);
        }

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) { break; }

        if (last_type == RANK_N) { continue; }

        // if the from_index is much ahead we skip the seeds that do not overlap
        if (seed_end < last_start) { continue; }        

        last_length = last_end - last_start;
        last_rlen  = last_rend - last_start;
        last_midx  = last_mlen - MINIMUM_SHIFT;

        // current seed and last seed have identical coordinates
        if (seed_start == last_start && seed_end == last_end) {
            //  if new seed is anchored and old seed is of higher rank ~ do not add new seed
            if (seed_type == RANK_A && last_type > RANK_A) {
                return tuple<int,int> { from_index_perfect, from_index_substut };
            }

            // if new seed is merged anchored and old seed is anchored ~ remove old seed
            else if (seed_type == RANK_C && last_type == RANK_A ) {
                seed_positions_anchored[i] = tuple<int, int, int, int> { last_start, last_end, last_mlen, RANK_N };
            }

            // if both seeds are of the same type
            else { identical.push_back(i); }
        }

        // current seed is nested in an existing seed
        else if (last_start <= seed_start && seed_end <= last_end) {

            //  if new (nested) seed type is of lesser rank than old (parent) seed
            if (last_type > seed_type) {
                return tuple<int,int> { from_index_perfect, from_index_substut };
            }
            else if (seed_type == RANK_C && last_type == RANK_A) { }

            else if ((seed_type == RANK_A && last_type == RANK_A) || (seed_type == RANK_C && last_type == RANK_C)) {
                // if new (nested) seed's motif length is a multiple of old (parent) seed's motif length ~ do not add new seed
                if (motif_length % last_mlen == 0 && (motif_length != 4)) {
                    return tuple<int,int> { from_index_perfect, from_index_substut };
                }

                // if new (nested) seed's motif length is a factor of old (parent) seed's motif length
                else if (last_mlen % motif_length == 0 && (last_mlen != 4)) {
                    if (seed_rlen >= last_mlen - 1 || seed_rlen >= last_length) {
                        seed_positions_anchored[i] = tuple<int, int, int, int> { last_start, last_end, last_mlen, RANK_N };
                        from_indices_new = addSeedToSeedPositionsAnchored(last_start, last_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                       seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                        return from_indices_new;
                    }
                    else { nestedin.push_back(i); continue; }
                }

                else {
                    bool retain = retainNestedSeedAnchored(motif_bsets, seed_start, seed_end, seed_midx, last_midx, bset_size);
                    if (!retain) { return tuple<int,int>{from_index_perfect, from_index_substut}; }
                    else { nestedin.push_back(i); continue; }
                }
            }
        }

        // current seed is parent in an existing seed
        else if (seed_start <= last_start && last_end <= seed_end) {

            if (last_type > seed_type) {
                if (motif_length % last_mlen == 0) {
                    if ((last_rlen >= motif_length - 2) || (last_rlen >= seed_length - 2)) {
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else if (last_type == RANK_S || last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_indices_new = addSeedToSeedPositionsAnchored(seed_start, seed_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                       seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, RANK_C);
                        return from_indices_new;
                        // return tuple<int,int>{ from_index_perfect, from_index_substut };
                    }

                    else {
                        parentof_subperf_factor.push_back(i);
                        parentof_subperf_factorsizes.push_back(last_mlen);
                        parentof_subperf_factortypes.push_back(last_type);
                    }
                }

                else if (last_mlen % motif_length == 0) {
                    if (last_mlen >= 4*motif_length || last_length >= 4*motif_length) {
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else if (last_type == RANK_S || last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_indices_new = addSeedToSeedPositionsAnchored(seed_start, seed_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                       seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, RANK_C);
                        return from_indices_new;
                        // return tuple<int,int>{ from_index_perfect, from_index_substut };
                    }

                    else {
                        parentof_subperf_multiple.push_back(i);
                        parentof_subperf_multipletypes.push_back(last_type);
                    }
                }

                else if (last_mlen > motif_length) {
                    if (last_mlen >= 4*motif_length || last_length >= 4*motif_length) {
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else if (last_type == RANK_S || last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_indices_new = addSeedToSeedPositionsAnchored(seed_start, seed_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                       seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, RANK_C);
                        return from_indices_new;
                        // return tuple<int,int>{ from_index_perfect, from_index_substut };
                    }
                }

                else {
                    parentof_subperf_nonfactor.push_back(i);
                    parentof_subperf_nonfactorsizes.push_back(last_mlen);
                    parentof_subperf_nonfactortypes.push_back(last_type);
                }
            }

            else if (seed_type == RANK_C && last_type == RANK_A ) {
                seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
            }

            else if ((seed_type == RANK_A && last_type == RANK_A) || (seed_type == RANK_C && last_type == RANK_C)) {
                if (last_mlen == motif_length) {
                    seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                }

                else {
                    bool retain = retainNestedSeedAnchored(motif_bsets, last_start, last_end, last_midx, seed_midx, bset_size);
                    if (!retain) { seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    else {
                        if (motif_length % last_mlen == 0) {
                            if ((last_rlen >= motif_length - 2) || (last_rlen >= seed_length - 2)) {
                                seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                                from_indices_new = addSeedToSeedPositionsAnchored(seed_start, seed_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                               seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                                return from_indices_new;                                
                            }
                            else { parentof_anchored_factor.push_back(i); }
                        }
                        else if (last_mlen % motif_length == 0) {
                            continue;
                        }
                        else { parentof_anchored_nonfactor.push_back(i); }
                    }
                }
            }
        }

        // overlap
        else {
            if (last_start < seed_start) {
                if (last_mlen <= motif_length) {
                    if (seed_end <= last_rend) overlap_length = seed_end - seed_start;
                    else overlap_length = last_rend - seed_start;
                }
                else {
                    if (seed_end <= last_end) overlap_length = seed_end - seed_start;
                    else overlap_length = last_end - seed_start;
                }
                merge_start = last_start; merge_end = seed_end;
            }

            else {
                if (motif_length <= last_mlen) {
                    if (last_end <= seed_rend) overlap_length = last_end - last_start;
                    else overlap_length = seed_rend - last_start;
                }
                else {
                    if (last_end <= seed_end) overlap_length = last_end - last_start;
                    else overlap_length = seed_end - last_start;
                }
                merge_start = seed_start; merge_end = last_end;
            }

            if (seed_type == RANK_A && last_type > RANK_C) {

                if (motif_length == last_mlen) {
                    if (overlap_length >= 4*motif_length) {
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else if (last_type == RANK_S || last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_indices_new = addSeedToSeedPositionsAnchored(merge_start, merge_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                       seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, RANK_C);
                        return from_indices_new;
                        // return tuple<int,int>{ from_index_perfect, from_index_substut };
                    }
                }

                if ((motif_length % last_mlen == 0) || (last_mlen % motif_length == 0)) { }

                else {
                    if ((overlap_length >= motif_length-1) || (overlap_length >= seed_length-1)) {
                        return tuple<int,int>{from_index_perfect, from_index_substut};
                    }
                }
            }

            else if ((seed_type == RANK_A && last_type == RANK_A) || (seed_type == RANK_C && last_type == RANK_C) || (seed_type == RANK_A && last_type == RANK_C) || (seed_type == RANK_C && last_type == RANK_A)) {
                if (motif_length == last_mlen) {
                    if (last_length >= seed_length) {
                        if ((seed_length >= 3*motif_length) && ((overlap_length >= 3*motif_length-1) || (overlap_length >= seed_length-1))) {
                            seed_type == (seed_type == RANK_C || last_type == RANK_C) ? RANK_C : RANK_A;
                            seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                            from_indices_new = addSeedToSeedPositionsAnchored(merge_start, merge_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                           seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                            return from_indices_new;
                            // return tuple<int,int>{from_index_perfect, from_index_substut};
                        }
                        else if ((seed_length < 3*motif_length) && ((overlap_length >= motif_length-1) || (overlap_length >= seed_length-1))) {
                            seed_type == (seed_type == RANK_C || last_type == RANK_C) ? RANK_C : RANK_A;
                            seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                            from_indices_new = addSeedToSeedPositionsAnchored(merge_start, merge_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                           seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                            return from_indices_new;
                            // return tuple<int,int>{from_index_perfect, from_index_substut};
                        }
                    }
                    else {
                        if ((last_length >= 3*last_mlen) && ((overlap_length >= 3*last_mlen-1) || (overlap_length >= last_length-1))) {
                            seed_type == (seed_type == RANK_C || last_type == RANK_C) ? RANK_C : RANK_A;
                            seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                            from_indices_new = addSeedToSeedPositionsAnchored(merge_start, merge_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                           seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                            return from_indices_new;
                            // return tuple<int,int>{from_index_perfect, from_index_substut};
                        }
                        else if ((seed_length < 3*last_mlen) && ((overlap_length >= last_mlen-1) || (overlap_length >= last_length-1))) {
                            seed_type == (seed_type == RANK_C || last_type == RANK_C) ? RANK_C : RANK_A;
                            seed_positions_anchored[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                            from_indices_new = addSeedToSeedPositionsAnchored(merge_start, merge_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                           seed_positions_anchored, seedlen_cutoffs, motif_bsets, bset_size, from_indices, seed_type);
                            return from_indices_new;
                            // return tuple<int,int>{from_index_perfect, from_index_substut};
                        }
                    }
                }
            }
        }
    }

    int nonfactor_coverage = 0, factor_coverage = 0, multiple_coverage = 0;
    uint32_t prev_start = -1;

    if (parentof_subperf_nonfactor.size() > 0) {
        for (int j=0; j < parentof_subperf_nonfactor.size(); j++) {
            int k = parentof_subperf_nonfactor[j];
            int ktype = parentof_subperf_nonfactortypes[j];
            if (ktype == RANK_P) {
                last_start = get<0> (seed_positions_perfect[j]);
                last_mlen  = get<2> (seed_positions_perfect[j]);
                last_end   = get<1> (seed_positions_perfect[j]);
                last_rend  = get<1> (seed_positions_perfect[j]) + last_mlen;
            }
            else if (ktype == RANK_S) {
                last_start = get<0> (seed_positions_substut[j]);
                last_mlen  = get<2> (seed_positions_substut[j]);
                last_end   = get<1> (seed_positions_substut[j]);
                last_rend  = get<1> (seed_positions_substut[j]) + last_mlen;
            }

            if (last_rend >= prev_start) { nonfactor_coverage += prev_start - last_start; }
            else if (last_rend < seed_end) { nonfactor_coverage += last_rend - last_start; }
            else { nonfactor_coverage += seed_end - last_start;}
            prev_start = last_start;
        }

        if (nonfactor_coverage > 0.5*seed_length) { return tuple<int,int>{ from_index_perfect, from_index_substut }; }
    }


    if (parentof_subperf_factor.size() > 0) {
        prev_start = -1;
        unordered_map<int,int> prev_starts;
        unordered_map<int,int> factor_coverages;
        for(int factorsize: parentof_subperf_factorsizes) {
            prev_starts[factorsize] = -1;
            factor_coverages[factorsize] = 0;
        }

        for (int j=0; j < parentof_subperf_factor.size(); j++) {
            int k = parentof_subperf_factor[j];
            int ktype = parentof_subperf_factortypes[j];
            if (ktype == RANK_P) {
                last_start = get<0> (seed_positions_perfect[j]);
                last_mlen  = get<2> (seed_positions_perfect[j]);
                last_end   = get<1> (seed_positions_perfect[j]);
                last_rend  = get<1> (seed_positions_perfect[j]) + last_mlen;
            }

            else if (ktype == RANK_S) {
                last_start = get<0> (seed_positions_substut[j]);
                last_mlen  = get<2> (seed_positions_substut[j]);
                last_end   = get<1> (seed_positions_substut[j]);
                last_rend  = get<1> (seed_positions_substut[j]) + last_mlen;
            }
            prev_start = prev_starts[last_mlen];

            if (last_rend >= prev_start) { factor_coverages[last_mlen] += prev_start - last_start; }
            else if (last_rend < seed_end) { factor_coverages[last_mlen] += last_rend - last_start; }
            else { factor_coverages[last_mlen] += seed_end - last_start;}
            prev_starts[last_mlen] = last_start;
        }

        vector<int> factors;
        for(auto it: factor_coverages) { factors.push_back(it.first); }
        std::sort(factors.begin(), factors.end());
        for (int factor: factors) {
            if (factor_coverages[factor] >= 0.8*seed_length) {
                motif_length = factor; seed_type = RANK_C;

                for (int j=0; j < parentof_subperf_factor.size(); j++) {
                    int k = parentof_subperf_factor[j];
                    int ktype = parentof_subperf_factortypes[j];
                    if (ktype == RANK_P) {
                        last_mlen  = get<2> (seed_positions_perfect[j]);
                        if (last_mlen==factor) { seed_positions_perfect[j] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    }                
                    else if (ktype == RANK_S) {
                        last_mlen  = get<2> (seed_positions_substut[j]);
                        if (last_mlen==factor) { seed_positions_substut[j] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    }
                }
                break;
            }
        }
    }

    // limiting the seeds to the edge
    if (seed_end > bset_size-motif_length) {
        seed_end = bset_size-motif_length;
    }
    seed_positions_anchored.push_back({seed_start, seed_end, motif_length, seed_type});
    return tuple<int,int> { from_index_perfect, from_index_substut };
}



vector<tuple<int,int,int,int>> processShiftXORsAnchored(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                        int &window_length, int &window_bitcount_threshold, vector<tuple<int, int, int, int>> &seed_positions_perfect,
                                                        vector<tuple<int, int, int, int>> &seed_positions_substut) {
    /*
     *  parsing the shift XORs of all shift sizes and picking seeds from each shift
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param N_bset N position bitset
     *  @param window_length length of the window to be scanned
     *  @param window_bitcount_threshold the threshold number of set bits in the window
     *  @param nshift number of shift XORs built
     *  @param min_shift motif size of the minimum shift generated
     *  @return vector<tuple<int, int, int>> vector of end position sorted seeds from all motif sizes
    */

    int bset_size = N_bset.size();          // size of the sequence
    int window_bitcount;        // stores window bitcount
    int valid_position = 0;     // position tracking valid bits in the window

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT, didx, motif_length;
    int last_starts[NMOTIFS];  // initialising a last record
    int last_ends[NMOTIFS];  // initialising a last record
    int current_starts[NMOTIFS];  // initialising a last record
    int seedlen_cutoffs[NMOTIFS];

    // initialising all to -1
    for (int _=0; _<NMOTIFS; _++) { last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1; seedlen_cutoffs[_] = 10;}

    tuple<int,int> from_indices = {0, 0};
    vector<tuple<int,int,int,int>> seed_positions_anchored;

    vector<boost::dynamic_bitset<>> window_bsets;
    for (int midx=0; midx < NMOTIFS; midx++) {
        boost::dynamic_bitset<> window_bset(window_length, 0ull);
        window_bsets.push_back(window_bset);   // initialised window bitset
        seedlen_cutoffs[midx] = ((midx+MINIMUM_MLEN) > 6) ? (midx+MINIMUM_MLEN) : 10;
        if (midx+MINIMUM_MLEN >= 10) { seedlen_cutoffs[midx] = 0.9 * (midx+MINIMUM_MLEN); }
    }

    int overlap_distance = 0;   // the allowed overlap distance between adjacent seeds

    int xor_idx = 0;
    int window_position = -1*window_length;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {
        window_position += 1;

        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < NMOTIFS+min_idx; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                if (current_starts[didx] != -1) {
                    // No seed is being tracked currently

                    // start position is the first nuc of the window; it is zero based
                    current_starts[didx] = window_position;

                    // if the last stored seed is beyond the overlapping distance
                    if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                        from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                      seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                      motif_bsets, bset_size, from_indices, RANK_A);
                        last_starts[didx] = -1; last_ends[didx] = -1;
                    }
                }
                window_bsets[didx] <<= window_length;
                current_starts[didx] = -1;
            }

            valid_position = 0;
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            valid_position += 1;

            for (int midx=min_idx; midx < NMOTIFS+min_idx; midx++) {
                didx = midx-min_idx;
                window_bsets[didx] <<= 1;
                window_bsets[didx][0] = motif_bsets[midx][xor_idx];
            }

            if (valid_position >= window_length) {
                for (int midx=min_idx; midx < NMOTIFS+min_idx; midx++) {
                    didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                    window_bitcount = window_bsets[didx].count();

                    if (window_bitcount >= window_bitcount_threshold) {
                        // window bitcount is above the threshold

                        if (current_starts[didx] == -1) {
                            // No seed is being tracked currently

                            // start position is the first nuc of the window; it is zero based
                            current_starts[didx] = window_position;

                            // if the last stored seed is beyond the overlapping distance
                            if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                                from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                              seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                              motif_bsets, bset_size, from_indices, RANK_A);
                                last_starts[didx] = -1; last_ends[didx] = -1;
                            }
                        }
                    }

                    else {
                        // window bitcount is not above the threshold

                        if (current_starts[didx] != -1) {
                            // if a seed is being tracked currently

                            if (last_starts[didx] == -1) {
                                // last seed is not recorded; save the current seed as the last seed
                                last_starts[didx] = current_starts[didx];
                                last_ends[didx] = window_position + window_length - 1; // end is exclusive
                            }

                            else {
                                // if the last seed is recorded it means that it is within the overlapping range
                                // hence we just update the end of the last record
                                last_ends[didx] = window_position + window_length - 1; // reassign end
                            }

                            current_starts[didx] = -1;
                        }

                        else {
                            // if there is no seed currently being tracked

                            // if the last stored seed is beyond the overlapping distance of current position
                            if (last_ends[didx] != -1 && last_ends[didx] < window_position - overlap_distance) {
                                from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                              seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                              motif_bsets, bset_size, from_indices, RANK_A);
                                // the last seed is reset
                                last_starts[didx] = -1; last_ends[didx] = -1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int midx=min_idx; midx < NMOTIFS+min_idx; midx++) {
        didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;

        // handling the records after the end of the sequence
        if (last_ends[didx] == -1) {
            if (current_starts[didx] != -1) {
                // presently not scanning through a passed window ~ save last record
                addSeedToSeedPositionsAnchored(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length, seed_positions_perfect,
                                                              seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                              motif_bsets, bset_size, from_indices, RANK_A);
            }
        }

        else {
            if (current_starts[didx] == -1) {
                // presently not scanning through a passed window ~ save last record
                addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                              seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                              motif_bsets, bset_size, from_indices, RANK_A);
            }

            else {
                if (last_ends[didx] >= current_starts[didx] - motif_length) {
                    // current passed window overlaps with last record ~ merge both and save
                    last_ends[didx] = bset_size - (xor_idx + 1); // reassign end
                    addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                  seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                  motif_bsets, bset_size, from_indices, RANK_A);
                }

                else {
                    // current passed window doesn't overlap with last record ~ save both separately
                    from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                  seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                  motif_bsets, bset_size, from_indices, RANK_A);

                    addSeedToSeedPositionsAnchored(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length, seed_positions_perfect,
                                                                  seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                  motif_bsets, bset_size, from_indices, RANK_A);
                }
            }
        }
    }

    return seed_positions_anchored;
}
