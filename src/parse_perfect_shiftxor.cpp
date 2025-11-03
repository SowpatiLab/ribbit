#include "parse_perfect_shiftxor.h"

using namespace std;


void addPerfectRepeatPositions(int seed_start, int seed_end, int motif_length, int bset_size,
                               vector<tuple<int,int,int,int,int,int,int>> &repeat_positions) {
    /*
     *  when the purity threshold is 1; this adds a identified perfect repeat locus to set of loci
     *  @param seed_start the start coordinate of the seed
     *  @param seed_end the end coordinate of the seed
     *  @param motif_length length of the motif of the repeat being added
     *  @param repeat_positions list of identified repeats
     *  @param bset_size size of the shift XOR bitset
     *  @returns none adds the perfect repeat locus to the list of identified perfect repeat loci
    */

    int last_start, last_end, last_rend, last_mlen;       // coordinate variables for existing seeds
    int seed_length = seed_end - seed_start;
    int seed_rend   = seed_end + motif_length;
    int seed_rlen   = seed_rend - seed_start;

    // indices for different shifts in motif_bsets
    int seed_midx = motif_length - MINIMUM_SHIFT;
    int last_slen, last_rlen, overlap_length;

    vector<int> remove_seeds;   // the indices of seeds that need to be removed

    for (int i=repeat_positions.size()-1; i>=0; i--) {
        // starting from the last seed and decrementing in indices
        last_start = get<0> (repeat_positions[i]);
        last_end   = get<1> (repeat_positions[i]);
        last_mlen  = get<2> (repeat_positions[i]);
        last_rend  = last_end + last_mlen;
        last_slen  = last_end - last_start;
        last_rlen  = last_slen + last_mlen;

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) break;

        // identical
        if (last_start == seed_start && last_rend == seed_rend) {
            if (last_mlen < motif_length) { return; }
            else { remove_seeds.push_back(i); }
        }

        // nested
        else if (last_start <= seed_start && last_rend >= seed_rend) {
            if (motif_length < last_mlen) {
                if (seed_rlen >= last_mlen || seed_rlen >= last_slen) {
                    remove_seeds.push_back(i);
                    for (int _=0; _<remove_seeds.size(); _++) repeat_positions.erase(repeat_positions.begin() + remove_seeds[_]);
                    addPerfectRepeatPositions(last_start, last_end, motif_length, bset_size, repeat_positions);
                    return;
                }
            }
            else { return; }
        }

        // parent
        else if (seed_start <= last_start && seed_rend >= last_rend) {
            if (last_mlen < motif_length) {
                if (last_rlen >= motif_length || last_rlen >= seed_length) {
                    remove_seeds.push_back(i);
                    for (int _=0; _<remove_seeds.size(); _++) repeat_positions.erase(repeat_positions.begin() + remove_seeds[_]);
                    addPerfectRepeatPositions(seed_start, seed_end, last_mlen, bset_size, repeat_positions);
                    return;
                }
            }
            else {
                remove_seeds.push_back(i);
            }
        }

        // overlap
        else {
            if (last_start < seed_start) { overlap_length = last_rend - seed_start; }
            else { overlap_length = seed_rend - last_start; }

            if (last_mlen < motif_length) {
                // if the overlap length is at least 1 less than the larger motif size
                //  the longer motif repeat with more than 3 units is retained
                if (motif_length - overlap_length <= 1 && seed_rlen/motif_length < 3) {
                    return;
                }
                else if (seed_rlen - overlap_length <= last_mlen) {
                    return;
                }
            }

            else if (motif_length < last_mlen) {
                // if the overlap length is at least 1 less than the larger motif size
                if (last_mlen - overlap_length <= 1 && last_rlen/last_mlen < 3) {
                    remove_seeds.push_back(i);
                }
                else if (last_rlen - overlap_length <= motif_length) {
                    remove_seeds.push_back(i);
                }
            }
        }
    }

    for (int i=0; i<remove_seeds.size(); i++) {
        // removing the redundant seeds
        // because the seeds are being removed in reverse order the index of the next
        // seed to be removed is not changed
        repeat_positions.erase(repeat_positions.begin() + remove_seeds[i]);
    }

    // limiting the seeds to the edge
    if (seed_end > bset_size-motif_length) {
        seed_end = bset_size-motif_length;
    }

    seed_length = seed_end - seed_start;
    repeat_positions.push_back(tuple<int,int,int,int,int,int,int> { seed_start, seed_end, motif_length, RANK_P, seed_length, 
                                                                    seed_length , seed_length });
}


void addSeedToSeedPositionsPerfect(int seed_start, int seed_end, int motif_length, int bset_size,
                                   vector<tuple<int,int,int,int,int,int,int>> &seed_positions) {
    /*
     *  adding a potential perfect repeat seed to seed positions
     *  @param seed_start the start coordinate of the seed
     *  @param seed_end the end coordinate of the seed
     *  @param motif_length length of the motif of the repeat being added
     *  @param seed_positions list of identified seeds
     *  @param bset_size size of the shift XOR bitset
     *  @returns none adds the perfect repeat locus to the list of identified perfect repeat loci
    */

    int last_start, last_end, last_rend, last_mlen;       // coordinate variables for existing seeds
    int seed_length = seed_end - seed_start;
    int seed_rend   = seed_end + motif_length;

    // indices for different shifts in motif_bsets
    int seed_midx = motif_length - MINIMUM_SHIFT;
    int overlap_length;

    vector<int> remove_seeds;   // the indices of seeds that need to be removed

    for (int i=seed_positions.size()-1; i>=0; i--) {
        // starting from the last seed and decrementing in indices
        last_start  = get<0> (seed_positions[i]);
        last_end    = get<1> (seed_positions[i]);
        last_mlen   = get<2> (seed_positions[i]);
        last_rend   = last_end + last_mlen;

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_rend < seed_start) break;

        int merge_start = 0, merge_end = 0;
        if (last_start < seed_start) {
            if (seed_rend <= last_rend) { overlap_length = seed_rend - seed_start; merge_start = last_start; merge_end = last_end; }
            else { overlap_length = last_rend - seed_start; merge_start = last_start; merge_end = seed_end; }
        }
        else {
            if (last_rend <= seed_rend) { overlap_length = last_rend - last_start; merge_start = seed_start; merge_end = seed_end; }
            else { overlap_length = seed_rend - last_start; merge_start = seed_start; merge_end = last_end; }
        }

        if (overlap_length <= 0) continue;

        if (last_mlen == motif_length && overlap_length >= motif_length) {
            remove_seeds.push_back(i);
            for (int _=0; _<remove_seeds.size(); _++) seed_positions.erase(seed_positions.begin() + remove_seeds[_]);
            addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, bset_size, seed_positions);
            return;
        }

        else if (last_mlen % motif_length == 0 && overlap_length >= last_mlen) {
            // if the overlap length is at least 1 less than the larger motif size
            // the longer motif repeat with more than 3 units is retained
            remove_seeds.push_back(i);
            for (int _=0; _<remove_seeds.size(); _++) seed_positions.erase(seed_positions.begin() + remove_seeds[_]);
            addSeedToSeedPositionsPerfect(merge_start, merge_end, motif_length, bset_size, seed_positions);
            return;
        }

        else if (motif_length % last_mlen == 0 && overlap_length >= motif_length) {
            // if the overlap length is at least 1 less than the larger motif size
            // the longer motif repeat with more than 3 units is retained
            remove_seeds.push_back(i);
            for (int _=0; _<remove_seeds.size(); _++) seed_positions.erase(seed_positions.begin() + remove_seeds[_]);
            addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, bset_size, seed_positions);
            return;
        }
    }

    for (int i=0; i<remove_seeds.size(); i++) {
        // removing the redundant seeds
        // because the seeds are being removed in reverse order the index of the next
        // seed to be removed is not changed
        seed_positions.erase(seed_positions.begin() + remove_seeds[i]);
    }

    // limiting the seeds to the edge
    if (seed_end > bset_size-motif_length) {
        seed_end = bset_size-motif_length;
    }

    seed_length = seed_end - seed_start;
    seed_positions.push_back(tuple<int,int,int,int,int,int,int> { seed_start, seed_end, motif_length, RANK_P, seed_length,
                                                                  seed_length, seed_length });
}


// function to identify windows based on the threshold of window bit counts
vector<tuple<int,int,int,int,int,int,int>> processShiftXORsPerfect(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset) {
    /*
     *  parsing the shift XORs of all shift sizes and picking seeds from each shift
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param N_bset N position bitset
     *  @return vector<tuple<int, int, int>> vector of end position sorted seeds from all motif sizes
    */

    int bset_size = N_bset.size();          // size of the sequence
    vector<tuple<int,int,int,int,int,int,int>> seed_positions;    // the vector of seed_positions // bool for perfect and imperfect

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT, didx, motif_length;

    int last_starts[NMLENS];     // stores the start of the previous seed
    int last_ends[NMLENS];       // stores the end of the previous seed
    int current_starts[NMLENS];  // stores the current seed start
    int seedlen_cutoffs[NMLENS];
    int valid_end;

    for (int _=0; _<NMLENS; _++) {
        last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1;

        if      (_+MINIMUM_MLEN <= 6) seedlen_cutoffs[_] = 12 - (_+MINIMUM_MLEN);
        else if (_+MINIMUM_MLEN < 20) seedlen_cutoffs[_] = (_+MINIMUM_MLEN);
        else seedlen_cutoffs[_] = 0.8*(_+MINIMUM_MLEN);
    }

    // min_idx - index of the motif_length in shift XOR bitsets
    // didx - index of the motif_length in the seed_positions
    int xor_idx = 0, window_position = 0;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {

        if (N_bset[xor_idx] == 1) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < min_idx+NMLENS; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                if (last_starts[didx] != -1) {
                    if (window_position - last_starts[didx] >= seedlen_cutoffs[motif_length-MINIMUM_MLEN]) {
                        if (PURITY_THRESHOLD == 1) {
                            addPerfectRepeatPositions(last_starts[didx], window_position - motif_length, motif_length, bset_size, seed_positions);
                        }
                        else {
                            addSeedToSeedPositionsPerfect(last_starts[didx], window_position - motif_length, motif_length, bset_size, seed_positions);
                        }
                    }
                    last_starts[didx] = -1;
                }
            }
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            for (int midx=min_idx; midx < min_idx+NMLENS; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                if (motif_bsets[midx][xor_idx] == 1) {
                    if (last_starts[didx] == -1) {
                        last_starts[didx] = window_position;
                    }
                }
                else {
                    if (last_starts[didx] != -1) {
                        if (window_position - last_starts[didx] >= seedlen_cutoffs[motif_length-MINIMUM_MLEN]) {
                            valid_end = (window_position < bset_size - motif_length) ? window_position : bset_size - motif_length;
                            adjustEndBasedonN(N_bset, valid_end, motif_length);
                            if (PURITY_THRESHOLD == 1) {
                                addPerfectRepeatPositions(last_starts[didx], valid_end, motif_length, bset_size, seed_positions);
                            }
                            else {
                                addSeedToSeedPositionsPerfect(last_starts[didx], valid_end, motif_length, bset_size, seed_positions);
                            }
                        }
                    }
                    last_starts[didx] = -1;
                }
            }
        }

        window_position += 1;
    }

    // handling the repeat at the end of the sequence
    for (int midx=min_idx; midx < min_idx+NMLENS; midx++) {
        didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
        if (last_starts[didx] != -1) {
            valid_end = (window_position < bset_size - motif_length) ? window_position : bset_size - motif_length;
            adjustEndBasedonN(N_bset, valid_end, motif_length);
            if (window_position - last_starts[didx] >= seedlen_cutoffs[motif_length-MINIMUM_MLEN]) {
                if (PURITY_THRESHOLD == 1) { addPerfectRepeatPositions(last_starts[didx], valid_end, motif_length, bset_size, seed_positions); }
                else {
                    addSeedToSeedPositionsPerfect(last_starts[didx], valid_end, motif_length, bset_size, seed_positions);
                }
            }
            last_starts[didx] = -1;
        }
    }

    return seed_positions;
}
