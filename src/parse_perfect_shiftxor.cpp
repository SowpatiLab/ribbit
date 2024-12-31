/*
 * Different methods for parsing shift XOR and identification of tandem repeats
*/


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

#include "global_variables.h"
#include "parse_perfect_shiftxor.h"

using namespace std;


bool retainNestedSeed(vector<boost::dynamic_bitset<>> &motif_bsets, int start, int end,
                      int nested_midx, int parent_midx, int bset_size) {

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

    int nested_count = 0, parent_count = 0;
    for(int i=start; i<end; i++) {
        if (motif_bsets[nested_midx][bset_size - 1 - i] == 1) nested_count += 1;
        if (motif_bsets[parent_midx][bset_size - 1 - i] == 1) parent_count += 1;
    }

    if (nested_count < parent_count) { return false; }
    else if (nested_count == parent_count) { return nested_midx < parent_midx; }
    else { return true; }
}



void addSeedToSeedPositionsPerfect(int seed_start, int seed_end, int motif_length,
                                   vector<tuple<int, int, int, int>> &seed_positions,
                                   vector<boost::dynamic_bitset<>> &motif_bsets, int bset_size) {

    int last_start, last_end, last_mlen;       // coordinate variables for existing seeds
    int seed_length = seed_end-seed_start, seed_rlen = seed_end - seed_start + motif_length;

    // indices for different shifts in motif_bsets
    int seed_midx = motif_length - MINIMUM_SHIFT;
    int last_length, last_rlen, overlap_length;

    vector<int> remove_seeds;   // the indeices of seeds that need to be removed

    for (int i=seed_positions.size()-1; i>=0; i--) {
        // starting from the last seed and decrementing in indices
        last_start = get<0> (seed_positions[i]);
        last_end   = get<1> (seed_positions[i]);
        last_mlen  = get<2> (seed_positions[i]);
        last_length = last_end - last_start;
        last_rlen = last_length + last_mlen;

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) break;

        // identical
        if (last_start == seed_start && last_end == seed_end) {
            if (last_mlen < motif_length) { return; }
            else { remove_seeds.push_back(i); }
        }

        // nested
        else if (last_start <= seed_start && last_end >= seed_end) {
            if (seed_rlen < last_mlen/3) { continue; }
            else { return; }
        }

        // parent
        else if (seed_start <= last_start && seed_end >= last_end) {
            if (last_rlen < motif_length/3) { continue; }
            else { remove_seeds.push_back(i); }
        }

        // overlap
        else {
            int merge_start = 0, merge_end = 0;
            if (last_start < seed_start) { overlap_length = last_end - seed_start + last_mlen; merge_start = last_start; merge_end = seed_end; }
            else { overlap_length = seed_end - last_start + motif_length; merge_start = seed_start; merge_end = last_end; }

            if (last_mlen == motif_length) {
                 addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, seed_positions,
                                               motif_bsets, bset_size);
                return;
            }

            else if (last_mlen < motif_length) {
                // if the overlap length is at least 1 less than the larger motif size
                //  the longer motif repeat with more than 3 units is retained
                if (motif_length - overlap_length <= 1 && seed_rlen/motif_length < 3) {
                    addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, seed_positions,
                                                  motif_bsets, bset_size);
                    return;
                }
                else if (seed_rlen - motif_length - overlap_length <= last_mlen) {
                    return;
                }
            }

            else if (motif_length < last_mlen) {
                // if the overlap length is at least 1 less than the larger motif size
                if (last_mlen - overlap_length <= 1 && last_rlen/last_mlen < 3) {
                    addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, seed_positions,
                                                  motif_bsets, bset_size);
                    return;
                }
                else if (last_rlen - last_mlen - overlap_length <= motif_length) {
                    remove_seeds.push_back(i);
                }
            }
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

    seed_positions.push_back(tuple<int, int, int, int> { seed_start, seed_end, motif_length, RANK_P});
}


// function to identify windows based on the threshold of window bit counts
vector<tuple<int, int, int, int>> processShiftXORsPerfect(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                          int &window_length, int &window_bitcount_threshold) {
    /*
     *  parsing the shift XORs of all shift sizes and picking seeds from each shift
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param N_bset N position bitset
     *  @param window_length length of the window to be scanned
     *  @param window_bitcount_threshold the threshold number of set bits in the window
     *  @return vector<tuple<int, int, int>> vector of end position sorted seeds from all motif sizes
    */

    int bset_size = N_bset.size();          // size of the sequence
    vector<tuple<int, int, int, int>> seed_positions;    // the vector of seed_positions // bool for perfect and imperfect

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT, cutoff, didx, motif_length;
    int last_starts[NMOTIFS] = {-1};  // initialising a last record
    int last_ends[NMOTIFS] = {-1};  // initialising a last record
    int current_starts[NMOTIFS] = {-1};  // initialising a last record

    vector<boost::dynamic_bitset<>> window_bsets;
    for (int midx=0; midx < NMOTIFS; midx++) {
        boost::dynamic_bitset<> window_bset(window_length, 0ull);
        window_bsets.push_back(window_bset);   // initialised window bitset
    }

    int xor_idx = 0;
    int window_position = 0;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {

        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < min_idx+NMOTIFS; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                cutoff = (motif_length <= 6) ? 12-motif_length : motif_length+midx;
                if (last_starts[didx] != -1) {
                    if (window_position - last_starts[didx] >= cutoff) {
                        addSeedToSeedPositionsPerfect(last_starts[didx], window_position, motif_length, seed_positions, motif_bsets, bset_size);
                    }
                    last_starts[didx] = -1;
                }
            }
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            for (int midx=min_idx; midx < min_idx+NMOTIFS; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                cutoff = (motif_length <= 6) ? 12-motif_length : motif_length;
                if (motif_bsets[midx][xor_idx]) {
                    if (last_starts[didx] == -1) {
                        last_starts[didx] = window_position;
                    }
                }
                else {
                    if (last_starts[didx] != -1) {
                        if (window_position - last_starts[didx] >= cutoff) {
                            addSeedToSeedPositionsPerfect(last_starts[didx], window_position, motif_length, seed_positions, motif_bsets, bset_size);
                        }
                    }
                    last_starts[didx] = -1;
                }
            }
        }

        window_position += 1;
    }

    window_position -= 1;
    for (int midx=min_idx; midx < min_idx+NMOTIFS; midx++) {
        didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
        cutoff = (motif_length <= 6) ? 12-motif_length : motif_length;
        if (last_starts[didx] != -1) {
            if (window_position - last_starts[didx] >= cutoff) {
                addSeedToSeedPositionsPerfect(last_starts[didx], window_position, motif_length, seed_positions, motif_bsets, bset_size);
            }
            last_starts[didx] = -1;
        }
    }

    return seed_positions;
}