#include "parse_anchored_shiftxor.h"

using namespace std;


vector<tuple<int, int>> getContinuousStretches(boost::dynamic_bitset<> &bset, int start_pos, int end_pos, int &motif_bitcount) {
    /*
     *  returns a list of positions and lengths of continuous stretches of 1s in a bitset
     *  @param bset the bitset to be processed
     *  @return a vector of tuples, each containing the start position and length of a stretch
    */

    vector<tuple<int, int>> stretches;
    int bset_size = bset.size();
    int start = -1, length = 0;
    for (int i = bset_size - 1 - start_pos; i >= bset_size - end_pos; i--) {
        if (bset[i] == 1) {
            motif_bitcount += 1;
            if (start == -1) { start = i; length = 1; }
            else { length += 1; }
        }
        else {
            if (start != -1) {
                if (length >= 5) { stretches.push_back(tuple<int, int> {start, length}); }
                start = -1; length = 0;
            }
        }
    }

    return stretches;
}


void getBitcount(boost::dynamic_bitset<> &bset, int start_pos, int end_pos, int &motif_bitcount) {
    /*
     *  calculates the number of 1s in a bitset between start and end positions
     *  @param bset the bitset to be processed
     *  @param start_pos the start position
     *  @param end_pos the end position
     *  @return motif_bitcount the number of 1s in the specified range
    */

    int bset_size = bset.size();
    motif_bitcount = 0;
    for (int i = bset_size - 1 - start_pos; i >= bset_size - end_pos; i--) {
        motif_bitcount += bset[i];
    }
}


void buildSupportCoverage(vector<int> last_indices, vector<int> last_types, int motif_length, int seed_start, int seed_end,
                          vector<tuple<int, int, int, int>> &seed_positions_perfect, vector<tuple<int, int, int, int>> &seed_positions_substut,
                          vector<tuple<int, int>> &support) {
    
    int last_start, last_end, last_rend, last_mlen;
    int last_length, last_rlen, last_type, last_midx;       // coordinate variables for existing seeds
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

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) { break; }

        if (last_type == RANK_N) { continue; }

        // if the from_index is much ahead we skip the seeds that do not overlap
        if (seed_end < last_start) { continue; }

        last_length = last_end - last_start;
        last_rlen  = last_rend - last_start;
        last_midx  = last_mlen - MINIMUM_SHIFT;

        // current seed is parent in an existing seed
        if (last_type >= RANK_S) {
            
            if ((seed_start <= last_start && last_start <= seed_end + motif_length) ||
                (seed_start <= last_rend && last_rend <= seed_end + motif_length)) {

                if (motif_length == last_mlen) {
                    support.push_back(tuple<int, int> {last_start, last_end});
                }
            }
        }
    }

    std::sort(support.begin(), support.end(), [](const tuple<int, int>& a, const tuple<int, int>& b) {
        return get<0>(a) < get<0>(b);
    });
}


boost::dynamic_bitset<> retainContinuous(boost::dynamic_bitset<> &bset, int x, int y) {


    int bset_size = bset.size();
    int start = -1;
    
    boost::dynamic_bitset<> anchor_bset(bset_size, 0ull);
    for (int idx = bset_size-1; idx >= 0; idx--) {

        if (bset[idx] == 1) {
            if (start == -1) start = idx;
        }

        else {
            if (start - idx >= x && start - idx <= y) {
                // the anchors to be retained should be at least of the minimum anchor size mentioned
                // and not more than twice of the motif length it being tagged in because this will retain all the
                // perfect repeats of that motif length
                anchor_bset.set(idx+1, start - idx, 1);
            } start = -1;
        }
    }

    return anchor_bset;
}


boost::dynamic_bitset<> retainContinuousBitOp(boost::dynamic_bitset<> &bset, int x, int y = -1) {


    int bset_size = bset.size();
    
    boost::dynamic_bitset<> anchor_bset(bset_size, 0ull);
    boost::dynamic_bitset<> atleast_anchor = bset;
    for (int i = 1; i < x; i++) {
        atleast_anchor &= (atleast_anchor >> 1);
    }
    for (int i = 1; i < x; i++) {
        atleast_anchor |= (atleast_anchor << 1);
    }

    if (y == -1) { return atleast_anchor; }

    boost::dynamic_bitset<> atleast_motif = bset;
    for (int i = 1; i < y + 1; i++) {
        atleast_motif &= (atleast_motif >> 1);
    }
    for (int i = 1; i < y + 1; i++) {
        atleast_motif |= (atleast_motif << 1);
    }

    anchor_bset = atleast_anchor & ~atleast_motif;

    return anchor_bset;
}


void generateAnchorShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                             vector<boost::dynamic_bitset<>> &lsxor_anchor_bsets, int anchor_size) {

    /*
     *  generates shift XOR bitsets only retaining the anchors
     *  @param lshift_xor_bsets vector of left shift XOR bitsets of all shifts
     *  @param N_bset bitset with information of N positions
     *  @param lsxor_anchor_bsets vector of the left shift anchor bitsets
     *  @param anchor_size the length of the anchor size
     *  @return void
    */

    int bset_size = N_bset.size();
    int anchor_start = -1;
    int motif_length;
    for (int lsxor_idx=0; lsxor_idx < NSHIFTS; lsxor_idx++) {
        motif_length = MINIMUM_SHIFT + lsxor_idx;
        if (motif_length < anchor_size) { motif_length = anchor_size; } 
        // boost::dynamic_bitset<> anchor_bset = retainContinuous(lshift_xor_bsets[lsxor_idx], anchor_size, motif_length);
        boost::dynamic_bitset<> anchor_bset = retainContinuousBitOp(lshift_xor_bsets[lsxor_idx], anchor_size, motif_length);
        lsxor_anchor_bsets.push_back(anchor_bset);
        anchor_start = -1;
    }
}


void generatePerfectShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                              vector<boost::dynamic_bitset<>> &lsxor_perfect_bsets, int anchor_size) {

    /*
     *  generates shift XOR bitsets only retaining the perfect repeats
     *  @param lshift_xor_bsets vector of left shift XOR bitsets of all shifts
     *  @param N_bset bitset with information of N positions
     *  @param lsxor_perfect_bsets vector of the left shift perfect bitsets
     *  @return void
    */

    int bset_size = N_bset.size();
    for (int lsxor_idx=0; lsxor_idx < NSHIFTS; lsxor_idx++) {
        boost::dynamic_bitset<> perfect_bset = retainContinuousBitOp(lshift_xor_bsets[lsxor_idx], anchor_size);
        lsxor_perfect_bsets.push_back(perfect_bset);
    }
}


tuple<int,int> addSeedToSeedPositionsAnchored(int seed_start, int seed_end, int motif_length, vector<tuple<int, int, int, int>> &seed_positions_perfect,
                                              vector<tuple<int, int, int, int>> &seed_positions_substut, vector<tuple<int, int, int, int>> &seed_positions_anchored,
                                              int* seedlen_cutoffs, vector<boost::dynamic_bitset<>> &motif_bsets,
                                              vector<boost::dynamic_bitset<>> &perfect_bsets,
                                              vector<boost::dynamic_bitset<>> &anchored_bsets,
                                              int bset_size, tuple<int,int> from_indices, int seed_type) {
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
    mergeAllLists(seed_positions_perfect, seed_positions_substut, from_index_perfect, from_index_substut, last_types, last_indices, seed_start);

    int seed_length = seed_end - seed_start;
    int seed_rlen   = seed_length + motif_length;

    // indices for different shifts in motif_bsets
    int merge_start = 0, merge_end = 0, overlap_length = 0;
    
    int motif_bitcount = 0, perfect_bitcount = 0, anchored_bitcount = 0;

    getBitcount(anchored_bsets[motif_length - MINIMUM_MLEN], seed_start, seed_end, anchored_bitcount);
    getBitcount(motif_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, motif_bitcount);
    getBitcount(perfect_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, perfect_bitcount);

    if (perfect_bitcount <= 0) { return tuple<int,int> { from_index_perfect, from_index_substut }; }
    if (motif_bitcount < 0.3*seed_length) { return tuple<int,int> { from_index_perfect, from_index_substut }; }

    int minimumSuccesses = minimumNumberOfSuccesses(seed_length, 5, 0.9);

    bool check = (anchored_bitcount >= minimumSuccesses && motif_bitcount >= 0.3*seed_length); // ? "True" : "False";

    if (check) {
        cout << SEQUENCE_ID << "\t" << CHUNK_START + seed_start << "\t" << CHUNK_START + seed_end << "\t" << motif_length << "\t" << seed_length 
             << "\t" << anchored_bitcount  << "\t" << motif_bitcount << "\t" << perfect_bitcount << "\n";
        seed_positions_anchored.push_back(tuple<int, int, int, int> {seed_start, seed_end, motif_length, seed_type});
    }

    return tuple<int,int> { from_index_perfect, from_index_substut };
}


vector<tuple<int,int,int,int>> processShiftXORsAnchored(vector<boost::dynamic_bitset<>> &lshift_anchored_bsets,
                                                        vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                                                        vector<boost::dynamic_bitset<>> &lshift_perfect_bsets,
                                                        boost::dynamic_bitset<> &N_bset,
                                                        vector<tuple<int, int, int, int>> &seed_positions_perfect,
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

    int current_starts[NMLENS];  // initialising a last record
    int seedlen_cutoffs[NMLENS];

    int motif_length;
    int window_length = 12;
    int offset = 0;   // offset to the start of the window

    

    tuple<int,int> from_indices = {0, 0};
    vector<tuple<int,int,int,int>> seed_positions_anchored;

    // vector<boost::dynamic_bitset<>> window_bsets;
    int window_bcounts[NMLENS];
    int minimum_window_length = 2*MAXIMUM_MLEN; // used for tracking the minimum valid position to start checking for seeds

    for (int midx=0; midx < NMLENS; midx++) {
        motif_length  = MINIMUM_MLEN + midx;
        window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];
        
        if (window_length < minimum_window_length) { minimum_window_length = window_length; }
        // boost::dynamic_bitset<> window_bset(window_length, 0ull);
        
        // window_bsets.push_back(window_bset);   // initialised window bitset
        window_bcounts[midx] = 0;
        seedlen_cutoffs[midx] = (motif_length > SMALL_MLEN_LIMIT) ? motif_length : 10;
        if (motif_length > SMALL_MLEN_LIMIT) { seedlen_cutoffs[midx] = 0.9 * (motif_length); }
    }

    int overlap_distance = 0;   // the allowed overlap distance between adjacent seeds

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT; // the minimum index in motif_bsets to be considered
    int didx; // tracks the index of current starts and last starts; 0 for MINIMUM_MLEN
    
    int last_starts[NMLENS], last_ends[NMLENS];  // initialising a last record
    // initialising all to -1
    for (int _=0; _<NMLENS; _++) { last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1; seedlen_cutoffs[_] = 10;}

    int xor_idx = 0;
    int window_position = -1*minimum_window_length;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {
        window_position += 1;

        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                motif_length = MINIMUM_SHIFT + midx; didx = midx-min_idx;
                if (current_starts[didx] != -1) {
                    // No seed is being tracked currently

                    // start position is the first nuc of the window; it is zero based
                    current_starts[didx] = window_position;

                    // if the last stored seed is beyond the overlapping distance
                    if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                        if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                            from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                          seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                          lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                        }
                        last_starts[didx] = -1; last_ends[didx] = -1;
                    }
                }
                // window_bsets[didx] <<= WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];
                current_starts[didx] = -1;
                window_bcounts[didx] = 0;
            }

            valid_position = 0;
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            valid_position += 1;

            for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                // window_bsets[didx] <<= 1;
                // window_bsets[didx][0] = lshift_anchored_bsets[didx][xor_idx];
                window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];
                window_bcounts[didx] += lshift_anchored_bsets[didx][xor_idx];
                if (valid_position > window_length) { window_bcounts[didx] -= lshift_anchored_bsets[didx][xor_idx + window_length]; }
            }

            if (valid_position >= minimum_window_length) {
                for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                    didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                    window_bitcount = window_bcounts[didx];
                    // window_bitcount = window_bsets[didx].count();
                    window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];

                    if (valid_position < window_length) continue;

                    if (window_bitcount >= WINDOW_THRESHOLDS[motif_length - MINIMUM_MLEN]) {
                        // window bitcount is above the threshold

                        if (current_starts[didx] == -1) {
                            // No seed is being tracked currently

                            // start position is the first nuc of the window; it is zero based
                            offset = 0;
                            // for (offset = window_length-1; offset >= 0; offset--) {
                            //     if (window_bsets[didx][offset] == 1) {
                            //         offset = (window_length - 1) - offset; break;
                            //     }
                            // }
                            int wstart = bset_size - 1 - (window_position - (window_length - minimum_window_length));
                            for (int x = wstart; x > wstart - window_length; x--) {
                                if (lshift_anchored_bsets[didx][x] == 1) {
                                    offset = wstart - x;
                                    break;
                                }
                            }
                            current_starts[didx] = window_position - (window_length - minimum_window_length) + offset; // start position is the first nuc of the window; it is zero based

                            // if the last stored seed is beyond the overlapping distance
                            if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                    from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                                  seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                                  lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                                }
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
                                offset = 0;
                                // for (offset = 0; offset <= window_length-1; offset++) {
                                //     if (window_bsets[didx][offset] == 1) break;
                                // }
                                // int other_offset = 0;
                                int wend = bset_size - 1 - (window_position + minimum_window_length - 1);
                                for (int x = wend; x < wend + window_length; x++) {
                                    if (lshift_anchored_bsets[didx][x] == 1) { break; }
                                    offset += 1;
                                }
                                last_ends[didx] = window_position + minimum_window_length - offset; // end is exclusive
                            }

                            else {
                                // if the last seed is recorded it means that it is within the overlapping range
                                // hence we just update the end of the last record
                                offset = 0;
                                // for (offset = 0; offset <= window_length-1; offset++) {
                                //     if (window_bsets[didx][offset] == 1) break;
                                // }
                                int wend = bset_size - 1 - (window_position + minimum_window_length - 1);
                                for (int x = wend; x < wend + window_length; x++) {
                                    if (lshift_anchored_bsets[didx][x] == 1) { break; }
                                    offset += 1;
                                }
                                last_ends[didx] = window_position + minimum_window_length - offset; // reassign end
                            }

                            current_starts[didx] = -1;
                        }

                        else {
                            // if there is no seed currently being tracked
                            // if the last stored seed is beyond the overlapping distance of current position
                            if (last_ends[didx] != -1 && last_ends[didx] < window_position - overlap_distance) {
                                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                    from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                                  seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                                  lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                                }
                                // the last seed is reset
                                last_starts[didx] = -1; last_ends[didx] = -1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
        didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;

        // handling the records after the end of the sequence
        if (last_ends[didx] == -1) {
            if (current_starts[didx] != -1) {
                // presently not scanning through a passed window ~ save last record
                if ((bset_size - (xor_idx + 1)) - current_starts[didx] >= seedlen_cutoffs[didx]) {
                    addSeedToSeedPositionsAnchored(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length, seed_positions_perfect,
                                                   seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                   lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                }
            }
        }

        else {
            if (current_starts[didx] == -1) {
                // presently not scanning through a passed window ~ save last record
                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                    addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                   seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                   lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                }
            }

            else {
                if (last_ends[didx] >= current_starts[didx] - overlap_distance) { 
                    // current passed window overlaps with last record ~ merge both and save
                    last_ends[didx] = bset_size - (xor_idx + 1); // reassign end
                    if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                        addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                       seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                       lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                    }
                }

                else {
                    // current passed window doesn't overlap with last record ~ save both separately
                    if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                        from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                      seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                      lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                    }

                    if ((bset_size - (xor_idx + 1)) - current_starts[didx] >= seedlen_cutoffs[didx]) {
                        addSeedToSeedPositionsAnchored(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length, seed_positions_perfect,
                                                       seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                       lshift_xor_bsets, lshift_perfect_bsets, lshift_anchored_bsets, bset_size, from_indices, RANK_A);
                    }
                }
            }
        }
    }

    return seed_positions_anchored;
}
