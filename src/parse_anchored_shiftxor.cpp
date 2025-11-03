#include "parse_anchored_shiftxor.h"

using namespace std;


tuple<int, int> adjustSeedPositions(boost::dynamic_bitset<> &anchored_bset, boost::dynamic_bitset<> &N_bset, int seed_start, int seed_end) {
    /*
     *  adjusts the seed positions to the actual positions in the anchored bitset
     *  @param anchor_bset the anchored bitset
     *  @param seed_start the start position of the seed in the original bitset
     *  @param seed_end the end position of the seed in the original bitset
     *  @return a tuple containing the adjusted start and end positions
     */

    int bset_size = N_bset.size();
    int adjusted_start = seed_start, adjusted_end = seed_end;

    int window_size = 5;
    int max_offset = 11;
    int xor_idx = 0;

    // Iteratively adjust start position extend the start position backwards
    // to include positions with 5 continuous 1s in anchored_bset
    bool start_found = false;
    int start_adjustment_count = 0;
    do {
        start_found = false;
        for (int offset = 5; (offset <= max_offset + window_size) && (seed_start - offset) >= 0; ++offset) {
            bool found = true;
            for (int j = 0; j < window_size; ++j) {
                xor_idx = bset_size - 1 - (seed_start - offset + j);
                if ((anchored_bset[xor_idx] == 0) || N_bset[xor_idx]) {
                    found = false; break;
                }
            }
            if (found) {
                adjusted_start = seed_start - offset;
                seed_start = adjusted_start; // update for next iteration
                start_found = true;
                start_adjustment_count++;
                break;
            }
        }
    } while (start_found && start_adjustment_count < 2);

    // Iteratively adjust end position
    bool end_found = false;
    int end_adjustment_count = 0;
    do {
        end_found = false;
        for (int offset = 1; offset <= max_offset && (seed_end + offset + window_size) < bset_size; ++offset) {
            bool found = true;
            for (int j = 0; j < window_size; ++j) {
                xor_idx = bset_size - 1 - (seed_end + offset + j);
                if ((anchored_bset[xor_idx] == 0) || N_bset[xor_idx]) {
                    found = false;
                    break;
                }
            }
            if (found) {
                adjusted_end = seed_end + offset + window_size;
                seed_end = adjusted_end; // update for next iteration
                end_found = true;
                end_adjustment_count++;
                break;
            }
        }
    } while (end_found && end_adjustment_count < 2);

    return make_tuple(adjusted_start, adjusted_end);
}


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


boost::dynamic_bitset<> retainContinuousBitOp(boost::dynamic_bitset<> &bset, int x, int y = -1) {

    /*
     *  Retains continuous bits in a bitset based on specified anchor and motif lengths.
     *  @param bset the original bitset
     *  @param x the anchor length
     *  @param y the motif length (optional)
     *  @return a new bitset with retained continuous bits
     */

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
        if (motif_length >= 10) { motif_length = -1; }
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


bool checkSeedValidity(int seed_start, int seed_end, int motif_length,
                       vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                       vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets,
                       int &anchored_bitcount, int &motif_bitcount, int &perfect_bitcount) {
    /*
     *  checks if the seed is valid based on the purity and probability thresholds
     *  @param seed_start start position of the seed
     *  @param seed_end end position of the seed
     *  @param motif_length motif length of the TR seed
     *  @param motif_bsets shift XOR bitsets of all the motif sizes
     *  @param N_bset bitset with information of N positions
     *  @param perfect_bsets bitsets with only perfect repeats
     *  @param anchored_bsets bitsets with anchored repeats
     *  @return bool indicating if the seed is valid
     */

    getBitCount(anchored_bsets[motif_length - MINIMUM_MLEN], seed_start, seed_end, anchored_bitcount);
    getBitCount(motif_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, motif_bitcount);
    getBitCount(perfect_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, perfect_bitcount);

    int seed_length = seed_end - seed_start;
    
    int minimumSuccesses = 0;
    long double probability_of_successes = 0.0;
    int anchor_size = 5;
    anchor_size = longestContinuousMatches(perfect_bsets[motif_length - MINIMUM_MLEN], seed_start, seed_end);
    if (anchor_size < 5) { anchor_size = 5; }

    if (perfect_bitcount <= 0) { return false; }
    if (motif_bitcount < 0.3*seed_length) { return false; }
    
    long double anchored_probability = 0.9;
    if (seed_length >= 500) anchored_probability = 0.8;
    minimumSuccesses = minimumNumberOfSuccesses(seed_length, anchor_size, anchored_probability);
    probability_of_successes = probabilityOfSuccesses(seed_length, anchor_size, anchored_probability, anchored_bitcount);

    bool check = (anchored_bitcount >= minimumSuccesses || (probability_of_successes >= pow(10.0, -1*log10(anchored_bitcount)) * 0.1));

    return check;
}


bool mergeWithPreviousSeed(int seed_start, int seed_end, int motif_length,
                          vector<tuple<int,int,int,int,int,int,int>> &seed_positions,
                          vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                          vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets) {
    /*
     *  merges the current seed with the previous seed if they are overlapping or adjacent
     *  @param seed_start start position of the current seed
     *  @param seed_end end position of the current seed
     *  @param motif_length motif length of the TR seed
     *  @param seed_positions vector of existing seed positions
     *  @return bool indicating if a merge occurred
    */

    int last_start, last_end, last_mlen, last_type;
    for (int i=seed_positions.size()-1; i>=0; i--) {
        last_start   = get<0> (seed_positions[i]);
        last_end     = get<1> (seed_positions[i]);
        last_mlen    = get<2> (seed_positions[i]);
        last_type    = get<3> (seed_positions[i]);

        if (last_type == RANK_N || last_mlen != motif_length) { continue; }

        if (last_end < seed_start - motif_length) { break; }
        // check if the current seed overlaps or is adjacent to the previous seed
        if (last_end >= seed_start - motif_length) {
            // merge the seeds by updating the end position and motif length
            int merge_start = (last_start < seed_start) ? last_start : seed_start;
            int merge_end   = (last_end > seed_end) ? last_end : seed_end;
            int anchored_bitcount = 0, motif_bitcount = 0, perfect_bitcount = 0;
            bool check_merge = checkSeedValidity(merge_start, merge_end, motif_length,
                                                 motif_bsets, N_bset, perfect_bsets, anchored_bsets,
                                                 anchored_bitcount, motif_bitcount, perfect_bitcount);
            if (!check_merge) { return false; }
            seed_positions[i] = tuple<int,int,int,int,int,int,int> { merge_start, merge_end, last_mlen, RANK_A, anchored_bitcount,
                                                                     motif_bitcount, perfect_bitcount };
            return true;
        }
    }

    return false;
}


tuple<int,int> addSeedToSeedPositionsAnchored(int seed_start, int seed_end, int motif_length, vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                                              vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut, vector<tuple<int,int,int,int,int,int,int>> &seed_positions_anchored,
                                              int* seedlen_cutoffs, vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                              vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets,
                                              int bset_size, tuple<int,int> from_indices, int seed_type, int chunk_start,
                                              unordered_map<int, int> &threshold_bits) {
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

    int last_start, last_end;

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
    for (int i=0; i < last_indices.size(); i++) {
        tuple<int,int,int,int,int,int,int> last_seed;
        if (last_types[i] == RANK_P) { last_seed = seed_positions_perfect[last_indices[i]]; }
        else if (last_types[i] == RANK_S) { last_seed = seed_positions_substut[last_indices[i]]; }
    }

    tuple<int, int> adjusted_positions = adjustSeedPositions(anchored_bsets[motif_length - MINIMUM_MLEN], N_bset, seed_start, seed_end);
    int adjusted_start = get<0>(adjusted_positions);
    int adjusted_end   = get<1>(adjusted_positions);
    
    int seed_length = seed_end - seed_start;
    int seed_rlen   = seed_length + motif_length;
    
    int motif_bitcount = 0, perfect_bitcount = 0, anchored_bitcount = 0;
    getBitCount(anchored_bsets[motif_length - MINIMUM_MLEN], seed_start, seed_end, anchored_bitcount);
    getBitCount(motif_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, motif_bitcount);
    getBitCount(perfect_bsets[motif_length - MINIMUM_SHIFT], seed_start, seed_end, perfect_bitcount);

    int minimumSuccesses = 0;
    long double probability_of_successes = 0.0;
    int anchor_size = 5;
    anchor_size = longestContinuousMatches(perfect_bsets[motif_length - MINIMUM_MLEN], seed_start, seed_end);
    if (anchor_size < 5) { anchor_size = 5; }

    if (perfect_bitcount <= 0) { return tuple<int,int> {from_index_perfect, from_index_substut }; }
    if (motif_bitcount < 0.3*seed_length) { return tuple<int,int> {from_index_perfect, from_index_substut }; }

    long double anchored_probability = 0.9;
    if (seed_length >= 500) anchored_probability = 0.8;

    minimumSuccesses = minimumNumberOfSuccesses(seed_length, anchor_size, anchored_probability);
    probability_of_successes = probabilityOfSuccesses(seed_length, anchor_size, anchored_probability, anchored_bitcount);

    bool check = (anchored_bitcount >= minimumSuccesses || (probability_of_successes >= pow(10.0, -1*log10(anchored_bitcount)) * 0.1));
    
    if (check) {
        bool merged = mergeWithPreviousSeed(seed_start, seed_end, motif_length, seed_positions_anchored,
                                            motif_bsets, N_bset, perfect_bsets, anchored_bsets);
        
        if (merged) { return tuple<int,int> { from_index_perfect, from_index_substut}; }
        seed_positions_anchored.push_back(tuple<int,int,int,int,int,int,int> { adjusted_start, adjusted_end, motif_length, seed_type,
                                                                               anchored_bitcount, motif_bitcount, perfect_bitcount });
    }

    return tuple<int,int> { from_index_perfect, from_index_substut };
}


vector<tuple<int,int,int,int,int,int,int>> processShiftXORsAnchored(vector<boost::dynamic_bitset<>> &lsxor_anchored_bsets,
                                                                    vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                                                                    vector<boost::dynamic_bitset<>> &lsxor_perfect_bsets,
                                                                    boost::dynamic_bitset<> &N_bset,
                                                                    vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                                                                    vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut,
                                                                    int chunk_start) {
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

    unordered_map<int, int> threshold_bits; // map to store the threshold bits for different seed lengths

    tuple<int,int> from_indices = {0, 0};
    vector<tuple<int,int,int,int,int,int,int>> seed_positions_anchored;

    // vector<boost::dynamic_bitset<>> window_bsets;
    int window_bitcounts[NMLENS];
    int minimum_window_length = 2*MAXIMUM_MLEN; // used for tracking the minimum valid position to start checking for seeds

    for (int midx=0; midx < NMLENS; midx++) {
        motif_length  = MINIMUM_MLEN + midx;
        window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];
        if (window_length < minimum_window_length) { minimum_window_length = window_length; }
        window_bitcounts[midx] = 0;
        seedlen_cutoffs[midx] = (motif_length > SMALL_MLEN_LIMIT) ? 0.9*motif_length : 10;
    }

    int overlap_distance = 0;   // the allowed overlap distance between adjacent seeds

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT; // the minimum index in motif_bsets to be considered
    int didx; // tracks the index of current starts and last starts; 0 for MINIMUM_MLEN
    
    int last_starts[NMLENS], last_ends[NMLENS];  // initialising a last record
    // initialising all to -1
    for (int _=0; _<NMLENS; _++) { last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1;}

    int xor_idx = 0;
    int window_position = -1*minimum_window_length;
    int valid_end;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {
        window_position += 1;

        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                valid_end = (bset_size - (xor_idx + 1)) - motif_length;
                valid_end = (valid_end < bset_size - motif_length) ? valid_end : bset_size - motif_length;

                // handling the records after the end of the sequence
                if (last_ends[didx] == -1) {
                    if (current_starts[didx] != -1) {
                        // presently not scanning through a passed window ~ save last record
                        if (valid_end - current_starts[didx] >= seedlen_cutoffs[didx]) {
                            addSeedToSeedPositionsAnchored(current_starts[didx], valid_end, motif_length, seed_positions_perfect,
                                                           seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                           lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                           from_indices, RANK_A, chunk_start, threshold_bits);
                        }
                    }
                }

                else {
                    if (current_starts[didx] == -1) {
                        // presently not scanning through a passed window ~ save last record
                        if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                            addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                           seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                           lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                           from_indices, RANK_A, chunk_start, threshold_bits);
                        }
                    }

                    else {
                        overlap_distance = max({0.1*(last_ends[didx]-last_starts[didx]), 0.1*(valid_end - current_starts[didx])});
                        overlap_distance = min({overlap_distance, motif_length/2});
                        if (last_ends[didx] >= current_starts[didx] - overlap_distance) { 
                            // current passed window overlaps with last record ~ merge both and save
                            adjustEndBasedonN(N_bset, valid_end, motif_length);
                            last_ends[didx] = valid_end; // reassign end
                            if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                               seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                               lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                               from_indices, RANK_A, chunk_start, threshold_bits);
                            }
                        }

                        else {
                            // current passed window doesn't overlap with last record ~ save both separately
                            if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                              seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                              lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                                              from_indices, RANK_A, chunk_start, threshold_bits);
                            }

                            if (valid_end - current_starts[didx] >= seedlen_cutoffs[didx]) {
                                addSeedToSeedPositionsAnchored(current_starts[didx], valid_end, motif_length, seed_positions_perfect,
                                                               seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                               lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                               from_indices, RANK_A, chunk_start, threshold_bits);
                            }
                        }
                    }
                }
                last_starts[didx] = -1; last_ends[didx] = -1; current_starts[didx] = -1;
                window_bitcounts[didx] = 0;
            }

            valid_position = 0;
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            valid_position += 1;

            for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];
                window_bitcounts[didx] += lsxor_anchored_bsets[didx][xor_idx];
                if (valid_position > window_length) { window_bitcounts[didx] -= lsxor_anchored_bsets[didx][xor_idx + window_length]; }
            }

            if (valid_position >= minimum_window_length) {
                for (int midx=min_idx; midx < NMLENS+min_idx; midx++) {
                    didx = midx-min_idx; motif_length = MINIMUM_SHIFT + midx;
                    window_bitcount = window_bitcounts[didx];
                    window_length = WINDOW_LENGTHS[motif_length - MINIMUM_MLEN];

                    if (valid_position < window_length) continue;

                    if (window_bitcount >= WINDOW_THRESHOLDS[motif_length - MINIMUM_MLEN]) {
                        // window bitcount is above the threshold

                        if (current_starts[didx] == -1) {
                            // No seed is being tracked currently

                            // start position is the first nuc of the window; it is zero based
                            offset = 0;
                            int wstart = bset_size - 1 - (window_position - (window_length - minimum_window_length));
                            for (int x = wstart; x > wstart - window_length; x--) {
                                if (lsxor_anchored_bsets[didx][x] == 1) {
                                    offset = wstart - x; break;
                                }
                            }
                            current_starts[didx] = window_position - (window_length - minimum_window_length) + offset; // start position is the first nuc of the window; it is zero based

                            // if the last stored seed is beyond the overlapping distance
                            overlap_distance = max({ 0.1*(last_ends[didx]-last_starts[didx]), 0.1*((bset_size - (xor_idx + 1)) - current_starts[didx]) });
                            overlap_distance = min({ overlap_distance, motif_length });
                            if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                    from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length,
                                                                                  seed_positions_perfect, seed_positions_substut,
                                                                                  seed_positions_anchored, seedlen_cutoffs,
                                                                                  lshift_xor_bsets, N_bset, lsxor_perfect_bsets,
                                                                                  lsxor_anchored_bsets, bset_size, from_indices,
                                                                                  RANK_A, chunk_start, threshold_bits);
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
                                int wend = bset_size - 1 - (window_position + minimum_window_length - 1);
                                for (int x = wend; x < wend + window_length; x++) {
                                    if (lsxor_anchored_bsets[didx][x] == 1) { break; }
                                    offset += 1;
                                }
                                valid_end = (window_position + minimum_window_length - offset < bset_size - motif_length) ? window_position + minimum_window_length - offset : bset_size - motif_length;
                                adjustEndBasedonN(N_bset, valid_end, motif_length);
                                last_ends[didx] = valid_end; // end is exclusive
                            }

                            else {
                                // if the last seed is recorded it means that it is within the overlapping range
                                // hence we just update the end of the last record
                                offset = 0;
                                int wend = bset_size - 1 - (window_position + minimum_window_length - 1);
                                for (int x = wend; x < wend + window_length; x++) {
                                    if (lsxor_anchored_bsets[didx][x] == 1) { break; }
                                    offset += 1;
                                }
                                valid_end = (window_position + minimum_window_length - offset < bset_size - motif_length) ? window_position + minimum_window_length - offset : bset_size - motif_length;
                                adjustEndBasedonN(N_bset, valid_end, motif_length);
                                last_ends[didx] = valid_end; // end is exclusive
                            }

                            current_starts[didx] = -1;
                        }

                        else {
                            // if there is no seed currently being tracked
                            // if the last stored seed is beyond the overlapping distance of current position
                            overlap_distance = max({ 0.1*(last_ends[didx]-last_starts[didx]), 0.1*((bset_size - (xor_idx + 1)) - current_starts[didx]) });
                            overlap_distance = min({ overlap_distance, motif_length });
                            if (last_ends[didx] != -1 && last_ends[didx] < window_position - overlap_distance) {
                                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                                    from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length,
                                                                                  seed_positions_perfect, seed_positions_substut,
                                                                                  seed_positions_anchored, seedlen_cutoffs,
                                                                                  lshift_xor_bsets, N_bset, lsxor_perfect_bsets,
                                                                                  lsxor_anchored_bsets, bset_size, from_indices,
                                                                                  RANK_A, chunk_start, threshold_bits);
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
                    valid_end = ((bset_size - (xor_idx + 1)) < bset_size - motif_length) ? (bset_size - (xor_idx + 1)) : bset_size - motif_length;
                    adjustEndBasedonN(N_bset, valid_end, motif_length);
                    addSeedToSeedPositionsAnchored(current_starts[didx], valid_end, motif_length, seed_positions_perfect,
                                                   seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                   lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                   from_indices, RANK_A, chunk_start, threshold_bits);
                }
            }
        }

        else {
            if (current_starts[didx] == -1) {
                // presently not scanning through a passed window ~ save last record
                if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                    addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                   seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                   lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                   from_indices, RANK_A, chunk_start, threshold_bits);
                }
            }

            else {
                overlap_distance = max({ 0.1*(last_ends[didx]-last_starts[didx]), 0.1*((bset_size - (xor_idx + 1)) - current_starts[didx]) });
                overlap_distance = min({ overlap_distance, motif_length });
                if (last_ends[didx] >= current_starts[didx] - overlap_distance) {
                    // current passed window overlaps with last record ~ merge both and save
                    valid_end = ((bset_size - (xor_idx + 1)) < bset_size - motif_length) ? (bset_size - (xor_idx + 1)) : bset_size - motif_length;
                    adjustEndBasedonN(N_bset, valid_end, motif_length);
                    last_ends[didx] = valid_end; // reassign end
                    if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                        addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                       seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                       lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                       from_indices, RANK_A, chunk_start, threshold_bits);
                    }
                }

                else {
                    // current passed window doesn't overlap with last record ~ save both separately
                    if (last_ends[didx] - last_starts[didx] >= seedlen_cutoffs[didx]) {
                        from_indices = addSeedToSeedPositionsAnchored(last_starts[didx], last_ends[didx], motif_length, seed_positions_perfect,
                                                                      seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                                      lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                                      from_indices, RANK_A, chunk_start, threshold_bits);
                    }

                    if ((bset_size - (xor_idx + 1)) - current_starts[didx] >= seedlen_cutoffs[didx]) {
                        valid_end = ((bset_size - (xor_idx + 1)) < bset_size - motif_length) ? (bset_size - (xor_idx + 1)) : bset_size - motif_length;
                        adjustEndBasedonN(N_bset, valid_end, motif_length);
                        addSeedToSeedPositionsAnchored(current_starts[didx], valid_end, motif_length, seed_positions_perfect,
                                                       seed_positions_substut, seed_positions_anchored, seedlen_cutoffs,
                                                       lshift_xor_bsets, N_bset, lsxor_perfect_bsets, lsxor_anchored_bsets, bset_size,
                                                       from_indices, RANK_A, chunk_start, threshold_bits);
                    }
                }
            }
        }
    }

    return seed_positions_anchored;
}
