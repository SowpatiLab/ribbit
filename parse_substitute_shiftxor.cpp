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
#include "parse_substitute_shiftxor.h"

using namespace std;

int addSeedToSeedPositionsSubstitutions(int seed_start, int seed_end, int motif_length, vector<tuple<int, int, int, int>> &seed_positions_perfect,
                                        vector<tuple<int, int, int, int>> &seed_positions_substut, int* seedlen_cutoff,
                                        vector<boost::dynamic_bitset<>> &motif_bsets, int bset_size, int from_index, int seed_type) {
    /*
     *  add seed to the existing seed positions list
     *  @param seed_start start position of the seed
     *  @param seed_end end position of the seed
     *  @param motif_length motif length of the TR seed
     *  @param seed_positions vector of existing seed positions
     *  @param motif_bsets shift XOR bitsets of all the motif sizes
     *  @param bset_size the total size of a shift XOR bitset
     *  @return none add the seed to seed_position
    */
    
    int last_start, last_end, last_rend, last_mlen;
    int last_length, last_rlen, last_type;       // coordinate variables for existing seeds
    for (int i=from_index; i<seed_positions_perfect.size(); i++) {
        last_start   = get<0> (seed_positions_perfect[i]);

        // go to the point where the start of the last seed is beyond the current seed
        // this logic will leave us with the last seed that is alteast overlapping at least by 1 base at the end
        if (last_start > seed_end) { break; }
        else if (from_index == seed_positions_perfect.size() - 1) { break; }
        else { from_index += 1;  }
    }

    if (seed_end-seed_start < seedlen_cutoff[motif_length-MINIMUM_MLEN]) { return from_index; }

    
    // merging the perfect and substitute seeds into one vector
    vector<int> last_types, last_indices;       // storing the type and indices of seeds that are to be compared
    bool mvnext_perfect = (seed_positions_perfect.size() == 0) ? false : true;
    bool mvnext_substut = (seed_positions_substut.size() == 0) ? false : true;
    int perfect_index = from_index;
    int substut_index = seed_positions_substut.size()-1;
    int perfect_end, substut_end;
    int perfect_type, substut_type;
    
    // while either we can move further in perfect or substitute seeds
    while (mvnext_perfect || mvnext_substut) {

        // if the substitute seeds are exhausted
        if (!mvnext_substut) {            
            while (mvnext_perfect) {    // add only perfect seeds
                perfect_end = get<1> (seed_positions_perfect[perfect_index]);
                perfect_type = get<3> (seed_positions_perfect[perfect_index]);
                if (perfect_end >= seed_start) {
                    if (perfect_type != RANK_N) {
                        last_types.push_back(RANK_P);
                        last_indices.push_back(perfect_index);
                    } perfect_index -= 1;
                }
                // checking if perfect seeds are exhausted
                if (perfect_index < 0 || perfect_end < seed_start) mvnext_perfect = false;
            }
        }

        // if the perfect seeds are exhausted
        else if (!mvnext_perfect) {
            while (mvnext_substut) {    // add only substitute seeds
                substut_end = get<1> (seed_positions_substut[substut_index]);
                substut_type = get<3> (seed_positions_substut[substut_index]);
                if (substut_end >= seed_start) {
                    if (substut_type != RANK_N) {
                        last_types.push_back(RANK_S);
                        last_indices.push_back(substut_index);
                    } substut_index -= 1;
                }
                // checking if substitute seeds are exhausted
                if (substut_index < 0 || substut_end < seed_start) mvnext_substut = false;
            }
        }

        // if neither are exhausted
        else {
            perfect_end = get<1> (seed_positions_perfect[perfect_index]);
            perfect_type = get<3> (seed_positions_perfect[perfect_index]);
            substut_end = get<1> (seed_positions_substut[substut_index]);
            substut_type = get<3> (seed_positions_substut[substut_index]);

            // adding the seed which has the greater end            
            if (substut_end > perfect_end) {
                if (substut_type != RANK_N) {
                    last_types.push_back(RANK_S);
                    last_indices.push_back(substut_index);
                } substut_index -= 1;
            }

            else if (substut_end <= perfect_end) {
                if (perfect_type != RANK_N) {
                    last_types.push_back(RANK_P);
                    last_indices.push_back(perfect_index);
                } perfect_index -= 1;
            }
            
            if (perfect_index < 0 || perfect_end < seed_start) mvnext_perfect = false;
            if (substut_index < 0 || substut_end < seed_start) mvnext_substut = false;
        }
    }
    
    int seed_rend   = seed_end + motif_length;
    int seed_length = seed_end - seed_start;
    int seed_rlen   = seed_length + motif_length;

    // indices for different shifts in motif_bsets
    int  seed_midx = motif_length - MINIMUM_SHIFT;
    int  last_midx = 0, new_type;
    int merge_start = 0, merge_end = 0, overlap_length = 0;

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
        last_length  = last_end - last_start;
        last_rlen    = last_rend - last_start;
        last_midx    = last_mlen - MINIMUM_SHIFT;

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) { break; }

        if (last_type == RANK_N) { continue; }
        
        // if the from_index is much ahead we skip the seeds that do not overlap
        if (seed_end < last_start) { continue; }

        // current seed and last seed have identical coordinates
        if (seed_start == last_start && seed_end == last_end) {
            //  last seed is a pefect seed ~ do not add the current seed
            if (seed_type == RANK_S && (last_type == RANK_P || last_type == RANK_Q)) { return from_index; }
            
            //  current seed is merged and last seed is perfect ~ do not add the current seed
            else if (seed_type == RANK_Q && last_type == RANK_P) { return from_index; }
            
            //  current seed is merged and last seed is substitute ~ remove last seed
            else if (seed_type == RANK_Q && last_type == RANK_S) { 
                seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
            }

            // current and last seeds are substitute or a merged type
            else if ((seed_type == RANK_Q && last_type == RANK_Q) || (seed_type == RANK_S && last_type == RANK_S)) {
                //  existing seed motif length is factor of the new seed motif length
                if (motif_length % last_mlen == 0) { return from_index; }

                //  existing seed's motif length is multiple of the new seed motif length
                else if (last_mlen % motif_length == 0) {
                    //  new motif length is factor of the existing seed motif length
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; //here seed_mlen is shorter.
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                                     seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                    return from_index;
                }
                
                // motif lengths are unequal
                else {
                    bool retain = retainIdenticalSeeds(motif_bsets, seed_start, seed_end, seed_midx, last_midx, bset_size);
                    if (!retain) { return from_index; }
                    else { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; break; }
                }
            }
        }

        // current seed is nested in an existing seed
        else if (last_start <= seed_start && seed_end <= last_end) {
            //  current seed is a subsitute type ~ last seed is either pefect or megred ~ do not add new seed
            if (seed_type == RANK_S && (last_type == RANK_P || last_type == RANK_Q)) { return from_index; }

            // current seed is merged type and last seed is perfect ~ do not add seed
            else if (seed_type == RANK_Q && last_type == RANK_P) { return from_index; }

            // if the nested seed is either equal or mutiple motif length we do not add the seed
            else if ((seed_type == RANK_Q && last_type == RANK_S) || (seed_type == RANK_Q && last_type == RANK_Q) || (seed_type == RANK_S && last_type == RANK_S)) {
                new_type = (seed_type == RANK_S && last_type == RANK_S) ? RANK_S : RANK_Q;
                
                // if the new (nested) and old (parent) seed have the same motif length
                if (motif_length == last_mlen) {
                    // update seed with the type as merged
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, new_type};
                    return from_index;
                }

                // if the new (nested) seed's motif length is a multiple of the old (parent) seed's motif length ~ do not add seed
                else if (motif_length % last_mlen == 0) { return from_index; }

                // if the new (nested) seed's motif length is a factor of or less than the old (parent) seed's motif length
                else if (last_mlen % motif_length == 0 || last_mlen < motif_length) {
                    // merge the seeds only if the new seed's repeat is covering at least 1bp less than the motif size
                    // or at least 1bp less than the seed length
                    if (seed_rlen >= last_mlen-1 || seed_rlen >= last_length - 1) {
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, new_type};
                        return from_index;
                    }
                    // else add the seed separately
                }

                // if the new (nested) seed's motif length is greater than the old (parent) seed's motif length
                else {
                    bool retain = retainNestedSeed(motif_bsets, seed_start, seed_end, seed_midx, last_midx, bset_size);
                    if (!retain) { return from_index; }
                }
            }
        }

        // current seed is parent to an existing seed
        else if (seed_start <= last_start && last_end <= seed_end) {
            if ((seed_type == RANK_S && (last_type == RANK_P || last_type == RANK_Q)) || (seed_type == RANK_Q && last_type == RANK_P)) {
                
                // if new (parent) seed's motif length a factor of old (nested) seed's motif length
                if (last_mlen % motif_length == 0) {
                    // tag existing seed as inactive ~ add merged seed
                    if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    else                     { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                                     seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                    return from_index;
                }

                // if new (parent) seed's motif length a multiple of or greater than old (nested) seed's motif length
                else if (motif_length % last_mlen == 0 || last_mlen < motif_length) {
                    // merge the seeds only if the old seed's repeat is covering at least 1bp less than the motif size
                    // or at least 1bp less than the seed length
                    if (last_rlen >= motif_length - 1 || last_rlen >= seed_length - 1) {
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else                     { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                                         seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                    // else add the seed separately
                }

                else if (motif_length < last_mlen) {
                    // retain both seeds separetly
                }

            }

            else if (seed_type == RANK_Q && last_type == RANK_S) {
                seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                break;
            }

            else if ((seed_type == RANK_Q && last_type == RANK_Q) || (seed_type == RANK_S && last_type == RANK_S)) {

                // if old (nested) seed's motif length is a multiple of new (parent) seed's motif length ~ remove existing seed
                if (last_mlen % motif_length == 0) {
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                }

                // if old (nested) seed's motif length is a factor of new (parent) seed's motif length
                else if ( (motif_length % last_mlen == 0) || (motif_length > last_mlen)) {
                    if (last_rlen >= motif_length - 1 || last_rlen >= seed_length - 1) {
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                                         seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }

                    else {
                        bool retain = retainNestedSeed(motif_bsets, last_start, last_end, last_midx, seed_midx, bset_size);
                        if (retain) { continue; }
                        else { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    }
                }

                else if (last_mlen > motif_length) {
                    bool retain = retainNestedSeed(motif_bsets, last_start, last_end, last_midx, seed_midx, bset_size);
                    if (retain) { continue; }
                    else {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                                         seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }
                }
            }
        }

        // current seed is overlapping with existing seed
        else {
            if (last_start < seed_start) {
                overlap_length = last_rend - seed_start;
                merge_start = last_start; merge_end = seed_end;
            }
            
            else {
                overlap_length = seed_rend - last_start;
                merge_start = seed_start; merge_end = last_end;
            }

            if (((last_mlen % motif_length == 0) || last_mlen > motif_length) && ((overlap_length >= last_mlen - 1) || (overlap_length >= last_length - 1))) {
                if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                else                     { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, motif_length, seed_positions_perfect, seed_positions_substut,
                                                                    seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                return from_index;
            }

            else if (((motif_length % last_mlen == 0) || motif_length > last_mlen) && ((overlap_length >= motif_length - 1) || (overlap_length >= seed_length - 1))) {
                if (last_type != RANK_P) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, last_mlen, seed_positions_perfect, seed_positions_substut,
                                                                    seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                return from_index;
            }
        }
    }

    // limiting the seeds to the edge
    if (seed_end > bset_size-motif_length) {
        seed_end = bset_size-motif_length;
    }
    
    seed_positions_substut.push_back(tuple<int, int, int, int> { seed_start, seed_end, motif_length, seed_type});
    return from_index;
}


vector<tuple<int, int, int, int>> processShiftXORswithSubstitutions(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                                    int &window_length, int &window_bitcount_threshold,
                                                                    vector<tuple<int, int, int, int>> &seed_positions_perfect) {
    /*
     *  parsing the shift XORs of all shift sizes and picking seeds from each shift
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param N_bset N position bitset
     *  @param window_length length of the window to be scanned
     *  @param window_bitcount_threshold the threshold number of set bits in the window
     *  @return vector<tuple<int, int, int>> vector of end position sorted seeds from all motif sizes
    */

    int bset_size = N_bset.size();          // size of the sequence
    int window_bitcount;        // stores window bitcount
    int valid_position = 0;     // position tracking valid bits in the window

    int min_idx = MINIMUM_MLEN-MINIMUM_SHIFT, didx, motif_length;
    int last_starts[NMOTIFS];     // stores the start of the previous seed
    int last_ends[NMOTIFS];       // stores the end of the previous seed
    int current_starts[NMOTIFS];  // stores the current seed start
    int seedlen_cutoffs[NMOTIFS];

    // initialising all to -1
    for (int _=0; _<NMOTIFS; _++) { last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1; seedlen_cutoffs[_] = 10;}

    int from_index = 0;
    vector<tuple<int, int, int, int>> seed_positions_substut;

    vector<boost::dynamic_bitset<>> window_bsets;
    for (int midx=0; midx < NMOTIFS; midx++) {
        boost::dynamic_bitset<> window_bset(window_length, 0ull);
        window_bsets.push_back(window_bset);   // initialised window bitset
        seedlen_cutoffs[midx] = ((midx+MINIMUM_MLEN) > 30) ? (midx+MINIMUM_MLEN)/3 : 10;
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
                        from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length,
                                                                         seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                         motif_bsets, bset_size, from_index, RANK_S);
                        
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
                                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length, 
                                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                                 motif_bsets, bset_size, from_index, RANK_S);
                                
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
                                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length,
                                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                                 motif_bsets, bset_size, from_index, RANK_S);

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
                from_index = addSeedToSeedPositionsSubstitutions(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length,
                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                 motif_bsets, bset_size, from_index, RANK_S);
            }
        }
        
        else {
            if (current_starts[didx] == -1) {
                // presently not scanning through a passed window ~ save last record
                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length,
                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                 motif_bsets, bset_size, from_index, RANK_S);
            }

            else {
                if (last_ends[didx] >= current_starts[didx] - motif_length) {
                    // current passed window overlaps with last record ~ merge both and save
                    last_ends[didx] = bset_size - (xor_idx + 1); // reassign end
                    from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                     motif_bsets, bset_size, from_index, RANK_S);
                }

                else {
                    // current passed window doesn't overlap with last record ~ save both separately
                    from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], motif_length,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                     motif_bsets, bset_size, from_index, RANK_S);
                    from_index = addSeedToSeedPositionsSubstitutions(current_starts[didx], (bset_size - (xor_idx + 1)), motif_length,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs,
                                                                     motif_bsets, bset_size, from_index, RANK_S);
                }
            }
        }
    }

    return seed_positions_substut;
}
