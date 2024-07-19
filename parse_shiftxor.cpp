/*
 * Different methods for parsing shift XOR and identification of tandem repeats
*/


#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>

#include "parse_shiftxor.h"
#include "parse_seed.h"
#include "global_variables.h"

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



void addSeedToSeedPositionsPerfect(int seed_start, int seed_end, int motif_length, int min_shift,
                                   vector<tuple<int, int, int, int>> &seed_positions,
                                   vector<boost::dynamic_bitset<>> &motif_bsets, int bset_size) {

    int  last_start, last_end, last_mlen;       // coordinate variables for existing seeds
    int seed_length = seed_end-seed_start, seed_rlen = seed_end - seed_start + motif_length;
    
    // indices for different shifts in motif_bsets
    int seed_midx = motif_length - min_shift;
    int last_length, last_rlen, overlap_length;
    
    vector<int> remove_seeds;   // the indeices of seeds that need to be removed

    // cerr << seed_start << "\t" << seed_end << "\t" << motif_length << "\n";

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
                 addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, min_shift, seed_positions,
                                               motif_bsets, bset_size);
                return;
            }
            
            else if (last_mlen < motif_length) {
                // if the overlap length is at least 1 less than the larger motif size
                //  the longer motif repeat with more than 3 units is retained
                if (motif_length - overlap_length <= 1 && seed_rlen/motif_length < 3) {
                    addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, min_shift, seed_positions,
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
                    addSeedToSeedPositionsPerfect(merge_start, merge_end, last_mlen, min_shift, seed_positions,
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


int addSeedToSeedPositionsSubstitutions(int seed_start, int seed_end, int motif_length, int min_shift, int min_motif_size,
                                         vector<tuple<int, int, int, int>> &seed_positions_perfect, vector<tuple<int, int, int, int>> &seed_positions_substut,
                                         int* seedlen_cutoff, vector<boost::dynamic_bitset<>> &motif_bsets, int bset_size, int from_index, int seed_type) {
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
    
    
    // cout << "\nSeed input: " << seed_start << "\t" << seed_end << "\t" << motif_length << "\t" << from_index << "\t";
    int  last_start, last_end, last_rend, last_mlen, last_type, last_rlen;       // coordinate variables for existing seeds
    for (int i=from_index; i<seed_positions_perfect.size(); i++) {
        last_start   = get<0> (seed_positions_perfect[i]);

        // go to the point where the start of the last seed is beyond the current seed
        // this logic will leave us with the last seed that is alteast overlapping at least by 1 base at the end
        if (last_start > seed_end) { break; }
        else if (from_index == seed_positions_perfect.size() - 1) { break; }
        else { from_index += 1;  }
    }
    // cout << from_index << "\n";

    vector<int> last_types, last_indices;
    int perfect_start_bool = false, substut_start_bool = false;
    int perfect_index = from_index, substut_index = seed_positions_substut.size()-1;
    int perfect_end, substut_end;
    int perfect_type, substut_type;
    // cout << "Seed list sizes: " << seed_positions_perfect.size() << "\t" << seed_positions_substut.size() << "\n";

    // when substut seeds are empty in the start
    if (seed_positions_substut.size() == 0) {
        substut_start_bool = true;
        if (seed_positions_perfect.size() == 0) {
            perfect_start_bool = true;
        }
        while (perfect_index > 0 || !perfect_start_bool) {
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

    else {
        if (seed_positions_perfect.size() == 0) {
            perfect_start_bool = true;
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
                perfect_type = get<3> (seed_positions_perfect[perfect_index]);
                substut_end = get<1> (seed_positions_substut[substut_index]);
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

    if (seed_end-seed_start < seedlen_cutoff[motif_length-min_motif_size]) { return from_index; }
    
    int seed_rend = seed_end + motif_length;
    int seed_length = seed_end - seed_start;
    int seed_rlen = seed_length + motif_length;

    // indices for different shifts in motif_bsets
    int  seed_midx = motif_length - min_shift;
    int  last_midx = 0;


    // cout << "Last Types:  "; 
    // for (int _=0; _<last_indices.size(); _++) { cout << last_indices[_] << ":" << last_types[_] << "  ";  }
    // cout << "\n";

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

        if (last_type == RANK_N) { continue; }
        
        last_rlen  = last_rend - last_start;
        last_midx  = last_mlen - min_shift;

        // seed positions are sorted based on the end position
        // once we encounter a seed that is beyond the start of the current seed
        if (last_end < seed_start) { break; }
        
        // if the from_index is much ahead we skip the seeds that do not overlap
        if (seed_end < last_start) { continue; }
        // if (last_types[_] == RANK_S && i >= seed_positions_substut.size()-1) { continue; }
        // cout << "Indices: " << last_types[_] << "\t" << i << "\t" << seed_positions_perfect.size() << "\t" << seed_positions_substut.size() << "\n";
        // cout << "chr1\t" << seed_start << "\t" << seed_end << "\t" << seed_rend << "\t" << motif_length << "\t" << seed_length << "\t" << seed_rlen << "\t" << seed_type << "\n";
        // cout << "chr1\t" << last_start << "\t" << last_end << "\t" << last_rend << "\t" << last_mlen << "\t" << last_end-last_start << "\t" << last_rlen << "\t" << last_type << "\n";

        // identical
        if (seed_start == last_start && seed_end == last_end) {
            //  if the last seed is a pefect seed we do not add the current seed
            if (seed_type == RANK_S && (last_type == RANK_P || last_type == RANK_Q)) { return from_index; }
            
            else if (seed_type == RANK_Q && last_type == RANK_S) { 
                // cout << "Last seed is tagged Negative!\n";
                seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                break;
            }

            else if (seed_type == RANK_Q && last_type == RANK_P) {
                return from_index;
            }

            else if (seed_type == RANK_Q && last_type == RANK_Q) {
                bool retain = retainIdenticalSeeds(motif_bsets, seed_start, seed_end, motif_length-min_shift, last_mlen-min_shift, bset_size);
                if (!retain) { return from_index; }
                else { 
                    // cout << "Last seed is tagged Negative!\n";
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; break; }
            }
            
            else {
                if (last_mlen % motif_length == 0) {
                    //  new motif length is factor of the existing seed motif length
                    // when removing a older seed decrease the from_index by 1
                    // from_index -= 1; remove_seeds.push_back(i);
                    // cout << "Last seed is tagged Negative!\n";
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; //here seed_mlen is shorter.
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, min_shift, min_motif_size,
                                                        seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                    return from_index;
                }
                
                //  existing seed motif length is factor of the new seed motif length
                else if (motif_length % last_mlen == 0) { return from_index; } // last_mlen is shorter, keeping it same

                // new motif length is greater and not a multiple
                // this case is mostly never occuring at this idenitical coordinates
                else if (motif_length >= last_mlen) { return from_index; } // prefering shorter 
                
                else {
                    // when removing an older seed decrease the from_index by 1
                    // from_index -= 1; remove_seeds.push_back(i);
                    // cout << "Last seed is tagged Negative!\n";
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, min_shift, min_motif_size,
                                                        seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                    return from_index;
                }
            }
        }

        // nested
        else if (last_start <= seed_start && seed_end <= last_end) {
            //  if the last seed is a pefect seed we do not add the current seed
            if (seed_type == RANK_S && (last_type == RANK_P || last_type == RANK_Q)) { return from_index; } //  perfect , substituted or merged (P+S)
            
            else if (seed_type == RANK_Q && last_type == RANK_P) { return from_index; } // vice versa
            
            else if (seed_type == RANK_Q && last_type == RANK_S) {
                // if the nested seed is either equal or mutiple motif length we do not add the seed
                if (motif_length == last_mlen) {
                    // cout << "Last seed is tagged Negative!\n";
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, RANK_Q}; //with motif size of new seed
                    return from_index;
                } // Instead of cheacking for the motif size can we just retain the RANK Q.

                else if (motif_length % last_mlen == 0) { return from_index; } // last_mlen smaller than seed_mlen (no change)

                // if the new seed motif length is smaller and factor of existing seed
                else if (last_mlen % motif_length == 0) {
                 // last_mlen - 1
                        // if new seed repeat covers half of the existing seed motif
                        // reassing the existing seed motif length to new seed motif length
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, RANK_Q};
                        return from_index;
                    // }
                    // else { continue; }
                }

                else if (motif_length < last_mlen) {
                    if (seed_rlen >= last_mlen) { 
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, RANK_Q};
                        return from_index;
                    }
                    else { continue; }
                }

                // if the nested seed is not the multiple motif length but greater than the parent motif length
                else {
                    bool retain = retainNestedSeed(motif_bsets, seed_start, seed_end, motif_length-min_shift, last_mlen-min_shift, bset_size);
                    if (!retain) { return from_index; }
                }

            } 

            else {
                // if the nested seed is either equal or mutiple motif length we do not add the seed
                if (motif_length % last_mlen == 0) {
                    return from_index; 
                }

                // if the new seed motif length is smaller and factor of existing seed
                else if (last_mlen % motif_length == 0) {
                    if (seed_rlen >= last_mlen - 1) {
                        // if new seed repeat covers half of the existing seed motif
                        // reassing the existing seed motif length to new seed motif length
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, motif_length, RANK_S};
                        return from_index;
                    }
                    else { continue; }
                }

                // if the nested seed is not the multiple motif length but greater than the parent motif length
                else {
                    bool retain = retainNestedSeed(motif_bsets, seed_start, seed_end, motif_length-min_shift, last_mlen-min_shift, bset_size);
                    if (!retain) { return from_index; }
                }
            }
        }

        // parent
        else if (seed_start <= last_start && last_end <= seed_end) {
            
            if (last_type == RANK_P || last_type == RANK_Q) {
                
                if (last_mlen % motif_length == 0) {
                    // tagging existing seed as inactive; and adding merged seed
                    if (last_types[_] == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                    else { 
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                    }
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, min_shift, min_motif_size,
                                                        seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                    return from_index;
                }
                
                // merging the seeds which have the same motif length or are factors
                else if ( motif_length % last_mlen == 0 ) {
                    // if the perfect repeat is only 1/3 of the repeat motif of the parent
                    // it makes sense to retain it
                    if (last_rlen <= motif_length -1 ) { continue; }  // new_seed motif length /3
                    
                    else {
                        // tagging existing seed as inactive; and adding merged seed
                        if (last_types[_] == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else { 
                            // cout << "Last seed is tagged Negative!\n";
                            seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        }
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                // same as the above condition
                else if (motif_length < last_mlen) {
                    if (last_rlen <= motif_length -1 ) { continue; } // 3

                    else {
                        // tagging existing seed as inactive; and adding merged seed
                        if (last_types[_] == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else { 
                            // cout << "Last seed is tagged Negative!\n";
                            seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        }
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                else if (motif_length > last_mlen) {

                    // if the parent's seed or the motif is covered by the nested repeat then it is skipped
                    // old seed is extended to match the motif
                    if (last_rlen >= motif_length || last_rlen >= seed_length) {
                        if (last_types[_] == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else { 
                            // cout << "Last seed is tagged Negative!\n";
                            seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        }
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }

                    else if (last_start-seed_start <= motif_length && seed_end-last_end <= motif_length) {
                        if (last_types[_] == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        else { 
                            // cout << "Last seed is tagged Negative!\n";
                            seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        }
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }

                    else {
                        continue;
                    }
                }
            }

            else if (seed_type == RANK_Q && last_type == RANK_S) {
                // cout << "Last seed is tagged Negative!\n";
                seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                break;
            }

            else {

                // merging the seeds which have the same motif length or are factors
                if ( last_mlen % motif_length == 0) {
                    // cout << "Last seed is tagged Negative!\n";
                    seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                    from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, min_shift, min_motif_size,
                                                        seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                    return from_index;
                }

                // merging the seeds which have the same motif length or are factors
                else if ( motif_length % last_mlen == 0 ) {
                    // if the perfect repeat is only 1/3 of the repeat motif of the parent
                    // it makes sense to retain it
                    if (last_rlen <= motif_length - 1) { continue; }
                    
                    else {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }
                }   

                else if (motif_length < last_mlen) {
                    bool retain = retainNestedSeed(motif_bsets, last_start, last_end, last_mlen-min_shift, motif_length-min_shift, bset_size);
                    if (retain) { continue; }
                    else {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, motif_length, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }
                }

                else if (motif_length > last_mlen) {
                    // if the parent's seed or the motif is covered by the nested repeat then it is skipped
                    if (last_rlen >= motif_length || last_rlen >= seed_length) {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }

                    else if (last_start-seed_start <= motif_length && ((seed_end)-last_end <= motif_length)) {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(seed_start, seed_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, seed_type);
                        return from_index;
                    }

                    else {
                        bool retain = retainNestedSeed(motif_bsets, last_start, last_end, last_mlen-min_shift, motif_length-min_shift, bset_size);
                        if (retain) { continue; }
                        else {
                            // remove the nested seed which is not having enough 1s
                            // cout << "Last seed is tagged Negative!\n";
                            seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        }
                    }
                }
            }
        }

        // overlap
        else {
            int merge_start = 0, merge_end = 0, overlap_length = 0;
            if (last_start < seed_start) {
                overlap_length = last_rend - seed_start;
                merge_start = last_start;
                merge_end = seed_end;
            }
            
            else {
                overlap_length = seed_rend - last_start;
                merge_start = seed_start;
                merge_end = last_end;
            }
            
            if (last_type == RANK_P || last_type == RANK_Q) {

                if (last_mlen % motif_length == 0) {
                    if (overlap_length >= last_mlen - 1 ) { // should be more than motif length - 1 
                        // cout << "Last seed is tagged Negative!\n";
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        if (last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, motif_length, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                else if (motif_length % last_mlen == 0) {
                    if (overlap_length >= motif_length -1 ) {
                        // cout << "Last seed is tagged Negative!\n";
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        if (last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                else if (motif_length > last_mlen) {
                    if (overlap_length >= motif_length) {
                        // cout << "Last seed is tagged Negative!\n";
                        if (last_type == RANK_P) { seed_positions_perfect[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        if (last_type == RANK_Q) { seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N}; }
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, last_mlen, min_shift, min_motif_size,
                                                                         seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }
            }

            else if (seed_type == RANK_Q && last_type == RANK_S) {
                
                if (last_mlen % motif_length == 0) {
                    if (overlap_length >= last_mlen - 1) {
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, motif_length, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                else if (motif_length % last_mlen == 0) {
                    if (overlap_length >= motif_length - 1) {                        
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, last_mlen, min_shift, min_motif_size,
                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets, bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }

                else if (motif_length < last_mlen) {
                    if (overlap_length >= last_mlen) {                        
                        // cout << "Last seed is tagged Negative!\n";
                        seed_positions_substut[i] = tuple<int, int, int, int> {last_start, last_end, last_mlen, RANK_N};
                        from_index = addSeedToSeedPositionsSubstitutions(merge_start, merge_end, motif_length, min_shift, min_motif_size,
                                                                         seed_positions_perfect, seed_positions_substut, seedlen_cutoff, motif_bsets,
                                                                         bset_size, from_index, RANK_Q);
                        return from_index;
                    }
                }
            }
        }
    }

    // limiting the seeds to the edge
    if (seed_end > bset_size-motif_length) {
        seed_end = bset_size-motif_length;
    }
    
    // cout << "Added: chr1\t" << seed_start << "\t" << seed_rend << "\t" << seed_end << "\t" << motif_length << "\t" << seed_length << "\t" << seed_rlen << "\t" << seed_type << "\n";
    seed_positions_substut.push_back(tuple<int, int, int, int> { seed_start, seed_end, motif_length, seed_type});
    return from_index;
}





// function to identify windows based on the threshold of window bit counts
vector<tuple<int, int, int, int>> processShiftXORsPerfect(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                          int &window_length, int &window_bitcount_threshold, int nshifts, int &min_motif_size,
                                                          int &min_shift) {
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
    vector<tuple<int, int, int, int>> seed_positions;    // the vector of seed_positions // bool for perfect and imperfect

    int min_idx = min_motif_size-min_shift, cutoff, didx;
    int last_starts[nshifts] = {-1};  // initialising a last record
    int last_ends[nshifts] = {-1};  // initialising a last record
    int current_starts[nshifts] = {-1};  // initialising a last record

    vector<boost::dynamic_bitset<>> window_bsets;
    for (int midx=0; midx < nshifts; midx++) {
        boost::dynamic_bitset<> window_bset(window_length, 0ull);
        window_bsets.push_back(window_bset);   // initialised window bitset
    }

    int xor_idx = 0;
    int window_position = 0;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {

        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < min_idx+nshifts; midx++) {
                didx = midx-min_idx;
                cutoff = (min_shift+midx <= 6) ? 12-(min_shift+midx) : min_shift+midx;
                if (last_starts[didx] != -1) {
                    if (window_position - last_starts[didx] >= cutoff) {
                        addSeedToSeedPositionsPerfect(last_starts[didx], window_position, midx+min_shift, min_shift, seed_positions, motif_bsets, bset_size);
                    }
                    last_starts[didx] = -1;
                }
            }
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            for (int midx=min_idx; midx < min_idx+nshifts; midx++) {
                didx = midx-min_idx;
                cutoff = (min_shift+midx <= 6) ? 12-(min_shift+midx) : min_shift+midx;
                if (motif_bsets[midx][xor_idx]) {
                    if (last_starts[didx] == -1) {
                        last_starts[didx] = window_position;
                    }
                }
                else {
                    if (last_starts[didx] != -1) {
                        if (window_position - last_starts[didx] >= cutoff) {
                            addSeedToSeedPositionsPerfect(last_starts[didx], window_position, midx + min_shift, min_shift, seed_positions, motif_bsets, bset_size);
                        }
                    }
                    last_starts[didx] = -1;
                }
            }
        }

        window_position += 1;
    }

    window_position -= 1;
    for (int midx=min_idx; midx < min_idx+nshifts; midx++) {
        didx = midx-min_idx;
        cutoff = (min_shift+midx <= 6) ? 12-(min_shift+midx) : min_shift+midx;
        if (last_starts[didx] != -1) {
            if (window_position - last_starts[didx] >= cutoff) {
                addSeedToSeedPositionsPerfect(last_starts[didx], window_position, midx+min_shift, min_shift, seed_positions, motif_bsets, bset_size);
            }
            last_starts[didx] = -1;
        }
    }

    return seed_positions;
}


int clearSeeds (int from_index, vector<tuple<int, int, int, int>> &seed_positions, vector<int> seed_indices, int seed_start) {

    vector<int> remove_seeds;
    int last_end, last_type;
    
    for (int _=from_index; _>=0; _--) {
        int i = seed_indices[_];
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


vector<tuple<int, int, int, int>> processShiftXORswithSubstitutions(vector<boost::dynamic_bitset<>> &motif_bsets, boost::dynamic_bitset<> &N_bset,
                                                                    int &window_length, int &window_bitcount_threshold, int nshifts, int &min_motif_size,
                                                                    int &min_shift, vector<tuple<int, int, int, int>> &seed_positions_perfect) {
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

    int min_idx = min_motif_size-min_shift; int didx;
    int last_starts[nshifts];  // initialising a last record
    int last_ends[nshifts];  // initialising a last record
    int current_starts[nshifts];  // initialising a last record
    int seedlen_cutoffs[nshifts];

    for (int _=0; _<nshifts; _++) { last_starts[_] = -1; last_ends[_] = -1; current_starts[_] = -1; seedlen_cutoffs[_] = 10;}

    vector<tuple<int, int, int, int>> seed_positions_substut;
    int from_index = 0;

    vector<boost::dynamic_bitset<>> window_bsets;
    for (int midx=0; midx < nshifts; midx++) {
        boost::dynamic_bitset<> window_bset(window_length, 0ull);
        window_bsets.push_back(window_bset);   // initialised window bitset
        seedlen_cutoffs[midx] = ((midx+min_motif_size) > 30) ? (midx+min_motif_size)/3 : 10;
    }

    int overlap_distance = 0;   // the allowed overlap distance between adjacent seeds

    int xor_idx = 0;
    int window_position = -1*window_length;
    for (xor_idx = bset_size-1; xor_idx >= 0; xor_idx--) {
        window_position += 1;
        if (N_bset[xor_idx]) {
            // N is present at this position reset the window
            for (int midx=min_idx; midx < nshifts+min_idx; midx++) {
                didx = midx-min_idx;
                if (current_starts[didx] != -1) {
                    // No seed is being tracked currently

                    // start position is the first nuc of the window; it is zero based
                    current_starts[didx] = window_position;

                    // if the last stored seed is beyond the overlapping distance
                    if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                        from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                            seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets,
                                                                            bset_size, from_index, RANK_S);
                        
                        last_starts[didx] = -1; last_ends[didx] = -1;
                    }
                }
                window_bsets[midx-min_idx] <<= window_length;
                current_starts[midx-min_idx] = -1;
            }

            valid_position = 0;
            /* Should either accept the N into the seed | Print out the passed seed */
        }

        else {
            valid_position += 1;

            for (int midx=min_idx; midx < nshifts+min_idx; midx++) {
                didx = midx-min_idx;
                window_bsets[didx] <<= 1;
                window_bsets[didx][0] = motif_bsets[midx][xor_idx];
            }

            if (valid_position >= window_length) {
                for (int midx=min_idx; midx < nshifts+min_idx; midx++) {
                    didx = midx-min_idx;
                    window_bitcount = window_bsets[didx].count();

                    if (window_bitcount >= window_bitcount_threshold) {
                        // window bitcount is above the threshold

                        if (current_starts[didx] == -1) {
                            // No seed is being tracked currently

                            // start position is the first nuc of the window; it is zero based
                            current_starts[didx] = window_position;

                            // if the last stored seed is beyond the overlapping distance
                            if (last_ends[didx] != -1 && last_ends[didx] < current_starts[didx] - overlap_distance) {
                                // from_index = clearSeeds(from_index, seed_positions, last_starts[didx]);
                                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets,
                                                                                 bset_size, from_index, RANK_S);
                                
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
                                // from_index = clearSeeds(from_index, seed_positions, last_starts[didx]);
                                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);

                                // the last seed is reset
                                last_starts[didx] = -1; last_ends[didx] = -1;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int midx=min_idx; midx < nshifts+min_idx; midx++) {
        didx = midx-min_idx;
        // handling the records after the end of the sequence
        if (last_ends[didx] == -1) {
            if (current_starts[didx] != -1) {
                // presently not scanning through a passed window ~ save last record
                // from_index = clearSeeds(from_index, seed_positions, current_starts[didx]);
                from_index = addSeedToSeedPositionsSubstitutions(current_starts[didx], (bset_size - (xor_idx + 1)), midx + min_shift, min_shift, min_motif_size,
                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);
            }
        }
        
        else {
            if (current_starts[didx] == -1) {
                // presently not scanning through a passed window ~ save last record
                // from_index = clearSeeds(from_index, seed_positions, last_starts[didx]);
                from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                 seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);
            }

            else {
                if (last_ends[didx] >= current_starts[didx] - (midx + min_shift)) {
                    // current passed window overlaps with last record ~ merge both and save
                    last_ends[didx] = bset_size - (xor_idx + 1); // reassign end
                    // from_index = clearSeeds(from_index, seed_positions, last_starts[didx]);
                    from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);
                }

                else {
                    // current passed window doesn't overlap with last record ~ save both separately
                    // from_index = clearSeeds(from_index, seed_positions, last_starts[didx]);
                    from_index = addSeedToSeedPositionsSubstitutions(last_starts[didx], last_ends[didx], midx + min_shift, min_shift, min_motif_size,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);
                    // from_index = clearSeeds(from_index, seed_positions, current_starts[didx]);
                    from_index = addSeedToSeedPositionsSubstitutions(current_starts[didx], (bset_size - (xor_idx + 1)), midx + min_shift, min_shift, min_motif_size,
                                                                     seed_positions_perfect, seed_positions_substut, seedlen_cutoffs, motif_bsets, bset_size, from_index, RANK_S);
                }
            }
        }
    }

    return seed_positions_substut;
}
