#include "seed_utils.h"

using namespace std;
using namespace boost;

void filterPerfectSeeds(vector<tuple<int,int,int,int,int,int,int>> &seed_positions_perfect,
                        vector<tuple<int,int,int,int,int,int,int>> &seed_positions_substut) {
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
                    seed_positions_perfect[i] = tuple<int,int,int,int,int,int,int>{perfect_start, perfect_end, perfect_mlen, RANK_N, 0,0,0};
                }
            }

            if (substut_end > perfect_end + perfect_mlen) break;
        }
    }
}


void filterShortSeeds(vector<tuple<int,int,int,int,int,int,int>> &seeds) {
    /*
     *  marks seeds which are shorter than the cutoff length as invalid
     *  @param seeds vector of seed positions
     *  @return none
     */

    int seed_start, seed_end, seed_mlen, seed_type;
    tuple<int,int,int,int,int,int,int> seed;
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
            seeds[i] = tuple<int,int,int,int,int,int,int>{seed_start, seed_end, seed_mlen, RANK_N, 0,0,0};
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


void getBitCount(boost::dynamic_bitset<> &bset, int start_pos, int end_pos, int &motif_bitcount) {
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


int longestContinuousMatches(boost::dynamic_bitset<> &bset) {
    /*
     * calculates the longest continuous stretch of 1s in a bitset
     * @param bset input bitset
     * @return int length of the longest continuous stretch of 1s
    */

    int nseq = bset.size(), l = 0, maxl = 0;
    for (int j=nseq-1; j >= 0; j--) {
        if (bset[j] == 1) l += 1;
        else {
            if (l > maxl) { maxl = l; }
            l = 0;
        }
    }
    if (l > maxl) { maxl = l; }

    return maxl;
}


int longestContinuousMatches(boost::dynamic_bitset<> &bset, int start_pos, int end_pos) {
    /*
     * calculates the longest continuous stretch of 1s in a bitset between start and end positions
     * @param bset input bitset
     * @param start_pos the start position
     * @param end_pos the end position
     * @return int length of the longest continuous stretch of 1s
    */

    int nseq = bset.size(), l = 0, maxl = 0;
    for (int j=nseq-1-start_pos; j >= nseq-end_pos; j--) {
        if (bset[j] == 1) l += 1;
        else {
            if (l > maxl) { maxl = l; }
            l = 0;
        }
    }
    if (l > maxl) { maxl = l; }

    return maxl;
}


void processOverlappingSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                             vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets, int bset_size,
                             vector<set<int>> &skip_atomicity) {
    
    vector<int> sorted_idx(overlapping_seeds.size());
    
    sort(overlapping_seeds.begin(), overlapping_seeds.end(), [](const tuple<int,int,int,int,int,int,int> &a, const tuple<int,int,int,int,int,int,int> &b) {
        if (get<0>(a) == get<0>(b))
            return get<1>(a) < get<1>(b);
        return get<0>(a) < get<0>(b);
    });
    
    skip_atomicity.resize(overlapping_seeds.size());
    for (size_t i = 0; i < overlapping_seeds.size(); ++i) {
        skip_atomicity[i] = set<int>();
    }
    
    int istart, iend, imlen, islen, itype, imcount, ipcount, iacount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jmcount, jpcount, jacount, jajump;

    int ostart, oend, olength;

    for (int i=0; i<overlapping_seeds.size(); i++) {
        istart = get<0>(overlapping_seeds[i]);
        iend   = get<1>(overlapping_seeds[i]);
        islen = iend - istart;
        imlen  = get<2>(overlapping_seeds[i]);
        itype  = get<3>(overlapping_seeds[i]);
        imcount = get<4>(overlapping_seeds[i]);
        ipcount = get<5>(overlapping_seeds[i]);
        iacount = get<6>(overlapping_seeds[i]);
        if (imlen <= 6) { iajump = 1; }
        else { iajump = ((2 * imlen / 10) > 1) ? (2 * imlen / 10) : 2; }
        
        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            jstart = get<0>(overlapping_seeds[j]);
            if (jstart > iend) break;
            jend   = get<1>(overlapping_seeds[j]);
            jslen  = jend - jstart;
            jmlen  = get<2>(overlapping_seeds[j]);
            jtype  = get<3>(overlapping_seeds[j]);
            jmcount = get<4>(overlapping_seeds[j]);
            jpcount = get<5>(overlapping_seeds[j]);
            jacount = get<6>(overlapping_seeds[j]);
            if (jmlen <= 6) { jajump = 1; }
            else { jajump = ((2 * jmlen / 10) > 1) ? (2 * jmlen / 10) : 2; }

            if (itype == RANK_N || jtype == RANK_N) continue;

            if (istart == jstart && iend == jend && imlen == jmlen) {
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                continue;
            }

            if (jstart <= iend) { ostart = jstart; }

            if (jend > iend) { oend = iend; }
            else { oend = jend; }
            olength = oend - ostart;

            if (imlen == jmlen && imlen <= SMALL_MLEN_LIMIT) {
                if (islen - olength < imlen && itype == RANK_A && jtype > RANK_A) {
                    overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0,0,0};
                }
                else if (jslen - olength < jmlen && jtype == RANK_A && itype > RANK_A) {
                    overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                }
                continue;
            }

            if ((jmlen < imlen && islen - olength < imlen) && jtype != RANK_N) {
                // do not need to processes if atomicity is identified
                skip_atomicity[i].insert(jmlen);
            }
            else if ((jmlen > imlen && jslen - olength < jmlen) && jtype != RANK_N) {
                // do not need to processes if atomicity is identified
                skip_atomicity[j].insert(imlen);
            }

            // skip the comparisons for seeds of motif length that differ only by the jump size
            // because these are result of the algorithm and the region could possibly be a repeat of those motif sizes
            if (imlen <= jmlen + jajump && imlen >= jmlen - jajump) {
                continue;
            }
            else if (jmlen <= imlen + iajump && jmlen >= imlen - iajump) {
                continue;
            }

            // cout << "Comparing seeds: " << istart << "-" << iend << " (" << imlen << "," << itype << ") and "
            //      << jstart << "-" << jend << " (" << jmlen << "," << jtype << ")\t";
            
            // things to consider now. The relationship of the motif lengths of the seeds
            // Only compare seeds which substantially overlap with each other
            if (jmlen < imlen) {
                // cout << "Overlap length: " << olength << "\t" << islen - olength << "\t";
                if (itype != RANK_P) {
                    if (islen - olength < imlen) {
                        int lmatches = longestContinuousMatches(perfect_bsets[jmlen - MINIMUM_SHIFT], ostart, oend) + imlen;
                        // cout << "Longest matches: " << lmatches << "\t";
                        if ((imlen <= 6 && lmatches >= 12) ||
                            ((imlen > 6 && imlen >= 2*jmlen) && ((imlen < 30 && lmatches >= imlen - 1) || (imlen >= 30 && lmatches >= 0.9*imlen))) ||
                            ((imlen > 6 && imlen < 2*jmlen) && ((imlen < 30 && lmatches >= 2*imlen - 1) || (imlen >= 30 && lmatches >= 0.9*2*imlen)))) {
                            overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0,0,0};
                            // cout << "Marking seed i as invalid";
                        }
                    }
                }
                else if (itype == RANK_P) {
                    if (islen - olength < imlen) {
                        int lmatches = longestContinuousMatches(perfect_bsets[jmlen - MINIMUM_SHIFT], ostart, oend) + imlen;
                        // cout << "Longest matches: " << lmatches << "\t";
                        if ((imlen <= 6 && lmatches >= 12) ||
                            ((imlen > 6 && imlen >= 2*jmlen) && (lmatches >= imlen - 1)) ||
                            ((imlen > 6 && imlen < 2*jmlen) && (lmatches >= 2*imlen - 1))) {
                            overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0,0,0};
                            // cout << "Marking seed i as invalid";
                        }
                    }
                }
            }
            else if (jmlen > imlen) {
                // cout << "Overlap length: " << olength << "\t" << jslen - olength << "\t";
                if (jtype != RANK_P) {
                    if (jslen - olength < jmlen) {
                        int lmatches = longestContinuousMatches(perfect_bsets[imlen - MINIMUM_SHIFT], ostart, oend) + jmlen;
                        // cout << "Longest matches: " << lmatches << "\t";
                        if ((jmlen <= 6 && lmatches >= 12) ||
                            (jmlen > 6 && jmlen >= 2*imlen) && ((jmlen < 30 && lmatches >= jmlen - 1) || (jmlen >= 30 && lmatches >= 0.9*jmlen)) ||
                            (jmlen > 6 && jmlen < 2*imlen) && ((jmlen < 30 && lmatches >= 2*jmlen - 1) || (jmlen >= 30 && lmatches >= 0.9*2*jmlen))) {
                            overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                            // cout << "Marking seed j as invalid";
                        }
                    }
                }
                else if (jtype == RANK_P) {
                    // cout << "Overlap length: " << olength << "\t";
                    if (jslen - olength < jmlen) {
                        int lmatches = longestContinuousMatches(perfect_bsets[imlen - MINIMUM_SHIFT], ostart, oend) + jmlen;
                        // cout << "Longest matches: " << lmatches << "\t";
                        if ((jmlen <= 6 && lmatches >= 12) ||
                            ((jmlen > 6 && jmlen >= 2*imlen) && (lmatches >= jmlen - 1)) ||
                            ((jmlen > 6 && jmlen < 2*imlen) && (lmatches >= 2*jmlen - 1))) {
                            overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                            // cout << "Marking seed j as invalid";
                        }
                    }
                }
                
            }
            // cout << "\n";
        }
    }
}