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
     *  calculates the longest continuous stretch of 1s in a bitset
     *  @param bset input bitset
     *  @return int length of the longest continuous stretch of 1s
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
     *  calculates the longest continuous stretch of 1s in a bitset between start and end positions
     *  @param bset input bitset
     *  @param start_pos the start position
     *  @param end_pos the end position
     *  @return int length of the longest continuous stretch of 1s
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


void previouslyIdentifiedMotif(int seed_start, int seed_end, int motif_length, int chunk_start,
                               vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                               string &motif) {

    int distance = 2000;
    int last_start, last_end, overlap_start, overlap_end, overlap_length;
    if (motif_length > SMALL_MLEN_LIMIT && repeat_loci.size() > 0) {
        for (int _=repeat_loci.size()-1; _ >= 0; _--) {
            
            if (_ >= repeat_loci.size()) { _ = repeat_loci.size()-1; }
            
            if (get<2>(repeat_loci[_]) < seed_start + chunk_start - distance) return;
            
            if (motif_length != get<6>(repeat_loci[_])) continue;
            
            // start comparing repeats from the end
            last_start  = get<1> (repeat_loci[_]);
            last_end    = get<2> (repeat_loci[_]);

            if ((seed_end + chunk_start < last_start) || (seed_start + chunk_start > last_end)) {
                // continue if repeat doesn't overlap with the repeat
                continue;
            }

            overlap_start = (seed_start + chunk_start < last_start) ? last_start : seed_start + chunk_start;
            overlap_end   = (seed_end + motif_length + chunk_start > last_end) ? last_end : seed_end + motif_length + chunk_start;
            overlap_length = overlap_end - overlap_start;
            if (overlap_length >= 0.8 * (seed_end + motif_length - seed_start)) {
                motif = get<3>(repeat_loci[_]); return;
            }
        }
    }

    return;
}


void sortSeedPositions(vector<tuple<int,int,int,int,int,int,int>> &seed_positions) {
    /*
     *  sorts the seed positions based on start position, end position and motif length
     *  @param seed_positions vector of seed positions
     *  @return none
     */

    sort(seed_positions.begin(), seed_positions.end(), [](const tuple<int,int,int,int,int,int,int> &a, const tuple<int,int,int,int,int,int,int> &b) {
        if (get<0>(a) == get<0>(b))
            return get<1>(a) > get<1>(b);
        return get<0>(a) < get<0>(b);
    });
}


void checkAtomicity(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                    vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets) {
    /*
     *  filters out the nested seeds based on the number of matches found in anchored bset and motif bset
     *  @param overlapping_seeds vector of overlapping seed positions
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param perfect_bsets perfect motif bsets of all shift sizes
     *  @param anchored_bsets anchored motif bsets of all shift sizes
     *  @return none
     */
    
    sortSeedPositions(overlapping_seeds);

    int istart, iend, imlen, islen, itype, iacount, imcount, ipcount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jacount, jmcount, jpcount, jajump;

    int ostart, oend, olength;
    int threshold;

    for (int i=0; i<overlapping_seeds.size(); i++) {
        tie(istart, iend, imlen, itype, iacount, imcount, ipcount) = overlapping_seeds[i];
        islen = iend - istart;

        if (imlen <= 6) { iajump = 1; }
        else { iajump = ((2 * imlen / 10) > 1) ? (2 * imlen / 10) : 2; }

        if (itype == RANK_N) continue;

        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            tie(jstart, jend, jmlen, jtype, jacount, jmcount, jpcount) = overlapping_seeds[j];
            jslen = jend - jstart;

            // as the seeds are sorted by start position, if the later seed start exceed the current seed end, we break
            if (jstart > iend) break;
            if (jtype == RANK_N) continue;

            // if the seeds are identical, we mark the later one as invalid
            if (istart == jstart && iend == jend && imlen == jmlen) {
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                continue;
            }

            if (jend <= iend && jmlen != imlen && jslen >= 3*imlen) {
                // seed j is nested and length of j seed is thrice the motif length of imlen
                
                int nia_count = 0, nim_count = 0;
                getBitCount(anchored_bsets[imlen - MINIMUM_MLEN], jstart, jend, nia_count);
                getBitCount(motif_bsets[imlen - MINIMUM_SHIFT],   jstart, jend, nim_count);

                int sd = sqrt(jslen * 0.9 * 0.1);
                threshold = 5 * sd;
                if (jacount - nia_count < threshold) { continue; } // checks if nested seed qualifies

                else if (imlen % jmlen == 0) {
                    int pja_count = 0, pjm_count = 0, pjp_count = 0;

                    // bitcounts of the motif length of nested seed in the span of the parent seed
                    getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], istart, iend, pja_count);
                    getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT],   istart, iend, pjm_count);
                    getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], istart, iend, pjp_count);

                    sd = sqrt(islen * 0.9 * 0.1);
                    threshold = 5 * sd;
                    // if the bitcount for jmlen in the span of i seed is greater in the motif_bset or greater in than the threshold
                    // in the anchored bitset then reassign atomicity for the seed i
                    if (pjm_count >= imcount || pja_count >= threshold + iacount) {
                        overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{ istart, iend, jmlen, RANK_A, pja_count, pjm_count, pjp_count};
                        overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ jstart, jend, jmlen, RANK_N, 0, 0, 0};
                        i -= 1; break;
                    }                    
                }
            }
        }
    }
}


void filterLowerMatchOverlapSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds,
                                  vector<boost::dynamic_bitset<>> &motif_bsets, vector<boost::dynamic_bitset<>> &perfect_bsets,
                                  vector<boost::dynamic_bitset<>> &anchored_bsets) {

    /*
     *  filters out the nested seeds based on the number of matches found in anchored bset and motif bset
     *  @param overlapping_seeds vector of overlapping seed positions
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param perfect_bsets perfect motif bsets of all shift sizes
     *  @param anchored_bsets anchored motif bsets of all shift sizes
     *  @return none
    */

    sortSeedPositions(overlapping_seeds);

    int istart, iend, imlen, islen, itype, iacount, imcount, ipcount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jacount, jmcount, jpcount, jajump;
    int ostart, oend, olength;
    int sd, threshold;

    int i = 0;
    while ( i < overlapping_seeds.size()) {
        tie(istart, iend, imlen, itype, iacount, imcount, ipcount) = overlapping_seeds[i];
        islen = iend - istart;
        if (itype == RANK_N) continue;

        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            tie(jstart, jend, jmlen, jtype, jacount, jmcount, jpcount) = overlapping_seeds[j];
            jslen = jend - jstart;

            // as the seeds are sorted by start position, if the later seed start exceed the current seed end, we break
            if (jstart >= iend) break;
            if (jtype == RANK_N) continue;

            // if the seeds are identical, we mark the later one as invalid
            if (istart == jstart && iend == jend && imlen == jmlen) {
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                continue;
            }

            if (jend <= iend && jmlen != imlen && jslen >= imlen ) {
                // seed j is nested in seed i and length of j seed is thrice the motif length of imlen
                int nia_count = 0, nim_count = 0, nip_count = 0;
                getBitCount(anchored_bsets[imlen - MINIMUM_MLEN], jstart, jend, nia_count);
                getBitCount(motif_bsets[imlen - MINIMUM_SHIFT],   jstart, jend, nim_count);
                getBitCount(perfect_bsets[imlen - MINIMUM_SHIFT], jstart, jend, nip_count);
                sd = sqrt(jslen * 0.9 * 0.1);
                threshold = 5*sd;
                if (jacount < (nia_count - threshold) && jmcount < nim_count) {
                    overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ jstart, jend, jmlen, RANK_N, 0, 0, 0};
                }
            }

            else if (jend > iend && jmlen != imlen) {
                // the seeds overlap partially
                olength = iend - jstart;
                bool icheck = false, jcheck = false;

                // In the overlapping condition
                // If any of the seed overlaps at least 80 % of its length and at least 3 motif lengths, we consider it for trimming
                if (olength >= 0.8*islen && olength >= 3*jmlen) { icheck = true; }
                if (olength >= 0.8*jslen && olength >= 3*imlen) { jcheck = true; }

                int oia_count = 0, oim_count = 0;
                int oja_count = 0, ojm_count = 0;
                getBitCount(anchored_bsets[imlen - MINIMUM_MLEN], jstart, iend, oia_count);
                getBitCount(motif_bsets[imlen - MINIMUM_SHIFT],   jstart, iend, oim_count);
                getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], jstart, iend, oja_count);
                getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT],   jstart, iend, ojm_count);
                
                sd = sqrt(olength * 0.9 * 0.1);
                threshold = 5*sd;
                if (icheck && !jcheck) {
                    if (oia_count < (oja_count - threshold) && imlen > SMALL_MLEN_LIMIT) {
                        // trim seed i to the non-overlapping part
                        int acount = 0, mcount = 0, pcount = 0;
                        getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], istart, jstart, acount);
                        getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT], istart, jstart, mcount);
                        getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], istart, jstart, pcount);
                        overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{ istart, jstart, imlen, RANK_A, acount, mcount, pcount};
                        i -= 1; break;
                    }
                }
                else if (!icheck && jcheck) {
                    if (oja_count < (oia_count - threshold) && jmlen > SMALL_MLEN_LIMIT) {
                        // trim seed j to the non-overlapping part
                        int acount = 0, mcount = 0, pcount = 0;
                        getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], iend, jend, acount);
                        getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT], iend, jend, mcount);
                        getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], iend, jend, pcount);
                        overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ iend, jend, jmlen, RANK_A, acount, mcount, pcount};
                    }
                }
                else if (icheck && jcheck) {
                    // both seeds engulf each other
                    int fia_count = 0, fim_count = 0, fip_count = 0;
                    int fja_count = 0, fjm_count = 0, fjp_count = 0;
                    getBitCount(anchored_bsets[imlen - MINIMUM_MLEN], istart, jend, fia_count);
                    getBitCount(motif_bsets[imlen - MINIMUM_SHIFT],   istart, jend, fim_count);
                    getBitCount(perfect_bsets[imlen - MINIMUM_SHIFT], istart, jend, fip_count);
                    getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], istart, jend, fja_count);
                    getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT],   istart, jend, fjm_count);
                    getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], istart, jend, fjp_count);

                    sd = sqrt((iend - jstart) * 0.9 * 0.1);
                    threshold = 5*sd;
                    if (fia_count - fja_count > threshold && fim_count >= fjm_count) {
                        // trim seed i to the non-overlapping part
                        overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{ istart, jend, imlen, RANK_A, fia_count, fim_count, fip_count};
                        overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ jstart, jend, jmlen, RANK_N, 0, 0, 0};
                        i = -1; break;
                    }
                    else if (fja_count - fia_count > threshold && fjm_count >= fim_count) {
                        // trim seed j to the non-overlapping part
                        overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{ istart, jend, jmlen, RANK_A, fja_count, fjm_count, fjp_count};
                        overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ jstart, jend, jmlen, RANK_N, 0, 0, 0};
                        i = -1; break;
                    }
                }
            }
        }

        if (i == -1) { sortSeedPositions(overlapping_seeds); }
        i += 1;
    }
}


void filterNearAtomicSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                           vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets) {
    /*
     *  filters out the seeds which are not strictly atomic but are very close to atomicity
     *  @param overlapping_seeds vector of overlapping seed positions
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param perfect_bsets perfect motif bsets of all shift sizes
     *  @param anchored_bsets anchored motif bsets of all shift sizes
     *  @return none
    */
    
    sortSeedPositions(overlapping_seeds);

    int istart, iend, imlen, islen, itype, iacount, imcount, ipcount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jacount, jmcount, jpcount, jajump;

    int ostart, oend, olength;

    for (int i=0; i<overlapping_seeds.size(); i++) {
        tie(istart, iend, imlen, itype, iacount, imcount, ipcount) = overlapping_seeds[i];
        islen = iend - istart;
        if (itype == RANK_N) continue;

        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            tie(jstart, jend, jmlen, jtype, jacount, jmcount, jpcount) = overlapping_seeds[j];
            jslen = jend - jstart;

            // as the seeds are sorted by start position, if the later seed start exceed the current seed end, we break
            if (jstart > iend) break;
            if (jtype == RANK_N) continue;

            // if the seeds are identical, we mark the later one as invalid
            if (istart == jstart && iend == jend && imlen == jmlen) {
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                continue;
            }

            ostart = jstart;
            if (jend > iend) { oend = iend; }
            else { oend = jend; }
            olength = oend - ostart;

            // things to consider now. The relationship of the motif lengths of the seeds
            // Only compare seeds which substantially overlap with each other
            if (jmlen < imlen) {
                if (islen - olength < imlen) {
                    int perfect_jlen = longestContinuousMatches(perfect_bsets[jmlen - MINIMUM_SHIFT], ostart, oend) + jmlen;
                    if (itype != RANK_P) {
                        bool check = (imlen <= 6 && perfect_jlen >= 12);  // if it's a short motif then the perfect match should be at least 12
                        check = check || ((imlen > 6 && imlen >= 2*jmlen) && ((imlen < 20 && perfect_jlen >= imlen - 1) || (imlen >= 20 && perfect_jlen >= 0.9*imlen)));
                        check = check || ((imlen > 6 && imlen <  2*jmlen) && ((imlen < 20 && perfect_jlen >= 2*imlen - 1) || (imlen >= 20 && perfect_jlen >= 0.9*2*imlen)));
                        if (check) {
                            overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0,0,0};
                            break;
                        }
                    }
                    else if (itype == RANK_P) {
                        bool check = (imlen <= 6 && perfect_jlen >= 12);  // if it's a short motif then the perfect match should be at least 12
                        check = check || ((imlen > 6 && imlen >= 2*jmlen) && (perfect_jlen >= imlen - 1));
                        check = check || ((imlen > 6 && imlen < 2*jmlen) && (perfect_jlen >= 2*imlen - 1));
                        if (check) {
                            overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0,0,0};
                            break;
                        }
                    }
                }
            }

            else if (jmlen > imlen) {
                if (jslen - olength < jmlen) {
                    int perfect_ilen = longestContinuousMatches(perfect_bsets[imlen - MINIMUM_SHIFT], ostart, oend) + imlen;
                    if (jtype != RANK_P) {
                        bool check = (jmlen <= 6 && perfect_ilen >= 12);  // if it's a short motif then the perfect match should be at least 12
                        check = check || ( (jmlen > 6 && jmlen >= 2*imlen) && ((jmlen < 20 && perfect_ilen >= jmlen - 1) || (jmlen >= 20 && perfect_ilen >= 0.9*jmlen)) );
                        check = check || ( (jmlen > 6 && jmlen < 2*imlen) && ((jmlen < 20 && perfect_ilen >= 2*jmlen - 1) || (jmlen >= 20 && perfect_ilen >= 0.9*2*jmlen)) );
                        if (check) {
                            overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                            continue;
                        }
                    }
                    else if (jtype == RANK_P) {
                        bool check = (jmlen <= 6 && perfect_ilen >= 12);  // if it's a short motif then the perfect match should be at least 12
                        check = check || ((jmlen > 6 && jmlen >= 2*imlen) && (perfect_ilen >= jmlen - 1));
                        check = check || ((jmlen > 6 && jmlen < 2*imlen) && (perfect_ilen >= 2*jmlen - 1));
                        if (check) {
                            overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0,0,0};
                            continue;
                        }
                    }
                }                
            }
        }
    }
}


void processOverlappingSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                             vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets, int bset_size,
                             vector<set<int>> &skip_atomicity) {
    
    /*
     *  processes the overlapping seeds to remove redundant seeds
     *  @param overlapping_seeds vector of overlapping seed positions
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param perfect_bsets perfect motif bsets of all shift sizes
     *  @param anchored_bsets anchored motif bsets of all shift sizes
     *  @param bset_size size of the shift XOR bitsets
     *  @param skip_atomicity vector of sets of motif lengths to skip atomicity checks
     *  @return none
    */

    sortSeedPositions(overlapping_seeds);

    skip_atomicity.resize(overlapping_seeds.size());
    for (size_t i = 0; i < overlapping_seeds.size(); ++i) {
        skip_atomicity[i] = set<int>();
    }
    
    int istart, iend, imlen, islen, itype, imcount, ipcount, iacount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jmcount, jpcount, jacount, jajump;
    int ostart, oend, olength;

    int i = 0;
    while ( i < overlapping_seeds.size()) {
        tie(istart, iend, imlen, itype, iacount, imcount, ipcount) = overlapping_seeds[i];
        islen = iend - istart;
        if (itype == RANK_N) continue;

        if (imlen <= 6) { iajump = 1; }
        else { iajump = ((2 * imlen / 10) > 1) ? (2 * imlen / 10) : 2; }

        vector<tuple<int, int, int>> unsupportive_nested_seeds;  // to store nested seeds information
        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            tie(jstart, jend, jmlen, jtype, jacount, jmcount, jpcount) = overlapping_seeds[j];
            jslen = jend - jstart;
            if (jmlen <= 6) { jajump = 1; }
            else { jajump = ((2 * jmlen / 10) > 1) ? (2 * jmlen / 10) : 2; }

            // as the seeds are sorted by start position, if the later seed start exceed the current seed end, we break
            if (jstart > iend) break;
            if (jtype == RANK_N) continue;

            // if the seeds are identical, we mark the later one as invalid
            if (istart == jstart && iend == jend && imlen == jmlen) {
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                continue;
            }

            ostart = jstart;
            if (jend > iend) { oend = iend; }
            else { oend = jend; }
            olength = oend - ostart;

            if (imlen == jmlen && imlen <= SMALL_MLEN_LIMIT) {
                if (islen - olength < imlen && itype == RANK_A && jtype > RANK_A) {
                    overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{istart, iend, imlen, RANK_N, 0, 0, 0};
                    break;
                }
                else if (jslen - olength < jmlen && jtype == RANK_A && itype > RANK_A) {
                    overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                    continue;
                }
            }

            if (jstart <= istart + 10 && imlen > SMALL_MLEN_LIMIT && jmlen < imlen - iajump && imlen >= 3* jmlen &&
                jtype >= RANK_Q && itype < RANK_Q && (olength + jmlen) >= imlen) {
                // seed j is nested in seed i and length of j seed is atleast the motif length of imlen
                int tia_count = 0, tim_count = 0, tip_count = 0;

                // bitcounts of the motif length of nested seed in the span of the parent seed
                getBitCount(anchored_bsets[imlen - MINIMUM_MLEN], jend, iend, tia_count);
                getBitCount(motif_bsets[imlen - MINIMUM_SHIFT],   jend, iend, tim_count);
                getBitCount(perfect_bsets[imlen - MINIMUM_SHIFT], jend, iend, tip_count);
                overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{ jend, iend, imlen, itype, tia_count, tim_count, tip_count};
                i = -1; break;
            }

            else if (jmlen > SMALL_MLEN_LIMIT && imlen < jmlen - jajump && jmlen >= 3* imlen &&
                     itype >= RANK_Q && jtype < RANK_Q && (olength + imlen) >= jmlen) {
                // seed i is nested in seed j and length of i seed is atleast the motif length of jmlen
                int tja_count = 0, tjm_count = 0, tjp_count = 0;

                // bitcounts of the motif length of nested seed in the span of the parent seed
                getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], iend, jend, tja_count);
                getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT],   iend, jend, tjm_count);
                getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], iend, jend, tjp_count);
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{ iend, jend, jmlen, jtype, tja_count, tjm_count, tjp_count};
                i = -1; break;
            }

            if ((jmlen < imlen && islen - olength < imlen) && jtype != RANK_N) {
                // do not need to processes if atomicity is identified
                skip_atomicity[i].insert(jmlen);
            }
            else if ((jmlen > imlen && jslen - olength < jmlen) && jtype != RANK_N) {
                // do not need to processes if atomicity is identified
                skip_atomicity[j].insert(imlen);
            }
        }

        if (i == -1) { sortSeedPositions(overlapping_seeds); }
        i += 1;
    }
}


void MergeIdenticalMotifSeeds(vector<tuple<int,int,int,int,int,int,int>> &overlapping_seeds, vector<boost::dynamic_bitset<>> &motif_bsets,
                              vector<boost::dynamic_bitset<>> &perfect_bsets, vector<boost::dynamic_bitset<>> &anchored_bsets) {
    /*
     *  merges identical motif seeds into a single seed with updated positions and bitcounts
     *  @param overlapping_seeds vector of overlapping seed positions
     *  @param motif_bsets shift XOR bsets of all shift sizes
     *  @param perfect_bsets perfect motif bsets of all shift sizes
     *  @param anchored_bsets anchored motif bsets of all shift sizes
     *  @return none
    */
    
    
    sortSeedPositions(overlapping_seeds);
    
    int istart, iend, imlen, islen, itype, imcount, ipcount, iacount, iajump;
    int jstart, jend, jmlen, jslen, jtype, jmcount, jpcount, jacount, jajump;

    int ostart, oend, olength;

    for (int i=0; i<overlapping_seeds.size(); i++) {
        tie(istart, iend, imlen, itype, iacount, imcount, ipcount) = overlapping_seeds[i];
        islen = iend - istart;
        if (itype == RANK_N) continue;
        
        // marking all the overlapping seeds as invalid
        for (int j=i+1; j<overlapping_seeds.size(); j++) {
            tie(jstart, jend, jmlen, jtype, jacount, jmcount, jpcount) = overlapping_seeds[j];
            jslen = jend - jstart;

            // as the seeds are sorted by start position, if the later seed start exceed the current seed end, we break
            if (jstart > iend) break;
            if (jtype == RANK_N) continue;

            // if the seeds are identical, we mark the later one as invalid
            if (imlen == jmlen && jmlen > SMALL_MLEN_LIMIT) {
                // merge the two seeds if they overlap and are of same long motif length
                int merge_start = (istart < jstart) ? istart : jstart;
                int merge_end   = (iend > jend) ? iend : jend;
                int acount = 0, mcount = 0, pcount = 0;
                getBitCount(anchored_bsets[jmlen - MINIMUM_MLEN], merge_start, merge_end, acount);
                getBitCount(motif_bsets[jmlen - MINIMUM_SHIFT], merge_start, merge_end, mcount);
                getBitCount(perfect_bsets[jmlen - MINIMUM_SHIFT], merge_start, merge_end, pcount);
                overlapping_seeds[i] = tuple<int,int,int,int,int,int,int>{merge_start, merge_end, jmlen, RANK_A, acount, mcount, pcount};
                overlapping_seeds[j] = tuple<int,int,int,int,int,int,int>{jstart, jend, jmlen, RANK_N, 0, 0, 0};
                i -= 1; break;
            }
        }
    }
}
