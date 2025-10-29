#include "parse_seed.h"

using namespace std;
using namespace boost;
using namespace boost::multiprecision;


int countBitsDiagonal(int row, int col, int diagonal_length, vector<boost::dynamic_bitset<>*> &MATRIX,
                      int &sequence_length, int neighbor_flank=0) {
    /*
     *  count the number of ones across a diagonal
     *  @param row index of the row in the sequence matrix
     *  @param col index of the column in the sequence matrix
     *  @param size length of the diagonal
     *  @param MATRIX dot matrix of the chromosome sequence
     *  @param sequence_length length of the chromosome sequence
     *  @return int position within the sequence that is the possible motif for the seed
     */

    int count = 0;
    bool neighbor_pass = 0;
    for (int s = 0; s < diagonal_length; s++) {
        if ((*MATRIX[row+s])[col-s] == 1) {
            // checking for one within the diagonal
            count += 1;
        }
        else if (neighbor_flank != 0) {
            neighbor_pass = 0;
            for (int n=1; n<=neighbor_flank; n++) {
                // if the right neighbor is within the length of the matrix and the position is 1
                if (col-s-n > 0 && (*MATRIX[row+s])[col-s-1] == 1) neighbor_pass = 1;

                // if the left neighbor is within the length of the matrix and the position is 1
                else if (col-s+n < sequence_length && (*MATRIX[row+s])[col-s-1] == 1) neighbor_pass = 1;
            }
            if (neighbor_pass) count += 1;
        }
    }

    return count;
}


int highestMatchRegion(boost::dynamic_bitset<> &lshift_xor_bset, boost::dynamic_bitset<> &lshift_anchored_bset,
                       int &o_start, int &o_end, int &o_mlen, int &sequence_length, boost::dynamic_bitset<> &N_bset) {
    /*
     *  picks the best region for searching for the motif in large motif seeds
     *  @param lshift_xor_bset the left shift XOR bitset
     *  @param lshift_anchored_bset the left shift anchored bitset
     *  @param o_start the start position of the region
     *  @param o_end the end position of the region
     *  @param o_mlen the length of the motif
     *  @param sequence_length the length of the sequence
     *  @param N_bset the N bitset
     *  @return int the best region for large motif
     */

    int best_region = -1;
    int best_score = -1;
    int window_length = 10000;
    int window_score = -1;
    int xor_idx;

    for (int i = o_start; i <= o_end - window_length; i++) {
        if (i == o_start) {
            window_score = 0;
            for (int j = 0; j < window_length; j++) {
                xor_idx = sequence_length - 1 - (i + j);
                if (lshift_anchored_bset[xor_idx] == 1 && N_bset[xor_idx] == 0) {
                    window_score += 1;
                }
            }
        }
        else {
            xor_idx = sequence_length - 1 - (i - 1);
            if (lshift_anchored_bset[xor_idx] == 1 && N_bset[xor_idx] == 0) {
                window_score -= 1;
            }
            xor_idx = sequence_length - 1 - (i - 1) + window_length;
            if (lshift_anchored_bset[xor_idx] == 1 && N_bset[xor_idx] == 0) {
                window_score += 1;
            }
        }
        if (window_score > best_score) {
            best_score = window_score;
            best_region = i;
        }
    }

    return best_region;
}


uint256_t mostFrequentLongMotif(boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset, int &seed_start,
                                int &seed_sequence_length, int &motif_length, int &sequence_length, vector<boost::dynamic_bitset<>*> &MATRIX) {
    /*
     *  return the index in the seed which is possible motif of the repeat
     *  @param left_bset the bitset of left bit of a sequence
     *  @param right_bset the bitset of right bit of a sequence
     *  @param seed_start starting position of seed in the sequence
     *  @param seed_sequence_length length of the seed sequence
     *  @param motif_length length of the motif
     *  @param sequence_length length of the chromosome sequence
     *  @param MATRIX dot matrix of the chromosome sequence
     *  @return int position within the sequence that is the possible motif for the seed
     */

    int seed_end = seed_start + seed_sequence_length;

    int max_motif_index = 0, max_count = 0, row_count = 0;
    int row, row_offset, cix;

    int dstream_index, ustream_index;
    int max_dindex, max_dcount;
    int dcount;

    int initial_lastrow, prefix_rows, pcindex;

    int doffset = 0;
    if (motif_length <= 20) { doffset = 1; }
    else { doffset = 2; }
    set<uint256_t> motif_units;
    uint256_t motif_unit, ONE = 1;
    int anchor_jump = ((2 * motif_length / 10) > 1) ? (2 * motif_length / 10) : 2;

    vector<int> row_starts;     // indices of all the rows where motif searches will be initiated
    // if seed_sequence length is lesser than 4 times the motif length then row_start is incremented by 1
    if (seed_sequence_length <= 5 * motif_length) {
        for (int rs = seed_start; rs <= seed_end - motif_length; rs += 1) { row_starts.push_back(rs); }
    }
    else {
        for (int rs = seed_start; rs <= seed_end - motif_length; rs += motif_length) {
            for (int _=0; _ < anchor_jump; _++) {
                if (rs + _ >= seed_end - motif_length + 1) break;
                row_starts.push_back(rs + _);
            }
        }
    }

    // outer loop for rows
    for (int row_start : row_starts) {
        if (row_start >= seed_end - motif_length + 1) continue;
        row_count = 0;
        int iterations = 0;

        motif_unit = 0;
        for (int j = row_start; j < row_start+motif_length; j++) {
            motif_unit <<= 1;
            if (left_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;
            motif_unit <<= 1;
            if (right_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;
        }

        // if the motif is encountered before, skip
        if (motif_units.find(motif_unit) != motif_units.end()) continue;
        motif_units.insert(motif_unit);

        dstream_index = row_start + motif_length;
        while (dstream_index < seed_end) {
            max_dindex = -2, max_dcount = 0;

            for (int x=0; x < doffset+1; x++) {
                for (int k: {-1*x, x}) {
                    dcount = 0;
                    for (int i=0; i<motif_length; i++) {
                        if (dstream_index + k + i >= seed_end) break;
                        if ((*MATRIX[row_start + i])[(sequence_length-1) - (dstream_index + k + i)] == 1) dcount += 1;
                    }

                    if (dcount > max_dcount) { max_dcount = dcount; max_dindex = k; }
                    if (k == 0) break;
                }
            }

            row_count += max_dcount;
            dstream_index += max_dindex;
            dstream_index += motif_length;
        }

        ustream_index = row_start - motif_length;
        while (ustream_index > seed_start) {
            max_dindex = -2, max_dcount = 0;

            for (int x=0; x < doffset+1; x++) {
                for (int k: {-1*x, x}) {
                    dcount = 0;
                    for (int i=0; i<motif_length; i++) {
                        if (ustream_index + k + i < 0) break;
                        if ((*MATRIX[row_start + i])[(sequence_length - 1) - (ustream_index + k + i)] == 1) dcount += 1;
                    }
                    if (dcount > max_dcount) { max_dcount = dcount; max_dindex = k; }
                    if (k == 0) break;
                }
            }

            row_count += max_dcount;
            ustream_index += max_dindex;
            ustream_index -= motif_length;
        }

        if (ustream_index < seed_start && abs(ustream_index-seed_start) < motif_length) {

            initial_lastrow = row_start + motif_length - 1;
            pcindex = seed_start + ((motif_length + (ustream_index-seed_start)) - 1);
            prefix_rows = motif_length + (ustream_index-seed_start);

            max_dindex = -2, max_dcount = 0;
            for (int x=0; x < doffset+1; x++) {
                for (int k: {-1*x, x}) {
                    dcount = 0;
                    for (int i=0; i < prefix_rows; i++) {
                        if (pcindex + k - i < seed_start || pcindex + k - i >= seed_end) break;
                        iterations += 1;
                        if ((*MATRIX[initial_lastrow - i])[(sequence_length - 1) - (pcindex + k - i)] == 1) dcount += 1;
                    }
                    if (dcount > max_dcount) { max_dcount = dcount; max_dindex = k; }
                    if (k == 0) break;
                }
            }

            row_count += max_dcount;
        }

        if (row_count > max_count) {
            max_count = row_count; max_motif_index = row_start;
        }
    }

    motif_unit = 0;
    for (int j = max_motif_index; j < max_motif_index+motif_length; j++) {
        motif_unit <<= 1;
        if (left_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;

        motif_unit <<= 1;
        if (right_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;
    }

    return motif_unit;
}


void processLargeMotifSeed(tuple<int, int> seed_position, int chunk_start, int &motif_length, int &seed_type, string &sequence_id,
                           string &sequence, int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset,
                           boost::dynamic_bitset<> &right_bset, boost::dynamic_bitset<> &N_bset, ostream *out, vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                           vector<boost::dynamic_bitset<>> &lshift_anchored_bsets, vector<boost::dynamic_bitset<>*> &MATRIX,
                           StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment,
                           vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, set<int> &skip_atomicity,
                           string identified_motif) {
    /*
     *  processes the seed and finds all the repeats in the sequence
     *  @param seed_position tuple with start and end position of the seed
     *  @param chunk_start the start of the seed sequence
     *  @param motif_length length of the motif
     *  @param seed_type type of the seed (perfect, substitution, anchored)
     *  @param sequence_id the name of the sequence
     *  @param sequence the complete sequence string
     *  @param sequence_length length of the complete sequence
     *  @param xor_bset the XOR conversion of the bitset (also includes anchor bitset)
     *  @param left_bset the dynamic bitset of the left bit of the sequence
     *  @param right_bset the dynamic bitset of the right bit of the sequence
     *  @param N_bset the dynamic bitset of the N positions in the sequence
     *  @param lshift_xor_bsets vector of dynamic bitsets of the shift XORs
     *  @param lshift_anchored_bsets vector of dynamic bitsets of the shift anchored
     *  @param MATRIX dot matrix of the chromosome sequence
     *  @param aligner the aligner object
     *  @param filter the filter object
     *  @param alignment the resultant alignment object
     *  @param repeat_loci set of identified repeat loci to be compared for redundancy
     *  @param skip_atomicity set of atomicities to be skipped
     *  @returns none prints out the repeat locations to the output file
     */

    int seed_start     = get<0> (seed_position);
    int seed_end       = get<1> (seed_position);
    int seed_bset_size = seed_end - seed_start;
    int seed_sequence_length = seed_bset_size + motif_length;

    for (int s=seed_start; s<seed_end+motif_length; s++) {
        if (N_bset[sequence_length-1-s] == 1) {
            seed_sequence_length = s - seed_start;
            break;
        }
    }
    string seed_sequence = sequence.substr(seed_start, seed_sequence_length);
    // the shift xor bitset of the complete repeat sequence
    boost::dynamic_bitset<> seed_bset(seed_bset_size, 0ull);
    for (int j = seed_start; j < seed_end; j++) {
        seed_bset[seed_end - 1 - j] = xor_bset[sequence_length - 1 - j];
    }

    // if the length of the seed is shorter than the motif size
    if (seed_end - seed_start < SEEDLEN_CUTOFF[motif_length - MINIMUM_MLEN]) return;

    int continuous_threshold = 3; // threshold for continuous matches
    // if the longest continuous stretch of 1s in the seed is lesser than threshold
    int longest_stretch = longestContinuousMatches(seed_bset);
    if (longest_stretch < continuous_threshold) { return; }

    vector<pair<int, int>> seed_repeat_loci;
    string perfect_repeat, motif;
    vector<int> cigar_values;

    int repeat_start, repeat_end, match_nucs, mismatch_nucs, match_units;
    int repeat_length, repeat_units;
    int alignment_length, interruptions, atomicity;
    double purity = 0, motifwise_purity = 0; int motifwise_indels = 0;
    string cigar_string, repeat_representation;
    int left_flank, right_flank;

    if (seed_type == RANK_P) {
        atomicity = motif_length;
        motif = seed_sequence.substr(0, atomicity);
        repeat_start = seed_start; repeat_end = seed_end + motif_length;
        repeat_length = repeat_end - repeat_start;
        if (repeat_length >= MINIMUM_LENGTH[atomicity]) {
            repeat_units = repeat_length/atomicity;
            match_units = repeat_units; purity = 1.0; cigar_string = to_string(repeat_length) + "M";
            purity = 1.0; motifwise_purity = 1.0; motifwise_indels = 0;

            if (((atomicity < 10  && (match_units >= PERFECT_UNITS[atomicity] || (purity > 0.9 && purity*repeat_length >= 2*atomicity))) 
                || ((atomicity >= 10) && (((purity * repeat_length) >= 3*atomicity) || (purity > 0.9 && purity*repeat_length >= 2*atomicity))))
                && atomicity >= MINIMUM_MLEN && atomicity <= MAXIMUM_MLEN
                && repeat_length >= MINIMUM_LENGTH[atomicity]
                && purity >= PURITY_THRESHOLD
                && motifwise_purity >= MOTIFPURITY_THRESHOLD) {

                repeat_start += chunk_start; repeat_end += chunk_start;
                int recursion_level = 0;
                vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
                addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                                 atomicity, repeat_length, repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
            }
        }
        return;
    }

    // length of the pseudo perfect sequence
    int ppr_length = seed_sequence_length + (2*motif_length) + ((1-PURITY_THRESHOLD)*seed_sequence_length);
    uint256_t motif_unit;

    if (identified_motif == "") {
        int motif_search_start = seed_start;
        int large_seed_threshold = 10000;
        if (seed_bset_size > large_seed_threshold) {
            motif_search_start = highestMatchRegion(lshift_xor_bsets[motif_length - MINIMUM_SHIFT],
                                                             lshift_anchored_bsets[motif_length - MINIMUM_MLEN],
                                                             seed_start, seed_end, motif_length, sequence_length, N_bset);
            motif_unit = mostFrequentLongMotif(left_bset, right_bset, motif_search_start, large_seed_threshold, motif_length, sequence_length, MATRIX);
        }
        else {
            motif_unit = mostFrequentLongMotif(left_bset, right_bset, seed_start, seed_sequence_length, motif_length, sequence_length, MATRIX);
        }

        if (THREADS > 1) MTX.lock();
        atomicity = calculateAtomicityLongMotif(motif_unit, motif_length);
        if (skip_atomicity.find(atomicity) != skip_atomicity.end()) {
            if (THREADS > 1) MTX.unlock();
            return;
        }
        if (THREADS > 1) MTX.unlock();
    
        if (atomicity <= SMALL_MLEN_LIMIT) {
            processSmallMotifSeed(tuple<int, int> { seed_start, seed_end }, chunk_start, atomicity, seed_type, sequence_id, sequence,
                                  sequence_length, lshift_xor_bsets[atomicity-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                  out, aligner, filter, alignment, repeat_loci);
            return;
        }
    
        if (motif_length % atomicity != 0) { return; }
    
        // the repeat should be treated based on the atomicity
        if (THREADS > 1) MTX.lock();
        motif = calculateMotif(motif_unit, motif_length);
        if (THREADS > 1) MTX.unlock();
    }

    else { motif = identified_motif; atomicity = motif.length(); }

    motif = motif.substr(0, atomicity);
    motif_unit >>= 2*(motif_length - atomicity);

    perfect_repeat = "";
    while(perfect_repeat.length() <= ppr_length) perfect_repeat += motif;
   
    string pruned_cigar = "";
    if (seed_sequence_length < 5000) {
        aligner.Align(seed_sequence.c_str(), perfect_repeat.c_str(), ppr_length, filter, &alignment, 15);
        pruned_cigar = alignment.cigar_string;
    }
    else {
        pruned_cigar = alignLargeSequence(seed_sequence, motif, atomicity, aligner, filter, alignment);
    }

    bool seed_trim = true;
    if (seed_sequence_length > 10000) {
        repeat_start = seed_start; repeat_end = seed_end + motif_length;
        vector<tuple<int, int, string, double>> largecigar_repeats = processLargeCigar(repeat_start, repeat_end, atomicity,
                                                                                       pruned_cigar, seed_repeat_loci);

        for (auto &lcr : largecigar_repeats) {
            repeat_start = get<0>(lcr); repeat_end = get<1>(lcr);
            cigar_string = get<2>(lcr); purity = get<3>(lcr); 
            repeat_length = repeat_end - repeat_start;
            repeat_units = repeat_length/atomicity;
            repeat_start += chunk_start; repeat_end += chunk_start;
            int recursion_level = 0;
            vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
            addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                             atomicity, repeat_length, repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
        }
    }

    else {

        processCIGARWithPruning(seed_start, seed_sequence_length, pruned_cigar, seed_sequence, atomicity,
                                repeat_start, repeat_end, alignment_length, match_units, cigar_string, purity,
                                motifwise_purity, motifwise_indels, seed_trim);
        trimCigarMotifPurity(cigar_string, atomicity, repeat_start, repeat_end, alignment_length, purity, motifwise_purity);
    
        if (seed_repeat_loci.size()==0) {
            seed_repeat_loci.push_back(pair<int, int> { repeat_start, repeat_end - atomicity });
        }
        else {
            bool inserted = false;
            for (int i=0; i<seed_repeat_loci.size(); i++) {
                if (seed_repeat_loci[i].first >= repeat_start) {
                    seed_repeat_loci.insert(seed_repeat_loci.begin()+i, pair<int, int> { repeat_start, repeat_end - atomicity } );
                    inserted = true;
                    break;
                }
            }
            if (!inserted) { seed_repeat_loci.push_back(pair<int, int> { repeat_start, repeat_end - atomicity }); }
        }
    
        if (alignment_length >= MINIMUM_LENGTH[atomicity]) {
            repeat_length = repeat_end - repeat_start;
            repeat_units = repeat_length/atomicity;
    
            if (((atomicity < 10  && (match_units >= PERFECT_UNITS[atomicity] || (purity > 0.9 && purity*repeat_length >= 2*atomicity))) 
                 || ((atomicity >= 10) && (((purity * repeat_length) >= 3*atomicity) || (purity > 0.9 && purity*repeat_length >= 2*atomicity))))
                && atomicity >= MINIMUM_MLEN && atomicity <= MAXIMUM_MLEN
                && repeat_length >= MINIMUM_LENGTH[atomicity]
                && purity >= PURITY_THRESHOLD
                && motifwise_purity >= MOTIFPURITY_THRESHOLD) {
    
                repeat_start += chunk_start; repeat_end += chunk_start;
                int recursion_level = 0;
                vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
                addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                                 atomicity, repeat_length, repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
            }
        }
    }

    if (seed_repeat_loci.size()==0) { return; }

    int flank_start = seed_start;
    for (int i=0; i<seed_repeat_loci.size(); i++) {

        if (flank_start >= seed_repeat_loci[i].first) { flank_start = seed_repeat_loci[i].second; continue;  }
        
        if (seed_repeat_loci[i].first - flank_start >= MINIMUM_LENGTH[atomicity]) {
            if (flank_start < seed_start) { flank_start = seed_start; }
            if (seed_repeat_loci[i].first > seed_end) { seed_repeat_loci[i].first = seed_end; }
            if (!((flank_start == seed_start) && (seed_repeat_loci[i].first == seed_end))) {
                skip_atomicity.clear();
                processLargeMotifSeed(tuple<int, int> { flank_start, seed_repeat_loci[i].first }, chunk_start, motif_length,
                                      seed_type, sequence_id, sequence, sequence_length, xor_bset, left_bset, right_bset,
                                      N_bset, out, lshift_xor_bsets, lshift_anchored_bsets, MATRIX, aligner, filter, alignment,
                                      repeat_loci, skip_atomicity);
            }
        }

        flank_start = seed_repeat_loci[i].second;
    }

    if (seed_end - flank_start >= MINIMUM_LENGTH[atomicity]) {
        if (flank_start < seed_start) { flank_start = seed_start; }
        if (flank_start != seed_start) {
            skip_atomicity.clear();
            processLargeMotifSeed(tuple<int, int> { flank_start, seed_end }, chunk_start, motif_length, seed_type, sequence_id, sequence,
                                  sequence_length, xor_bset, left_bset, right_bset, N_bset, out, lshift_xor_bsets, lshift_anchored_bsets,
                                  MATRIX, aligner, filter, alignment, repeat_loci, skip_atomicity);
        }
    }
}
