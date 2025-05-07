#include "parse_seed.h"

using namespace std;
using namespace boost;

using namespace boost::multiprecision;


int countBitsDiagonal(int row, int col, int diagonal_length, vector<boost::dynamic_bitset<>*> &MATRIX,
                      int &sequence_length, int neighbor_flank=0) {
    /*
     * count the number of ones across a diagonal
     * @param row index of the row in the sequence matrix
     * @param col index of the column in the sequence matrix
     * @param size length of the diagonal
     * @param MATRIX dot matrix of the chromosome sequence
     * @param sequence_length length of the chromosome sequence
     * @return int position within the sequence that is the possible motif for the seed
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


double calculateArrayMean(int arr[], int size) {
    /*
     * calculating the mean for array of integers
     * @param arr array of integers
     * @param size size of the array
     * @return double the average of all the integer values
    */
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += arr[i];
    }
    return static_cast<double>(sum) / size; // Convert sum to double before division
}


uint256_t mostFrequentLongMotif(boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset, int &seed_start,
                                  int &seed_sequence_length, int &motif_length, int &sequence_length, vector<boost::dynamic_bitset<>*> &MATRIX) {
    /*
     * return the index in the seed which is possible motif of the repeat
     * @param left_bset the bitset of left bit of a sequence
     * @param right_bset the bitset of right bit of a sequence
     * @param seed_start starting position of seed in the sequence
     * @param seed_sequence_length length of the seed sequence
     * @param motif_length length of the motif
     * @param sequence_length length of the chromosome sequence
     * @param MATRIX dot matrix of the chromosome sequence
     * @return int position within the sequence that is the possible motif for the seed
    */

    int seed_end = seed_start + seed_sequence_length;

    int mmotif_index = 0, max_count = 0, row_count = 0;
    int row, row_offset, cix;

    int dstream_index, ustream_index;
    int max_dindex, max_dcount;
    int dcount;

    int initial_lastrow, prefix_rows, pcindex;

    // outer loop for rows
    for (int row_start = seed_start; row_start < seed_end - motif_length + 1; row_start++) {
        row_count = 0;
        int iterations = 0;

        dstream_index = row_start + motif_length;
        while (dstream_index < seed_end) {
            max_dindex = -2, max_dcount = 0;

            for (int x=-2; x < 3; x++) {
                dcount = 0;
                for (int i=0; i<motif_length; i++) {
                    if (dstream_index + x + i >= seed_end) break;
                    if ((*MATRIX[row_start + i])[(sequence_length-1) - (dstream_index + x + i)] == 1) dcount += 1;
                }

                if (dcount > max_dcount) { max_dcount = dcount; max_dindex = x; }
            }

            row_count += max_dcount;
            dstream_index += max_dindex;
            dstream_index += motif_length;
        }

        ustream_index = row_start - motif_length;
        while (ustream_index > seed_start) {
            max_dindex = -2, max_dcount = 0;

            for (int x=-2; x < 3; x++) {
                dcount = 0;
                for (int i=0; i<motif_length; i++) {
                    if (ustream_index + x + i < 0) break;
                    if ((*MATRIX[row_start + i])[(sequence_length - 1) - (ustream_index + x + i)] == 1) dcount += 1;
                }
                if (dcount > max_dcount) { max_dcount = dcount; max_dindex = x; }
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
            for (int x=-2; x < 3; x++) {
                dcount = 0;
                for (int i=0; i < prefix_rows; i++) {
                    if (pcindex + x - i >= seed_end || pcindex + x - i < seed_start) break;
                    iterations += 1;
                    if ((*MATRIX[initial_lastrow - i])[(sequence_length - 1) - (pcindex + x - i)] == 1) dcount += 1;
                }
                if (dcount > max_dcount) { max_dcount = dcount; max_dindex = x; }
            }

            row_count += max_dcount;
        }

        if (row_count > max_count) {
            max_count = row_count; mmotif_index = row_start;
        }
    }

    // return mmotif_index;
    uint256_t motif_unit, ONE = 1;
    for (int j = mmotif_index; j < mmotif_index+motif_length; j++) {
        motif_unit <<= 1;
        if (left_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;

        motif_unit <<= 1;
        if (right_bset[sequence_length -1 -j] == 1) motif_unit |= ONE;
    }

    return motif_unit;
}


void processSeed(tuple<int, int> seed_position, int seq_start, int &motif_length, int &seed_type, string &sequence_id, string &sequence,
                 int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                 boost::dynamic_bitset<> &N_bset, int &continuous_threshold, ostream &out, vector<boost::dynamic_bitset<>> &lshift_xor_bsets,
                 vector<boost::dynamic_bitset<>*> &MATRIX, StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter,
                 StripedSmithWaterman::Alignment &alignment, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     * processes the seed and finds all the repeats in the sequence
     * @param seed_position tuple with start and end position of the seed
     * @param motif_length length of the motif
     * @param sequence_id name of the sequence
     * @param sequence_length length of the complete sequence
     * @param xor_bset the XOR conversion of the bitset (also includes anchor bitset)
     * @param left_bset the dynamic bitset of the left bit of the sequence
     * @param right_bset the dynamic bitset of the right bit of the sequence
     * @param purity_cutoff_mode if the impurity levels are given by purity threshold or mismatches
     * @param minimum_length the minimum length cutoff for different motif sizes
     * @param purity_threshold the allowed minimum purity
     * @param mismatches_threshold the allowed maximum mismatches
     * @param perfect_units the minimum number of perfect units for different motif sizes
     * @param continuous_threshold minimum length of continuous stretch of 1s in the seed
     * @param out outfile
     * @param aligner the aligner object
     * @param filter the filter object
     * @param alignment the resultant alignment object
     * @param repeat_loci set of identified repeat loci to be compared for redundancy
     * @returns none prints out the repeat locations to the output file
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

    // length of the pseudo perfect sequence
    int ppr_length = seed_sequence_length + motif_length + ((1-PURITY_THRESHOLD)*seed_sequence_length);
    uint256_t motif_unit;

    motif_unit = mostFrequentLongMotif(left_bset, right_bset, seed_start, seed_sequence_length,
                                            motif_length, sequence_length, MATRIX);
    atomicity = calculateAtomicityLongMotif(motif_unit, motif_length);

    if (atomicity <= SMALL_MLEN_LIMIT) {
        processSeedMotifWise(tuple<int, int> { seed_start, seed_end }, seq_start, atomicity, seed_type, sequence_id, sequence,
                                sequence_length, lshift_xor_bsets[atomicity-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                continuous_threshold, out, aligner, filter, alignment, repeat_loci);
        return;
    }

    if (motif_length % atomicity != 0) { return; }

    // the repeat should be treated based on the atomicity
    motif = calculateMotif(motif_unit, motif_length);
    motif = motif.substr(0, atomicity);
    motif_unit >>= 2*(motif_length - atomicity);

    perfect_repeat = "";
    while(perfect_repeat.length() <= ppr_length) perfect_repeat += motif;
   
    aligner.Align(seed_sequence.c_str(), perfect_repeat.c_str(), ppr_length, filter, &alignment, 15);
    processCIGARWithPruning(seed_start, seed_sequence_length, alignment.cigar_string, seed_sequence, atomicity,
                            repeat_start, repeat_end, alignment_length, match_units, cigar_string, purity,
                            motifwise_purity, motifwise_indels);

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

        if ((   (atomicity < 10  && (match_units >= PERFECT_UNITS[atomicity] || (purity > 0.9 && purity*repeat_length >= 2*atomicity))) 
             || ((atomicity >= 10) && (((purity * repeat_length) >= 3*atomicity) || (purity > 0.9 && purity*repeat_length >= 2*atomicity))))
            && atomicity >= MINIMUM_MLEN && atomicity <= MAXIMUM_MLEN
            && repeat_length >= MINIMUM_LENGTH[atomicity] 
            && purity >= PURITY_THRESHOLD 
            && motifwise_purity >= MOTIFPURITY_THRESHOLD) {

            addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                             atomicity, repeat_length, repeat_units, out, repeat_loci);
        }
    }

    if (seed_repeat_loci.size()==0) { return; }

    int flank_start = seed_start;
    for (int i=0; i<seed_repeat_loci.size(); i++) {

        if (flank_start >= seed_repeat_loci[i].first) { flank_start = seed_repeat_loci[i].second; continue;  }

        if (seed_repeat_loci[i].first - flank_start >= MINIMUM_LENGTH[motif_length]) {
            if (flank_start < seed_start) { flank_start = seed_start; }
            if (seed_repeat_loci[i].first > seed_end) { seed_repeat_loci[i].first = seed_end; }
            if (!((flank_start == seed_start) && (seed_repeat_loci[i].first == seed_end))) {
                processSeed(tuple<int, int> { flank_start, seed_repeat_loci[i].first }, seq_start, motif_length, seed_type, sequence_id,
                            sequence, sequence_length, xor_bset, left_bset, right_bset, N_bset, continuous_threshold, out,
                            lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci);
            }
        }

        flank_start = seed_repeat_loci[i].second;
    }

    if (seed_end - flank_start >= MINIMUM_LENGTH[motif_length]) {
        if (flank_start < seed_start) { flank_start = seed_start; }
        if (flank_start != seed_start) {
            processSeed(tuple<int, int> { flank_start, seed_end }, seq_start, motif_length, seed_type, sequence_id, sequence,
                        sequence_length, xor_bset, left_bset, right_bset, N_bset, continuous_threshold, out, lshift_xor_bsets,
                        MATRIX, aligner, filter, alignment, repeat_loci);
        }
    }
}
