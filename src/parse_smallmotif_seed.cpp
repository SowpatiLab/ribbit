/*
 * Different methods for parsing motif_length XOR and identification of tandem repeats
*/

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <mutex>
#include <boost/dynamic_bitset.hpp>

#include <cstdint>
#include <numeric>
#include <iomanip>
#include <boost/multiprecision/cpp_int.hpp>

#include "ssw_cpp.h"

#include "global_variables.h"
#include "bitseq_utils.h"
#include "process_cigar.h"
#include "output_utils.h"
#include "parse_seed.h"
#include "parse_smallmotif_seed.h"

using namespace std;
using namespace boost;


int calculateMotifUnits(boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset, int &start,
                        int &length, int &motif_length, int &sequence_length, uint32_t motif_unit) {
    /*
     * calculates the number of perfect motif units in a repeat sequence
     * @param left_bset the dynamic bitset of the left bit of the sequence
     * @param right_bset the dynamic bitset of the right bit of the sequence
     * @param start start coordinate of the sequence
     * @param length length of the sequence to be looked at
     * @param motif_length length of the motif
     * @param sequence_length total length of the sequence
     * @param motif_unit the motif unit that is being repeated
     * @returns int number of perfect motif units repeated
    */

    unordered_map<uint32_t, int> motif_position, motif_units;
    unordered_map<uint32_t, int> maxfrequency_motifs;
    uint32_t motif;
    int seed_end = start + length;
    if (seed_end > sequence_length - 1) { seed_end = sequence_length - 1; }

    boost::dynamic_bitset<> window(2*motif_length, 0ull); // window to track the motif
    for (int j = start; j < seed_end; j++) {
        window[0] = right_bset[sequence_length -1 -j]; window[1] = left_bset[sequence_length -1 -j];

        if (j-start >= (0.9*motif_length)-1) {   // window is atleast the size of motif length
            motif = calculateRepeatClass(window, motif_length);

            if (motif_position.find(motif) == motif_position.end()) {
                // if the motif is not tracked for its position
                // motif position is position of the first nucleotide in the motif
                motif_position[motif] = j - (motif_length - 1);
                motif_units[motif] = 1;
            }
            else {
                if ((j - (motif_length - 1)) - motif_position[motif] >= motif_length) {
                    // if the present position of motif occurence is atleast motif length away from previous occurrence
                    motif_position[motif] = j - (motif_length - 1);     // update motif position
                    motif_units[motif] += 1;                            // increase the number of units for motif
                }
            }
        }

        window <<= 2;
    }

    return motif_units[motif_unit];
}


void possibleMotifs(boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset, int &seed_start,
                    int &seed_sequence_length, int &motif_length, int &sequence_length,
                    vector<uint32_t> &motifs, vector<int> &starts, vector<int> &ends, string &sequence) {
    /*
     * finding the most repeating motif without converting the seed to string; applying the KMP algorithm
     * @param left_bset the dynamic bitset of the left bit of the sequence
     * @param right_bset the dynamic bitset of the right bit of the sequence
     * @param seed_start start of the seed sequence
     * @param seed_sequence_length length of the seed sequence
     * @param motif_length length of the motif
     * @param sequence_length total length of the sequence
     * @param motifs the vector of identified possible motifs; passed as reference; updated
     * @param starts the vector of starts of the identified motifs
     * @param ends the vector of the ends of the identified motifs
     * @param sequence the nucleotide sequence of the seed
     * @returns none updates motifs, starts and ends of the identified motifs
    */

    unordered_map<uint32_t, int> new_motif_start;

    uint32_t motif; int wstart, wend;
    int seed_end = seed_start + seed_sequence_length;
    if (seed_end > sequence_length - 1) { seed_end = sequence_length - 1; }
    int lc_motif_length = motif_length;

    boost::dynamic_bitset<> window(2*motif_length, 0ull); // window to track the motif
    for (int j = seed_start; j < seed_end; j++) {
        window[0] = right_bset[sequence_length -1 -j]; window[1] = left_bset[sequence_length -1 -j];
        motif = calculateRepeatClass(window, motif_length);
        wstart = j - (motif_length - 1);
        wend = j + 1;

        if (j-seed_start >= (0.9*motif_length)-1) {   // window is atleast the size of motif length

            if (new_motif_start.find(motif) == new_motif_start.end()) {
                // if the motif is not tracked for its position
                // motif position is position of the first nucleotide in the motif
                new_motif_start[motif] = wstart;
                MOTIF_START[motif] = wstart;
                MOTIF_END[motif] = wend;
                MOTIF_UNITS[motif] = 1;
                MOTIF_GAPS[motif] = 0;
                MOTIF_GAPSIZE[motif] = 0;
                MOTIF_NEXT[motif] = ((window << 2) | (window >> (motif_length-1)*2)).to_ulong();

            }

            else {
                if (wstart - MOTIF_END[motif] > 3*motif_length) {
                    // if the new position of the motif is beyond three motif lengths of the old
                    if (MOTIF_END[motif] - MOTIF_START[motif] >= MINIMUM_LENGTH[motif_length] && MOTIF_UNITS[motif] >= PERFECT_UNITS[motif_length]) {
                        // check if the previous repeat is of valid length and
                        // if (MOTIF_GAPS[motif] < (MOTIF_UNITS[motif]/2 + 1)) {
                            motifs.push_back(motif);
                            starts.push_back(MOTIF_START[motif]);
                            ends.push_back(MOTIF_END[motif]);
                        // }
                    }

                    // reinitialise all the values
                    MOTIF_START[motif] = wstart;
                    MOTIF_END[motif] = wend;
                    MOTIF_UNITS[motif] = 1;
                    MOTIF_GAPS[motif] = 0;
                    MOTIF_GAPSIZE[motif] = 0;
                    MOTIF_NEXT[motif] = ((window << 2) | (window >> (motif_length-1)*2)).to_ulong();
                    new_motif_start[motif] = wstart;
                }

                else {
                    // if the motif is not occurring consecutively
                    if (MOTIF_END[motif] < j) {
                        if (j - MOTIF_END[motif] < motif_length) {
                            MOTIF_GAPS[motif] += 1;
                            MOTIF_GAPSIZE[motif] += 1;
                        }
                        else if ((j - MOTIF_END[motif]) % motif_length > 0) {
                            MOTIF_GAPS[motif] += ((j - MOTIF_END[motif]) / motif_length) + 1;
                            MOTIF_GAPSIZE[motif] += (j - MOTIF_END[motif]) + 1;
                        }
                        else {
                            MOTIF_GAPS[motif] += ((j - MOTIF_END[motif]) / motif_length);
                            MOTIF_GAPSIZE[motif] += (j - MOTIF_END[motif]);
                        }
                    }
                    else if (MOTIF_END[motif] == j && MOTIF_NEXT[motif] != window.to_ulong()) {
                        MOTIF_GAPS[motif] += 1;
                        MOTIF_GAPSIZE[motif] += 1;
                    }

                    if (wstart - new_motif_start[motif] >= motif_length) {
                        // if the present position of motif occurence is atleast motif length away from previous occurrence
                        new_motif_start[motif] = wstart;     // update motif position
                        MOTIF_UNITS[motif] += 1;             // increase the number of units for motif
                    }
                    MOTIF_END[motif] = wend;
                    MOTIF_NEXT[motif] = ((window << 2) | (window >> (motif_length-1)*2)).to_ulong();
                }
            }
        }

        window <<= 2;
    }

    for (auto& it: new_motif_start) {
        // reiterate through all the left over motifs and record them
        motif = it.first;
        if (MOTIF_END[motif] - MOTIF_START[motif] >= MINIMUM_LENGTH[motif_length] && MOTIF_UNITS[motif] >= PERFECT_UNITS[motif_length]) {
            motifs.push_back(motif);
            starts.push_back(MOTIF_START[motif]);
            ends.push_back(MOTIF_END[motif]);
        }
    }
}


void processSeedMotifWise(tuple<int, int> seed_position, int seq_start, int &motif_length, int &seed_type, string &sequence_id, string &sequence,
                          int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                          boost::dynamic_bitset<> &N_bset, int &continuous_threshold, ostream &out,
                          StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment,
                          vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     * processes the seed and finds all the repeats in the sequence
     * @param seed_position tuple with start and end position of the seed
     * @param seq_start the start of the seed sequence
     * @param motif_length length of the motif
     * @param seed_type type of the seed
     * @param sequence_id name of the sequence
     * @param sequence nucleotide sequence of the seed
     * @param sequence_length length of the complete sequence
     * @param xor_bset shift XOR bitset of the motif size
     * @param left_bset the dynamic bitset of the left bit of the sequence
     * @param right_bset the dynamic bitset of the right bit of the sequence
     * @param N_bset the bitset indicating the presence of Ns at a position
     * @param continuous_threshold minimum length of continuous stretch of 1s in the seed
     * @param out output file name
     * @param aligner the aligner object pf ssw alignment
     * @param filter the filter object of ssw alignment
     * @param alignment the resultant alignment object
     * @param repeat_loci the list of identified repeat loci
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
        seed_bset[seed_end - 1 - j] = xor_bset[sequence_length -j - 1];
    }

    int longest_stretch = longestContinuousMatches(seed_bset);
    if (longest_stretch < continuous_threshold) { return; }
    if (THREADS > 1) MTX.lock();
    vector<uint32_t> motifs; vector<int> starts, ends;
    possibleMotifs(left_bset, right_bset, seed_start, seed_sequence_length,
                   motif_length, sequence_length, motifs, starts, ends, sequence);

    if (THREADS > 1) MTX.unlock();
    if (motifs.size() == 0) return;

    string perfect_repeat, motif;
    vector<int> cigar_values;
    int ppr_length;

    int repeat_start, repeat_end, match_nucs, mismatch_nucs, match_units;
    int repeat_length, repeat_units;
    int alignment_length, interruptions, atomicity;
    int motif_seed_length, motifwise_indels, avg_matchlen;
    double purity = 0, motifwise_purity = 0;
    string cigar_string, motif_seed_sequence;


    int motif_idx; uint32_t motif_unit;
    for(motif_idx=0; motif_idx < motifs.size(); motif_idx++) {
        motif_unit = motifs[motif_idx];
        atomicity = calculateAtomicity(motif_unit, motif_length);

        // the repeat should be treated based on the atomicity
        motif = calculateMotif(motif_unit, motif_length);
        motif = motif.substr(0, atomicity);

        motif_unit >>= 2*(motif_length - atomicity);
        // seed sequence limiting to the coordinates where full motif alignment matches are found 
        motif_seed_sequence = sequence.substr(starts[motif_idx], ends[motif_idx] - starts[motif_idx]);
        motif_seed_length = ends[motif_idx] - starts[motif_idx];

        ppr_length = ends[motif_idx] - starts[motif_idx] + motif_length + ((1-PURITY_THRESHOLD)*(ends[motif_idx] - starts[motif_idx]));
        perfect_repeat = "";
        while(perfect_repeat.length() <= ppr_length) perfect_repeat += motif;

        aligner.Align(seed_sequence.c_str(), perfect_repeat.c_str(), ppr_length, filter, &alignment, 15);
        processCIGARMotifWise(starts[motif_idx], motif_seed_length, alignment.cigar_string, motif_seed_sequence, atomicity,
                              repeat_start, repeat_end, alignment_length, cigar_string, purity, motifwise_purity, motifwise_indels, avg_matchlen);
        
        if (THREADS > 1) MTX.lock();
        match_units = calculateMotifUnits(left_bset, right_bset, repeat_start, repeat_length, atomicity, sequence_length, motif_unit);
        if (THREADS > 1) MTX.unlock();

        repeat_length = repeat_end - repeat_start;
        repeat_units = repeat_length/atomicity;

        if ((match_units >= PERFECT_UNITS[atomicity] && match_units >= (0.7*repeat_units))
            && (motifwise_purity >= MOTIFPURITY_THRESHOLD || avg_matchlen > 2*atomicity)
            && repeat_length >= MINIMUM_LENGTH[atomicity]
            && atomicity >= MINIMUM_MLEN && atomicity <= MAXIMUM_MLEN
            && (!(repeat_units < 3 && purity < PURITY_THRESHOLD))) {
            // a small motif seed is considered valid based on a set of criteria
            // - match units are more than threshold AND 70% of the total units are perfect
            // - average motif purity is 80% OR the average continuous match length twice the atomicity

            addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                             atomicity, repeat_length, repeat_units, out, repeat_loci);
        }
    }
}
