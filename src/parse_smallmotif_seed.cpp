#include "parse_smallmotif_seed.h"

using namespace std;
using namespace boost;


tuple<string,int> getMostCommonKmer(const string &sequence, int k) {
    /*
     *  finds the most common kmer in a sequence
     *  @param sequence the input DNA sequence
     *  @param k length of the kmer
     *  @return a tuple containing the most common kmer and its count
     */
    if (k <= 0 || sequence.size() < (size_t)k) return make_tuple(string(""), 0);

    unordered_map<string, int> count;
    unordered_map<string, size_t> last_pos;  // last accepted start position

    for (size_t i = 0; i + k <= sequence.size(); ++i) {
        string kmer = sequence.substr(i, k);

        // Check if previous accepted occurrence overlaps
        if (last_pos.find(kmer) == last_pos.end() || i >= last_pos[kmer] + k) {
            count[kmer]++;
            last_pos[kmer] = i;  // mark this position
        }
    }

    int max_count = 0;
    string max_kmer;
    for (auto &p : count) {
        if (p.second > max_count) {
            max_count = p.second;
            max_kmer = p.first;
        }
    }
    return make_tuple(max_kmer, max_count);
}


bool qualifyShortMotifRepeat(string &sequence, int kmer) {

    /*
     *  checks if the short motif repeat qualifies the minimum criteria to be reported
     *  @param sequence the sequence of the repeat
     *  @param kmer length of the motif
     *  @return bool if the repeat qualifies the minimum criteria to be reported
     */

    int seq_length = sequence.length();

    tuple<string,int> most_common_kmer = getMostCommonKmer(sequence, kmer);
    string motif = get<0>(most_common_kmer);
    int count = get<1>(most_common_kmer);

    if (count < PERFECT_UNITS[kmer]) return false;

    return true;
}


int calculateMotifUnits(boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset, int &seed_start,
                        int &length, int &motif_length, int &sequence_length, uint32_t motif_unit) {
    /*
     *  calculates the number of perfect motif units in a repeat sequence
     *  @param left_bset the dynamic bitset of the left bit of the sequence
     *  @param right_bset the dynamic bitset of the right bit of the sequence
     *  @param start start coordinate of the sequence
     *  @param length length of the sequence to be looked at
     *  @param motif_length length of the motif
     *  @param sequence_length total length of the sequence
     *  @param motif_unit the motif unit that is being repeated
     *  @returns int number of perfect motif units repeated
     */

    unordered_map<uint32_t, int> motif_position, motif_units;
    unordered_map<uint32_t, int> maxfrequency_motifs;
    uint32_t motif = 0ull;
    int seed_end = seed_start + length;
    if (seed_end > sequence_length - 1) { seed_end = sequence_length - 1; }

    boost::dynamic_bitset<> window(2*motif_length, 0ull); // window to track the motif
    for (int j = seed_start; j < seed_end; j++) {
        window[0] = right_bset[sequence_length -1 -j]; window[1] = left_bset[sequence_length -1 -j];

        if (j-seed_start >= motif_length-1) {   // window is atleast the size of motif length
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
     *  finding the most repeating motif without converting the seed to string; applying the KMP algorithm
     *  @param left_bset the dynamic bitset of the left bit of the sequence
     *  @param right_bset the dynamic bitset of the right bit of the sequence
     *  @param seed_start start of the seed sequence
     *  @param seed_sequence_length length of the seed sequence
     *  @param motif_length length of the motif
     *  @param sequence_length total length of the sequence
     *  @param motifs the vector of identified possible motifs; passed as reference; updated
     *  @param starts the vector of starts of the identified motifs
     *  @param ends the vector of the ends of the identified motifs
     *  @param sequence the nucleotide sequence of the seed
     *  @returns none updates motifs, starts and ends of the identified motifs
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

        if (j-seed_start >= motif_length-1) {   // window is atleast the size of motif length

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
                    if (MOTIF_END[motif] - MOTIF_START[motif] >= MINIMUM_LENGTH[motif_length]
                        && MOTIF_UNITS[motif] >= PERFECT_UNITS[motif_length]) {
                        // check if the previous repeat is of valid length and
                            motifs.push_back(motif);
                            starts.push_back(MOTIF_START[motif]);
                            ends.push_back(MOTIF_END[motif]);
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
                    // if the motif end is j the motif is coming bookended
                    // if the motif is not occurring consecutively
                    if (MOTIF_END[motif] < j) {
                        if (j - MOTIF_END[motif] < motif_length) {
                            MOTIF_GAPS[motif] += 1;
                            MOTIF_GAPSIZE[motif] += 1;
                            if (motif_length > 2 && (motif_length - (j - MOTIF_END[motif])) == 1) {
                                // for motifs longer than 2, if the gap is less than motif length
                                MOTIF_UNITS[motif] += 1;
                            }
                        }
                        else if ((wstart - MOTIF_END[motif]) % motif_length > 0) {
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
        if (MOTIF_END[motif] - MOTIF_START[motif] >= MINIMUM_LENGTH[motif_length]
            && MOTIF_UNITS[motif] >= PERFECT_UNITS[motif_length]) {
            motifs.push_back(motif);
            starts.push_back(MOTIF_START[motif]);
            ends.push_back(MOTIF_END[motif]);
        }
    }
}


void orderMotifs(vector<uint32_t> &motifs, vector<int> &starts, vector<int> &ends, int seed_start, int seed_end, int motif_length) {
    /*
     *  orders the motifs and adjusts the starts and ends such that the whole seed is covered
     *  @param motifs the vector of identified possible motifs; passed as reference; updated
     *  @param starts the vector of starts of the identified motifs
     *  @param ends the vector of the ends of the identified motifs
     *  @param seed_start start coordinate of the seed sequence
     *  @param seed_end end coordinate of the seed sequence
     *  @param motif_length length of the motif
     *  @returns none updates motifs, starts and ends of the identified motifs
     */

    seed_end = seed_end + motif_length; // extend the seed end by motif length to cover the last motif
    int motif_start = 0, motif_end = 0;
    int upstream_start = 0, upstream_end = 0;
    int downstream_start = 0, downstream_end = 0;
    vector<int> new_starts(motifs.size()), new_ends(motifs.size()); int d,u;
    vector<int> sorted_indices(ends.size());
    iota(sorted_indices.begin(), sorted_indices.end(), 0);
    sort(sorted_indices.begin(), sorted_indices.end(), [&ends](int i, int j) { return ends[i] < ends[j]; });
    int idx; int min_start = seed_end, min_idx = 0;
    for (int _=sorted_indices.size()-1; _ >= 0; _--) {
        idx = sorted_indices[_];
        motif_start = starts[idx];
        motif_end = ends[idx];
        // if the motif is the last one in the list or the motif ends one motif length away from the seed end
        if (_ == sorted_indices.size()-1 || (seed_end - motif_end <= motif_length)) { motif_end = seed_end; }

        // if the motif is the first one in the list or the motif starts one motif length away from the seed start
        if (motif_start - seed_start <= motif_length) { motif_start = seed_start; }

        if (motif_start < min_start)  { min_start = motif_start; min_idx = idx; }

        d = 1;
        while ( _ + d < sorted_indices.size()) {
            downstream_start = starts[sorted_indices[_ + d]]; downstream_end = ends[sorted_indices[_ + d]];
            if (motif_start > downstream_start && motif_end <= downstream_end) {
                break;
            }
            else if (downstream_start >= motif_start && downstream_end <= motif_end) {
                d = d + 1;
            }
            else {
                if (motif_end - downstream_start < motif_length) {
                    motif_end = downstream_start + motif_length;
                    if (motif_end > seed_end) motif_end = seed_end;
                }
                break;
            }
        }

        u = 1;
        while ( _ - u >= 0) {
            upstream_start = starts[sorted_indices[_ - u]]; upstream_end = ends[sorted_indices[_ - u]];
            if (motif_start > upstream_start && motif_end <= upstream_end) {
                break;
            }
            else if (upstream_start >= motif_start && upstream_end <= motif_end) {
                u = u + 1;
            }
            else {
                if (upstream_end - motif_start < motif_length) {
                    motif_start = upstream_end - motif_length;
                    if (motif_start < seed_start) motif_start = seed_start;
                }
                break;
            }
        }
        new_starts[idx] = motif_start;
        new_ends[idx] = motif_end;
    }

    for (int _=0; _ < motifs.size(); _++) {
        starts[_] = new_starts[_];
        ends[_] = new_ends[_];
    }
    starts[min_idx] = seed_start;
}


void processSmallMotifSeed(tuple<int, int> seed_position, int chunk_start, int &motif_length, int &seed_type, string &sequence_id,
                           string &sequence, int &sequence_length, boost::dynamic_bitset<> &xor_bset, boost::dynamic_bitset<> &left_bset,
                           boost::dynamic_bitset<> &right_bset, boost::dynamic_bitset<> &N_bset, ostream* out,
                           StripedSmithWaterman::Aligner &aligner, StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment,
                           vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     *  processes the seed and finds all the repeats in the sequence
     *  @param seed_position tuple with start and end position of the seed
     *  @param chunk_start the start of the seed sequence
     *  @param motif_length length of the motif
     *  @param seed_type type of the seed
     *  @param sequence_id name of the sequence
     *  @param sequence nucleotide sequence of the seed
     *  @param sequence_length length of the complete sequence
     *  @param xor_bset shift XOR bitset of the motif size
     *  @param left_bset the dynamic bitset of the left bit of the sequence
     *  @param right_bset the dynamic bitset of the right bit of the sequence
     *  @param N_bset the bitset indicating the presence of Ns at a position
     *  @param out output file name
     *  @param aligner the aligner object pf ssw alignment
     *  @param filter the filter object of ssw alignment
     *  @param alignment the resultant alignment object
     *  @param repeat_loci the list of identified repeat loci
     *  @returns none prints out the repeat locations to the output file
     */

    int seed_start     = get<0> (seed_position);
    int seed_end       = get<1> (seed_position);

    // check if there are Ns in the seed sequence
    // if yes, break the seed around the N position and consider the parts
    for (int s=seed_start; s < seed_end+motif_length; s++) {
        if (N_bset[sequence_length-1-s] == 1) {
            if (s - seed_start >= motif_length) {
                processSmallMotifSeed(tuple<int,int>{seed_start, s - motif_length}, chunk_start, motif_length, seed_type,
                                      sequence_id, sequence, sequence_length, xor_bset, left_bset, right_bset, N_bset, out,
                                      aligner, filter, alignment, repeat_loci);
            }
            seed_start = s - motif_length + 1;
        }
    }

    int seed_bset_size = seed_end - seed_start;
    int seed_sequence_length = seed_bset_size + motif_length;
    string repeat_sequence = "";

    // the shift xor bitset of the complete repeat sequence
    boost::dynamic_bitset<> seed_bset(seed_bset_size, 0ull);
    for (int j = seed_start; j < seed_end; j++) {
        seed_bset[seed_end - 1 - j] = xor_bset[sequence_length -j - 1];
    }

    int continuous_threshold = 3; // threshold for continuous matches
    int longest_stretch = longestContinuousMatches(seed_bset);
    if (longest_stretch < continuous_threshold) { return; }

    if (THREADS > 1) MTX.lock();
    vector<uint32_t> motifs; vector<int> starts, ends;

    if (seed_type == RANK_P) {
        uint32_t motif = 0ull;
        for (int _=0; _ < motif_length; _++) {
            motif <<= 1; motif |= left_bset[sequence_length - 1 - (seed_start + _)];
            motif <<= 1; motif |= right_bset[sequence_length - 1 - (seed_start + _)];
        }
        dynamic_bitset<> window(2*motif_length, motif); // window to track the motif
        motifs.push_back(calculateRepeatClass(window, motif_length));
        starts.push_back(seed_start); ends.push_back(seed_end + motif_length);
    }
    else {
        possibleMotifs(left_bset, right_bset, seed_start, seed_sequence_length,
                       motif_length, sequence_length, motifs, starts, ends, sequence);
    }

    if (THREADS > 1) MTX.unlock();

    // if no motifs are found in the seed. No processing the seed further
    if (motifs.size() == 0) return;

    // if only one motif is found, extend the start and end for the motif to cover the whole seed
    if (motifs.size() == 1) {
        starts[0] = seed_start;ends[0] = seed_end + motif_length;
    }
    else if (motifs.size() > 1) {
        // if multiple motifs are found in the seed, adjust that starts and ends of the motifs
        // to cover the whole seed
        orderMotifs(motifs, starts, ends, seed_start, seed_end, motif_length);
    }

    string perfect_repeat, motif;
    vector<int> cigar_values;
    int ppr_length;

    int repeat_start, repeat_end, match_nucs, mismatch_nucs, match_units;
    int repeat_length, repeat_units;
    int alignment_length, substitutions, indels, atomicity;
    int motif_seed_length, motifwise_indels, avg_matchlen;
    double purity = 0, motifwise_purity = 0;
    string cigar_string, motif_seed_sequence;

    int motif_idx; uint32_t motif_unit;
    for(motif_idx=0; motif_idx < motifs.size(); motif_idx++) {
        motif_unit = motifs[motif_idx];

        if (THREADS > 1) MTX.lock();
        atomicity = calculateAtomicity(motif_unit, motif_length);
        motif = calculateMotif(motif_unit, motif_length);
        if (THREADS > 1) MTX.unlock();

        motif = motif.substr(0, atomicity);
        motif_unit >>= 2*(motif_length - atomicity);

        // seed sequence limiting to the coordinates where full motif alignment matches are found 
        motif_seed_sequence = sequence.substr(starts[motif_idx], ends[motif_idx] - starts[motif_idx]);
        motif_seed_length   = ends[motif_idx] - starts[motif_idx];

        ppr_length = motif_seed_length + motif_length + ((1-PURITY_THRESHOLD)* motif_seed_length);
        perfect_repeat = "";
        while(perfect_repeat.length() <= ppr_length) perfect_repeat += motif;

        string intermediate_cigar = "";
        if (motif_seed_sequence.length() < 5000) {
            aligner.Align(motif_seed_sequence.c_str(), perfect_repeat.c_str(), ppr_length, filter, &alignment, 15);
            intermediate_cigar = alignment.cigar_string;
        }
        else {
            intermediate_cigar = alignLargeSequence(motif_seed_sequence, motif, motif_length, aligner, filter, alignment);
        }
        processCIGARMotifWise(starts[motif_idx], motif_seed_length, intermediate_cigar, motif_seed_sequence, atomicity,
                              repeat_start, repeat_end, alignment_length, cigar_string, purity, substitutions, indels,
                              motifwise_purity, motifwise_indels, avg_matchlen);
        repeat_length = repeat_end - repeat_start;

        if (THREADS > 1) MTX.lock();
        match_units = calculateMotifUnits(left_bset, right_bset, repeat_start, repeat_length, atomicity, sequence_length, motif_unit);
        if (THREADS > 1) MTX.unlock();

        repeat_units = repeat_length/atomicity;

        // conditions used to cover some edge cases
        // if match units are more than 10 and the number of interruptions is less than 80% of the match units
        if (match_units > 10 && (indels > 0.8*match_units)) { continue; }
        if (repeat_length < 3*atomicity && purity < 1) { continue; }

        if (atomicity >= MINIMUM_MLEN && atomicity <= MAXIMUM_MLEN
            && (match_units >= PERFECT_UNITS[atomicity])
            && (repeat_length >= MINIMUM_LENGTH[atomicity])
            && (motifwise_purity >= MOTIFPURITY_THRESHOLD || avg_matchlen >= 2*atomicity)
            && (purity >= PURITY_THRESHOLD)) {
            // a small motif seed is considered valid based on a set of criteria
            // - match units are more than threshold AND 70% of the total units are perfect
            // - average motif purity is 80% OR the average continuous match length twice the atomicity

            repeat_sequence = sequence.substr(repeat_start, repeat_length);
            if (!qualifyShortMotifRepeat(repeat_sequence, atomicity)) { continue; }

            repeat_start += chunk_start; repeat_end += chunk_start;
            int recursion_level = 0;
            vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
            if (out) {
                addLocusToOutput(sequence_id, repeat_start, repeat_end, motif.substr(0, atomicity), purity, cigar_string,
                                 atomicity, repeat_length, repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
            }
        }
    }
}
