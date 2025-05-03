#include "process_cigar.h"

using namespace std;
using namespace boost;


double calculateMedian(std::vector<double>& numbers) {
    /*
     *  calculates the median of list of double numbers
     *  @param numbers vector of double numbers
     *  @return double the median value of the numbers
    */
    int size = numbers.size();
    std::sort(numbers.begin(), numbers.end());

    if (size % 2 == 0) {
        return (numbers[size / 2 - 1] + numbers[size / 2]) / 2.0;
    } else {
        return numbers[size / 2];
    }
}


tuple<vector<int>, vector<char>> cigarSplit(string cigar) {
    /*
     *  splits a CIGAR string into consecutive operations and lengths
     *  @param cigar character pointer of the CIGAR string
     *  @return tuple<vector<int>, vector<char>> tuple of two vectors cigar lengths and cigar types
    */

    string length = "";
    vector<int> clens; vector<char> ctypes;

    for (int i = 0; i<cigar.length(); i++) {
        if (isdigit(cigar[i])) {
            length += cigar[i];
        }
        else {
            clens.push_back(stoi(length));
            if (cigar[i] == '=') { ctypes.push_back('M'); }
            else { ctypes.push_back(cigar[i]); }
            length = "";
        }
    }

    return tuple<vector<int>, vector<char>> {clens, ctypes};
}


pair<int, int> calculateTrimEdges(double &purity_threshold, double &purity, vector<int> &ccigar_lengths,
                                  int &alignment_length, int &motif_length) {
    /*
     *  calculates the lengths of trims on either sides based on threshold number of mismatches and indels
     *  @param purity_threshold threshold value for the total repeat purity
     *  @param purity the repeat purity value; passed as reference
     *  @param ccigar_lengths the lengths of cigar operations; has negative lengths for mismatches
     *  @param alignment_length the total length of alignment with repeat and the perfect repeat; passed as reference
     *  @param motif_length the length of the repeating motif
     *  return the pair of trim lengths from left and right
    */
    int trim_length = 0;        // length of the trim
    pair<int, int> trim_edges;  // the final pair of trim lengths

    // parameters for a certain pair of trim lengths
    int pair_match, pair_alignment; double pair_purity;
    // keeping track of the maximum purity and alignment length for all the trim combinations
    double max_purity = 0; int max_alength = 0;


    // iteratively trim till the purity threshold is reached
    while (purity < purity_threshold) {
        trim_length += 1;

        // clear everything before moving further to next trim lengths
        max_purity = 0; max_alength = 0;

        // building all trim combinations based on the trim length
        for (int i=0; i<=trim_length; i++) {
            pair_match = 0, pair_alignment = 0;

            // clength vector has the lengths of all alternative match and mismatch operations
            // starting with a match, hence we only consider even indices for matches
            for (int j=2*i; j <= (ccigar_lengths.size()-1)-(2*(trim_length-i)); j++) {
                if (j%2 == 0) pair_match += ccigar_lengths[j];    // for even indices increase the repeat length
                pair_alignment += ccigar_lengths[j];
            }
            pair_purity = (double)pair_match/(double)pair_alignment;

            // among all the trim combinations that pass the purity threshold
            // we take the one with the highest alignment length
            if (pair_purity >= purity_threshold) {
                if (max_alength < pair_alignment) {
                    max_purity = pair_purity; max_alength = pair_alignment;
                    trim_edges = {i, trim_length-i};
                }
            }
        }

        // consider the trim only if the purity has increased from previous
        // otherwise move to the next trim with increased length
        if (max_purity > purity) {
            purity = max_purity; alignment_length = max_alength;
        }

        // break more iterations if the alignment length is less than minimum length
        if (alignment_length < MINIMUM_LENGTH[motif_length]) break;
    }

    return trim_edges;
}


void motifwiseParametersInTrim(vector<int> &clens, vector<char> &ctypes, int cidx_start, int cidx_end,
                               int motif_length, double &avg_motifpurity, int &alignment_length) {
    /*
     *  calculates the motif wise parameters including motifwise_purity and motifwise_indels when trim values are given
     *  @param clens the lengths of consecutive CIGAR operations
     *  @param ctypes the consecutive operations in CIGAR
     *  @param cidx_start the start of the CIGAR index to be considered after trimming
     *  @param cidx_end the send of the CIGAR index to be considered after trimming
     *  @param motif_length the length of the repeating motif
     *  @param avg_motifpurity the average motif wise purity across the repeat; passed as reference
     *  @param alignment_length length of the total alignment; passed as reference and updated
     *  return void
    */

    char ctype; int clength;
    int  cidx = cidx_start;
    alignment_length = 0;

    bool mismatch_continue = false;         // to check if there are contiguous mismatches

    int motif_covered = 0, motif_matches = 0, motif_mismatches = 0, motif_indels = 0;
    int excess = 0;
    vector<double> motifwise_matchpercent;

    for ( ; cidx <= cidx_end; cidx++) {
        clength = clens[cidx]; ctype = ctypes[cidx];

        switch (ctype) {
            case 'S':
                // soft clip: edit the repeat start and end
                break;

            case 'X':
                alignment_length += clength; mismatch_continue = true;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_mismatches += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motif_covered = excess, motif_matches = 0, motif_mismatches = excess, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_mismatches += clength;
                }
                break;
            case 'I':
                alignment_length += clength; mismatch_continue = true;

                motif_indels += clength;
                break;
            case 'D':
                mismatch_continue = true;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_indels += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motif_covered = excess, motif_matches = 0, motif_mismatches = 0, motif_indels = excess;
                }
                else {
                    motif_covered += clength;
                    motif_indels += clength;
                }
                break;
            case '=': case 'M':
                alignment_length += clength; mismatch_continue = false;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_matches += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motif_covered = excess, motif_matches = excess, motif_mismatches = 0, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_matches += clength;
                }
                break;
            default: break;
        }
    }

    if (motif_covered > 0) {
        motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
    }

    avg_motifpurity = 0;
    if (motifwise_matchpercent.size() > 0) {
        for (int _=0; _<motifwise_matchpercent.size(); _++) { avg_motifpurity += motifwise_matchpercent[_]; }
        avg_motifpurity = avg_motifpurity / ((double) (motifwise_matchpercent.size()));
    }
}


pair<int, int> calculateTrimEdgesMotifPurity(double &motifpurity_threshold, double &avg_motifpurity, vector<int> &clens,
                                             vector<char> &ctypes, int &alignment_length, int &motif_length) {
    /*
     *  calculates trim edges for a short motif repeat based on average motif purity
     *  @param motifpurity_threshold threshold value of the motifwise purity
     *  @param avg_motifpurity average motif wise purity of the repeat; passed as reference; updated
     *  @param clens the lengths of the consecutive CIGAR operations
     *  @param ctypes the types of the consecutive CIGAR operations
     *  @param alignment_length length of the total alignment; passed as reference and updated
     *  @param motif_length the length of the repeating motif
     *  return pair<int, int> left and right trim lengths of the CIGAR
    */
    int trim_length = 0;        // length of the trim
    int cidx_start, cidx_end;
    pair<int, int> trim_edges;  // the final pair of trim lengths

    // parameters for a certain pair of trim lengths
    int pair_match, pair_alignment; double pair_motifpurity;
    // keeping track of the maximum purity and alignment length for all the trim combinations
    double max_purity = 0; int max_alength = 0;

    // iteratively trim till the purity threshold is reached
    while (avg_motifpurity < motifpurity_threshold) {
        trim_length += 1;

        // clear everything before moving further to next trim lengths
        max_purity = 0; max_alength = 0;

        // building all trim combinations based on the trim length
        for (int i=0; i<=trim_length; i++) {

            cidx_start = -1;
            int _ = 0;
            while(_ <= i) {
                cidx_start += 1;
                if (ctypes[cidx_start] == 'M' || ctypes[cidx_start] == '=') { _ += 1; }
            }
            cidx_end = ctypes.size(); _ = 0;
            while (_ <= trim_length-i) {
                cidx_end -= 1;
                if (ctypes[cidx_end] == 'M' || ctypes[cidx_end] == '=') { _ += 1; }
            }

            motifwiseParametersInTrim(clens, ctypes, cidx_start, cidx_end, motif_length, pair_motifpurity, pair_alignment);

            // among all the trim combinations that pass the purity threshold
            // we take the one with the highest alignment length
            if (pair_motifpurity >= motifpurity_threshold) {
                if (max_alength < pair_alignment) {
                    max_purity = pair_motifpurity; max_alength = pair_alignment;
                    trim_edges = {cidx_start, cidx_end};
                }
            }
        }

        // consider the trim only if the purity has increased from previous
        // otherwise move to the next trim with increased length
        if (max_purity > avg_motifpurity) {
            avg_motifpurity = max_purity; alignment_length = max_alength;
        }

        // break more iterations if the alignment length is less than minimum length
        if (alignment_length < MINIMUM_LENGTH[motif_length]) break;
    }

    return trim_edges;
}


pair<int, int> calculateTrimEdges(int &mismatches_threshold, double &purity, int &mismatches,
                                  vector<int> &ccigar_lengths, int &alignment_length, int &motif_length) {
    /*
     *  calculates the lengths of trims on either sides based on threshold number of mismatches and indels
     *  @param mismatches_threshold allowed number of mismatches across the repeat
     *  @param purity the repeat purity value; passed as reference
     *  @param mismatches the number of mismatches in the repeat; passed as reference
     *  @param ccigar_lengths the lengths of cigar operations; has negative lengths for mismatches
     *  @param alignment_length the total length of alignment with repeat and the perfect repeat; passed as reference
     *  @param motif_length the length of the repeating motif
     *  return the pair of trim lengths from left and right
    */
    int trim_length = 0;        // length of the trim
    pair<int, int> trim_edges;  // the final pair of trim lengths

    // parameters for a certain pair of trim lengths
    int pair_matches, pair_mismatches, pair_alignment;
    // keeping track of the minimum mismatches and alignment length for all the trim combinations
    int min_mismatches = alignment_length; int max_alength = 0;
    int rtrim = 0, ltrim = 0;
    int max_stretch = (ccigar_lengths.size() >= (2*mismatches_threshold)+1) ? (2*mismatches_threshold)+1 : ccigar_lengths.size();

    for (int s = 1; s <= max_stretch; s = s+2) {
        for (ltrim=0; ltrim<=ccigar_lengths.size()-s; ltrim=ltrim+2) {
            rtrim =  ccigar_lengths.size() - ltrim - s;
            pair_matches = 0; pair_alignment = 0;

            for (int j=ltrim; j < ccigar_lengths.size()-rtrim; j++) {
                if (j%2 == 0) pair_matches += ccigar_lengths[j];    // for even indices increase the repeat length
                pair_alignment += ccigar_lengths[j];
            }
            pair_mismatches = pair_alignment - pair_matches;

            if (pair_mismatches <= mismatches_threshold) {
                if (max_alength < pair_alignment) {
                    min_mismatches = pair_mismatches; max_alength = pair_alignment;
                    trim_edges = {ltrim/2, rtrim/2};
                }
            }
        }
    }

    return trim_edges;
}


void motifwiseParameters(string &cigar, int motif_length, double &avg_motifpurity, int &avg_motifindels) {
    /*
     *  calculates the motif wise parameters including motifwise_purity and motifwise_indels
     *  @param cigar the cigar of alignment with perfect repeat
     *  @param avg_motifpurity the average motif wise purity across the repeat; passed as reference
     *  @param avg_motifindels the average indels seen motif wise; passed as reference
     *  return void
    */
    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int>  clens  = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);

    char ctype; int clength, cidx = 0;
    int qpos = 0;

    bool mismatch_continue = false;         // to check if there are contiguous mismatches

    int motif_covered = 0, motif_matches = 0, motif_mismatches = 0, motif_indels = 0;
    int excess = 0;
    vector<double> motifwise_matchpercent;
    vector<int> motifwise_indels;

    for ( ; cidx < clens.size(); cidx++) {
        clength = clens[cidx]; ctype = ctypes[cidx];

        switch (ctype) {
            case 'S':
                // soft clip: edit the repeat start and end
                qpos += clength;
                break;

            case 'X':
                qpos += clength;
                mismatch_continue = true;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_mismatches += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = 0, motif_mismatches = excess, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_mismatches += clength;
                }
                break;
            case 'I':
                qpos += clength;
                mismatch_continue = true;

                motif_indels += clength;
                break;
            case 'D':
                mismatch_continue = true;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_indels += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = 0, motif_mismatches = 0, motif_indels = excess;
                }
                else {
                    motif_covered += clength;
                    motif_indels += clength;
                }
                break;
            case '=': case 'M':
                qpos += clength;
                mismatch_continue = false;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_matches += clength - excess;
                    if (clength > motif_length) {
                        for(int _=0; _< clength/motif_length; _++) motifwise_matchpercent.push_back(1.0);
                    }
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = excess, motif_mismatches = 0, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_matches += clength;
                }
                break;
            default: break;
        }
    }

    if (motif_covered > 0) {
        motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
        motifwise_indels.push_back(motif_indels);
    }

    avg_motifpurity = 0;
    avg_motifindels = 0;

    if (motifwise_matchpercent.size() > 0) {
        for (int _=0; _<motifwise_matchpercent.size(); _++) { avg_motifpurity += motifwise_matchpercent[_]; }
        avg_motifpurity = avg_motifpurity / ((double) (motifwise_matchpercent.size()));
        for (int _=0; _<motifwise_indels.size(); _++) { avg_motifindels += motifwise_indels[_]; }
        avg_motifindels = avg_motifindels / motifwise_indels.size();
    }
}


void processCIGARWithPruning(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length,
                             int &repeat_start, int &repeat_end, int &alignment_length, int &match_units, string &new_cigar,
                             double &purity, double &avg_motifpurity, int &avg_motifindels) {
    /*
     * processes the CIGAR string and returns the repeat based on the purity threshold
     * @param seed_start position of the start of the seed sequence
     * @param seed_sequence_length length of the seed sequence
     * @param cigar CIGAR string of the alignment of the repeat with a perfect repeat
     * @param seed_sequence sequence of the seed
     * @param motif_length length of the repeating unit motif
     * @param repeat_start start of the repeat sequence; passed as reference; updated
     * @param repeat_end end of the repeat sequence; passed as reference; updated
     * @param alignment_length total alignment length between repeat sequence and a perfect repeat; passed as reference updated
     * @param match_units number of complete units of motif found in the repeat; passed as reference; updated
     * @param new_cigar new CIGAR string after trimming
     * @param purity purity of the complete repeat stretch; passed as reference; updated
     * @param avg_motifpurity average motif wise purity with a complete motif
     * @param avg_motifindels average number of indels identified per motif
     * @returns vector having corrected attributes of the repeat sequence
    */
    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int>  clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);

    // initialise the repeat coordinates to the seed coordinates
    repeat_start = seed_start; repeat_end = seed_start + seed_sequence_length;
    alignment_length = 0;

    char ctype; int clength, cidx = 0;
    int qpos = 0;
    int matches = 0; match_units = 0;

    // ccigar attributes denote the vectors regarding the compressed cigar notation
    vector<int> ccigar_indices;         // index mapping from actual cigar to compressed cigar
    vector<int> ccigar_lengths;         // compressed cigar lengths

    bool mismatch_continue = false;         // to check if there are contiguous mismatches
    int  start_soft_clip = 0;
    pair<int, int> trim_edges = {0, 0};
    new_cigar = "";

    for ( ; cidx < clens.size(); cidx++) {
        clength = clens[cidx]; ctype = ctypes[cidx];

        switch (ctype) {
            case 'S':
                // soft clip: edit the repeat start and end
                if (cidx == 0) { repeat_start += clength; start_soft_clip = clength;}
                else { repeat_end -= clength; }
                qpos += clength;
                break;

            case 'X':
                qpos += clength; alignment_length += clength;

                if (mismatch_continue) ccigar_lengths[ccigar_lengths.size() - 1] += clength;
                else { ccigar_lengths.push_back(clength); }
                ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = true; new_cigar += to_string(clength) + ctype;
                break;
            case 'I':
                qpos += clength; alignment_length += clength;

                if (mismatch_continue) ccigar_lengths[ccigar_lengths.size() - 1] += clength;
                else { ccigar_lengths.push_back(clength); }
                ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = true; new_cigar += to_string(clength) + ctype;
                break;
            case 'D':
                alignment_length += clength;

                if (mismatch_continue) ccigar_lengths[ccigar_lengths.size() - 1] += clength;
                else { ccigar_lengths.push_back(clength); }
                ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = true; new_cigar += to_string(clength) + ctype;
                break;
            case '=': case 'M':
                qpos += clength; alignment_length += clength;
                matches += clength; match_units += clength/motif_length;

                ccigar_lengths.push_back(clength); ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = false; new_cigar += to_string(clength) + ctype;
                break;
            default: break;
        }
    }

    purity = (double)matches / (double)alignment_length;
    int mismatches = alignment_length - matches;

    bool trim = false;
    if (purity < PURITY_THRESHOLD) trim = true;

    if (trim) {

        trim_edges = calculateTrimEdges(PURITY_THRESHOLD, purity, ccigar_lengths, alignment_length, motif_length);

        // based on the trim edges we adjust all the repeat parameters
        new_cigar = ""; matches = 0; match_units = 0; qpos = start_soft_clip;

        for (int i=0; i<ccigar_indices.size(); i++) {

            int ccidx = ccigar_indices[i];
            if (start_soft_clip) { clength = clens[i+1]; ctype = ctypes[i+1]; }
            else { clength = clens[i]; ctype = ctypes[i]; }

            if (ccidx < 2*trim_edges.first) {
                if (ctype != 'D') {
                    repeat_start += clength;
                    qpos += clength;
                }
            }

            else if (ccidx >= 2*trim_edges.first && ccidx <= ccigar_lengths.size()-1-(2*trim_edges.second)) {
                new_cigar += to_string(clength) + ctype;
                switch(ctype) {
                    case 'M': case '=':
                        matches += clength; qpos += clength;
                        match_units += clength/motif_length;
                        break;
                    case 'X': case 'I':
                        qpos += clength;
                        break;
                    default: break;
                }
            }

            else {
                if (ctype != 'D') repeat_end -= clength;
            }
        }
    }

    motifwiseParameters(new_cigar, motif_length, avg_motifpurity, avg_motifindels);
}


void processCIGARMotifWise(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length,
                           int &repeat_start, int &repeat_end, int &alignment_length, string &new_cigar, double &purity,
                           int &substitutions, int &indels, double &avg_motifpurity, int &avg_motifindels, int &avg_matchlen) {
    /*
     * processes the CIGAR string and returns the repeat based on the purity threshold
     * @param seed_start position of the start of the seed sequence
     * @param seed_sequence_length length of the seed sequence
     * @param cigar CIGAR string of the alignment
     * @param seed_sequence nucleotide sequence of the seed
     * @param motif_length periodicity of the repeat seed
     * @param repeat_start start of the repeat passed as reference
     * @param repeat_end end of the repeat passed as reference
     * @param alignment_length alignment length of the seed sequence and perfect repeat passed as reference
     * @param new_cigar cigar after processing the locus passed as reference
     * @param purity purity of the repeat passed as reference
     * @param avg_motifpurity motif wise purity of the repeat passed as reference
     * @param avg_motifindels average number of indels observed for repeat passed as reference
     * @returns vector having corrected attributes of the repeat sequence
    */

    // should add condition to trim based on motif purity

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int>  clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);

    if (ctypes.size()==1 && ctypes[0]=='S') { return; }
    alignment_length = 0;
    substitutions = 0, indels = 0;

    // initialise the repeat coordinates to the seed coordinates
    repeat_start = seed_start, repeat_end = seed_start + seed_sequence_length;

    char ctype; int clength, cidx = 0;
    int qpos = 0;
    int matches = 0, match_units = 0;

    int  start_soft_clip = 0;
    pair<int, int> trim_edges = {0, 0};
    new_cigar = "";

    int motif_covered = 0, motif_matches = 0, motif_mismatches = 0, motif_indels = 0;
    int excess = 0;
    vector<double> motifwise_matchpercent;
    vector<int> motifwise_indels;
    vector<int>  match_lens; int match_length = 0;

    for ( ; cidx < clens.size(); cidx++) {
        clength = clens[cidx]; ctype = ctypes[cidx];

        switch (ctype) {
            case 'S':
                // soft clip: edit the repeat start and end
                if (cidx == 0) { repeat_start += clength; start_soft_clip = clength;}
                else { repeat_end -= clength; }
                qpos += clength;
                break;

            case 'X':
                qpos += clength; alignment_length += clength;
                match_length += clength; substitutions += 1;

                new_cigar += to_string(clength) + ctype;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_mismatches += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = 0, motif_mismatches = excess, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_mismatches += clength;
                }
                break;
            case 'I':
                qpos += clength; alignment_length += clength;
                match_lens.push_back(match_length); match_length = 0; indels += 1;

                new_cigar += to_string(clength) + ctype;

                motif_indels += clength;
                break;
            case 'D':
                alignment_length += clength;
                match_lens.push_back(match_length); match_length = 0; indels += 1;

                new_cigar += to_string(clength) + ctype;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_indels += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = 0, motif_mismatches = 0, motif_indels = excess;
                }
                else {
                    motif_covered += clength;
                    motif_indels += clength;
                }
                break;
            case '=': case 'M':
                qpos += clength; alignment_length += clength;
                matches += clength; match_units += clength/motif_length;
                match_length += clength;

                new_cigar += to_string(clength) + ctype;

                if ((motif_covered+clength) >= motif_length) {
                    excess = (motif_covered+clength) % motif_length;
                    motif_covered += clength - excess;
                    motif_matches += clength - excess;
                    motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
                    motifwise_indels.push_back(motif_indels);
                    motif_covered = excess, motif_matches = excess, motif_mismatches = 0, motif_indels = 0;
                }
                else {
                    motif_covered += clength;
                    motif_matches += clength;
                }
                break;
            default: break;
        }
    }

    if (match_length > 0) match_lens.push_back(match_length);

    if (motif_covered > 0) {
        motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
        motifwise_indels.push_back(motif_indels);
    }

    purity = double(matches)/double(alignment_length);

    avg_motifpurity = 0;
    avg_motifindels = 0;
    if (motifwise_matchpercent.size() > 0) {
        for (int _=0; _<motifwise_matchpercent.size(); _++) { avg_motifpurity += motifwise_matchpercent[_]; }
        avg_motifpurity = avg_motifpurity / ((double) (motifwise_matchpercent.size()));
        for (int _=0; _<motifwise_indels.size(); _++) { avg_motifindels += motifwise_indels[_]; }
        avg_motifindels = avg_motifindels / motifwise_indels.size();
    }

    avg_matchlen = 0;
    if (match_lens.size() > 0) {
        for (int _=0; _<match_lens.size(); _++) { avg_matchlen += match_lens[_]; }
        avg_matchlen = (int) (((double) avg_matchlen) / ((double) (match_lens.size())));
    }

    bool trim = false;
    // do not trim if the average perfect length stretch is atleast twice as motif length
    if (avg_matchlen > 2*motif_length) { trim = false; }

    else if (avg_motifpurity < MOTIFPURITY_THRESHOLD
             && match_units > PERFECT_UNITS[motif_length]
             && alignment_length > MINIMUM_LENGTH[motif_length]) {
        // trims only when the motif wise purity is lower than the threshold
        // and match units and alignment length are more than threshold
        trim = true;
    }

    // if (repeat_start == 3866 && repeat_end == 4134) {
    //     cout << repeat_start << "\t" << repeat_end << "\t" << cigar << "\n";
    //     cout << avg_motifpurity << "\t" << trim << "\n";
    // }


    if (trim) {
        trim_edges = calculateTrimEdgesMotifPurity(MOTIFPURITY_THRESHOLD, avg_motifpurity, clens, ctypes, alignment_length, motif_length);
        // based on the trim edges we adjust all the repeat parameters
        for (int i=0; i<get<0>(trim_edges); i++) {
            if (ctypes[i] == 'M' || ctypes[i] == '=' || ctypes[i] == 'X' || ctypes[i] == 'I') {
                repeat_start += clens[i];
            }
        }
        for (int i=clens.size()-1; i>get<1>(trim_edges); i--) {
            if (ctypes[i] == 'M' || ctypes[i] == '=' || ctypes[i] == 'X' || ctypes[i] == 'I') {
                repeat_end -= clens[i];
            }
        }

        cigar = "";
        for (int i=get<0> (trim_edges); i<=get<1> (trim_edges); i++) {
            clength = clens[i]; ctype = ctypes[i];
            cigar += to_string(clength) + ctype;
        }
        processCIGARMotifWise(repeat_start, repeat_end-repeat_start, cigar, seed_sequence, motif_length,
                              repeat_start, repeat_end, alignment_length, new_cigar, purity,
                              substitutions, indels, avg_motifpurity, avg_motifindels, avg_matchlen);
        return;
    }
}
