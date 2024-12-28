#include <iostream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <numeric>
#include <utility>

#include "global_variables.h"
#include "ssw_cpp.h"

using namespace std;
using namespace boost;


tuple<vector<int>, vector<char>> cigarSplit(const char* cigar) {
    const char * t; // first copy the pointer to not change the original
    string length = "";
    vector<int> clens; vector<char> ctypes;

    for (t = cigar; *t != '\0'; t++) {
        if (isdigit(*t)) {
            length += *t;
        }
        else {
            clens.push_back(stoi(length));
            ctypes.push_back(*t);
            length = "";
        }
    }

    return tuple<vector<int>, vector<char>> {clens, ctypes};
}

double calculateMedian(std::vector<double>& numbers) {
    int size = numbers.size();
    std::sort(numbers.begin(), numbers.end());
    
    if (size % 2 == 0) {
        return (numbers[size / 2 - 1] + numbers[size / 2]) / 2.0;
    } else {
        return numbers[size / 2];
    }
}


pair<int, int> calculateTrimEdges(double &purity_threshold, double &purity, vector<int> &ccigar_lengths,
                                  int &alignment_length, int &motif_length) {

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


pair<int, int> calculateTrimEdges(int &mismatches_threshold, double &purity, int &mismatches,
                                  vector<int> &ccigar_lengths, int &alignment_length, int &motif_length) {

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
    
    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar.c_str());
    vector<int>  clens = get<0> (csplit);
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
        avg_motifpurity = (double)(std::reduce(motifwise_matchpercent.begin(), motifwise_matchpercent.end(), 0)) / (double)motifwise_matchpercent.size();
        // avg_motifpurity = calculateMedian(motifwise_matchpercent);
        avg_motifindels = std::reduce(motifwise_indels.begin(), motifwise_indels.end()) / motifwise_indels.size();
    }
}


void processCIGARWithPruning(int seed_start, int seed_sequence_length, string &cigar, string &seed_sequence, int motif_length, 
                             int &repeat_start, int &repeat_end, int &alignment_length, int &match_units, string &new_cigar,
                             double &purity, double &avg_motifpurity, int &avg_motifindels) {
    /*
     * processes the CIGAR string and returns the repeat based on the purity threshold
     * @param seed_start position of the start of the seed sequence
     * @param seed_sequence_length length of the seed sequence
     * @param cigar      CIGAR string of the alignment
     * @param seed_sequence sequence of the seed
     * @returns vector having corrected attributes of the repeat sequence
    */
    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar.c_str());
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
                           int &repeat_start, int&repeat_end, int &alignment_length, string &new_cigar, double &purity,
                           double &avg_motifpurity, int &avg_motifindels) {
    /*
     * processes the CIGAR string and returns the repeat based on the purity threshold
     * @param seed_start position of the start of the seed sequence
     * @param seed_sequence_length length of the seed sequence
     * @param cigar      CIGAR string of the alignment
     * @param seed_sequence sequence of the seed
     * @returns vector having corrected attributes of the repeat sequence
    */
    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar.c_str());
    vector<int>  clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);    

    // initialise the repeat coordinates to the seed coordinates
    repeat_start = seed_start, repeat_end = seed_start + seed_sequence_length;
    alignment_length = 0;

    char ctype; int clength, cidx = 0;
    int qpos = 0;
    int matches = 0, match_units = 0;    

    // ccigar attributes denote the vectors regarding the compressed cigar notation
    vector<int> ccigar_indices;         // index mapping from actual cigar to compressed cigar
    vector<int> ccigar_lengths;         // compressed cigar lengths

    bool mismatch_continue = false;         // to check if there are contiguous mismatches
    int  start_soft_clip = 0;
    pair<int, int> trim_edges = {0, 0};
    new_cigar = "";

    int motif_covered = 0, motif_matches = 0, motif_mismatches = 0, motif_indels = 0;
    int excess = 0;
    vector<double> motifwise_matchpercent;
    vector<int> motifwise_indels;

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

                if (mismatch_continue) ccigar_lengths[ccigar_lengths.size() - 1] += clength;
                else { ccigar_lengths.push_back(clength); } 
                ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = true; new_cigar += to_string(clength) + ctype;
                
                motif_indels += clength;
                break;
            case 'D':
                alignment_length += clength;

                if (mismatch_continue) ccigar_lengths[ccigar_lengths.size() - 1] += clength;
                else { ccigar_lengths.push_back(clength); } 
                ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = true; new_cigar += to_string(clength) + ctype;
                
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

                ccigar_lengths.push_back(clength); ccigar_indices.push_back(ccigar_lengths.size() - 1);
                mismatch_continue = false; new_cigar += to_string(clength) + ctype;

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

    if (motif_covered > 0) {
        motifwise_matchpercent.push_back((double) motif_matches/ (double) (motif_matches + motif_mismatches + motif_indels));
        motifwise_indels.push_back(motif_indels);
    }

    purity = double(matches)/double(alignment_length);
    
    avg_motifpurity = 0;
    avg_motifindels = 0;
    if (motifwise_matchpercent.size() > 0) {
        avg_motifpurity = (double)(std::reduce(motifwise_matchpercent.begin(), motifwise_matchpercent.end())) / (double)motifwise_matchpercent.size();
        // avg_motifpurity = calculateMedian(motifwise_matchpercent);
        avg_motifindels = std::reduce(motifwise_indels.begin(), motifwise_indels.end()) / motifwise_indels.size();
    }
}
