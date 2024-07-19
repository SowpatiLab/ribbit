#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>

#include "global_variables.h"
#include "bitseq_utils.h"

using namespace std;
using namespace boost;


unordered_map<int, unordered_map<uint256_t, string>> motifs_map; // Store motif to int32
string calculateMotif(uint256_t motif_unit, int motif_length) {
    /*
     * converts a integer representation of motif to string
     * @param motif the 32-bit integer representation of motif
     * @param motif_length length of the motif
     * @returns string string representation of the motif
    */

    if (motifs_map[motif_length].find(motif_unit) != motifs_map[motif_length].end()) {
        return motifs_map[motif_length][motif_unit];
    }

    uint256_t mask = 3;
    string motif = "";
    for (int i=0; i<motif_length; i++) {
        uint256_t val = (motif_unit >> (2*(motif_length-1-i))) & mask;
        if (val == 0) motif += 'A';
        else if (val == 1) motif += 'C';
        else if (val == 2) motif += 'G';
        else if (val == 3) motif += 'T';
    }

    motifs_map[motif_length][motif_unit] = motif;
    return motif;
}


unordered_map<int, unordered_map<uint256_t, int>> atomicitybits_map; // store atomicity of motifs
int calculateAtomicity(dynamic_bitset<> &window, int &motif_length) {
    /*
     * finding the atomicity of a motif represented in dynamic bits
     * @param window the dynamic bitset of the motif sequence; length = 2*motif_length
     * @param motif_length length of the motif
     * @returns int atomicity of the motif
    */
    uint256_t motif = window.to_ulong();
    if (atomicitybits_map[motif_length].find(motif) != atomicitybits_map[motif_length].end()) {
        // if motif already existing in the atomicity map
        return atomicitybits_map[motif_length][motif];
    }

    vector<int> motif_factors;
    for (int f = 1; f <= motif_length/2; f++) {
        if (motif_length % f == 0) { motif_factors.push_back(f); }
    }

    uint256_t shift, original;
    uint256_t mask;
    uint256_t cycle;
    for (int f: motif_factors) {
        mask = 0;
        for (int i=0; i<2*(motif_length-f); i++) { mask <<= 1; mask |= 1; }
        shift = motif >> 2*f;
        original = mask & motif;
        if (shift == original) {
            cycle = motif;
            atomicitybits_map[motif_length][cycle] = f;
            for (int i = 0; i < motif_length-1; i++) {
                cycle = ((window >> (2*(motif_length-(i+1)))) | (window << (2*(i+1)))).to_ulong();
                atomicitybits_map[motif_length][cycle] = f;
            }
            return f;
        }
    }

    cycle = motif;
    atomicitybits_map[motif_length][cycle] = motif_length;
    for (int i = 0; i < motif_length-1; i++) {
        cycle = ((window >> (2*(motif_length-(i+1)))) | (window << (2*(i+1)))).to_ulong();
        atomicitybits_map[motif_length][cycle] = motif_length;
    }
    return motif_length;
}

int calculateAtomicity(uint256_t &motif, int &motif_length) {
    /*
     * finding the atomicity of a motif represented in dynamic bits
     * @param window the dynamic bitset of the motif sequence; length = 2*motif_length
     * @param motif_length length of the motif
     * @returns int atomicity of the motif
    */

    vector<int> motif_factors;
    for (int f = 1; f <= motif_length/2; f++) {
        if (motif_length % f == 0) { motif_factors.push_back(f); }
    }

    uint256_t shift, original;
    uint256_t mask;
    for (int f: motif_factors) {
        mask = 0;
        for (int i=0; i<2*(motif_length-f); i++) { mask <<= 1; mask |= 1; }
        shift = motif >> 2*f;
        original = mask & motif;
        if (shift == original) {
            return f;
        }
    }

    return motif_length;
}

int calculateAtomicityLongMotif(uint256_t &motif, int &motif_length) {
    /*
     * finding the atomicity of a motif represented in dynamic bits
     * @param window the dynamic bitset of the motif sequence; length = 2*motif_length
     * @param motif_length length of the motif
     * @returns int atomicity of the motif
    */

    uint256_t shift, original;
    uint256_t mask;
    for (int f=1; f<motif_length - motif_length/3; f++) {
        mask = 0;
        for (int i=0; i<2*(motif_length-f); i++) { mask <<= 1; mask |= 1; }
        shift = motif >> 2*f;
        original = mask & motif;
        if (shift == original) {
            return f;
        }
    }

    return motif_length;
}

int calculateAtomicity(uint32_t &motif, int &motif_length) {
    /*
     * finding the atomicity of a motif represented in dynamic bits
     * @param window the dynamic bitset of the motif sequence; length = 2*motif_length
     * @param motif_length length of the motif
     * @returns int atomicity of the motif
    */

    if (atomicitybits_map[motif_length].find(motif) != atomicitybits_map[motif_length].end()) {
        // if motif already existing in the atomicity map
        return atomicitybits_map[motif_length][motif];
    }

    vector<int> motif_factors;
    for (int f = 1; f <= motif_length/2; f++) {
        if (motif_length % f == 0) { motif_factors.push_back(f); }
    }

    uint32_t shift, original;
    uint32_t mask;
    uint32_t cycle;
    for (int f: motif_factors) {
        mask = 0;
        for (int i=0; i<2*(motif_length-f); i++) { mask <<= 1; mask |= 1; }
        shift = motif >> 2*f;
        original = mask & motif;
        if (shift == original) {
            cycle = motif;
            atomicitybits_map[motif_length][cycle] = f;
            for (int i = 0; i < motif_length-1; i++) {
                cycle = (motif >> (2*(motif_length-(i+1)))) | (motif << (2*(i+1)));
                atomicitybits_map[motif_length][cycle] = f;
            }
            return f;
        }
    }

    cycle = motif;
    atomicitybits_map[motif_length][cycle] = motif_length;
    for (int i = 0; i < motif_length-1; i++) {
        cycle = (motif >> (2*(motif_length-(i+1)))) | (motif << (2*(i+1)));
        atomicitybits_map[motif_length][cycle] = motif_length;
    }
    return motif_length;
}

uint32_t calculateRepeatClass(boost::dynamic_bitset<> &window, int &motif_length) {
    /*
     * calculates repeat class of the motif within the given window
     * @param window the dynamic bitset which has the motif in bits
     *               length of the window bset is same as the motif size
     * @param motif_length length of the motif
     * @returns uint32_t the repeat class represented as bits
    */

    uint32_t motif = window.to_ulong();
    uint32_t repeat_class = REPEAT_CLASSES[motif_length-1][motif];
    if (repeat_class != NUM_MOTIFS) {
        // if motif already exists in the repeat class map
        return repeat_class;
    }

    repeat_class = motif; // initialize repeat_class as the motif
    uint32_t cycle = motif;
    uint32_t cycles[motif_length];
    cycles[0] = cycle;
    for (int i = 0; i < motif_length-1; i++) {
        // shifting the window to the right by motif_length-(i+1) | shifting the window to the left by i+1
        // Example GATA
        // 1. 000G | ATA0   2. 00GA | TA00  3. 0GAT | A000
        cycle = ((window >> (2*(motif_length-(i+1)))) | (window << (2*(i+1)))).to_ulong();
        // if cycle is smaller than repeat class replace repeat class
        if (cycle < repeat_class)  repeat_class = cycle;
        cycles[i+1] = cycle;
    }

    // Store the repeat class for all cyclic permutations of the motif
    for (int i=0; i < motif_length; i++) {
        REPEAT_CLASSES[motif_length-1][cycles[i]] = repeat_class;
    }

    return repeat_class;
}
