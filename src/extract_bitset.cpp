#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>

using namespace std;

int MINIMUM_MLEN = 2; int MAXIMUM_MLEN = 100;
int MINIMUM_SHIFT = 1; int MAXIMUM_SHIFT = 105;
int NSHIFTS = MAXIMUM_SHIFT - MINIMUM_SHIFT + 1;

void generateAnchoredShiftXORs(vector<boost::dynamic_bitset<>> &lshift_xor_bsets, boost::dynamic_bitset<> &N_bset,
                               vector<boost::dynamic_bitset<>> &lsxor_anchor_bsets, int anchor_size) {

    /*
     *  generates shift XOR bitsets only retaining the anchors
     *  @param lshift_xor_bsets vector of left shift XOR bitsets of all shifts
     *  @param N_bset bitset with information of N positions
     *  @param lsxor_anchor_bsets vector of the left shift anchor bitsets
     *  @param anchor_size the length of the anchor size
     *  @return void
    */

    int bset_size = N_bset.size();
    int anchor_start = -1;
    int motif_length;
    for (int lsxor_idx=0; lsxor_idx < NSHIFTS; lsxor_idx++) {
        motif_length = MINIMUM_SHIFT + lsxor_idx;
        boost::dynamic_bitset<> anchor_bset(bset_size, 0ull);
        for (int xor_idx = bset_size-1; xor_idx >= lsxor_idx + MINIMUM_SHIFT; xor_idx--) {

            if (lshift_xor_bsets[lsxor_idx][xor_idx] == 1) {
                if (anchor_start == -1) anchor_start = xor_idx;
            }

            else {
                if (anchor_start - xor_idx >= anchor_size && anchor_start - xor_idx <= motif_length) {
                    // the anchors to be retained should be at least of the minimum anchor size mentioned
                    // and not more than twice of the motif length it being tagged in because this will retain all the
                    // perfect repeats of that motif length
                    anchor_bset.set(xor_idx+1, anchor_start - xor_idx, 1);
                } anchor_start = -1;
            }
        }

        lsxor_anchor_bsets.push_back(anchor_bset);
        anchor_start = -1;
    }
}

void buildBitDatastructures(string &sequence, int sequence_length, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                            boost::dynamic_bitset<> &N_bset, boost::dynamic_bitset<> &A, boost::dynamic_bitset<> &T, boost::dynamic_bitset<> &G,
                            boost::dynamic_bitset<> &C, vector<boost::dynamic_bitset<>*> &MATRIX) {
    /*
     * builds the bit datastructures for the sequence
     * @param left_bset dynamic bitset for left side of the sequence
     * @param right_bset dynamic bitset for right side of the sequence
     * @param N_bset dynamic bitset for N nucleotides
     * @param A, T, G, C dynamic bitsets for A, T, G, C nucleotides respectively
     * @param MATRIX vector of pointers to the nucleotide bitsets
    */

    int seq_idx = 0, bidx; char nuc;
    for (; seq_idx < sequence_length; seq_idx++) {
        nuc = sequence[seq_idx];                // nucleotide
        bidx = (sequence_length-1) - seq_idx;   // bit indexing starts from right
        switch (nuc) {
            case 'A': case 'a': // 00
                left_bset[bidx] = 0; right_bset[bidx] = 0;
                MATRIX.push_back(&A);
                A[bidx] = 1; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'C': case 'c': // 01
                left_bset[bidx] = 0; right_bset[bidx] = 1;
                MATRIX.push_back(&C);
                C[bidx] = 1; A[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'G': case 'g': // 10
                left_bset[bidx] = 1; right_bset[bidx] = 0;
                MATRIX.push_back(&G);
                G[bidx] = 1; A[bidx] = 0; C[bidx] = 0; T[bidx] = 0; break;
            case 'T': case 't': // 11
                left_bset[bidx] = 1; right_bset[bidx] = 1;
                MATRIX.push_back(&T);
                T[bidx] = 1; A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; break;
            default: // no match probably N or any other nuc
                N_bset[bidx] = 1;
                MATRIX.push_back(NULL);
                A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
        }
    }
}


void processSequence(string sequence_id, string &sequence, vector<boost::dynamic_bitset<>> &lshift_xor_bsets) {
    /*
     *  processes each sequence from 2-bit conversion to identifying repeats
     *  @param sequence_id name of the sequence from fasta
     *  @param sequence string of the fasta sequence
     *  @return void generates the dynamic bitsets of shift XOR matches and
     *               proceeds to identifying repeats
    */

    double seconds_since_start;

    // converting the sequencing to bitsets
    int sequence_length = sequence.length();
    boost::dynamic_bitset<> left_bset(sequence_length, 0ull);
    boost::dynamic_bitset<> right_bset(sequence_length, 0ull);
    boost::dynamic_bitset<> N_bset(sequence_length, 0ull);

    // bitsets of A,C,G,T to build the matrix
    boost::dynamic_bitset<> A(sequence_length, 0ull);
    boost::dynamic_bitset<> T(sequence_length, 0ull);
    boost::dynamic_bitset<> G(sequence_length, 0ull);
    boost::dynamic_bitset<> C(sequence_length, 0ull);
    vector<boost::dynamic_bitset<>*> MATRIX;

    // builds the bit datastructures that are necessary
    buildBitDatastructures(sequence, sequence_length, left_bset, right_bset, N_bset, A, T, G, C, MATRIX);

    lshift_xor_bsets.clear();
    // vector<boost::dynamic_bitset<>> lshift_xor_bsets;       // vector of dynamic bitsets for each shift XOR
    // generating the shift XORs from minimum shift size to maximum shift size
    for (int i = MINIMUM_SHIFT; i <= MAXIMUM_SHIFT; i++) {
        lshift_xor_bsets.push_back( ~(left_bset ^ (left_bset<<(i))) & ~(right_bset ^ (right_bset<<(i))) );
    }

     vector<boost::dynamic_bitset<>> lsxor_anchor_bsets; 
    int anchor_jump   = 1;
    int anchor_length = 5;
    generateAnchoredShiftXORs(lshift_xor_bsets, N_bset, lsxor_anchor_bsets, anchor_length);
    boost::dynamic_bitset<> anchor_bset(sequence_length, 0ull);
    int motif_length = MINIMUM_MLEN;
    for (; motif_length <= MAXIMUM_MLEN; motif_length++) {
        anchor_bset.reset();

        if (motif_length <= 6) { anchor_jump = 1; }
        else if (motif_length <= 20) { anchor_jump = (2*MAXIMUM_MLEN)/10; }
        else { anchor_jump = 4; }

        int i = motif_length - anchor_jump;
        if (i < MINIMUM_SHIFT) { i = MINIMUM_SHIFT; }
        for (; i <= motif_length + anchor_jump; i++) {
            int shift_idx = i - MINIMUM_SHIFT;
            // OR with actual shift XOR for same motif size
            if (i == motif_length) { anchor_bset |= lshift_xor_bsets[shift_idx]; }
            // OR with anchor bitset for neigboring shifts
            else { anchor_bset |= lsxor_anchor_bsets[shift_idx]; }
        }

        lshift_xor_bsets[motif_length-MINIMUM_SHIFT] = anchor_bset;
    }

    lsxor_anchor_bsets.clear();
    // for (int i=MINIMUM_MLEN; i <= MAXIMUM_MLEN; i++) {
    //     lshift_xor_bsets[i-MINIMUM_SHIFT] |= (lshift_xor_bsets[i-MINIMUM_SHIFT] >> i);
    // }
}


int main(int argc, char* argv[]) {
    /*
     * main function to process the fasta file and identify repeats
    */

    string infile = argv[1];
    string bed_file = argv[2];

    cout << "Input fasta file: " << infile << "\n";
    cout << "Input bed file: " << bed_file << "\n";
    ifstream fasta_file(infile);
    vector<boost::dynamic_bitset<>> lshift_xor_bsets;     // vector of dynamic bitsets for anchor bitsets

    string line, sequence_id = "", sequence = "";
    while (getline(fasta_file, line)) {
        if (line.empty()) { continue; }
        if (line[0] == '>') {
            // process the previous sequence
            if (!sequence.empty()) {
                processSequence(sequence_id, sequence, lshift_xor_bsets);
                sequence = "";
            }
            // new sequence
            sequence_id = line.substr(1); // remove '>'
        } else {
            // append to current sequence
            sequence += line;
        }
    }
    // process the last sequence
    if (!sequence.empty()) {
        processSequence(sequence_id, sequence, lshift_xor_bsets);
    }
    fasta_file.close();

    int sequence_length = sequence.length();
    string chrom, start_str, end_str, motif_size_str;
    int start = 0, end = 0, motif_size = 0;
    int length = 0;
    
    ifstream bed(bed_file);
    string bed_line;
    vector<tuple<string, int, int, int>> bed_entries;
    while (getline(bed, bed_line)) {
        if (bed_line.empty()) continue;
        istringstream iss(bed_line);

        // Use tab delimiter
        getline(iss, chrom, '\t');
        getline(iss, start_str, '\t');
        getline(iss, end_str, '\t');
        getline(iss, motif_size_str, '\t');

        start = stoi(start_str);
        end = stoi(end_str);
        motif_size = stoi(motif_size_str);
        
        cout << chrom << "\t" << start << "\t" << end << "\t" << motif_size << "\n";
        length = end - start;
        boost::dynamic_bitset<> region_bset(length, 0ull);

        for (int i = start; i < end; i++) {
            if (lshift_xor_bsets[motif_size - MINIMUM_SHIFT][sequence_length - 1 - i] == 1) {
                region_bset[length - 1 - (i-start)] = 1;
            }
        }
        cout << "Region bitset: " << region_bset << "\n\n";
    }
    bed.close();

    return 0;
}