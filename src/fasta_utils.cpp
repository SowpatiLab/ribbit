#include "fasta_utils.h"

using namespace std;


void parseFai(string infai, int &nseqs, unordered_map<string, int> &seq_lens) {
    /*
     * parses the fasta index file
     * @param infai file identifier of the fasta index
     * @param nseqs number of sequences to be updated
     * @param seq_lens unordered map of sequence lengths
    */

    ifstream ins(infai);
    if (!ins) { return; }
    string line, chrom;
    int seq_len;

    while (getline(ins, line)) {
        int delim_pos = line.find('\t');
        chrom = line.substr(0, delim_pos);
        seq_len = stoi(line.substr(delim_pos, line.length()-delim_pos));
        seq_lens[chrom] = seq_len; nseqs += 1;
    }

    ins.close();
}


int failedSeeds(vector<tuple<int, int, int, int>> &seed_positions) {
    /*
     * counts the number of failed seeds in a seed positions vector
     * @param seed_positions vector of seed_positions
    */
    int count = 0; int seed_type = 0;
    for (int seed_idx=seed_positions.size()-1; seed_idx>=0; seed_idx--) {
        seed_type  = get<3> (seed_positions[seed_idx]);
        if (seed_type == -1) { count += 1; }
    }
    return count;
}


int calculateBasesinSeeds(vector<tuple<int, int, int, int>> &seed_positions) {
    /*
     * calculates the number of bases in the seeds
     * @param seed_positions vector of seed positions
     * @return int number of bases in the seeds
    */

    int total_bases = 0;
    if (seed_positions.empty()) return 0;

    int prev_end = -1;
    for (const auto& seed : seed_positions) {
        int start = get<0>(seed);
        int end = get<1>(seed);
        if (start > prev_end) {
            total_bases += end - start;
            prev_end = end;
        } else if (end > prev_end) {
            total_bases += end - prev_end;
            prev_end = end;
        }
    }
    return total_bases;
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


void processSequence(string sequence_id, string &sequence, ostream &out, int chunk_start, int chunk_end,
                     vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     *  processes each sequence from 2-bit conversion to identifying repeats
     *  @param sequence_id name of the sequence from fasta
     *  @param sequence string of the fasta sequence
     *  @param out the output file to which the output has to be printed
     *  @return void generates the dynamic bitsets of shift XOR matches and
     *               proceeds to identifying repeats
    */

    START_TIME = time(0);
    double seconds_since_start;

    CHUNK_START = chunk_start;

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

    vector<boost::dynamic_bitset<>> lshift_xor_bsets;       // vector of dynamic bitsets for each shift XOR
    // generating the shift XORs from minimum shift size to maximum shift size
    for (int i = MINIMUM_SHIFT; i <= MAXIMUM_SHIFT; i++) {
        lshift_xor_bsets.push_back( ~(left_bset ^ (left_bset<<(i))) & ~(right_bset ^ (right_bset<<(i))) );
    }
    seconds_since_start = difftime( time(0), START_TIME);

    // generating seed positions; vector of tuple with start and end of the seeds
    vector<tuple<int, int, int, int>> seed_positions_perfect;
    vector<tuple<int, int, int, int>> seed_positions_substut;
    vector<tuple<int, int, int, int>> seed_positions_anchored;

    if (PURITY_THRESHOLD == 1) {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset);
    }

    else {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset);

        seed_positions_substut = processShiftXORswithSubstitutions(lshift_xor_bsets, N_bset, seed_positions_perfect);

        // filtering out the perfect seeds which are inside substituted seeds with short flanks
        filterPerfectSeeds(seed_positions_perfect, seed_positions_substut);

        // generating the anchor bitsets for all shift sizes
        vector<boost::dynamic_bitset<>> lsxor_anchor_bsets;     // vector of dynamic bitsets for anchor bitsets

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
        for (int i=MINIMUM_MLEN; i <= MAXIMUM_MLEN; i++) {
            lshift_xor_bsets[i-MINIMUM_SHIFT] |= (lshift_xor_bsets[i-MINIMUM_SHIFT] >> i);
        }
        seed_positions_anchored = processShiftXORsAnchored(lshift_xor_bsets, N_bset, seed_positions_perfect, seed_positions_substut);
        filterShortSeeds(seed_positions_perfect);
        filterShortSeeds(seed_positions_substut);
        filterShortSeeds(seed_positions_anchored);
    }

    // Objects used by complete striped smithwater algorithm
    StripedSmithWaterman::Aligner   aligner;
    StripedSmithWaterman::Filter    filter;
    StripedSmithWaterman::Alignment alignment;

    tuple<int,int,int,int> seed;
    int seed_start, seed_end, seed_mlen, seed_type, seed_bset_size;
    uint64_t smallest; int smallest_type = -1;
    int spidx_p=0, spidx_s=0, spidx_a=0;
    int processed_seeds = 0;

    while (spidx_p < seed_positions_perfect.size() || spidx_s < seed_positions_substut.size() || spidx_a < seed_positions_anchored.size()) {
        smallest = -1;
        smallest_type = RANK_N;

        while (spidx_p < seed_positions_perfect.size() && get<3> (seed_positions_perfect[spidx_p]) == RANK_N) { spidx_p += 1; }
        while (spidx_s < seed_positions_substut.size() && get<3> (seed_positions_substut[spidx_s]) == RANK_N) { spidx_s += 1; }
        while (spidx_a < seed_positions_anchored.size() && get<3> (seed_positions_anchored[spidx_a]) == RANK_N) { spidx_a += 1; }

        // Find the smallest element among the current elements of the three vectors
        if (spidx_p < seed_positions_perfect.size() && (smallest > get<0> (seed_positions_perfect[spidx_p]))) {
            smallest = get<0> (seed_positions_perfect[spidx_p]); smallest_type = RANK_P;
        }
        if (spidx_s < seed_positions_substut.size() && (smallest > get<0> (seed_positions_substut[spidx_s]))) {
            smallest = get<0> (seed_positions_substut[spidx_s]); smallest_type = RANK_S;
        }
        if (spidx_a < seed_positions_anchored.size() && (smallest > get<0> (seed_positions_anchored[spidx_a]))) {
            smallest = get<0> (seed_positions_anchored[spidx_a]); smallest_type = RANK_A;
        }
        // Print the smallest element and move the corresponding pointer
        if (smallest_type == RANK_P) {
            seed = seed_positions_perfect[spidx_p]; spidx_p += 1;
        }
        else if (smallest_type == RANK_S) {
            seed = seed_positions_substut[spidx_s]; spidx_s += 1;
        }
        else if (smallest_type == RANK_A) {
            seed = seed_positions_anchored[spidx_a]; spidx_a += 1;
        }

        if (smallest_type == RANK_N) { continue; } // skip if no seeds left

        seed_type  = get<3> (seed);
        if (seed_type == RANK_N) { continue; }
        seed_start = get<0> (seed);
        seed_end   = get<1> (seed);
        seed_mlen  = get<2> (seed);
        if (seed_end + seed_mlen > sequence_length) { seed_end = sequence_length - seed_mlen; }

        seed_bset_size = seed_end - seed_start;
        if (seed_bset_size < SEEDLEN_CUTOFF[seed_mlen - MINIMUM_MLEN]) { continue; }

        boost::dynamic_bitset<> seed_bset(seed_bset_size, 0ull);
        for (int j = seed_start; j < seed_end; j++) {
            seed_bset[seed_end - 1 - j] = lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT][sequence_length - 1 - j];
        }

        cout << sequence_id << "\t" << seed_start + CHUNK_START << "\t" << seed_end + CHUNK_START << "\t" << seed_mlen << "\t" << seed_end - seed_start
             << "\t" << seed_type << "\n";
        continue;

        // process seed if it is alteast the size of the motif length
        processed_seeds += 1;
        int slice_length = 20000 - 2*seed_mlen;
        if (seed_bset_size > slice_length) {
            int slice_start = 0, slice_end = 0;
            while (slice_end < seed_bset_size) {
                if (slice_start + slice_length > seed_bset_size) { slice_end = seed_bset_size; }
                else { slice_end = slice_start + slice_length; }

                if (seed_mlen <= SMALL_MLEN_LIMIT) {
                    processSeedMotifWise(tuple<int, int> { seed_start + slice_start, seed_start + slice_end }, 0, seed_mlen,
                                         seed_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT],
                                         left_bset, right_bset, N_bset, out, aligner, filter, alignment, repeat_loci);
                }

                else {
                    processSeed(tuple<int, int> { seed_start + slice_start, seed_start + slice_end }, 0, seed_mlen,
                                seed_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT],
                                left_bset, right_bset, N_bset, out, lshift_xor_bsets, MATRIX,
                                aligner, filter, alignment, repeat_loci);
                }

                slice_start += slice_length - 500;
            }
        }

        else {

            if (seed_mlen <= SMALL_MLEN_LIMIT) {
                processSeedMotifWise(tuple<int, int> { seed_start, seed_end }, 0, seed_mlen, seed_type, sequence_id, sequence,
                                        sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                        out, aligner, filter, alignment, repeat_loci);
            }

            else {
                processSeed(tuple<int, int> { seed_start, seed_end }, 0, seed_mlen, seed_type, sequence_id, sequence, sequence_length,
                            lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                            out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci);
            }
        }
    }

    if (out) {
        int remove_loci_index = 0;
        if (repeat_loci.size() > 0) {
            for (int i=0; i<repeat_loci.size(); i++) {
                if (get<2> (repeat_loci[i]) >= chunk_end - SPLIT_OVERLAP) { continue; }
                if (get<1> (repeat_loci[i]) >= chunk_end - SPLIT_OVERLAP) { break; }
                out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) << "\t" << get<2> (repeat_loci[i]) << "\t"
                    << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<6> (repeat_loci[i]) << "\t" 
                    << get<7> (repeat_loci[i]) << "\t" << get<8> (repeat_loci[i]);
                if (CIGAROUTPUT) { out << "\t" << get<5> (repeat_loci[i]); }
                out << "\n";
                remove_loci_index += 1;
            }
        }
    
        if (remove_loci_index > 0) {
            repeat_loci.erase(repeat_loci.begin(), repeat_loci.begin() + remove_loci_index);
        }
    }

    seed_positions_perfect.clear();
    seed_positions_substut.clear();
    seed_positions_anchored.clear();
}


vector<tuple<size_t, size_t, string>> splitSequenceIntoBins(string& sequence, size_t bin_size = 5000000, size_t overlap = 50000) {
    /*
     * Splits a sequence into bins of specified size with specified overlap.
     * @param sequence: The DNA sequence to split.
     * @param bin_size: Size of each bin (default 5,000,000).
     * @param overlap: Number of bases each bin overlaps with the next (default 50,000).
     * @return vector of tuples: (start, end, bin_sequence) where start is inclusive, end is exclusive.
     */
    vector<tuple<size_t, size_t, string>> bins;
    size_t sequence_length = sequence.length();
    if (sequence_length == 0 || bin_size == 0) return bins;

    size_t start = 0;
    while (start < sequence_length) {
        size_t end = start + bin_size;
        if (end-start > sequence.length()) end = start + sequence.length();
        bins.emplace_back(start, end, sequence.substr(0, end - start));
        sequence.erase(0, end - start);
        if (end == sequence_length) break;
        start += (bin_size > overlap) ? (bin_size - overlap) : bin_size;
    }
    return bins;
}


void splitProcessSequence(const string& sequence_id, string& sequence, ostream &out) {
    /*
     * Splits a sequence into bins and writes them to the output stream.
     * @param sequence_id: ID of the sequence.
     * @param sequence: The DNA sequence to split.
     * @param out: Output stream to write the bins.
    */

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;
    int start = 0, end = 0;

    vector<tuple<size_t, size_t, string>> bins = splitSequenceIntoBins(sequence, SPLIT_LENGTH, SPLIT_OVERLAP);

    for (const auto& bin : bins) {
        start = get<0>(bin);
        end = get<1>(bin);
        string bin_sequence = get<2>(bin);

        if (THREADS > 1) {
            // create an output file for each thread and print the output of the thread to an output file

        }

        else {
            // process the chunk of the sequence and prints the output to the same file retaining the overlappint repeat loci
            processSequence(sequence_id, bin_sequence, out, start, end, repeat_loci);
        }
    }
    if (repeat_loci.size() > 0) {
        for (int i=0; i<repeat_loci.size(); i++) {
            out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) + start << "\t" << get<2> (repeat_loci[i]) + start << "\t"
                << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<6> (repeat_loci[i]) << "\t" 
                << get<7> (repeat_loci[i]) << "\t" << get<8> (repeat_loci[i]);
            if (CIGAROUTPUT) { out << "\t" << get<5> (repeat_loci[i]); }
            out << "\n";
        }
    }
}


void parseFasta(string fasta_file, string out_file) {
    /*
     * parses a fasta file input
     * @param fasta_file input fasta file name
     * @param out_file output file name
    */

    vector<string> sequence_ids;
    ifstream fastain(fasta_file);
    string line;

    // assigns the output to either a file or standard output
    streambuf * buf; ofstream outstream;

    // if the output file is not given by default: input file + ".ribbit"
    if (out_file == "") { out_file = fasta_file + ".ribbit"; }
    outstream.open(out_file);
    buf = outstream.rdbuf();    // output file buffer is created
    ostream out(buf);

    // adding header to the output file
    out << "#Chrom\t" << "Start\t" << "Stop\t" << "Motif\t" << "Purity\t"  << "Strand\t" << "Motif length\t"
        << "Repeat length\t" << "Repeat Units";
    if (CIGAROUTPUT) {out << "\tCigar"; } out << "\n";

    while (getline(fastain, line)) {
        if (line[0] == '>') {
            if (SEQUENCE != "") {
                std::cerr << "Processing " << SEQUENCE_ID << "\n";
                std::cerr << "Length of the sequence: " << SEQUENCE.length() << "\n";
                splitProcessSequence(SEQUENCE_ID, SEQUENCE, out);
            }
            SEQUENCE_ID = line.substr(1, line.find(' ') - 1);
            sequence_ids.push_back(SEQUENCE_ID);
            SEQUENCE = "";
        }
        else { SEQUENCE += line; }
    }

    if (SEQUENCE != "") {
        std::cerr << "Processing " << SEQUENCE_ID << "\n";
        std::cerr << "Length of the sequence: " << SEQUENCE.length() << "\n";
        splitProcessSequence(SEQUENCE_ID, SEQUENCE, out);
    }
    // processSequence(SEQUENCE_ID, SEQUENCE, out);
    fastain.close(); outstream.close();
}
