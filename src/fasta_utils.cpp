#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>
#include <numeric>
#include <algorithm>
#include <thread>
#include <mutex>

// local imports
#include "global_variables.h"
#include "fasta_utils.h"
#include "bitseq_utils.h"
#include "concatenate_output.h"
#include "parse_perfect_shiftxor.h"
#include "parse_substitute_shiftxor.h"
#include "parse_anchored_shiftxor.h"
#include "parse_seed.h"
#include "parse_smallmotif_seed.h"
#include "ssw_cpp.h"

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


void parseFasta(string fasta_file, int window_length, int window_bitcount_threshold,
                int anchor_length, int continuous_ones_threshold, string out_file) {
    /*
     * parses a fasta file input
     * @param fasta_file input fasta file name
     * @param window_length length of the window to be considered for picking seeds
     * @param window_bitcount_threshold threshold number of set bits in window
     * @param anchor_length length of continuous ones in anchor shift XORs
     * @param continuous_ones_threshold minimum continuous set bits in the window
     * @param out_file output file name
    */
    vector<string> seq_names;
    ifstream fastain(fasta_file);
    string line, seq_name, sequence="";

    if (THREADS == 1) {
        // assigns the output to either a file or standard output
        streambuf * buf; ofstream outstream;

        // if the output file is not given by default: input file + ".ribbit"
        if (out_file == "") { out_file = fasta_file + ".ribbit"; }
        outstream.open(out_file);
        buf = outstream.rdbuf();    // output file buffer is created
        ostream out(buf);

        // adding header to the output file
        out << "#Chrom\t" << "Start\t" << "Stop\t" << "Motif\t" << "Purity\t"  << "Strand\t" << "Cigar\t" << "Motif length\t"
            << "Repeat length\t" << "Repeat Units\n";

        while (getline(fastain, line)) {
            if (line[0] == '>') {
                if (sequence != "") {
                    // cerr << "\nProcessing sequence " << seq_name << "\n";
                    processSequence(seq_name, sequence, window_length, window_bitcount_threshold, anchor_length,
                                    continuous_ones_threshold, out);
                }
                seq_name = line.substr(1, line.find(' ') - 1);
                seq_names.push_back(seq_name);
                sequence = "";
            }
            else { sequence += line; }
        }
        processSequence(seq_name, sequence, window_length, window_bitcount_threshold, anchor_length,
                        continuous_ones_threshold, out);
        fastain.close(); outstream.close();
    }

    else {
        unordered_map<string, int> seq_lens; int nseqs = 0;
        parseFai(fasta_file+".fai", nseqs, seq_lens);
        if (nseqs == 0) {
            // cerr << "ERROR: Index for fasta file missing! Required when running in threads mode.\n";
            // cerr << "NOTE: Index can generated using `samtools faidx [fasta_file]`\n";
            return;
        }
        int chunk_size = 0;
        int tnum = 1;
        int toverlap = 100000;
        int seq_start = 0;

        std::vector<std::thread> threads;
        string output_name = "";

        while (getline(fastain, line)) {
            if (line[0] == '>') {
                if (sequence != "") {
                    threads.clear();
                    // cerr << "Processing sequence " << seq_name << "\n";
                    output_name = out_file + "_" + seq_name + "_" + to_string(tnum);
                    threads.emplace_back(processSequenceThread, seq_name, sequence, seq_start, window_length, window_bitcount_threshold,
                                        anchor_length, continuous_ones_threshold, tnum, output_name);
                    for (int _=0; _<THREADS; _++) { threads[_].join(); }
                }
                tnum = 1; seq_start = 0;
                seq_name = line.substr(1, line.find(' ') - 1);
                seq_names.push_back(seq_name);
                chunk_size = ((seq_lens[seq_name] + (THREADS*toverlap))/THREADS) + 10;
                sequence = "";
            }
            else {
                sequence += line;
                if (sequence.length() >= chunk_size) {
                    output_name = out_file + "_" + seq_name + "_" + to_string(tnum);
                    // cerr << "Writing thread "<< tnum << " output to " << output_name << "\n";
                    threads.emplace_back(processSequenceThread, seq_name, sequence, seq_start, window_length, window_bitcount_threshold,
                                         anchor_length, continuous_ones_threshold, tnum, output_name);
                    seq_start += sequence.length() - toverlap;
                    sequence = "" + sequence.substr(sequence.length() - (toverlap + 1), toverlap);
                    tnum += 1;
                }
            }
        }
        if (sequence != "") {
            output_name = out_file + "_" + seq_name + "_" + to_string(tnum);
            // cerr << "Writing thread "<< tnum << " output to " << output_name << "\n";
            threads.emplace_back(processSequenceThread, seq_name, sequence, seq_start, window_length, window_bitcount_threshold,
                                 anchor_length, continuous_ones_threshold, tnum, output_name);
            for (int _=0; _<THREADS; _++) { threads[_].join(); }
        }

        START_TIME = time(0);
        double seconds_since_start;
        if (THREADS > 1) {
            concatenateOutputs(out_file, seq_names, THREADS);
        }
        seconds_since_start = difftime(time(0), START_TIME);
        // std::cerr << "Concatenated all outputs.\t Time elapsed: " << seconds_since_start << "secs\n";
    }
}


void processSequence(string sequence_id, string sequence, int window_length, int window_bitcount_threshold,
                     int anchor_length, int continuous_ones_threshold, ostream &out) {
    /*
     *  processes each sequence from 2-bit conversion to identifying repeats
     *  @param sequence_id name of the sequence from fasta
     *  @param sequence string of the fasta sequence
     *  @param window_length length of the window to be scanned from the shift XOR
     *  @param window_bitcount_threshold threshold bitcount in the bit window
     *  @param anchor_length the length of the continuous stretch of 1s to be considered as anchor
     *  @param continuous_ones_threshold the minimum number of continuous stretch of 1s in seed
     *  @param out the output file to which the output has to be printed
     *  @return void generates the dynamic bitsets of shift XOR matches and
     *               proceeds to identifying repeats
    */

    START_TIME = time(0);
    double seconds_since_start;

    std::cerr << "Processing " << sequence_id << "\n";
    std::cerr << "Length of the sequence: " << sequence.length() << "\n";

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

    int seq_idx = 0, bidx; char nuc;
    for (; seq_idx < sequence_length; seq_idx++) {
        nuc = sequence[seq_idx];                // nucleotide
        bidx = (sequence_length-1) - seq_idx;   // bit indexing starts from right
        switch (nuc) {
            case 'A': case 'a': // 00
                left_bset[bidx] = 0; right_bset[bidx] = 0;
                MATRIX.push_back(&A); A[bidx] = 1;
                C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'C': case 'c': // 01
                left_bset[bidx] = 0; right_bset[bidx] = 1;
                MATRIX.push_back(&C); C[bidx] = 1;
                A[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'G': case 'g': // 10
                left_bset[bidx] = 1; right_bset[bidx] = 0;
                MATRIX.push_back(&G); G[bidx] = 1;
                A[bidx] = 0; C[bidx] = 0; T[bidx] = 0; break;
            case 'T': case 't': // 11
                left_bset[bidx] = 1; right_bset[bidx] = 1;
                MATRIX.push_back(&T); T[bidx] = 1;
                A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; break;
            default: // no match probably N or any other nuc
                N_bset[bidx] = 1; MATRIX.push_back(NULL);
                A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
        }
    }

    vector<boost::dynamic_bitset<>> lshift_xor_bsets;       // vector of dynamic bitsets for each shift XOR

    // generating the shift XORs from minimum shift size to maximum shift size
    for (int i = MINIMUM_SHIFT; i <= MAXIMUM_SHIFT; i++) {
        lshift_xor_bsets.push_back( ~(left_bset ^ (left_bset<<(i))) & ~(right_bset ^ (right_bset<<(i))) );
    }
    seconds_since_start = difftime( time(0), START_TIME);
    // std::cerr << "Generated shift XORs!\t Time elapsed:" << seconds_since_start << "secs\n";


    // generating seed positions; vector of tuple with start and end of the seeds
    vector<tuple<int, int, int, int>> seed_positions_perfect;
    vector<tuple<int, int, int, int>> seed_positions_substut;
    vector<tuple<int, int, int, int>> seed_positions_anchored;
    int failed_seeds = 0;

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;

    if (PURITY_THRESHOLD == 1) {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset, window_length);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Total number of perfect seeds: " << seed_positions_perfect.size() << "\t Time elapsed: " << seconds_since_start << "secs\n";
    }

    else {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset, window_length);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Total number of perfect seeds: " << seed_positions_perfect.size() << "\t Time elapsed: " << seconds_since_start << "secs\n";

        seed_positions_substut = processShiftXORswithSubstitutions(lshift_xor_bsets, N_bset, window_length, window_bitcount_threshold, seed_positions_perfect);
        failed_seeds = failedSeeds(seed_positions_perfect); failed_seeds += failedSeeds(seed_positions_substut);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Total number of seeds considering substitutions: " << seed_positions_perfect.size() + seed_positions_substut.size() - failed_seeds
                //   << "\t Time elapsed: " << seconds_since_start << "secs\n";


        // generating the anchor bitsets for all shift sizes
        vector<boost::dynamic_bitset<>> lsxor_anchor_bsets;     // vector of dynamic bitsets for anchor bitsets
        generateAnchoredShiftXORs(lshift_xor_bsets, N_bset, lsxor_anchor_bsets, anchor_length);
        boost::dynamic_bitset<> anchor_bset(sequence_length, 0ull);
        int motif_length = MINIMUM_MLEN;
        for (; motif_length <= MAXIMUM_MLEN; motif_length++) {
            anchor_bset.reset();

            int i = (motif_length > 2) ? motif_length - 2 : 1;
            for (; i <= motif_length + 2; i++) {
                int shift_idx = i - MINIMUM_SHIFT;
                // OR with actual shift XOR for same motif size
                if (i == motif_length) { anchor_bset |= lshift_xor_bsets[shift_idx]; }
                // OR with anchor bitset for neigboring shifts
                else { anchor_bset |= lsxor_anchor_bsets[shift_idx]; }
            }

            lshift_xor_bsets[motif_length-MINIMUM_SHIFT] = anchor_bset;
        }
        lsxor_anchor_bsets.clear();
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Generated anchored shift XORs!\t Time elapsed: " << seconds_since_start << "secs\n";

        window_bitcount_threshold = 6;  // threshold selected for identifying repeats with indels
        seed_positions_anchored = processShiftXORsAnchored(lshift_xor_bsets, N_bset, window_length, window_bitcount_threshold,
                                                           seed_positions_perfect, seed_positions_substut);
        seconds_since_start = difftime( time(0), START_TIME);
        failed_seeds = failedSeeds(seed_positions_perfect); failed_seeds += failedSeeds(seed_positions_substut); failed_seeds += failedSeeds(seed_positions_anchored);
        // std::cerr << "Total number of seeds considering indels: " << seed_positions_perfect.size() + seed_positions_substut.size() + seed_positions_anchored.size() - failed_seeds
                //   << "\t Time elapsed: " << seconds_since_start << "secs\n";
    }


    // Objects used by complete striped smithwater algorithm
    StripedSmithWaterman::Aligner   aligner;
    StripedSmithWaterman::Filter    filter;
    StripedSmithWaterman::Alignment alignment;

    // shift XORs for desired motif sizes; combination of shift XOR and anchor XOR

    tuple<int,int,int,int> seed;
    int seed_start, seed_end, seed_mlen, seed_type, seed_bset_size;
    uint64_t smallest; int smallest_type = -1;
    int spidx_p=0, spidx_s=0, spidx_a=0;
    int processed_seeds = 0;

    while (spidx_p < seed_positions_perfect.size() || spidx_s < seed_positions_substut.size() || spidx_a < seed_positions_anchored.size()) {
        smallest = -1;

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
            seed = seed_positions_perfect[spidx_p]; ++spidx_p;
        }
        else if (smallest_type == RANK_S) {
            seed = seed_positions_substut[spidx_s]; ++spidx_s;
        }
        else if (smallest_type == RANK_A) {
            seed = seed_positions_anchored[spidx_a]; ++spidx_a;
        }

        seed_type  = get<3> (seed);
        if (seed_type == -1) { continue; }
        seed_start = get<0> (seed);
        seed_end   = get<1> (seed);
        seed_mlen  = get<2> (seed);

        seed_bset_size = seed_end - seed_start;
        boost::dynamic_bitset<> seed_bset(seed_bset_size, 0ull);
        for (int j = seed_start; j < seed_end; j++) {
            seed_bset[seed_end - 1 - j] = lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT][sequence_length - 1 - j];
        }

        if (seed_end - seed_start >= 0.9*seed_mlen) {
            // process seed if it is alteast the size of the motif length
            processed_seeds += 1;

            if (seed_bset_size > 5000) {
                int slice_length = 5000, slice_overlap = ((int)(500/seed_mlen)) * seed_mlen;
                int slices = seed_bset_size / (slice_length - slice_overlap);
                if (seed_bset_size % (slice_length - slice_overlap) != 0) slices += 1;
                int slice_start, slice_end;
                for (int i=0; i<slices; i++) {
                    slice_start = seed_start + (i*(slice_length - slice_overlap));
                    if (i == slices - 1) slice_end = seed_end;
                    else slice_end = slice_start + slice_length;
                    
                    if (seed_mlen <= SMALL_MLEN_LIMIT) {
                        processSeedMotifWise(tuple<int, int> { slice_start, slice_end }, 0, seed_mlen, seed_type, sequence_id, sequence,
                                             sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                             continuous_ones_threshold, out, aligner, filter, alignment, repeat_loci);
                    }

                    else {
                        processSeed(tuple<int, int> { slice_start, slice_end }, 0, seed_mlen, seed_type, sequence_id, sequence, 
                                    sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                    continuous_ones_threshold, out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci);
                    }
                }
            }

            else {
                if (seed_mlen <= SMALL_MLEN_LIMIT) {
                    processSeedMotifWise(tuple<int, int> { seed_start, seed_end }, 0, seed_mlen, seed_type, sequence_id, sequence,
                                        sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                        continuous_ones_threshold, out, aligner, filter, alignment, repeat_loci);
                }

                else {
                    processSeed(tuple<int, int> { seed_start, seed_end }, 0, seed_mlen, seed_type, sequence_id, sequence, sequence_length,
                                lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset, continuous_ones_threshold,
                                out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci);
                }
            }
        }
    }

    if (repeat_loci.size() > 0) {
        for (int i=0; i<repeat_loci.size(); i++) {
            out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) << "\t"    << get<2> (repeat_loci[i]) << "\t"
                << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<5> (repeat_loci[i]) << "\t"
                << get<6> (repeat_loci[i]) << "\t" << get<7> (repeat_loci[i]) << "\t"    << get<8> (repeat_loci[i]) << "\n";
        }
    }

    seed_positions_perfect.clear();
    seed_positions_substut.clear();
    seed_positions_anchored.clear();

    seconds_since_start = difftime( time(0), START_TIME);
    // std::cerr << "Total number of seeds that are processed for alignment: " << processed_seeds << "\t Time elapsed: " << seconds_since_start << "secs\n\n";
    std::cerr << "Done!" << " Time elapsed: " << seconds_since_start << "secs\n\n";
}


void processSequenceThread(string sequence_id, string sequence, int seq_start, int window_length, int window_bitcount_threshold,
                           int anchor_length, int continuous_ones_threshold, int tnum, string out_file) {
    /*
     *  processes each sequence from 2-bit conversion to identifying repeats
     *  @param sequence_id name of the sequence from fasta
     *  @param sequence string of the fasta sequence
     *  @param seq_start the start coordinate of the sequence to be processed in the thread
     *  @param window_length length of the window to be scanned from the shift XOR
     *  @param window_bitcount_threshold threshold bitcount in the bit window
     *  @param anchor_length the length of the continuous stretch of 1s to be considered as anchor
     *  @param continuous_ones_threshold the minimum number of continuous stretch of 1s in seed
     *  @param out the output file to which the output has to be printed
     *  @return void generates the dynamic bitsets of shift XOR matches and
     *               proceeds to identifying repeats
    */

    streambuf * buf; ofstream outstream;
    outstream.open(out_file);
    buf = outstream.rdbuf();
    ostream out(buf);
    if (THREADS > 1) MTX.lock();
    START_TIME = time(0);
    if (THREADS > 1) MTX.unlock();
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

    int seq_idx = 0, bidx; char nuc;
    for (; seq_idx < sequence_length; seq_idx++) {
        nuc = sequence[seq_idx];                // nucleotide
        bidx = (sequence_length-1) - seq_idx;   // bit indexing starts from right
        switch (nuc) {
            case 'A': case 'a': // 00
                left_bset[bidx] = 0; right_bset[bidx] = 0;
                MATRIX.push_back(&A); A[bidx] = 1;
                C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'C': case 'c': // 01
                left_bset[bidx] = 0; right_bset[bidx] = 1;
                MATRIX.push_back(&C); C[bidx] = 1;
                A[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
            case 'G': case 'g': // 10
                left_bset[bidx] = 1; right_bset[bidx] = 0;
                MATRIX.push_back(&G); G[bidx] = 1;
                A[bidx] = 0; C[bidx] = 0; T[bidx] = 0; break;
            case 'T': case 't': // 11
                left_bset[bidx] = 1; right_bset[bidx] = 1;
                MATRIX.push_back(&T); T[bidx] = 1;
                A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; break;
            default: // no match probably N or any other nuc
                N_bset[bidx] = 1; MATRIX.push_back(NULL);
                A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0; break;
        }
    }

    vector<boost::dynamic_bitset<>> lshift_xor_bsets;       // vector of dynamic bitsets for each shift XOR

    // generating the shift XORs from minimum shift size to maximum shift size
    for (int i = MINIMUM_SHIFT; i <= MAXIMUM_SHIFT; i++) {
        lshift_xor_bsets.push_back( ~(left_bset ^ (left_bset<<(i))) & ~(right_bset ^ (right_bset<<(i))) );
    }
    seconds_since_start = difftime( time(0), START_TIME);
    // std::cerr << "Thread " << tnum << ": Generated shift XORs!\t Time elapsed:" << seconds_since_start << "secs\n";

    // generating seed positions; vector of tuple with start and end of the seeds
    vector<tuple<int, int, int, int>> seed_positions_perfect;
    vector<tuple<int, int, int, int>> seed_positions_substut;
    vector<tuple<int, int, int, int>> seed_positions_anchored;
    int failed_seeds = 0;

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;

    if (PURITY_THRESHOLD == 1) {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset, window_length);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Thread " << tnum << ": Total number of perfect seeds: " << seed_positions_perfect.size() << "\t Time elapsed: " << seconds_since_start << "secs\n";
    }

    else {
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset, window_length);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Thread " << tnum << ": Total number of perfect seeds: " << seed_positions_perfect.size() << "\t Time elapsed: " << seconds_since_start << "secs\n";

        seed_positions_substut = processShiftXORswithSubstitutions(lshift_xor_bsets, N_bset, window_length, window_bitcount_threshold, seed_positions_perfect);
        failed_seeds = failedSeeds(seed_positions_perfect); failed_seeds += failedSeeds(seed_positions_substut);
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Thread " << tnum << ": Total number of seeds considering substitutions: " << seed_positions_perfect.size() + seed_positions_substut.size() - failed_seeds
                // << "\t Time elapsed: " << seconds_since_start << "secs\n";


        // generating the anchor bitsets for all shift sizes
        vector<boost::dynamic_bitset<>> lsxor_anchor_bsets;     // vector of dynamic bitsets for anchor bitsets
        generateAnchoredShiftXORs(lshift_xor_bsets, N_bset, lsxor_anchor_bsets, anchor_length);
        boost::dynamic_bitset<> anchor_bset(sequence_length, 0ull);
        int motif_length = MINIMUM_MLEN;
        for (; motif_length <= MAXIMUM_MLEN; motif_length++) {
            anchor_bset.reset();

            int i = (motif_length > 2) ? motif_length - 2 : 1;
            for (; i <= motif_length + 2; i++) {
                int shift_idx = i - MINIMUM_SHIFT;
                // OR with actual shift XOR for same motif size
                if (i == motif_length) { anchor_bset |= lshift_xor_bsets[shift_idx]; }
                // OR with anchor bitset for neigboring shifts
                else { anchor_bset |= lsxor_anchor_bsets[shift_idx]; }
            }

            lshift_xor_bsets[motif_length-MINIMUM_SHIFT] = anchor_bset;
        }
        lsxor_anchor_bsets.clear();
        seconds_since_start = difftime( time(0), START_TIME);
        // std::cerr << "Thread " << tnum << ": Generated anchored shift XORs!\t Time elapsed: " << seconds_since_start << "secs\n";

        window_bitcount_threshold = 6;
        seed_positions_anchored = processShiftXORsAnchored(lshift_xor_bsets, N_bset, window_length, window_bitcount_threshold,
                                                        seed_positions_perfect, seed_positions_substut);
        seconds_since_start = difftime( time(0), START_TIME);
        failed_seeds = failedSeeds(seed_positions_perfect); failed_seeds += failedSeeds(seed_positions_substut); failed_seeds += failedSeeds(seed_positions_anchored);
        // std::cerr << "Thread " << tnum << ": Total number of seeds considering indels: "
                //   << seed_positions_perfect.size() + seed_positions_substut.size() + seed_positions_anchored.size() - failed_seeds
                //   << "\t Time elapsed: " << seconds_since_start << "secs\n";
    }


    // Objects used by complete striped smithwater algorithm
    StripedSmithWaterman::Aligner   aligner;
    StripedSmithWaterman::Filter    filter;
    StripedSmithWaterman::Alignment alignment;

    // shift XORs for desired motif sizes; combination of shift XOR and anchor XOR
    tuple<int,int,int,int> seed;
    int seed_start, seed_end, seed_mlen, seed_type, seed_bset_size;
    uint64_t smallest; int smallest_type = -1;
    int spidx_p=0, spidx_s=0, spidx_a=0;
    int processed_seeds = 0;

    // MTX.lock();

    while (spidx_p < seed_positions_perfect.size() || spidx_s < seed_positions_substut.size() || spidx_a < seed_positions_anchored.size()) {
        smallest = -1;

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
            seed = seed_positions_perfect[spidx_p]; ++spidx_p;
        }
        else if (smallest_type == RANK_S) {
            seed = seed_positions_substut[spidx_s]; ++spidx_s;
        }
        else if (smallest_type == RANK_A) {
            seed = seed_positions_anchored[spidx_a]; ++spidx_a;
        }

        seed_type  = get<3> (seed);
        if (seed_type == -1) { continue; }
        seed_start = get<0> (seed);
        seed_end   = get<1> (seed);
        seed_mlen  = get<2> (seed);

        seed_bset_size = seed_end - seed_start;
        boost::dynamic_bitset<> seed_bset(seed_bset_size, 0ull);
        for (int j = seed_start; j < seed_end; j++) {
            seed_bset[seed_end - 1 - j] = lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT][sequence_length - 1 - j];
        }

        if (seed_end - seed_start >= 0.9*seed_mlen) {

            processed_seeds += 1;

            if (seed_mlen <= SMALL_MLEN_LIMIT) {
                processSeedMotifWise(tuple<int, int> { seed_start, seed_end }, seq_start, seed_mlen, seed_type, sequence_id, sequence,
                                     sequence_length, lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                     continuous_ones_threshold, out, aligner, filter, alignment, repeat_loci);
            }

            else {
                processSeed(tuple<int, int> { seed_start, seed_end }, seq_start, seed_mlen, seed_type, sequence_id, sequence, sequence_length,
                            lshift_xor_bsets[seed_mlen-MINIMUM_SHIFT], left_bset, right_bset, N_bset, continuous_ones_threshold,
                            out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci);
            }
        }
    }
    // MTX.unlock();

    if (repeat_loci.size() > 0) {
        for (int i=0; i<repeat_loci.size(); i++) {
            out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) << "\t"    << get<2> (repeat_loci[i]) << "\t"
                << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<5> (repeat_loci[i]) << "\t"
                << get<6> (repeat_loci[i]) << "\t" << get<7> (repeat_loci[i]) << "\t"    << get<8> (repeat_loci[i]) << "\n";
        }
    }

    seed_positions_perfect.clear();
    seed_positions_substut.clear();
    seed_positions_anchored.clear();

    outstream.close();
    seconds_since_start = difftime( time(0), START_TIME);
    // std::cerr << "Thread " << tnum << ": Total number of seeds that are processed for alignment: " << processed_seeds << "\t Time elapsed: " << seconds_since_start << "secs\n";
}
