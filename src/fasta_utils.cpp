#include "fasta_utils.h"

using namespace std;


void buildBitDatastructures(string &sequence, int sequence_length, boost::dynamic_bitset<> &left_bset, boost::dynamic_bitset<> &right_bset,
                            boost::dynamic_bitset<> &N_bset, boost::dynamic_bitset<> &A, boost::dynamic_bitset<> &T, boost::dynamic_bitset<> &G,
                            boost::dynamic_bitset<> &C, vector<boost::dynamic_bitset<> *> &MATRIX) {
    /*
     * builds the bit datastructures for the sequence
     * @param sequence the input DNA sequence
     * @param sequence_length length of the input sequence
     * @param left_bset dynamic bitset for left side of the sequence
     * @param right_bset dynamic bitset for right side of the sequence
     * @param N_bset dynamic bitset for N nucleotides
     * @param A, T, G, C dynamic bitsets for A, T, G, C nucleotides respectively
     * @param MATRIX vector of pointers to the nucleotide bitsets
    */

    int seq_idx = 0, bidx;
    char nuc;
    for (; seq_idx < sequence_length; seq_idx++) {
        nuc = sequence[seq_idx];                // nucleotide
        bidx = (sequence_length - 1) - seq_idx; // bit indexing starts from right
        switch (nuc) {
        case 'A': case 'a': // 00
            left_bset[bidx] = 0; right_bset[bidx] = 0;
            MATRIX.push_back(&A);
            A[bidx] = 1; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0;
            break;
        case 'C': case 'c': // 01
            left_bset[bidx] = 0; right_bset[bidx] = 1;
            MATRIX.push_back(&C);
            A[bidx] = 0; C[bidx] = 1; G[bidx] = 0; T[bidx] = 0;
            break;
        case 'G': case 'g': // 10
            left_bset[bidx] = 1; right_bset[bidx] = 0;
            MATRIX.push_back(&G);
            A[bidx] = 0; C[bidx] = 0; G[bidx] = 1; T[bidx] = 0;
            break;
        case 'T': case 't': // 11
            left_bset[bidx] = 1; right_bset[bidx] = 1;
            MATRIX.push_back(&T);
            A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; T[bidx] = 1;
            break;
        default: // no match probably N or any other nuc
            N_bset[bidx] = 1;
            MATRIX.push_back(NULL);
            A[bidx] = 0; C[bidx] = 0; G[bidx] = 0; T[bidx] = 0;
            break;
        }
    }
}


void processSequence(string sequence_id, string sequence, ofstream &out, ofstream &seeds_out, int chunk_start, int chunk_end,
                     vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     * processes each sequence from 2-bit conversion to identifying repeats
     * @param sequence_id name of the sequence from fasta
     * @param sequence string of the fasta sequence
     * @param out the output file to which the output has to be printed
     * @param chunk_start the start position of the chunk being processed
     * @param chunk_end the end position of the chunk being processed
     * @param repeat_loci vector to store the repeat loci information
     * @return void generates the dynamic bitsets of shift XOR matches and
     * proceeds to identifying repeats
     */

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
    vector<boost::dynamic_bitset<> *> MATRIX;

    // builds the bit datastructures that are necessary
    buildBitDatastructures(sequence, sequence_length, left_bset, right_bset, N_bset, A, T, G, C, MATRIX);

    vector<boost::dynamic_bitset<>> lshift_xor_bsets; // vector of dynamic bitsets for each shift XOR

    // generating the anchor bitsets for all shift sizes
    vector<boost::dynamic_bitset<>> lsxor_anchor_bsets;   // vector of dynamic bitsets for anchor bitsets
    vector<boost::dynamic_bitset<>> lsxor_perfect_bsets;  // vector of dynamic bitsets for perfect bitsets
    vector<boost::dynamic_bitset<>> lsxor_anchored_bsets; // vector of dynamic bitsets for each shift XOR with added anchors

    // generating the shift XORs from minimum shift size to maximum shift size
    for (int i = MINIMUM_SHIFT; i <= MAXIMUM_SHIFT; i++) {
        lshift_xor_bsets.push_back(~(left_bset ^ (left_bset << (i))) & ~(right_bset ^ (right_bset << (i))));
    }

    // generating seed positions; vector of tuple with start and end of the seeds
    vector<tuple<int, int, int, int, int, int, int>> seed_positions_perfect;
    vector<tuple<int, int, int, int, int, int, int>> seed_positions_substut;
    vector<tuple<int, int, int, int, int, int, int>> seed_positions_anchored;

    if (PURITY_THRESHOLD == 1) {
        // Only identify perfect repeat sequences if the purity threshold is set to 1
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset);
    }

    else {
        
        int anchor_length = 5;
        generateAnchorShiftXORs(lshift_xor_bsets, N_bset, lsxor_anchor_bsets, anchor_length);
        generatePerfectShiftXORs(lshift_xor_bsets, N_bset, lsxor_perfect_bsets, anchor_length);
        
        int anchor_jump = 1;
        boost::dynamic_bitset<> anchor_bset(sequence_length, 0ull);
        int motif_length = MINIMUM_MLEN;
        for (; motif_length <= MAXIMUM_MLEN; motif_length++) {
            anchor_bset.reset();
            
            if (motif_length <= 6) { anchor_jump = 1; }
            else {
                anchor_jump = ((2 * motif_length / 10) > 1) ? (2 * motif_length / 10) : 2;
            }
            
            int i = motif_length - anchor_jump;
            if (i < MINIMUM_SHIFT) {
                i = MINIMUM_SHIFT;
            }
            for (; i <= motif_length + anchor_jump; i++) {
                int shift_idx = i - MINIMUM_SHIFT;
                // OR with actual shift XOR for same motif size
                if (i == motif_length) {
                    anchor_bset |= lshift_xor_bsets[shift_idx];
                    anchor_bset |= (lsxor_perfect_bsets[shift_idx] << motif_length);
                }
                // OR with anchor bitset for neigboring shifts
                else {
                    anchor_bset |= lsxor_anchor_bsets[shift_idx];
                }
            }
            
            lsxor_anchored_bsets.push_back(anchor_bset);
        }
        
        // freeing up memory and keeping only the final anchored bitsets
        lsxor_anchor_bsets.clear();

        // Identifying perfect repeat seeds
        seed_positions_perfect = processShiftXORsPerfect(lshift_xor_bsets, N_bset);

        // Identifying repeat seeds with only allowed substitutions or mismatches between motifs
        seed_positions_substut = processShiftXORswithSubstitutions(lshift_xor_bsets, lsxor_perfect_bsets, lsxor_anchored_bsets, N_bset, seed_positions_perfect);
        
        // filtering out the perfect seeds which are inside substituted seeds with short flanks
        // filters out the perfect seeds which are contained within substitute seeds of almost same length
        filterPerfectSeeds(seed_positions_perfect, seed_positions_substut);


        seed_positions_anchored = processShiftXORsAnchored(lsxor_anchored_bsets, lshift_xor_bsets, lsxor_perfect_bsets,
                                                            N_bset, seed_positions_perfect, seed_positions_substut, seeds_out, chunk_start);
        filterShortSeeds(seed_positions_perfect);
        filterShortSeeds(seed_positions_substut);
        filterShortSeeds(seed_positions_anchored);
        // return;
    }

    // Objects used by complete striped smithwater algorithm
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    tuple<int, int, int, int, int, int, int> seed;
    int seed_start, seed_end, seed_mlen, seed_type, seed_bset_size;
    uint64_t smallest;
    int smallest_type = -1;
    int spidx_p = 0, spidx_s = 0, spidx_a = 0;
    int processed_seeds = 0;
    int seedlen_cutoff = 0;

    vector<tuple<int, int, int, int, int, int, int>> overlapping_seeds;
    vector<set<int>> skip_atomicity;
    int overlap_end = -1;

    while (spidx_p < seed_positions_perfect.size() || spidx_s < seed_positions_substut.size() || spidx_a < seed_positions_anchored.size()) {
        smallest = -1;
        smallest_type = RANK_N;

        // considers the seed based on the start position of the seed
        while (spidx_p < seed_positions_perfect.size() && get<3>(seed_positions_perfect[spidx_p]) == RANK_N) {
            spidx_p += 1;
        }
        while (spidx_s < seed_positions_substut.size() && get<3>(seed_positions_substut[spidx_s]) == RANK_N) {
            spidx_s += 1;
        }
        while (spidx_a < seed_positions_anchored.size() && get<3>(seed_positions_anchored[spidx_a]) == RANK_N) {
            spidx_a += 1;
        }

        // Find the smallest element among the current elements of the three vectors
        if (spidx_p < seed_positions_perfect.size() && (smallest > get<0>(seed_positions_perfect[spidx_p]))) {
            smallest = get<0>(seed_positions_perfect[spidx_p]);
            smallest_type = RANK_P;
        }
        if (spidx_s < seed_positions_substut.size() && (smallest > get<0>(seed_positions_substut[spidx_s]))) {
            smallest = get<0>(seed_positions_substut[spidx_s]);
            smallest_type = RANK_S;
        }
        if (spidx_a < seed_positions_anchored.size() && (smallest > get<0>(seed_positions_anchored[spidx_a]))) {
            smallest = get<0>(seed_positions_anchored[spidx_a]);
            smallest_type = RANK_A;
        }
        // Print the smallest element and move the corresponding pointer
        if (smallest_type == RANK_P) {
            seed = seed_positions_perfect[spidx_p];
            spidx_p += 1;
        }
        else if (smallest_type == RANK_S) {
            seed = seed_positions_substut[spidx_s];
            spidx_s += 1;
        }
        else if (smallest_type == RANK_A) {
            seed = seed_positions_anchored[spidx_a];
            spidx_a += 1;
        }

        if (smallest_type == RANK_N) { continue; } // skip if no seeds left

        seed_type = get<3>(seed);
        if (seed_type == RANK_N) { continue; }
        
        seed_start = get<0>(seed);
        seed_end = get<1>(seed);
        if (seed_start < 0 || seed_end < 0) { continue;}
        if (seed_end + seed_mlen > sequence_length) {
            seed_end = sequence_length - seed_mlen;
        }
        seed_mlen = get<2>(seed);
        if (overlap_end == -1) { overlapping_seeds.push_back(seed); overlap_end = seed_end + seed_mlen; }
        else {
            if (seed_start <= overlap_end) {
                overlapping_seeds.push_back(seed);
                if (seed_end + seed_mlen > overlap_end) overlap_end = seed_end + seed_mlen;
            }
            
            else {
                // process the overlapping seeds to remove redundant seeds
                processOverlappingSeeds(overlapping_seeds, lshift_xor_bsets, lsxor_perfect_bsets, lsxor_anchored_bsets, sequence_length, skip_atomicity);

                for (int j = 0; j < overlapping_seeds.size(); j++) {
                    if (get<3>(overlapping_seeds[j]) == RANK_N) { continue; }
                    
                    seeds_out << sequence_id << "\t" << get<0>(overlapping_seeds[j]) + chunk_start << "\t"
                              << get<1>(overlapping_seeds[j]) + chunk_start << "\t"
                              << get<2>(overlapping_seeds[j]) << "\t" << get<3>(overlapping_seeds[j]) << "\t" << get<4>(overlapping_seeds[j]) << "\t"
                              << get<5>(overlapping_seeds[j]) << "\t" << get<6>(overlapping_seeds[j]) << "\t";
                    for (int s=0; s<skip_atomicity[j].size(); s++) {
                        if (s != 0) seeds_out << ",";
                        seeds_out << *(next(skip_atomicity[j].begin(), s));
                    }
                    seeds_out << "\n";

                    int o_start = get<0>(overlapping_seeds[j]);
                    int o_end = get<1>(overlapping_seeds[j]);
                    int o_mlen = get<2>(overlapping_seeds[j]);
                    if (o_start < 0 || o_end < 0) { continue; }
                    if (o_end + o_mlen > sequence_length) {
                        o_end = sequence_length - o_mlen;
                    }
                    int o_type = get<3>(overlapping_seeds[j]);

                    int o_bset_size = o_end - o_start;
                    if (THREADS > 1) MTX.lock();
                    seedlen_cutoff = SEEDLEN_CUTOFF[o_mlen - MINIMUM_MLEN];
                    if (THREADS > 1) MTX.unlock();
                    if (o_bset_size < seedlen_cutoff) { continue; }
    
                    // process seed if it is alteast the size of the motif length
                    processed_seeds += 1;
                    int slice_length = 20000 - 2 * o_mlen;
                    if (o_bset_size > slice_length) {
                        int slice_start = 0, slice_end = 0;
                        while (slice_end < o_bset_size) {
                            if (slice_start + slice_length > o_bset_size)
                                slice_end = o_bset_size;
                            else
                                slice_end = slice_start + slice_length;

                            if (o_mlen <= SMALL_MLEN_LIMIT) {
                                processSeedMotifWise(tuple<int, int>{o_start + slice_start, o_start + slice_end}, chunk_start, o_mlen,
                                                    o_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT],
                                                    left_bset, right_bset, N_bset, out, aligner, filter, alignment, repeat_loci);
                            }
    
                            else {
                                processSeed(tuple<int, int>{o_start + slice_start, o_start + slice_end}, chunk_start, o_mlen,
                                            o_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT],
                                            left_bset, right_bset, N_bset, out, lshift_xor_bsets, MATRIX,
                                            aligner, filter, alignment, repeat_loci, skip_atomicity[j]);
                            }
    
                            slice_start += slice_length - 500;
                        }
                    }
    
                    else {

                        if (o_mlen <= SMALL_MLEN_LIMIT) {
                            processSeedMotifWise(tuple<int, int>{o_start, o_end}, chunk_start, o_mlen, o_type, sequence_id,
                                                sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT], left_bset, right_bset,
                                                N_bset, out, aligner, filter, alignment, repeat_loci);
                        }
    
                        else {
                            processSeed(tuple<int, int>{o_start, o_end}, chunk_start, o_mlen, o_type, sequence_id, sequence,
                                        sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                        out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci, skip_atomicity[j]);
                        }
                    }
                }

                overlapping_seeds.clear();
                skip_atomicity.clear();
                overlapping_seeds.push_back(seed); overlap_end = seed_end + seed_mlen;
            }
        }
    }

    if (overlapping_seeds.size() > 0) {
        processOverlappingSeeds(overlapping_seeds, lshift_xor_bsets, lsxor_perfect_bsets, lsxor_anchored_bsets, sequence_length, skip_atomicity);

        for (int j = 0; j < overlapping_seeds.size(); j++) {
            if (get<3>(overlapping_seeds[j]) == RANK_N) { continue; }
            
            seeds_out << sequence_id << "\t" << get<0>(overlapping_seeds[j]) + chunk_start << "\t"
                        << get<1>(overlapping_seeds[j]) + chunk_start << "\t"
                        << get<2>(overlapping_seeds[j]) << "\t" << get<3>(overlapping_seeds[j]) << "\t" << get<4>(overlapping_seeds[j]) << "\t"
                        << get<5>(overlapping_seeds[j]) << "\t" << get<6>(overlapping_seeds[j]) << "\t";
            for (int s=0; s<skip_atomicity[j].size(); s++) {
                if (s != 0) seeds_out << ",";
                seeds_out << *(next(skip_atomicity[j].begin(), s));
            }
            seeds_out << "\n";

            int o_start = get<0>(overlapping_seeds[j]);
            int o_end = get<1>(overlapping_seeds[j]);
            int o_mlen = get<2>(overlapping_seeds[j]);
            if (o_start < 0 || o_end < 0) { continue; }
            if (o_end + o_mlen > sequence_length) {
                o_end = sequence_length - o_mlen;
            }
            int o_type = get<3>(overlapping_seeds[j]);

            int o_bset_size = o_end - o_start;
            if (THREADS > 1) MTX.lock();
            seedlen_cutoff = SEEDLEN_CUTOFF[o_mlen - MINIMUM_MLEN];
            if (THREADS > 1) MTX.unlock();
            if (o_bset_size < seedlen_cutoff) { continue; }

            // process seed if it is alteast the size of the motif length
            processed_seeds += 1;
            int slice_length = 20000 - 2 * o_mlen;
            if (o_bset_size > slice_length) {
                int slice_start = 0, slice_end = 0;
                while (slice_end < o_bset_size) {
                    if (slice_start + slice_length > o_bset_size)
                        slice_end = o_bset_size;
                    else
                        slice_end = slice_start + slice_length;

                    if (o_mlen <= SMALL_MLEN_LIMIT) {
                        processSeedMotifWise(tuple<int, int>{o_start + slice_start, o_start + slice_end}, chunk_start, o_mlen,
                                            o_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT],
                                            left_bset, right_bset, N_bset, out, aligner, filter, alignment, repeat_loci);
                    }

                    else {
                        processSeed(tuple<int, int>{o_start + slice_start, o_start + slice_end}, chunk_start, o_mlen,
                                    o_type, sequence_id, sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT],
                                    left_bset, right_bset, N_bset, out, lshift_xor_bsets, MATRIX,
                                    aligner, filter, alignment, repeat_loci, skip_atomicity[j]);
                    }

                    slice_start += slice_length - 500;
                }
            }

            else {

                if (o_mlen <= SMALL_MLEN_LIMIT) {
                    processSeedMotifWise(tuple<int, int>{o_start, o_end}, chunk_start, o_mlen, o_type, sequence_id,
                                        sequence, sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT], left_bset, right_bset,
                                        N_bset, out, aligner, filter, alignment, repeat_loci);
                }

                else {
                    processSeed(tuple<int, int>{o_start, o_end}, chunk_start, o_mlen, o_type, sequence_id, sequence,
                                sequence_length, lshift_xor_bsets[o_mlen - MINIMUM_SHIFT], left_bset, right_bset, N_bset,
                                out, lshift_xor_bsets, MATRIX, aligner, filter, alignment, repeat_loci, skip_atomicity[j]);
                }
            }
        }

    }

    if (THREADS > 1) {
        // if multi threaded print till the end of the chunk into the temp file
        if (repeat_loci.size() > 0) {
            printRepeatsToOutput(out, repeat_loci, repeat_loci.size()-1);
        }
    }

    else {
        if (repeat_loci.size() > 0) printRepeatsToOutput(out, repeat_loci, repeat_loci.size()-1);
    }

    seed_positions_perfect.clear();
    seed_positions_substut.clear();
    seed_positions_anchored.clear();
}


vector<tuple<size_t, size_t>> splitSequenceIntoBins(string &sequence, size_t bin_size = 5000000, size_t overlap = 50000) {
    /*
     * Splits a sequence into bins of specified size with specified overlap.
     * @param sequence: The DNA sequence to split.
     * @param bin_size: Size of each bin (default 5,000,000).
     * @param overlap: Number of bases each bin overlaps with the next (default 50,000).
     * @return vector of tuples: (start, end, bin_sequence) where start is inclusive, end is exclusive.
     */
    vector<tuple<size_t, size_t>> bins;
    size_t sequence_length = sequence.length();
    if (sequence_length == 0 || bin_size == 0)
        return bins;

    size_t start = 0;
    while (start < sequence_length) {
        size_t end = start + bin_size;
        if (end > sequence_length)
            end = sequence_length;
        bins.emplace_back(start, end);
        if (end == sequence_length)
            break;
        start += (bin_size > overlap) ? (bin_size - overlap) : bin_size;
    }
    return bins;
}


void splitProcessSequence(const string &sequence_id, string &sequence, ofstream &out, ofstream &seeds_out, string output_file) {
    /*
     * Splits a sequence into bins and writes them to the output stream.
     * @param sequence_id: ID of the sequence.
     * @param sequence: The DNA sequence to split.
     * @param out: Output stream to write the bins.
     * @param seeds_out: Output stream to write the seed regions.
     */

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;
    int start = 0, end = 0;

    vector<tuple<size_t, size_t>> bins = splitSequenceIntoBins(sequence, SPLIT_LENGTH, SPLIT_OVERLAP);

    if (THREADS > 1) {
        // create an output file for each thread and print the output of the thread to an output file

        int nbins = bins.size();
        std::vector<std::thread> threads;
        std::vector<std::string> temp_files(nbins);
        std::vector<std::ofstream> temp_streams(nbins);

        std::vector<std::string> temp_seed_files(nbins);
        std::vector<std::ofstream> temp_seed_streams(nbins);

        ofstream seq_out(output_file + '.' + sequence_id);
        ofstream seq_seeds_out(output_file + '.' + sequence_id + ".seeds");

        bool original_CIGAROUTPUT = CIGAROUTPUT;
        CIGAROUTPUT = true; // enable cigar output for temp files to merge them properly

        for (int i = 0; i < nbins; i++) {
            temp_files[i] = output_file + "." + sequence_id + ".thread" + std::to_string(i) + ".tmp";
            temp_streams[i].open(temp_files[i]);

            temp_seed_files[i] = output_file + "." + sequence_id + ".thread" + std::to_string(i) + ".seeds.tmp";
            temp_seed_streams[i].open(temp_seed_files[i]);

        }

        for (size_t i = 0; i < nbins; ++i) {
            start = get<0>(bins[i]);
            end = get<1>(bins[i]);
            // Each thread processes its bin and writes to its temp file
            threads.emplace_back([&, i, start, end]() {
                vector<tuple<string, int, int, string, double, string, int, int, int>> thread_repeat_loci;
                processSequence(sequence_id, sequence.substr(start, end - start), temp_streams[i], temp_seed_streams[i], start, end, thread_repeat_loci); });
        }

        // Wait for all threads to finish
        for (auto &t : threads) { t.join(); }

        for (auto &ts : temp_streams) { ts.close(); }

        CIGAROUTPUT = original_CIGAROUTPUT; // restore original CIGAROUTPUT setting

        concatenateThreadOutputs(temp_files, seq_out);

        string line;
        for (string thread_seeds_outfile: temp_seed_files) {
            ifstream thread_seeds_out(thread_seeds_outfile);
            while(getline(thread_seeds_out, line)) {
                seq_seeds_out << line << "\n";
            }
            thread_seeds_out.close();
            remove(thread_seeds_outfile.c_str());
        }
        seq_out.close();
        seq_seeds_out.close();
    }
    else {
        for (const auto &bin : bins) {
            start = get<0>(bin);
            end = get<1>(bin);
            // process the chunk of the sequence and prints the output to the same file retaining the overlapping repeat loci
            processSequence(sequence_id, sequence.substr(start, end - start), out, seeds_out, start, end, repeat_loci);
        }
        if (repeat_loci.size() > 0) {
            for (int i = 0; i < repeat_loci.size(); i++) {
                out << get<0>(repeat_loci[i]) << "\t" << get<1>(repeat_loci[i]) + start << "\t" << get<2>(repeat_loci[i]) + start << "\t"
                    << get<3>(repeat_loci[i]) << "\t" << std::setprecision(2) << (get<4>(repeat_loci[i])) << "\t" << get<6>(repeat_loci[i]) << "\t"
                    << get<7>(repeat_loci[i]) << "\t" << get<8>(repeat_loci[i]);
                if (CIGAROUTPUT) {
                    out << "\t" << get<5>(repeat_loci[i]);
                }
                out << "\n";
            }
        }
    }
}


void parseFasta(string fasta_file, string output_file) {
    /*
     * parses a fasta file input
     * @param fasta_file input fasta file name
     * @param output_file output file name
     */

    vector<string> sequence_ids;
    bool is_gzipped = false;
    if (fasta_file.size() > 3 && fasta_file.substr(fasta_file.size() - 3) == ".gz") {
        is_gzipped = true;
    }

    string line;
    // if the output file is not given by default: input file + ".ribbit"
    if (output_file == "") { output_file = fasta_file + ".ribbit"; }
    string seeds_output_file = output_file + ".seeds";
    ofstream out(output_file);
    ofstream seeds_out(seeds_output_file);

    // adding header to the output file
    out << "#chrom\t" << "start\t" << "stop\t" << "motif\t" << "purity\t" << "motif_length\t"
        << "repeat_length\t" << "repeat_units";
    if (CIGAROUTPUT) { out << "\tcigar"; }
    out << "\n";

    if (is_gzipped) {
        gzFile gzfin = gzopen(fasta_file.c_str(), "rb");
        if (!gzfin) {
            cerr << "Error opening gzipped fasta file: " << fasta_file << endl;
            return;
        }
        char buffer[65536];
        string gzline;
        while (gzgets(gzfin, buffer, sizeof(buffer))) {
            gzline = buffer;
            // Remove trailing newline
            if (!gzline.empty() && gzline.back() == '\n') gzline.pop_back();
            if (!gzline.empty() && gzline.back() == '\r') gzline.pop_back();
            if (gzline[0] == '>') {
                if (SEQUENCE != "") {
                    std::cerr << "Processing " << SEQUENCE_ID << "\n";
                    std::cerr << "Length of the sequence: " << SEQUENCE.length() << "\n";
                    splitProcessSequence(SEQUENCE_ID, SEQUENCE, out, seeds_out, output_file);
                }
                SEQUENCE_ID = gzline.substr(1, gzline.find(' ') - 1);
                sequence_ids.push_back(SEQUENCE_ID);
                SEQUENCE = "";
            } else {
                SEQUENCE += gzline;
            }
        }
        gzclose(gzfin);
    }
    
    else {
        ifstream fastain(fasta_file);
        if (!fastain) {
            cerr << "Error opening fasta file: " << fasta_file << endl;
            return;
        }
        while (getline(fastain, line)) {
            if (line[0] == '>') {
                if (SEQUENCE != "") {
                    std::cerr << "Processing " << SEQUENCE_ID << "\n";
                    std::cerr << "Length of the sequence: " << SEQUENCE.length() << "\n";
                    splitProcessSequence(SEQUENCE_ID, SEQUENCE, out, seeds_out, output_file);
                }
                SEQUENCE_ID = line.substr(1, line.find(' ') - 1);
                sequence_ids.push_back(SEQUENCE_ID);
                SEQUENCE = "";
            }
            else {
                SEQUENCE += line;
            }
        }
        fastain.close();
    }

    if (SEQUENCE != "") {
        std::cerr << "Processing " << SEQUENCE_ID << "\n";
        std::cerr << "Length of the sequence: " << SEQUENCE.length() << "\n";
        splitProcessSequence(SEQUENCE_ID, SEQUENCE, out, seeds_out, output_file);
    }

    if (THREADS > 1) {
        for (string sequence_id: sequence_ids) {
            string seq_outfile = output_file + '.' + sequence_id;
            ifstream seq_out(seq_outfile);
            while(getline(seq_out, line)) {
                out << line << "\n";
            }
            seq_out.close();
            remove(seq_outfile.c_str());

            string seq_seeds_outfile = output_file + '.' + sequence_id + ".seeds";
            ifstream seq_seeds_out(seq_seeds_outfile);
            while(getline(seq_seeds_out, line)) {
                seeds_out << line << "\n";
            }
            seq_seeds_out.close();
            remove(seq_seeds_outfile.c_str());
        }
    }

    out.close();
    seeds_out.close();
}
