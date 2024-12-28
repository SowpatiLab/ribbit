#include <iostream>
#include <bitset>
#include <string>
#include <unordered_map>
#include <mutex>

using namespace std;

// function to parse fasta index
void parseFai(string infai, int &nseqs, unordered_map<string, int> &seq_lens);

// // process sequence using dynamic bitset
void processSequence(string sequence_id, string sequence, int window_length, int window_bitcount_threshold, int anchor_size,
                     int continuous_threshold, ostream &out);

// // process sequence using dynamic bitset
void processSequenceThread(string sequence_id, string sequence, int seq_start, int window_length, int window_bitcount_threshold,
                           int anchor_size, int continuous_threshold, int tnum, string out_file);
