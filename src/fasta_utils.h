#include <thread>
#include <mutex>

#include "global_variables.h"
#include "bitseq_utils.h"
#include "concatenate_output.h"
#include "parse_perfect_shiftxor.h"
#include "parse_substitute_shiftxor.h"
#include "parse_anchored_shiftxor.h"
#include "parse_seed.h"
#include "parse_smallmotif_seed.h"

using namespace std;

// function to parse fasta index
void parseFai(string infai, int &nseqs, unordered_map<string, int> &seq_lens);

// function to parse the fasta file
void parseFasta(string fasta_file, string out_file);

void processSequence(string sequence_id, string &sequence, ostream &out, int chunk_start, int chunk_end,
                     vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci);
