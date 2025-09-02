#include <boost/program_options.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

namespace po = boost::program_options;

#include "global_variables.h"
#include "fasta_utils.h"
#include "concatenate_output.h"

using namespace std;


void parseFunctionArguments(int minimum_mlen = 2, int maximum_mlen = 100, double purity_threshold = 0.8,
                            double motif_purity_threshold = 0.8) {
    /*
     *  parsing input arguments for the program
     *  @param minimum_mlen minimum length of the motif
     *  @param maximum_mlen maximum length of the motif
     *  @param purity_threshold purity of the complete repeat
     *  @param motif_purity_threshold average match of each motif with consensus motif
     *  @return void
    */

    MINIMUM_MLEN = minimum_mlen;
    MAXIMUM_MLEN = maximum_mlen;
    PURITY_THRESHOLD = purity_threshold;
    MOTIFPURITY_THRESHOLD = motif_purity_threshold;
}


vector<tuple<string, int, int, string, double, string, int, int, int>> parseSequence(string sequence, string sequence_id="test", int minimum_mlen = 2,
                                                                                     int maximum_mlen = 100, double purity_threshold = 0.8, double motif_purity_threshold = 0.8) {
    /*
     *  parsing the sequence and identifying repeat loci
     *  @param sequence the sequence to be parsed
     *  @param sequence_id id of the sequence
     *  @param out output file stream to write results
     *  @param repeat_loci vector to store identified repeat loci
     *  @return void
    */

    int default_minimum_length = 12; // default minimum length for the motif

    parseFunctionArguments(minimum_mlen, maximum_mlen, purity_threshold, motif_purity_threshold);
    SEQUENCE = sequence;
    SEQUENCE_ID = sequence_id;

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;

    if (SEQUENCE.length() < 2*MINIMUM_MLEN) {
        cerr << "Sequence length is less than twice the minimum motif length. Skipping sequence: " << SEQUENCE_ID << "\n";
        return repeat_loci;
    }

    if (SEQUENCE.length() <  2*MAXIMUM_MLEN) {
        MAXIMUM_MLEN = SEQUENCE.length() / 2;
    }

    // uses minimum length of 12 as default if no input for minimum length or units are provided
    for (int key=MINIMUM_MLEN; key<=MAXIMUM_MLEN; key++) {            
        if (default_minimum_length < 2*key) {
            // if the minimum length is not atleast twice as the motif
            MINIMUM_LENGTH[key] = 2*key;
        }
        else MINIMUM_LENGTH[key] = default_minimum_length;
    }

    // for the motif sizes which are not in selected range but are factors of motif sizes
    // thresholds are set to the nearest selected motif size
    for (int m=MINIMUM_MLEN; m<=MAXIMUM_MLEN; m++) {
        vector<int> motif_factors;
        for (int f = 1; f <= m/2; f++) {
            // calculate all the factors of the motif size
            if (m % f == 0) { motif_factors.push_back(f); }
        }
        for (int f: motif_factors) {
            if (MINIMUM_LENGTH.find(f) == MINIMUM_LENGTH.end()) {
                // if threshold for factor not set; set the threshold
                MINIMUM_LENGTH[f] = MINIMUM_LENGTH[m];
            }
            if (PERFECT_UNITS.find(f) == PERFECT_UNITS.end()) {
                // if threshold for factor not set; set the threshold
                PERFECT_UNITS[f] = PERFECT_UNITS[m] * (m/f);
            }
        }
    }
    
    // minimum shift XOR to be generated; should be one less than the minimum motif size
    NMLENS = MAXIMUM_MLEN - MINIMUM_MLEN + 1;
    MINIMUM_SHIFT = (MINIMUM_MLEN > 2) ? MINIMUM_MLEN-2 : 1;

    // maximum anchor jump based on maximum motif size
    if      (MAXIMUM_MLEN <= 6)  { MAXIMUM_SHIFT = MAXIMUM_MLEN + 1; }
    else if (MAXIMUM_MLEN <= 20) { MAXIMUM_SHIFT = MAXIMUM_MLEN + (2 * MAXIMUM_MLEN) / 10; }
    else                         { MAXIMUM_SHIFT = MAXIMUM_MLEN + 4; }

    NSHIFTS = MAXIMUM_SHIFT - MINIMUM_SHIFT + 1;

    // Dynamically allocate memory for the matrix
    REPEAT_CLASSES = new uint32_t*[SMALL_MLEN_LIMIT];
    NUM_MOTIFS = pow(4, SMALL_MLEN_LIMIT);
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        REPEAT_CLASSES[i] = new uint32_t[NUM_MOTIFS];
    }
    MOTIF_START   = new int[NUM_MOTIFS];
    MOTIF_END     = new int[NUM_MOTIFS];
    MOTIF_NEXT    = new uint32_t[NUM_MOTIFS];
    MOTIF_UNITS   = new int[NUM_MOTIFS];
    MOTIF_GAPS    = new int[NUM_MOTIFS];
    MOTIF_GAPSIZE = new int[NUM_MOTIFS];

    SEEDLEN_CUTOFF = new int[NMLENS];
    for (int i = 0; i < NMLENS; ++i) {
        SEEDLEN_CUTOFF[i] = ((i+MINIMUM_MLEN) > SMALL_MLEN_LIMIT) ? 0.9*(i+MINIMUM_MLEN) : (MINIMUM_LENGTH[i+MINIMUM_MLEN]-(i+MINIMUM_MLEN));
    }

    // Initialize the matrix (optional)
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        for (int j = 0; j < NUM_MOTIFS; ++j) {
            REPEAT_CLASSES[i][j] = NUM_MOTIFS;
        }
    }
    
    ofstream out;
    processSequence(SEQUENCE_ID, SEQUENCE, out, 0, SEQUENCE.length(), repeat_loci);

    // Don't forget to free the memory when done
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        delete[] REPEAT_CLASSES[i];
    }
    delete[] REPEAT_CLASSES;
    delete[] MOTIF_START;
    delete[] MOTIF_END;
    delete[] MOTIF_NEXT;
    delete[] MOTIF_UNITS;
    delete[] MOTIF_GAPS;
    delete[] MOTIF_GAPSIZE;

    return repeat_loci;
}


PYBIND11_MODULE(ribbit, m) {
    /*
     *  Pybind11 module for ribbit
     *  @param m module name
    */
    m.doc() = "Ribbit: A tool to identify tandem repeats in DNA sequences";
    m.def("find_repeats", &parseSequence, "Parse a sequence and identify repeat loci",
          pybind11::arg("sequence"), pybind11::arg("sequence_id") = "test",
          pybind11::arg("minimum_mlen") = 2, pybind11::arg("maximum_mlen") = 100,
          pybind11::arg("purity_threshold") = 0.8, pybind11::arg("motif_purity_threshold") = 0.8);
}
