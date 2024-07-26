#include <stdio.h>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "global_variables.h"
#include "fasta_utils.h"

using namespace std;


bool isNumber(const string &s) {
    /*
     *  checks if a string is number or not
     *  @param s string to be checked if it is numeric
     *  @return bool if the string is numeric or not
    */ 
    return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}


bool parseDualtypeArgs(po::variables_map &args, const string &option, unordered_map<int, int> &cutoff,
                         int &minimum_motif_length, int &maximum_motif_length) {
    /*
     *  parsing input for arguments with either an integer or file options
     *  @param args arguments object from the program options object
     *  @param option name of the option
     *  @param cutoff unordered_map the option to be updated for each motif length
     *  @param minimum_motif_length minimum length of the motif
     *  @param maximum_motif_length maximum length of the motif
     *  @return bool for successful completion of the function
    */
    int key, value;    
    if (isNumber(args[option].as<string>())) {
        // if the input is just a number; set the same cutoff for all motif lengths
        value = stoi(args[option].as<string>());
        for (key = minimum_motif_length; key<=maximum_motif_length; key++) {
            cutoff[key] = value;
        }
    }

    else {
        // if the input is a file; take assigned inputs for each motif size
        ifstream infile;
        if (infile.fail()) {
            cerr << "Please provide a valid file for " << option << ".\n";
            return false;
        }
        infile.open(args[option].as<string>());
        string line; int delim_pos;
        while(getline(infile, line)) {
            delim_pos = line.find('\t');
            key = stoi(line.substr(0, delim_pos));
            value = stoi(line.substr(delim_pos+1, line.length()-(delim_pos+1)));
            cutoff[key] = value;
        }
        infile.close();
    }

    return true;
}



bool parseArguments(int &argc, char* argv[], string &fasta_file, string &out_file, int &window_length,
                     int &window_bitcount_threshold, int &anchor_length, int &continuous_ones_threshold) {
    /*
     *  parsing input arguments for the program
     *  @param argc number of commandline arguments
     *  @param argv list of commandline arguments
     *  @param fasta_file stores the name of the fasta file
     *  @param out_file stores the name of the output file
     *  @param window_length stores the length of the window
     *  @param window_bitcount_threshold bitcount threshold in the window; default: 4
     *  @param anchor_length minimum length of continuous ones to be considered in the neighboring shift; default 3
     *  @param continuous_ones_threshold minimum number of continuous set bits in the shift XOR
     *  @return bool for successful completion of the function
    */
    po::options_description argparser("Below are the running options for the tool.");
    argparser.add_options()
        ("help,h", "Ribbit tool identifies short tandem repeats with allowed levels of inpurity.")

        ("input-file,i", po::value<string>(), "File path for the input fasta file.")
        ("output-file,o", po::value<string>(), "File path for the input fasta file.")        

        ("min-motif-length,m", po::value<int>(), "The minimum length of the motif of the repeats to be identified. Default: 2")
        ("max-motif-length,M", po::value<int>(), "The maximum length of the motif of the repeats to be identified, Default: 100")

        ("purity,p", po::value<float>(), "Threshold value for cotinuous number of ones found in a seed. Default: 0.85")

        ("min-length,l", po::value<string>(), "The minimum length of the repeat. Default: 12")
        ("min-units", po::value<string>(), "The minimum number of units of the repeat. Can be a integer value, for cutoff across all motif sizes.\
                                            Tab separated file with two columns, first is the motif size and second unit cutoff. Default: 2")
        ("perfect-units", po::value<string>(), "The minimum number of complete units of the repeat. Can be a integer value, for cutoff across all motif sizes.\
                                                Tab separated file with two columns, first is the motif size and second unit cutoff. Default: 2")

        /*
          currently all are set to default parameters and non-accessible to the user
          ("anchor", "Run the identification in anchor mode.")        // should be on by default; making non-accessible to the user
          ("window-length,w", po::value<int>(), "The length of window to be considered during seed identification. Default: 8")
          ("window-threshold,t", po::value<int>(), "The threshold value for number of 1s in the window. Defaut: 4")        
          ("anchor-length,a", po::value<int>(), "If running in anchor mode the length of the anchor to be considered. Default: 3")
          ("cones-threshold,c", po::value<int>(), "Threshold value for cotinuous number of ones found in a seed. Default: 0")        
        */
    ;

    po::variables_map args;
    po::store(po::parse_command_line(argc, argv, argparser), args);
    po::notify(args);

    if (args.count("help")) {
        cerr << argparser << "\n";
        return 0;
    }

    int default_perfect_units = 2;
    int default_minimum_length = 12;

    if (args.count("input-file"))  fasta_file = args["input-file"].as<string>();
    else {
        cerr << "ERROR: Please specify an input fasta file!\n";
        return 0;
    }

    if (args.count("output-file")) out_file = args["output-file"].as<string>();

    if (args.count("min-motif-length")) { MINIMUM_MLEN = args["min-motif-length"].as<int>(); }
    if (args.count("max-motif-length")) { MAXIMUM_MLEN = args["max-motif-length"].as<int>(); }

    /*
      currently all are set to default parameters and non-accessible to the user
      if (args.count("anchor")) run_mode = "anchor";
      if (args.count("window-length")) window_length = args["window-length"].as<int>();
      if (args.count("window-threshold")) window_bitcount_threshold = args["window-threshold"].as<int>();
      if (args.count("anchor-length")) anchor_length = args["anchor-length"].as<int>();
      if (args.count("cones-threshold")) continuous_ones_threshold = args["cones-threshold"].as<int>();
    */


    if (args.count("min-length")) {
        // either take minimum length as the input or minimum units
        parseDualtypeArgs(args, "min-length", MINIMUM_LENGTH, MINIMUM_MLEN, MAXIMUM_MLEN);
    }
    else if (args.count("min-units")) {
        LENGTH_CUTOFF_MODE = false;
        parseDualtypeArgs(args, "min-units", MINIMUM_UNITS, MINIMUM_MLEN, MAXIMUM_MLEN);
    }
    else {
        // uses minimum length of 12 as default if no input for minimum length or units are provided
        for (int key=MINIMUM_MLEN; key<=MAXIMUM_MLEN; key++) {            
            if (default_minimum_length < 2*key) {
                // if the minimum length is not atleast twice as the motif
                MINIMUM_LENGTH[key] = 2*key;
            }
            else MINIMUM_LENGTH[key] = default_minimum_length;
        }
    }

    // input for minimum number of perfect units. Is set to 2 by default
    if (args.count("perfect-units")) {
        parseDualtypeArgs(args, "perfect-units", PERFECT_UNITS, MINIMUM_MLEN, MAXIMUM_MLEN);
    } else {
        for (int m=1; m<=MAXIMUM_MLEN; m++) {
            switch (m) {
                case 1:  PERFECT_UNITS[m] = 8; break;
                case 2:  PERFECT_UNITS[m] = 4; break;
                case 3:  PERFECT_UNITS[m] = 3; break;
                default: PERFECT_UNITS[m] = 2; break;
            }
        }
    }

    return 1;
}


int main(int argc, char *argv[]) {
    /*
     *  main entry point to the ribbit programme
     *  @param argc number of commandline arguments
     *  @param argv list of commandline arguments
    */

    // exception for handling missing fasta index files handling gzip inputs
    string fasta_file = "", out_file = "";

    // defaults which are not be changed
    int window_length = 8, window_bitcount_threshold = 7, anchor_length = 3, cones_threshold = 3;

    bool success = parseArguments(argc, argv, fasta_file, out_file, window_length,
                                   window_bitcount_threshold, anchor_length, cones_threshold);
    if (!success) exit(1);

    // assigns the output to either a file or standard output
    streambuf * buf; ofstream outstream;
    if (out_file != "") {
        outstream.open(out_file);
        buf = outstream.rdbuf();    // output file buffer is created
    }
    // if output file is not provided the default is set at stdout
    else { buf = std::cerr.rdbuf(); }
    ostream out(buf);

    ifstream fastain(fasta_file);
    string line, seq_name, sequence="";

    if (!LENGTH_CUTOFF_MODE) {
        // if length cutoff is mentioned as units we convert that into bases
        for(auto pair: MINIMUM_UNITS) {
            MINIMUM_LENGTH[pair.first] = pair.first * pair.second;
        }
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

    cerr << "Minimum motif:\t" << MINIMUM_MLEN << "\n";
    cerr << "Maximum motif:\t" << MAXIMUM_MLEN << "\n";
    // minimum shift XOR to be generated; should be one less than the minimum motif size
    NMOTIFS = MAXIMUM_MLEN - MINIMUM_MLEN + 1;
    MINIMUM_SHIFT = (MINIMUM_MLEN > 2) ? MINIMUM_MLEN-2 : 1;
    MAXIMUM_SHIFT = MAXIMUM_MLEN + 2;
    NSHIFTS = MAXIMUM_SHIFT - MINIMUM_SHIFT + 1;

    cerr << "Purity threshold: " << PURITY_THRESHOLD << "\n";

    // Dynamically allocate memory for the matrix
    int SMALL_MLEN_LIMIT = 10;    // only save repeat classes for smaller motif sizes
    REPEAT_CLASSES = new uint32_t*[SMALL_MLEN_LIMIT];
    NUM_MOTIFS = pow(4, SMALL_MLEN_LIMIT);
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        REPEAT_CLASSES[i] = new uint32_t[NUM_MOTIFS];
    }
    MOTIF_FREQUENCY = new int[NUM_MOTIFS];
    MOTIF_UNITS = new int[NUM_MOTIFS];
    MOTIF_START = new int[NUM_MOTIFS];
    MOTIF_END = new int[NUM_MOTIFS];
    MOTIF_GAPS = new int[NUM_MOTIFS];
    MOTIF_GAPSIZE = new int[NUM_MOTIFS];
    MOTIF_NEXT = new uint32_t[NUM_MOTIFS];

    // Initialize the matrix (optional)
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        for (int j = 0; j < NUM_MOTIFS; ++j) {
            REPEAT_CLASSES[i][j] = NUM_MOTIFS;
        }
    }

    while (getline(fastain, line)) {
        if (line[0] == '>') {
            if (sequence != "") {
                cerr << "Processing sequence " << seq_name << "\n";
                processSequence(seq_name, sequence, window_length, window_bitcount_threshold, anchor_length, cones_threshold, out);
            }
            seq_name = line.substr(1, line.find(' ') - 1);
            sequence = "";
        }
        else { sequence += line; }
    }
    processSequence(seq_name, sequence, window_length, window_bitcount_threshold, anchor_length, cones_threshold, out);

    // Don't forget to free the memory when done
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        delete[] REPEAT_CLASSES[i];
    }
    delete[] REPEAT_CLASSES;
    delete[] MOTIF_FREQUENCY;
    delete[] MOTIF_UNITS;
    delete[] MOTIF_START;
    delete[] MOTIF_END;
    delete[] MOTIF_GAPS;
    delete[] MOTIF_GAPSIZE;
    delete[] MOTIF_NEXT;

    fastain.close(); outstream.close();
    return 0;
}
