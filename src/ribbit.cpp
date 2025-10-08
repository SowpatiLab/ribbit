#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include "global_variables.h"
#include "fasta_utils.h"
#include "concatenate_output.h"
#include "binomial_thresholds.h"

using namespace std;


bool isNumber(const string &s) {
    /*
     *  checks if a string is number or not
     *  @param s string to be checked if it is numeric
     *  @return bool if the string is numeric or not
     */
    return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c)
                                      { return !std::isdigit(c); }) == s.end();
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
        for (key = minimum_motif_length; key <= maximum_motif_length; key++) {
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
        string line;
        int delim_pos;
        while (getline(infile, line)) {
            delim_pos = line.find('\t');
            key = stoi(line.substr(0, delim_pos));
            value = stoi(line.substr(delim_pos + 1, line.length() - (delim_pos + 1)));
            cutoff[key] = value;
        }
        infile.close();
    }

    return true;
}


bool parseArguments(int &argc, char *argv[], string &input_file, string &output_file) {
    /*
     *  parsing input arguments for the program
     *  @param argc number of commandline arguments
     *  @param argv list of commandline arguments
     *  @param input_file stores the name of the fasta file
     *  @param output_file stores the name of the output file
     *  @return bool for successful completion of the function
     */

    po::options_description argparser("Below are the running options for the tool.");
    argparser.add_options()("help,h", "Ribbit is designed to identify tandem repeats in DNA sequences with specific focus\
                                       on annotating complex TR loci.")

        ("input-file,i",  po::value<string>(), "File path for the input fasta file.")
        ("output-file,o", po::value<string>(), "File path for the input fasta file.\
                                                Default: adds a ribbit suffix to input file.")

        ("min-motif-length,m", po::value<int>(), "The minimum length of the motif of identified TR loci. Default: 2")
        ("max-motif-length,M", po::value<int>(), "The maximum length of the motif of identified TR loci. Default: 100")

        ("min-purity,p", po::value<double>(), "The minimum allowed purity of complete repeat. Default: 0.8")
        ("min-motif-purity,q", po::value<double>(), "Minimum match of each motif with consensus motif. Default: 0.8")

        ("min-length,l",  po::value<string>(), "The minimum length of the repeat. Default: 12")
        ("min-units",     po::value<string>(), "The minimum number of units of the repeat. Can be a integer value, for cutoff across all motif sizes.\
                                                Tab separated file with two columns, first is the motif size and second unit cutoff. Default: 2")
        ("perfect-units", po::value<string>(), "The minimum number of complete units of the repeat. Can be a integer value, for cutoff across all motif sizes.\
                                                Tab separated file with two columns, first is the motif size and second unit cutoff. Default: 2")

        ("cigar",     po::bool_switch()->default_value(false), "Include cigar string in the output. Default is off.")
        ("threads,t", po::value<int>(), "Number of threads to be used for running. default: 1");

    po::variables_map args;
    po::store(po::parse_command_line(argc, argv, argparser), args);
    po::notify(args);

    if (args.count("help")) {
        cerr << argparser << "\n";
        return 0;
    }

    int default_perfect_units = 2;
    int default_minimum_length = 12;

    if (args.count("input-file")) input_file = args["input-file"].as<string>();
    else {
        cerr << "ERROR: Please specify an input fasta file!\n";
        return 0;
    }

    if (args.count("output-file")) output_file = args["output-file"].as<string>();

    if (args.count("min-motif-length")) {
        MINIMUM_MLEN = args["min-motif-length"].as<int>();
    }
    if (args.count("max-motif-length")) {
        MAXIMUM_MLEN = args["max-motif-length"].as<int>();
    }
    if (args.count("purity")) {
        PURITY_THRESHOLD = args["purity"].as<double>();
    }
    if (args.count("motif-purity")) {
        MOTIFPURITY_THRESHOLD = args["motif-purity"].as<double>();
    }
    if (args.count("threads")) {
        THREADS = args["threads"].as<int>();
    }

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
        for (int key = MINIMUM_MLEN; key <= MAXIMUM_MLEN; key++)
        {
            if (default_minimum_length < 2 * key)
            {
                // if the minimum length is not atleast twice as the motif
                MINIMUM_LENGTH[key] = 2 * key;
            }
            else
                MINIMUM_LENGTH[key] = default_minimum_length;
        }
    }

    // input for minimum number of perfect units. Is set to 2 by default
    if (args.count("perfect-units")) {
        parseDualtypeArgs(args, "perfect-units", PERFECT_UNITS, MINIMUM_MLEN, MAXIMUM_MLEN);
    }
    else {
        for (int m = 1; m <= MAXIMUM_MLEN; m++) {
            switch (m) {
            case 1:
                PERFECT_UNITS[m] = 8; break;
            case 2:
                PERFECT_UNITS[m] = 4; break;
            case 3:
                PERFECT_UNITS[m] = 3; break;
            default:
                PERFECT_UNITS[m] = 2; break;
            }
        }
    }

    // cigar to be included in the output
    if (args.count("cigar")) CIGAROUTPUT = args["cigar"].as<bool>();

    return 1;
}


int main(int argc, char *argv[]) {
    /*
     *  main entry point to the ribbit programme
     *  @param argc number of commandline arguments
     *  @param argv list of commandline arguments
    */

    cerr << "\nRibbit: A fast and accurate tandem repeat finder.\n";
    cerr << "Version 1.0.0\n";

    // exception for handling missing fasta files handling gzip inputs
    string input_file = "", output_file = "";

    bool success = parseArguments(argc, argv, input_file, output_file);
    // Check if argument parsing was successful if not exit the program
    if (!success) exit(1);

    // check if the input fasta file exists
    fs::path input_path(input_file);
    if (!fs::exists(input_path)) {
        cerr << "ERROR: Input fasta file " << input_file << " does not exist!\n";
        exit(1);
    }

    if (!LENGTH_CUTOFF_MODE) {
        // if length cutoff is mentioned as units we convert that into bases
        for (auto pair : MINIMUM_UNITS) {
            MINIMUM_LENGTH[pair.first] = pair.first * pair.second;
        }
    }

    // for the motif sizes which are not in selected range but are factors of motif sizes
    // thresholds are set to the nearest selected motif size
    for (int m = MINIMUM_MLEN; m <= MAXIMUM_MLEN; m++) {
        vector<int> motif_factors;
        for (int f = 1; f <= m / 2; f++) {
            // calculate all the factors of the motif size
            if (m % f == 0) { motif_factors.push_back(f); }
        }
        for (int f : motif_factors) {
            if (MINIMUM_LENGTH.find(f) == MINIMUM_LENGTH.end()) {
                // if threshold for factor not set; set the threshold
                MINIMUM_LENGTH[f] = MINIMUM_LENGTH[m];
            }
            if (PERFECT_UNITS.find(f) == PERFECT_UNITS.end()) {
                // if threshold for factor not set; set the threshold
                PERFECT_UNITS[f] = PERFECT_UNITS[m] * (m / f);
            }
        }
    }

    cerr << "Minimum motif:\t" << MINIMUM_MLEN << "\n";
    cerr << "Maximum motif:\t" << MAXIMUM_MLEN << "\n";

    // minimum shift XOR to be generated; should be one less than the minimum motif size
    NMLENS = MAXIMUM_MLEN - MINIMUM_MLEN + 1;
    MINIMUM_SHIFT = (MINIMUM_MLEN > 2) ? MINIMUM_MLEN - 2 : 1;

    // maximum anchor jump based on maximum motif size
    // for motif sizes <= 6 the anchor jump is 1
    // for motif sizes > 6 the anchor jump is integer value of 20% of the motif size
    if (MAXIMUM_MLEN <= 6) { MAXIMUM_SHIFT = MAXIMUM_MLEN + 1; }
    else {
        MAXIMUM_SHIFT = ((2 * MAXIMUM_MLEN / 10) > 1) ? MAXIMUM_MLEN + (2 * MAXIMUM_MLEN / 10) : MAXIMUM_MLEN + 2;
    }

    NSHIFTS = MAXIMUM_SHIFT - MINIMUM_SHIFT + 1;

    cerr << "Purity threshold: " << PURITY_THRESHOLD << "\n";
    cerr << "Motif purity threshold: " << MOTIFPURITY_THRESHOLD << "\n\n";

    cerr << "NOTE: Purity threshold is a strict cutoff for VNTRs (motif length >= 7bp) but for STRs (motif length <= 6bp)\n\
             the cutoff is flexible as length of the perfect stretches is also considered for STRs.\n\n";

    // Dynamically allocate memory for the matrix
    REPEAT_CLASSES = new uint32_t *[SMALL_MLEN_LIMIT];
    NUM_MOTIFS = pow(4, SMALL_MLEN_LIMIT);
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        REPEAT_CLASSES[i] = new uint32_t[NUM_MOTIFS];
    }
    MOTIF_START = new int[NUM_MOTIFS];
    MOTIF_END = new int[NUM_MOTIFS];
    MOTIF_NEXT = new uint32_t[NUM_MOTIFS];
    MOTIF_UNITS = new int[NUM_MOTIFS];
    MOTIF_GAPS = new int[NUM_MOTIFS];
    MOTIF_GAPSIZE = new int[NUM_MOTIFS];

    WINDOW_LENGTHS    = new int[NMLENS];
    WINDOW_THRESHOLDS = new int[NMLENS];
    calculateWindowThresholds();    // calculates the window lengths and thresholds for all motif sizes 

    SEEDLEN_CUTOFF = new int[NMLENS];
    for (int i = 0; i < NMLENS; ++i) {
        SEEDLEN_CUTOFF[i] = ((i + MINIMUM_MLEN) > SMALL_MLEN_LIMIT) ? 0.9 * (i + MINIMUM_MLEN) : 
                                                                      (MINIMUM_LENGTH[i + MINIMUM_MLEN] - (i + MINIMUM_MLEN));
    }

    // Initialize the matrix (optional)
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i) {
        for (int j = 0; j < NUM_MOTIFS; ++j) {
            REPEAT_CLASSES[i][j] = NUM_MOTIFS;
        }
    }

    parseFasta(input_file, output_file);

    // Don't forget to free the memory when done
    for (int i = 0; i < SMALL_MLEN_LIMIT; ++i)  delete[] REPEAT_CLASSES[i];

    delete[] REPEAT_CLASSES;
    delete[] MOTIF_START;
    delete[] MOTIF_END;
    delete[] MOTIF_NEXT;
    delete[] MOTIF_UNITS;
    delete[] MOTIF_GAPS;
    delete[] MOTIF_GAPSIZE;

    double seconds_since_start = difftime(time(0), START_TIME);
    std::cerr << "Total time elapsed: " << seconds_since_start << "secs\n";

    return 0;
}
