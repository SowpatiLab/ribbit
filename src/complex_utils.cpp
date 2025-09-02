#include <boost/program_options.hpp>
#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <vector>
#include <tuple>

namespace po = boost::program_options;

using namespace std;


bool isNumber(const string &s) {
    /*
     *  checks if a string is number or not
     *  @param s string to be checked if it is numeric
     *  @return bool if the string is numeric or not
    */ 
    return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

string buildCigar(vector<int> &clens, vector<char> &ctypes) {
    /*
     * builds the cigar string from the lengths and types of cigar operations
     * @param clens vector of lengths of cigar operations
     * @param ctypes vector of types of cigar operations
     * @return string the cigar string
    */
    string cigar = "";
    for (int _=0; _<clens.size(); _++) {
        cigar += to_string(clens[_]); cigar += ctypes[_];
    }
    return cigar;
}


string buildCigar(tuple<vector<int>, vector<char>> &cigar_values) {
    /*
     * builds the cigar string from the lengths and types of cigar operations
     * @param cigar_values tuple of lengths and types of cigar operations
     * @return string the cigar string
    */
    vector<int> clens = get<0> (cigar_values);
    vector<char> ctypes = get<1> (cigar_values);
    string cigar = "";
    for (int _=0; _<clens.size(); _++) {
        cigar += to_string(clens[_]); cigar += ctypes[_];
    }
    return cigar;
}



tuple<vector<int>, vector<char>> cigarSplit(string cigar) {
    /*
     *  splits a CIGAR string into consecutive operations and lengths
     *  @param cigar character pointer of the CIGAR string
     *  @return tuple<vector<int>, vector<char>> tuple of two vectors cigar lengths and cigar types
    */

    string length = "";
    vector<int> clens; vector<char> ctypes;

    for (int i = 0; i<cigar.length(); i++) {
        if (isdigit(cigar[i])) {
            length += cigar[i];
        }
        else {
            clens.push_back(stoi(length));
            if (cigar[i] == '=') { ctypes.push_back('M'); }
            else { ctypes.push_back(cigar[i]); }
            length = "";
        }
    }

    return tuple<vector<int>, vector<char>> {clens, ctypes};
}


tuple<vector<int>, vector<char>> extractRegionCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end) {
    /*
     * extracts the cigar between two coordinates within a locus
     *  @param cigar tuple of lengths nd types of cigar operations
     *  @param start the start of the region for which cigar should be pulled
     *  @param end the end of the region for which cigar should be pulled
     *  @return the cigar of the the region to be extracted as tuple<vector<int>, vector<char>>
    */

    vector<int> clens = get<0> (cigar);
    vector<char> ctypes = get<1> (cigar);
    vector<int> new_clens; vector<char> new_ctypes;
    new_clens.clear(); new_ctypes.clear();
    if (start == end) {
        new_clens.push_back(0); new_ctypes.push_back('M');
        return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
    }

    int rpos = 0, i = 0;
    int clen; char ctype;
    bool region = false;
    for (i=0; i<clens.size(); i++) {
        clen = clens[i]; ctype = ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }

        if (!region) {
            if (rpos == start) { region = true; }
            else if (rpos > start) {
                region = true;
                new_clens.push_back(rpos - start);
                new_ctypes.push_back(ctype);
            }
            if (rpos == end) { return tuple<vector<int>, vector<char>> {new_clens, new_ctypes}; }
            else if (rpos > end) {
                new_clens[new_clens.size()-1] = new_clens[new_clens.size()-1] - (rpos - end);
                return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
            }
            continue;
        }

        if (region) {
            if (new_clens.size() == 0 && ctype == 'D') {}
            else {
                new_clens.push_back(clen);
                new_ctypes.push_back(ctype);
                if (rpos == end) break;
                else if (rpos > end) {
                    new_clens[new_clens.size()-1] = new_clens[new_clens.size()-1] - (rpos - end);
                    return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
                }
            }
        }
    }

    return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
}


void getMatches(string &cigar, int &matches, int &longest_match) {
    /*
     * calculates the number of matches from the cigar string
     *  @param cigar the cigar string
     *  @param matches the number of matches
     *  @param longest_match the length of the longest match
    */

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int> clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);
    matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
        }
    }
}


void getMatches(tuple<vector<int>, vector<char>> &cigar, int &matches, int &longest_match) {
    /*
     * calculates the number of matches from the lengths and types of cigar operations
     *  @param cigar tuple of lengths and types of cigar operations
     *  @param matches the number of matches
     *  @param longest_match the length of the longest match
    */
    vector<int> clens = get<0> (cigar);
    vector<char> ctypes = get<1> (cigar);
    matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
            if (clens[_] > longest_match) longest_match = clens[_];
        }
    }
}


void getMatches(vector<int> &clens, vector<char> &ctypes, int &matches, int &longest_match) {
    /*
     * calculates the number of matches from the lengths and types of cigar operations
     * @param clens vector of lengths of cigar operations
     * @param ctypes vector of types of cigar operations
     *  @param matches the number of matches
     *  @param longest_match the length of the longest match
    */
    matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
            if (clens[_] > longest_match) longest_match = clens[_];
        }
    }
}


int getRepeatLength(vector<int> &clens, vector<char> &ctypes) {
    /*
     * calculates the alignment length from the lengths and types of cigar operations
     * @param clens vector of lengths of cigar operations
     * @param ctypes vector of types of cigar operations
     * @return int the alignment length
    */
    int alignment_length = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M' || ctypes[_] == 'X' || ctypes[_] == 'I') {
            alignment_length += clens[_];
        }
    }
    return alignment_length;
}


int getRepeatLength(tuple<vector<int>, vector<char>> &cigar_values) {
    /*
     * calculates the alignment length from the lengths and types of cigar operations
     * @param cigar_values tuple of lengths and types of cigar operations
     * @return int the alignment length
    */
    vector<int> clens = get<0> (cigar_values);
    vector<char> ctypes = get<1> (cigar_values);
    return getRepeatLength(clens, ctypes);
}


tuple<vector<int>, vector<char>> extractDownCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end) {
    /*
     * extracts the downstream cigar from the start position
     * @param cigar tuple of lengths and types of cigar operations
     * @param start start position of the downstream locus
     * @param end end position of the downstream locus
    */
    vector<int> clens = get<0> (cigar);
    vector<char> ctypes = get<1> (cigar);
    vector<int> new_clens; vector<char> new_ctypes;
    int rpos = start, i = 0;
    int clen; char ctype;
    for (i=0; i<clens.size(); i++) {
        clen = clens[i]; ctype = ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }
        new_clens.push_back(clen);
        new_ctypes.push_back(ctype);

        if (rpos >= end) break;
    }
    if (rpos > end) new_clens[i] -= (rpos - end);
    return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
}


tuple<vector<int>, vector<char>> extractUpCigar(tuple<vector<int>, vector<char>> &cigar, int start, int end) {
    /*
     * extracts the upstream cigar from the end position
     * @param cigar tuple of lengths and types of cigar operations
     * @param start start position of the upstream locus
     * @param end end position of the upstream locus
    */

    vector<int> clens = get<0> (cigar);
    vector<char> ctypes = get<1> (cigar);
    vector<int> new_clens; vector<char> new_ctypes;
    int rpos = start, i = 0;
    int clen; char ctype;
    for (i=clens.size()-1; i>=0; i--) {
        clen = clens[i]; ctype = ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos -= clen;
        }

        if (rpos <= end) break;
    }

    new_clens  = vector<int>(clens.begin() + i, clens.end());
    new_ctypes = vector<char>(ctypes.begin() + i, ctypes.end());
    if (rpos < end) new_clens[0] -= (end - rpos);
    return tuple<vector<int>, vector<char>> {new_clens, new_ctypes};
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


bool parseArguments(int &argc, char* argv[], string &input_file, string &out_file, int &minimum_mlen, int &maximum_mlen,
                    double &purity_threshold, unordered_map<int, int> &minimum_length, unordered_map<int, int> &minimum_units,
                    int &distance, bool &variant_mode) {
    /*
     *  parsing input arguments for the program
     *  @param argc number of commandline arguments
     *  @param argv list of commandline arguments
     *  @param input_file stores the name of the fasta file
     *  @param out_file stores the name of the output file
     *  @param window_length stores the length of the window
     *  @param window_bitcount_threshold bitcount threshold in the window; default: 4
     *  @param anchor_length minimum length of continuous ones to be considered in the neighboring shift; default 3
     *  @param continuous_ones_threshold minimum number of continuous set bits in the shift XOR
     *  @return bool for successful completion of the function
    */
    po::options_description argparser("Below are the running options for the tool.");
    argparser.add_options()
        ("help,h", "This script processes the output of ribbit to give complex tandem repeat loci.\n\
                    Usage: ribbit_complex -i <input_file> -o <output_file> [options]\n\
                    Note: The ribbit output file should include the cigar column. Make sure to run ribbit with --cigar option.\n\
                    Note: The accurate minimum and maximum motif length of repeats needs to be provided otherwise the script fails.")

        ("input-file,i", po::value<string>(), "The ribbit output file path. Please make sure the file is sorted by start position of the repeat loci.")
        ("output-file,o", po::value<string>(), "File path for the input fasta file. Default: adds a complex suffix to input file.")

        ("min-motif-length,m", po::value<int>(), "The minimum length of the motif of identified TR loci. Default: 2")
        ("max-motif-length,M", po::value<int>(), "The maximum length of the motif of identified TR loci. Default: 100")
        ("purity,p", po::value<double>(), "Filters out any repeats below this purity. Default: 0.8")

        ("min-length,l", po::value<string>(), "Filters out any repeats below this length. Default: 12")
        ("min-units", po::value<string>(), "The minimum number of units of the repeat. Can be a integer value, for cutoff across all motif sizes.\
                                            Tab separated file with two columns, first is the motif size and second unit cutoff. Default: 2")
        
        ("distance,d", po::value<int>(), "The maximum distance between neighboring repeats that need to be merged into a complex repeat. Default merges \
                                          bookend or overlapping repeats. Default: 0")
        ("variant-mode,v", po::bool_switch()->default_value(false), "Include cigar string in the output. Default is off." )

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

    if (args.count("input-file"))  input_file = args["input-file"].as<string>();
    else {
        cerr << "ERROR: Please specify an input file!\n";
        return 0;
    }

    if (args.count("output-file")) out_file = args["output-file"].as<string>();

    if (args.count("min-motif-length")) { minimum_mlen = args["min-motif-length"].as<int>(); }
    if (args.count("max-motif-length")) { maximum_mlen = args["max-motif-length"].as<int>(); }
    if (args.count("purity")) { purity_threshold = args["purity"].as<double>(); }

    if (args.count("min-length")) {
        // either take minimum length as the input or minimum units
        parseDualtypeArgs(args, "min-length", minimum_length, minimum_mlen, maximum_mlen);
    }
    else if (args.count("min-units")) {
        parseDualtypeArgs(args, "min-units", minimum_units, minimum_mlen, maximum_mlen);
    }
    else {
        // uses minimum length of 12 as default if no input for minimum length or units are provided
        for (int key=minimum_mlen; key<=maximum_mlen; key++) {            
            if (default_minimum_length < 2*key) {
                // if the minimum length is not atleast twice as the motif
                minimum_length[key] = 2*key;
            }
            else minimum_length[key] = default_minimum_length;
        }
    }

    if (args.count("distance")) { distance = args["distance"].as<int>(); }

    if (args.count("variant-mode")) variant_mode = args["variant-mode"].as<bool>();

    return 1;
}


void splitLine(const string &line, string &chrom, int &start, int &stop, string &motif, double &purity,
                string &strand, string &cigar, int &motif_length, int &repeat_length, int &repeat_units) {
    /*
     *  splits a line from the input file into its components
     *  @param line the line to be split
     *  @param chrom chromosome name
     *  @param start start position of the repeat
     *  @param stop end position of the repeat
     *  @param motif the sequence of the repeat motif
     *  @param purity purity value of the repeat
     *  @param strand strand information of the repeat
     *  @param cigar CIGAR string for the repeat alignment
     *  @param motif_length length of the motif
     *  @param repeat_length total length of the repeat locus
     *  @param repeat_units number of units of the motif in the repeat
    */
    istringstream iss(line);
    string token;
    getline(iss, chrom, '\t');
    iss >> start >> stop >> motif >> purity >> strand >> motif_length >> repeat_length >> repeat_units >> cigar;
}



int getBoundaryOverlappingLoci(int upstart, int upend, string &upmotif, tuple<vector<int>, vector<char>> &up_cigarvalues, string &upcigar,
                               int dnstart, int dnend, string &dnmotif, tuple<vector<int>, vector<char>> &dn_cigarvalues, string &dncigar) {
    /*
     * gets the boundary of the overlapping loci
     *  @param upstart start position of the upstream locus
     *  @param upend end position of the upstream locus
     *  @param upmotif motif of the upstream locus
     *  @param up_cigarvalues the lengths and types of cigar operations of the upstream locus
     *  @param dnstart start position of the downstream locus
     *  @param dnend end position of the downstream locus
     *  @param dnmotif motif of the downstream locus
     *  @param dn_cigarvalues the lengths and types of cigar operations of the downstream locus
    */

    int overlap_length = upend - dnstart;
    int up_matches = 0, dn_matches = 0, longest_match = 0, shortmotif_length = 0;
    int max_matches = 0, max_longest_match = 0;
    int max_shortmotif_length = 0;
    tuple<vector<int>, vector<char>> up_olseg_cigarvalues;
    tuple<vector<int>, vector<char>> dn_olseg_cigarvalues;
    tuple<vector<int>, vector<char>> up_new_cigarvalues;
    tuple<vector<int>, vector<char>> dn_new_cigarvalues;

    int boundary = 0, uppos = 0, dnpos = 0;
    string new_upcigar = "", new_dncigar = "";
    for (int i=0; i <= overlap_length; i++) {
        up_matches = 0, dn_matches = 0;

        up_olseg_cigarvalues = extractRegionCigar(up_cigarvalues, dnstart-upstart, (dnstart+i)-upstart);
        dn_olseg_cigarvalues = extractRegionCigar(dn_cigarvalues, i, upend-dnstart);
            
        getMatches(up_olseg_cigarvalues, up_matches, longest_match);
        getMatches(dn_olseg_cigarvalues, dn_matches, longest_match);
            
        if (up_matches + dn_matches > max_matches) {
            max_matches = up_matches + dn_matches;
            max_longest_match = longest_match;
            max_shortmotif_length = (upmotif.length() <= dnmotif.length()) ? getRepeatLength(up_olseg_cigarvalues) : getRepeatLength(dn_olseg_cigarvalues);
            boundary = i;
            up_new_cigarvalues = extractDownCigar(up_cigarvalues, upstart, dnstart + i);
            dn_new_cigarvalues = extractUpCigar(dn_cigarvalues, dnend, dnstart + i);
            new_upcigar = buildCigar(up_new_cigarvalues);
            new_dncigar = buildCigar(dn_new_cigarvalues);
        }

        else if (up_matches + dn_matches == max_matches) {
            if (longest_match > max_longest_match) {
                max_longest_match = longest_match;
                max_shortmotif_length = (upmotif.length() <= dnmotif.length()) ? getRepeatLength(up_olseg_cigarvalues) : getRepeatLength(dn_olseg_cigarvalues);
                boundary = i;
                up_new_cigarvalues = extractDownCigar(up_cigarvalues, upstart, dnstart + i);
                dn_new_cigarvalues = extractUpCigar(dn_cigarvalues, dnend, dnstart + i);
                new_upcigar = buildCigar(up_new_cigarvalues);
                new_dncigar = buildCigar(dn_new_cigarvalues);
            }
            else if (longest_match == max_longest_match) {
                if (upmotif.length() <= dnmotif.length())  shortmotif_length = getRepeatLength(up_olseg_cigarvalues);
                else  shortmotif_length = getRepeatLength(dn_olseg_cigarvalues);
                
                if (shortmotif_length > max_shortmotif_length) {
                    max_shortmotif_length = shortmotif_length;
                    max_longest_match = longest_match;
                    boundary = i;
                    up_new_cigarvalues = extractDownCigar(up_cigarvalues, upstart, dnstart + i);
                    dn_new_cigarvalues = extractUpCigar(dn_cigarvalues, dnend, dnstart + i);
                    new_upcigar = buildCigar(up_new_cigarvalues);
                    new_dncigar = buildCigar(dn_new_cigarvalues);
                }
            }
        }
    }

    upcigar = new_upcigar;
    dncigar = new_dncigar;
    return boundary + dnstart;
}


void processOverlappingRepeats(vector<tuple<string, int, int, string, double, int, int, int, string>> &overlapping_repeats, int distance) {
    /*
     *  processes the overlapping repeats and merges them based on the distance
     *  @param overlapping_repeats vector of tuples containing repeat information
     *  @param distance maximum distance for merging repeats
     *  @return void
    */
    // This function would typically process the overlapping repeats and merge them based on the distance.
    // For now, it is just a placeholder.
    cout << "Processing overlapping repeats with distance: " << distance << "\n";

    string chrom, motif, cigar;
    int start, stop, motif_length, repeat_length, repeat_units;
    double purity;
    tuple<vector<int>, vector<char>> cigarvalues;

    string other_chrom, other_motif, other_cigar;
    int other_start, other_stop, other_motif_length, other_repeat_length, other_repeat_units;
    double other_purity;
    tuple<vector<int>, vector<char>> other_cigarvalues;

    int ridx = 0, ridy = 0;   
    for (const auto &repeat : overlapping_repeats) {
        // compare each repeat with all the ones it is overlapping with
        chrom = get<0>(repeat);
        start = get<1>(repeat);
        stop = get<2>(repeat);
        motif = get<3>(repeat);
        purity = get<4>(repeat);
        motif_length = get<5>(repeat);
        repeat_length = get<6>(repeat);
        repeat_units = get<7>(repeat);
        cigar = get<8>(repeat);
        

        ridy = 0;
        for (const auto &other_repeat : overlapping_repeats) {
            if (ridy <= ridx) { ridy++; continue; }  // skip self-comparison
            
            other_chrom  = get<0>(other_repeat);
            other_start  = get<1>(other_repeat);
            other_stop   = get<2>(other_repeat);
            other_motif  = get<3>(other_repeat);
            other_purity = get<4>(other_repeat);
            other_motif_length  = get<5>(other_repeat);
            other_repeat_length = get<6>(other_repeat);
            other_repeat_units  = get<7>(other_repeat);
            other_cigar = get<8>(other_repeat);

            if (other_start < stop) {

                if (other_stop <= stop) {
                    // other repeat is completely within the current repeat
                    bool update = false;
                    // trimNestedRepeat(other_chrom, other_start, other_stop, other_motif, other_purity,
                    //                 other_motif_length, other_repeat_length, other_repeat_units, other_cigar, update);
                    overlapping_repeats[ridy] = make_tuple(other_chrom, other_start, other_stop, other_motif, other_purity,
                                                           other_motif_length, other_repeat_length, other_repeat_units, other_cigar);
                }

                else if (other_start == start && other_stop == stop) {
                    // both positions are identical
                }

                else if (other_start == start) {
                    // current repeat is nested within the other repeat

                    bool update = false;
                    // trimNestedRepeat(chrom, start, stop, motif, purity, motif_length, repeat_length, repeat_units, cigar, update);
                    if (update) {
                        overlapping_repeats[ridx] = make_tuple(chrom, start, stop, motif, purity, motif_length, repeat_length, repeat_units, cigar);
                    }
                }

                else {
                    // both are overlapping but not nested
                    cigarvalues = cigarSplit(cigar);
                    other_cigarvalues = cigarSplit(other_cigar);
                    int boundary = getBoundaryOverlappingLoci(start, stop, motif, cigarvalues, cigar,
                                                              other_start, other_stop, other_motif, other_cigarvalues, other_cigar);
                    cout << chrom << "\t" << start << "\t" << stop << "\t" << motif << "\t" 
                         << purity << "\t" << motif_length << "\t" << repeat_length << "\t" 
                         << repeat_units << "\t" << cigar << "\n";
                    cout << other_chrom << "\t" << other_start << "\t" << other_stop << "\t"
                         << other_motif << "\t" << other_purity << "\t" << other_motif_length << "\t" 
                         << other_repeat_length << "\t" << other_repeat_units << "\t" << other_cigar << "\n";
                    cout << "Boundary: " << boundary << "\n\n";

                }

            }
            else break;
        }

        overlapping_repeats[ridx] = make_tuple(chrom, start, stop, motif, purity, motif_length, repeat_length, repeat_units, cigar);
        ridx++;
    }

}


void cleanCigar(string &cigar) {
    /*
     *  cleans the cigar string by removing any unnecessary operations
     *  @param cigar the cigar string to be cleaned
     *  @return string the cleaned cigar string
    */
    // This function would typically clean the cigar string by removing unnecessary operations.
    // For now, it is just a placeholder.
    tuple<vector<int>, vector<char>> cigar_values = cigarSplit(cigar);
    cigar = "";
    char ctype = 'N'; int clen = 0;
    char pctype = 'N'; int pclen = 0;
    for (int i = 0; i < get<0>(cigar_values).size(); i++) {
        ctype = get<1>(cigar_values)[i];
        clen = get<0>(cigar_values)[i];
        if (pctype != ctype) {
            if (pctype != 'N') {
                // if the previous type is not empty, add it to the cigar string
                cigar += to_string(pclen) + pctype;
            }
            pctype = ctype;  // update the previous type
            pclen = clen;  // update the previous length
        }
        else {
            pclen += clen;  // if the previous type is the same, add the length to the previous length
        }
    }
    if (pctype != 'N') {
        // if the last type is not empty, add it to the cigar string
        cigar += to_string(pclen) + pctype;
    }
}


void processFile(string infile, int distance) {
    /*
     *  processes the input file and generates the output
     *  @param infile input file path
     *  @return void
    */
    // This function would typically read the input file, process it, and write to the output file.
    // For now, it is just a placeholder.
    cout << "Processing file: " << infile << "\n";

    string line;
    ifstream file(infile);
    if (!file.is_open()) {
        cerr << "ERROR: Could not open file " << infile << "\n";
        return;
    }

    string chrom, motif, strand, cigar;
    int start, stop, motif_length, repeat_length, repeat_units;
    double purity;
    vector<tuple<string, int, int, string, double, int, int, int, string>> overlapping_repeats;
    int region_end = -1;
    while (getline(file, line)) {
        // Process each line of the file
        // This is where you would implement the logic to process the input file

        if (line[0] == '#') continue;  // Skip comment lines
        splitLine(line, chrom, start, stop, motif, purity, strand, cigar, motif_length, repeat_length, repeat_units);
        cleanCigar(cigar);

        if (region_end == -1 || start - distance <= region_end) {
            // If this is the first repeat or a new region, reset the region end
            overlapping_repeats.push_back(make_tuple(chrom, start, stop, motif, purity, motif_length, repeat_length, repeat_units, cigar));
            if (region_end < stop) region_end = stop;  // Update the end of the current region
        }
        else if (start - distance > region_end) {
            // If the repeat overlaps with the current region, update the end
            if (overlapping_repeats.size() > 1) {
                // Process the overlapping repeats before adding the new one
                for (auto &repeat : overlapping_repeats) {
                    chrom = get<0>(repeat);
                    start = get<1>(repeat);
                    stop = get<2>(repeat);
                    motif = get<3>(repeat);
                    purity = get<4>(repeat);
                    motif_length = get<5>(repeat);
                    repeat_length = get<6>(repeat);
                    repeat_units = get<7>(repeat);
                    cigar = get<8>(repeat);

                    cout << chrom << "\t" << start << "\t" << stop << "\t" << motif << "\t"
                         << purity << "\t" << motif_length << "\t" << repeat_length << "\t"
                         << repeat_units << "\t" << cigar << "\n";
                }
                processOverlappingRepeats(overlapping_repeats, distance);
            }
            overlapping_repeats.clear();
            overlapping_repeats.push_back(make_tuple(chrom, start, stop, motif, purity, motif_length, repeat_length, repeat_units, cigar));
            region_end = stop;  // Reset the region end
        }

    }
    file.close();
}

int main(int argc, char *argv[]) {

    string input_file, out_file;
    int minimum_mlen = 2, maximum_mlen = 100;
    double purity_threshold = 0.8;
    unordered_map<int, int> minimum_length, minimum_units;
    int distance = 0;
    bool variant_mode = false;

    bool success = parseArguments(argc, argv, input_file, out_file, minimum_mlen, maximum_mlen,
                                  purity_threshold, minimum_length, minimum_units, distance, variant_mode);
    if (!success) exit(1);

    cout << "Input file: " << input_file << "\n";
    cout << "Output file: " << out_file << "\n";
    cout << "Minimum motif length: " << minimum_mlen << "\n";
    cout << "Maximum motif length: " << maximum_mlen << "\n";
    cout << "Purity threshold: " << purity_threshold << "\n";
    cout << "Minimum length: \n";
    cout << "Distance for merging repeats: " << distance << "\n";

    processFile(input_file, distance);
    // Here you would typically call the function to process the input file and generate the output
    
}