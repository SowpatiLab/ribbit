#include "concatenate_output.h"

using namespace std;


struct BedRecord {
    /*
     *  a structure to store the bed records
     *  @param start the start of the record
     *  @param end the end of the record
     *  @param record the record
     *  @returns void
     */
    int start;
    int end;
    string record;

    // Define a comparator for sorting
    bool operator<(const BedRecord& other) const {
        if (start != other.start) return start < other.start;
        return end < other.end;
    }
};


void sortBedRecords(std::vector<BedRecord>& records) {
    /*
     *  sorts the bed records based on the start position
     *  @param records vector of bed records
     *  @returns void
     */

    std::sort(records.begin(), records.end());
}


vector<tuple<string, int, int, string, double, string, int, int, int>> splitInfoField(string &info, int nloci,
                                                                                      string &sequence_id) {
    /*
     *  splits the info field of the output into individual repeat loci
     *  @param info the info field string
     *  @param sequence_id the sequence id of the repeats
     *  @returns vector of tuples of repeat loci
     */

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;
    vector<tuple<int, int, int, double>> coords;
    vector<string> motifs;
    vector<string> cigars;
    size_t pos = 0, dash_pos = 0;
    string token;
    int repeat_start, repeat_end, motif_length;
    double purity;

    string coords_str = info.substr(0, info.find(':'));
    info = info.substr(info.find(':') + 1);
    while ((pos = coords_str.find(',')) != string::npos) {
        token = coords_str.substr(0, pos); dash_pos = token.find('-');
        repeat_start = stoi(token.substr(0, dash_pos));
        token = token.substr(dash_pos + 1); dash_pos = token.find('-');
        repeat_end = stoi(token.substr(0, dash_pos));
        token = token.substr(dash_pos + 1); dash_pos = token.find('-');
        motif_length = stoi(token.substr(0, dash_pos));
        token = token.substr(dash_pos + 1);
        purity = stod(token);
        coords.push_back(make_tuple(repeat_start, repeat_end, motif_length, purity));
        coords_str.erase(0, pos + 1);
    }

    // last coordinate
    token = coords_str;
    token = coords_str.substr(0, pos); dash_pos = token.find('-');
    repeat_start = stoi(token.substr(0, dash_pos));
    token = token.substr(dash_pos + 1); dash_pos = token.find('-');
    repeat_end = stoi(token.substr(0, dash_pos));
    token = token.substr(dash_pos + 1); dash_pos = token.find('-');
    motif_length = stoi(token.substr(0, dash_pos));
    token = token.substr(dash_pos + 1);
    purity = stod(token);
    coords.push_back(make_tuple(repeat_start, repeat_end, motif_length, purity));
    coords_str.erase(0, pos + 1);

    string motifs_str = info.substr(0, info.find(':'));
    info = info.substr(info.find(':') + 1);
    while ((pos = motifs_str.find(',')) != string::npos) {
        token = motifs_str.substr(0, pos);
        motifs.push_back(token);
        motifs_str.erase(0, pos + 1);
    }
    motifs.push_back(motifs_str); // last motif

    string cigars_str = info;
    while ((pos = cigars_str.find(',')) != string::npos) {
        token = cigars_str.substr(0, pos);
        cigars.push_back(token);
        cigars_str.erase(0, pos + 1);
    }
    cigars.push_back(cigars_str); // last cigar

    for (int i=0; i<nloci; i++) {
        tie(repeat_start, repeat_end, motif_length, purity) = coords[i];
        string motif = motifs[i];
        string cigar_string = cigars[i];
        int repeat_length = repeat_end - repeat_start;
        int repeat_units = repeat_length / motif_length;
        repeat_loci.push_back(make_tuple(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                         motif_length, repeat_length, repeat_units));
    }

    return repeat_loci;
}


void concatenateThreadOutputs(const vector<string> &temp_files, ofstream &out) {
    /*
     *  concatenate the output files from different threads into a single output file
     *  @param temp_files list of temporary files from different threads
     *  @param out output stream to write the final output
    */

    string line;
    string sequence_id;
    int start, end;
    string motif, orientation;
    double purity;
    string info;
    string cigar_string;
    int motif_length, repeat_length, repeat_units;

    ostream* out_ptr = &out;

    vector<tuple<string, int, int, string, double, string, int, int, int>> repeat_loci;
    for (const auto& temp_file : temp_files) {
        ifstream ts(temp_file);
        if (!ts.is_open()) {
            cerr << "Could not open temporary file: " << temp_file << "\n";
            continue;
        }
        vector<string> fields;
        while (getline(ts, line)) {
            size_t start = 0, end;
            while ((end = line.find('\t', start)) != string::npos) {
                fields.push_back(line.substr(start, end - start));
                start = end + 1;
            }
            fields.push_back(line.substr(start)); // last field

            sequence_id = fields[0];
            start = stoi(fields[1]);
            end = stoi(fields[2]);
            motif = fields[3];
            purity = stod(fields[4]);
            motif_length = stoi(fields[5]);
            repeat_length = stoi(fields[6]);
            repeat_units = stoi(fields[7]);
            info = fields[8];

            if (info[0] == 'I') {
                cigar_string = info.substr(2);
                int recursion_level = 0;
                vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
                addLocusToOutput(sequence_id, start, end, motif, purity, cigar_string, motif_length, repeat_length,
                                 repeat_units, out_ptr, repeat_loci, recursion_level, new_repeat_loci);
            }
            else if (info[0] == 'M') {
                info = info.substr(2);
                cigar_string = info.substr(0, info.find(':'));
                info = info.substr(info.find(':') + 1);
                int nloci = stoi(info.substr(0, info.find(':')));
                info = info.substr(info.find(':') + 1);
                vector<tuple<string, int, int, string, double, string, int, int, int>> nested_loci = splitInfoField(info, nloci, sequence_id);
                // add the merged locus
                for (auto locus: nested_loci) {
                    int recursion_level = 0;
                    vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
                    tie(sequence_id, start, end, motif, purity, cigar_string, motif_length, repeat_length, repeat_units) = locus;
                    addLocusToOutput(sequence_id, start, end, motif, purity, cigar_string, motif_length, repeat_length,
                                     repeat_units, out_ptr, repeat_loci, recursion_level, new_repeat_loci);
                }
            }

            fields.clear();
        }
        ts.close();
    }

    if (repeat_loci.size() > 0) {
        for (int i=0; i<repeat_loci.size(); i++) {
            out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) << "\t" << get<2> (repeat_loci[i]) << "\t"
                << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<6> (repeat_loci[i]) << "\t" 
                << get<7> (repeat_loci[i]) << "\t" << get<8> (repeat_loci[i]);
            if (CIGAROUTPUT) { out << "\t" << get<5> (repeat_loci[i]); }
            out << "\n";
        }
    }

    for (const auto& temp_file : temp_files) {
        std::remove(temp_file.c_str());
    }

    out_ptr = nullptr;
}
