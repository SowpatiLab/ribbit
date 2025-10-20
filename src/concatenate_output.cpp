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


void concatenateOutputs(string out_file, vector<string>seq_names, int THREADS) {
    /*
     *  sorts and concatenates all the outputs from different threads
     *  @param out_file the name of the output file
     *  @param seq_names vector of the sequence names
     *  @param THREADS number of threads used by the program
     *  @returns void
     */

    ofstream out(out_file);
    vector<int> starts;
    string line; vector<string> lines;
    string chrom; int start, end;

    std::vector<BedRecord> records;
    for (string seq_name: seq_names) {
        for (int tnum = 1; tnum <= THREADS; tnum++) {
            records.clear();
            string output_name = out_file + "_" + seq_name + "_" + to_string(tnum);
            ifstream chunk(output_name);
            cerr << "Concatenating " << output_name << "\n";
            while (getline(chunk, line)) {
                chrom = line.substr(0, line.find('\t')); line = line.substr(line.find('\t')+1, line.length()-(line.find('\t') +1));
                start = stoi(line.substr(0, line.find('\t'))); line = line.substr(line.find('\t')+1, line.length()-(line.find('\t') +1));
                end = stoi(line.substr(0, line.find('\t'))); line = line.substr(line.find('\t')+1, line.length()-(line.find('\t') +1));
                records.push_back({start, end, line});
            }
            sortBedRecords(records);

            for (BedRecord record: records) {
                out << seq_name << "\t" << record.start << "\t" << record.end << "\t" << record.record << "\n";
            }
            chunk.close();
            remove(output_name.c_str());
        }
    }
    out.close();
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
    string cigar_string;
    int motif_length, repeat_length, repeat_units;

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
            orientation = fields[5];
            motif_length = stoi(fields[6]);
            repeat_length = stoi(fields[7]);
            repeat_units = stoi(fields[8]);
            cigar_string = fields[9];

            int recursion_level = 0;
            vector<tuple<string, int, int, string, double, string, int, int, int>> new_repeat_loci;
            addLocusToOutput(sequence_id, start, end, motif, purity, cigar_string,
                             motif_length, repeat_length, repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);

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
}
