#include "concatenate_output.h"

using namespace std;


struct BedRecord {
    /*
     * A structure to store the bed records
     * @param start the start of the record
     * @param end the end of the record
     * @param record the record
     * @returns void
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
     * sorts the bed records based on the start position
     * @param records vector of bed records
     * @returns void
    */
    std::sort(records.begin(), records.end());
}


void concatenateOutputs(string out_file, vector<string>seq_names, int THREADS) {
    /*
     * sorts and concatenates all the outputs from different threads
     * @param out_file the name of the output file
     * @param seq_names vector of the sequence names
     * @param THREADS number of threads used by the program
     * @returns void
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
