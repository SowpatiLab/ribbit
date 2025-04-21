#include "consensus.h"


bool sortByVal(const pair<string, int> &a, const pair<string, int> &b) { 
    /*
     * sorts the variation count map
     * @param a first pair
     * @param b second pair
     * @return bool if the first pair is less than the second pair
    */
    return (a.second > b.second); 
} 


vector<pair<string, int>> sortVariationCount(unordered_map<string, int> var_count) {
    /*
     * sorts the variation count map
     * @param var_count the variation count map
     * @return unordered_map<string, int> the sorted variation count map
    */

    vector<pair<string, int>> vec;
    unordered_map<string, int> :: iterator it;
    for (it=var_count.begin(); it!=var_count.end(); it++) {
        vec.push_back(make_pair(it->first, it->second));
    }

    sort(vec.begin(), vec.end(), sortByVal);
    return vec;
}


string consensusMotif(string motif, string cigar) {
    /*
     * generates the consensus motif from the cigar string
     * @param motif the motif sequence
     * @param cigar the cigar string
     * @return string the consensus motif
    */

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int> clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);
    unordered_map<string, int> var_count;

    int motif_position = 0;
    int position = 0;
    string variation = "";
    cout << "CIGAR: " << cigar << "\n";

    for (int i=0; i<clens.size(); i++) {
        if (ctypes[i] == 'M' || ctypes[i] == '=') {
            position += clens[i];
            motif_position += clens[i];
        }
        else if (ctypes[i] == 'X') {
            variation = to_string(motif_position) + "X";
            if (var_count.find(variation) == var_count.end()) {
                var_count[variation] = 0;
            }
            var_count[variation] += 1;
            position += clens[i];
            motif_position += clens[i];
        }
        else if (ctypes[i] == 'I') {
            variation = to_string(motif_position) + "I";
            if (var_count.find(variation) == var_count.end()) {
                var_count[variation] = 0;
            }
            var_count[variation] += 1;
            position += clens[i];
        }
        else if (ctypes[i] == 'D') {
            variation = to_string(motif_position) + "D";
            if (var_count.find(variation) == var_count.end()) {
                var_count[variation] = 0;
            }
            var_count[variation] += 1;
            motif_position += clens[i];
        }

        if (motif_position >= motif.length()) {
            motif_position = motif_position % motif.length();
        }
    }

    vector<pair<string, int>> sorted_var_count = sortVariationCount(var_count);
    for (int i=0; i<sorted_var_count.size(); i++) {
        cout << sorted_var_count[i].first << "\t" << sorted_var_count[i].second << "\n";
    }

    return motif;
}