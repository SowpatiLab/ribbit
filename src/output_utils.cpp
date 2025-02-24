#include <iostream>
#include <fstream>
#include <unordered_map>
#include <mutex>
#include <boost/dynamic_bitset.hpp>

#include <cstdint>
#include <numeric>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <boost/multiprecision/cpp_int.hpp>

#include "global_variables.h"
#include "output_utils.h"

using namespace std;
using namespace boost;


int levenshteinDistance(string word1, string word2) {
    /*
     * calculates the levenshtein distance betwee two words
     * @param word1 string of word 1
     * @param word2 string of word 2
     * @return int levenshtein distance between the two words
    */
    int size1 = word1.size();
    int size2 = word2.size();
    int matrix[size1 + 1][size2 + 1]; // Verification matrix i.e. 2D array which will store the calculated distance.

    // If one of the words has zero length, the distance is equal to the size of the other word.
    if (size1 == 0) return size2;
    if (size2 == 0) return size1;

    // Sets the first row and the first column of the verification matrix with the numerical order from 0 to the length of each word.
    for (int i = 0; i <= size1; i++)
        matrix[i][0] = i;
    for (int j = 0; j <= size2; j++)
        matrix[0][j] = j;

    // Verification step / matrix filling.
    for (int i = 1; i <= size1; i++) {
        for (int j = 1; j <= size2; j++) {
            // Sets the modification cost.
            // 0 means no modification (i.e. equal letters) and 1 means that a modification is needed (i.e. unequal letters).
            int cost = (word2[j - 1] == word1[i - 1]) ? 0 : 1;

            // Sets the current position of the matrix as the minimum value between a (deletion), b (insertion) and c (substitution).
            // a = the upper adjacent value plus 1: matrix[i - 1][j] + 1
            // b = the left adjacent value plus 1: matrix[i][j - 1] + 1
            // c = the upper left adjacent value plus the modification cost: matrix[i - 1][j - 1] + cost
            matrix[i][j] = min(
                min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1),
                matrix[i - 1][j - 1] + cost
            );
        }
    }

    // The last position of the matrix will contain the Levenshtein distance.
    return matrix[size1][size2];
}


void printRepeatsToOutput(ostream &out, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                          int end_index) {
    /*
     * print repeats from the recorded repeat loci and clears the repeat loci
     * @param out outstream of the output file
     * @param repeat_loci vector of recorded repeat loci sorted based on the start
     * @param end_index index to which point repeats should be printed to the output
    */
    for (int i=0; i<=end_index; i++) {
        out << get<0> (repeat_loci[i]) << "\t" << get<1> (repeat_loci[i]) << "\t" << get<2> (repeat_loci[i]) << "\t"
            << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<5> (repeat_loci[i]) << "\t"
            << get<6> (repeat_loci[i]) << "\t" << get<7> (repeat_loci[i]) << "\t" << get<8> (repeat_loci[i]) << "\n";
    }

    repeat_loci.erase(repeat_loci.begin(), repeat_loci.begin() + end_index + 1);
}


bool expectedPurityDifference(string parent_motif, double parent_purity, string child_motif, int child_len, double child_purity) {
    /*
     * expected difference in purity between a parent and a child repeat
     * @param parent_motif motif sequence of the parent repeat
     * @param parent_puriy purity of the parent repeat
     * @param child_motif sequence of the child repeat
     * @param child_length length of the child repeat
     * @param child_puriy purity of the parent repeat
     * @return bool boolean value if the child repeat should be retained
    */

    if (parent_motif.length() < child_motif.length()) { // motif length of parent STR is shorter than motif length of nested STR
        double least_d = INFINITY;      // least edit distance
        int parent_units = 0;           // number of parent units

        // identifying length of perfect parent STR that has least edit distance with one motif of nested STR
        for (int i=1; i<parent_motif.length(); i++) {
            string extended_parent_motif = "";
            for (int _=0; _ < i; _++) extended_parent_motif += parent_motif;
            int edit_d = levenshteinDistance(extended_parent_motif, child_motif);
            if (edit_d < least_d) {
                least_d = edit_d;
                parent_units = i;
            }
        }
        for (int i=0; i<parent_units-1; i++) parent_motif += parent_motif;
    }

    else if (child_motif.length() < parent_motif.length()) { // motif length of the nested STR is shorter than motif length of parent STR
        double least_d = INFINITY;      // least edit distance
        int nested_units = 0;           // number of nested units

        // identifying the length of nested STR that has least edit distance with one motif of parent STR
        for (int i=1; i<child_motif.length(); i++) {
            string extended_child_motif = "";
            for (int _=0; _ < i; _++) extended_child_motif += child_motif;
            int edit_d = levenshteinDistance(extended_child_motif, parent_motif);
            if (edit_d < least_d) {
                least_d = edit_d;
                nested_units = i;
            }
        }
        for (int i=0; i<nested_units-1; i++) child_motif += child_motif;
    }


    // identifying the least edit distance between different cyclical variations of motifs of parent and nested STRs
    double edit_d = INFINITY;
    for (int i=0; i < parent_motif.length(); i++) {
        string pm = parent_motif.substr(i) + parent_motif.substr(0,i);
        for (int j=0; j < child_motif.length(); j++) {
            string cm = child_motif.substr(j) + child_motif.substr(0,j);
            if (levenshteinDistance(pm, cm) < edit_d) {
                edit_d = levenshteinDistance(pm, cm);
            }
        }
    }

    int parent_imperfections = ((1 - parent_purity) * child_len);       // number of imperfections in the parent STR based on the length of nested STR
    int total_edit_d = (child_len / child_motif.length()) * edit_d;     // possible total edit distance between parent STR and nested STR for length of nested STR

    // if the purity of nested STR is greater than the edit distance between parent and nested STR
    // retain the nested STR
    return !(child_purity > (1 - ( ((double)(abs(parent_imperfections - total_edit_d))) / ((double)(child_len)) )));
}
    

tuple<vector<int>, vector<char>> cigarSplit(string cigar){
    /*
     * cigar string is split to the operation and the length
     * @param cigar cigar string of the repeat
     * @return tuple<vector<int>, vector<char>> two vectors one with cigar lengths and the other with cigar operations
    */
    string length = "";
    vector<int> clens; vector<char> ctypes;

    for (char c: cigar) {
        if (isdigit(c)) length += c;
        else {
            clens.push_back(stoi(length));
            ctypes.push_back(c);
            length = "";
        }
    }

    return {clens, ctypes};
}


tuple<vector<int>, vector<char>> extractNonOverlapCigar(int a_end, int b_start, vector<int> b_clens, vector<char> b_ctypes) {
     /*
     * extracts the non-overlapping cigar from the downstream locus of two overlapping loci 
     * @param a_end     end position of the upstream locus
     * @param b_start   start position of the downstream locus
     * @param b_clens   the lengths of the cigar operations for the downstream locus
     * @param b_ctypes  the consecutive cigar operations of the downstream locus
     * @return tuple<vector<int>, vector<char>> the lengths of continuous cigar operations and continuous cigar operations
    */
    vector<int> nover_clens; vector<char> nover_ctypes;
    int rpos = b_start, i = 0;
    for (i=0; i<b_clens.size(); i++) {  // pass through the locus
        int clen = b_clens[i]; char ctype = b_ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }

        if (rpos == a_end) {
            // if reached the end of the upstream overlapping locus; record the remaining cigar
            nover_clens = vector<int>(b_clens.begin() + i + 1, b_clens.end());
            nover_ctypes = vector<char>(b_ctypes.begin() + i + 1, b_ctypes.end());
            return {nover_clens, nover_ctypes};
        }
        else if (rpos > a_end) {
            // if reached beyond the end of the upstream overlapping locus
            if (i < b_clens.size() - 1) {
                nover_clens.push_back(rpos-a_end);
                for (int _=i; _<b_clens.size(); _++) { nover_clens.push_back(b_clens[_]); }
                nover_ctypes = vector<char>(b_ctypes.begin() + i, b_ctypes.end());
                return {nover_clens, nover_ctypes};
            }
            else if (i == b_clens.size() - 1){
                nover_clens.push_back(rpos-a_end);
                nover_ctypes = vector<char>(b_ctypes.begin() + i, b_ctypes.end());
                return {nover_clens, nover_ctypes};
            }
        }
    }
    
    if (rpos >= a_end) {
        // if the end is reached only at the last cigar operation
        if (rpos == a_end) return {nover_clens, nover_ctypes};
        else if (rpos > a_end) {
            nover_clens.push_back(rpos-a_end);
            nover_ctypes.push_back(b_ctypes[i]);
            return {nover_clens, nover_ctypes};
        }
    }
}


tuple<string, double> mergeRepeats(int a_end, int b_start, string a_cigar, string b_cigar) {
    /*
     * extracts the non-overlapping cigar from the downstream locus of two overlapping loci 
     * @param a_end     end position of the upstream locus
     * @param b_start   start position of the downstream locus
     * @param b_clens   the lengths of the cigar operations for the downstream locus
     * @param b_ctypes  the consecutive cigar operations of the downstream locus
     * @return tuple<vector<int>, vector<char>> the lengths of continuous cigar operations and continuous cigar operations
    */
    tuple<vector<int>, vector<char>> a_splitcigar_values = cigarSplit(a_cigar);
    vector<int>  a_clens  = get<0> (a_splitcigar_values);
    vector<char> a_ctypes = get<1> (a_splitcigar_values);    
    tuple<vector<int>, vector<char>> b_splitcigar_values = cigarSplit(b_cigar);
    vector<int>  b_clens  = get<0> (b_splitcigar_values);
    vector<char> b_ctypes = get<1> (b_splitcigar_values);

    tuple<vector<int>, vector<char>> nover_splitcigar_values = extractNonOverlapCigar(a_end, b_start, b_clens, b_ctypes);
    vector<int>  nover_clens   = get<0> (nover_splitcigar_values);
    vector<char> nover_ctypes  = get<1> (nover_splitcigar_values);

    vector<int> merged_clens;
    for (int _=0; _< a_clens.size(); _++)     { merged_clens.push_back(a_clens[_]); }
    for (int _=0; _< nover_clens.size(); _++) { merged_clens.push_back(nover_clens[_]); }

    vector<char> merged_ctypes;
    for (int _=0; _< a_ctypes.size(); _++)     { merged_ctypes.push_back(a_ctypes[_]); }
    for (int _=0; _< nover_ctypes.size(); _++) { merged_ctypes.push_back(nover_ctypes[_]); }
    
    string merged_cigar = ""; //.join([f'{clens[i]}{ctypes[i]}' for i in range(len(clens))])

    int alignment_length = 0;
    for (int _=0; _< merged_clens.size(); _++) {
        alignment_length += merged_clens[_];
    }

    int matches = 0, mismatches = 0;
    int clen = 0; char ctype;
    for (int _=0; _<merged_ctypes.size(); _++) {
        clen = merged_clens[_]; ctype = merged_ctypes[_];
        
        merged_cigar += to_string(clen); merged_cigar += ctype;
        
        if (ctype == '=' || ctype == 'M') { matches += clen; }
        else { mismatches += clen; }
    }
    
    double merge_purity = ((double)matches)/((double)alignment_length);

    return { merged_cigar, merge_purity };
}
    

// Define a comparison function for tuples (e.g., comparing the first element)
bool compareRepeatLoci(const tuple<string, int, int, string, double, string, int, int, int> &a,
                       const tuple<string, int, int, string, double, string, int, int, int> &b) {
    
    return get<1> (a) < get<1> (b); // Compare based on the first element (int)
}


void addLocusToOutput(string &sequence_id, int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                      int atomicity, int repeat_length, int repeat_units, ostream &out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    

    string last_seqid = "";
    if (repeat_loci.size() > 0) { last_seqid = get<0> (repeat_loci[repeat_loci.size()-1]); }
    if (last_seqid != "" && sequence_id != last_seqid) {
        printRepeatsToOutput(out, repeat_loci, repeat_loci.size()-1);
        tuple<string, int, int, string, double, string, int, int, int> repeat_locus = {sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, atomicity, repeat_length, repeat_units};
        auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
        // Insert the new element at the found position
        repeat_loci.insert(it, repeat_locus);
        return;
    }

    int last_start, last_end, last_length;
    string last_motif, last_cigar; double last_purity;
    vector<int> remove_loci;
    int i = 0;
    for (i=repeat_loci.size()-1; i >= 0; i--) {
        last_start  = get<1> (repeat_loci[i]);
        last_end    = get<2> (repeat_loci[i]);
        last_motif  = get<3> (repeat_loci[i]);
        last_purity = get<4> (repeat_loci[i]);
        last_cigar  = get<5> (repeat_loci[i]);
        last_length = get<7> (repeat_loci[i]);

        if (repeat_start == last_start && repeat_end == last_end) {
            // new location is nested in previous location
            if (purity <= last_purity) return;
            else remove_loci.push_back(i);
        }

        else if (last_start <= repeat_start && repeat_end <= last_end) {
            // new location is nested in previous location
            if (purity <= last_purity || motif == last_motif || expectedPurityDifference(last_motif, last_purity, motif, repeat_length, purity)) {
                return;
            }
        }

        else if (repeat_start <= last_start && last_end <= repeat_end) {
            // new location is nested in previous location
            if (last_purity <= purity || motif == last_motif || expectedPurityDifference(motif, purity, last_motif, last_length, last_purity)) {
                remove_loci.push_back(i);
            }
        }

        else {
            if (last_start <= repeat_start && repeat_start <= last_end) {     // last-repeat is upstream of current repeat
                
                // STR-i and STR-j have the same motif ~ Could happen if they are identified from different motif shifts
                if (motif == last_motif) {
                    tuple <string, double>merge_values = mergeRepeats(last_end, repeat_start, last_cigar, cigar_string);
                    remove_loci.push_back(i);
                    repeat_start = last_start; cigar_string = get<0> (merge_values); purity = get<1> (merge_values);
                }
                
                // current repeat is shorter 
                else if ((repeat_length < last_length) && ( ((repeat_end - last_end) < atomicity) || ((repeat_end - last_end) < 3) )) {
                    return;
                }

                // last repeat is shorter
                else if ((repeat_length > last_length) && ( ((repeat_start - last_start) < last_motif.size()) || ((repeat_start - last_start) < 3) )) {
                    remove_loci.push_back(i);
                }
                
                // STR-i length equals STR-j length ~ choose repeat with higher purity
                else if (repeat_length == last_length) {
                    if ( (((repeat_start - last_start) < last_motif.size()) || ((repeat_start - last_start) < 3))  && (last_purity < purity)) { remove_loci.push_back(i); }
                    else if (( ((repeat_end - last_end) < atomicity) || ((repeat_end - last_end) < 3) ) && (purity < last_purity)) { return; }
                }
            }

            else if (repeat_start <= last_start && last_start <= repeat_end) {  // last-repeat is downstream of current repeat

                if (motif == last_motif) {
                    tuple <string, double>merge_values = mergeRepeats(repeat_end, last_start, cigar_string, last_cigar);
                    remove_loci.push_back(i);
                    repeat_end = last_end; cigar_string = get<0> (merge_values); purity = get<1> (merge_values);
                }

                // STR-i is shorter than STR-j ~ Unique length of STR-i is shorter than STR-i motif length or less than 3 bp
                else if ((repeat_length > last_length) && ( ((last_end - repeat_end) < last_motif.size()) || ((last_end - repeat_end) < 3) )) {
                    remove_loci.push_back(i);
                }
                
                // STR-i is shorter than STR-j ~ Unique length of STR-i is shorter than STR-i motif length or less than 3 bp
                else if ((repeat_length < last_length) && ( ((last_start - repeat_start) < atomicity) || ((last_start - repeat_start) < 3) )) {
                    return;
                }
                
                // STR-i length equals STR-j length ~ choose repeat with higher purity
                else if (repeat_length == last_length) {
                    if (( ((last_end - repeat_end) < last_motif.length()) || ((last_end - repeat_end) < 3) ) && (last_purity < purity)) { remove_loci.push_back(i); }
                    else if ( (((last_start - repeat_start) < atomicity) || ((last_start - repeat_start) < 3))  && (purity < last_purity)) { return; }
                }
            }
        }
        
        if (last_end < repeat_start) {
            break;
        }
    }

    if (remove_loci.size()>0) {
        for (int _=0; _<remove_loci.size(); _++) {
            int i = remove_loci[_];
            repeat_loci.erase(repeat_loci.begin() + remove_loci[_], repeat_loci.begin() + remove_loci[_] + 1);
        }
    }

    for (int j=i; j>=0; j--) {
        if (repeat_start - get<2> (repeat_loci[j]) > 10000) {
            printRepeatsToOutput(out, repeat_loci, j);
            break;
        }
    }

    tuple<string, int, int, string, double, string, int, int, int> repeat_locus = {sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, atomicity, repeat_length, repeat_units};
    auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
    // Insert the new element at the found position
    repeat_loci.insert(it, repeat_locus);
}
