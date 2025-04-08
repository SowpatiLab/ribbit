#include "cigar_utils.h"

using namespace std;
using namespace boost;


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


int getAlignmentLength(string &cigar) {
    /*
     * calculates the alignment length from the cigar string
     * @param cigar the cigar string
     * @return int the alignment length
    */

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int> clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);
    int alignment_length = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M' || ctypes[_] == 'X' || ctypes[_] == 'I') {
            alignment_length += clens[_];
        }
    }
    return alignment_length;
}


int getAlignmentLength(tuple<vector<int>, vector<char>> &cigar_values) {
    /*
     * calculates the alignment length from the lengths and types of cigar operations
     * @param cigar_values tuple of lengths and types of cigar operations
     * @return int the alignment length
    */
    vector<int> clens = get<0> (cigar_values);
    vector<char> ctypes = get<1> (cigar_values);
    int alignment_length = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M' || ctypes[_] == 'X' || ctypes[_] == 'I') {
            alignment_length += clens[_];
        }
    }
    return alignment_length;
}


int getAlignmentLength(vector<int> &clens, vector<char> &ctypes) {
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


int getAlignmentLengthWithDeletions(string &cigar) {
    /*
     * calculates the alignment length from the cigar string
     * @param cigar the cigar string
     * @return int the alignment length
    */

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int> clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);
    int alignment_length = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M' || ctypes[_] == 'X' || ctypes[_] == 'I' || ctypes[_] == 'D') {
            alignment_length += clens[_];
        }
    }
    return alignment_length;
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


int getMatches(string &cigar) {
    /*
     * calculates the number of matches from the cigar string
     *  @param cigar the cigar string
     *  @param matches the number of matches
     *  @param longest_match the length of the longest match
    */

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int> clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);
    int matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
        }
    }

    return matches;
}


int getMatches(tuple<vector<int>, vector<char>> &cigar) {
    /*
     * calculates the number of matches from the lengths and types of cigar operations
     *  @param cigar tuple of lengths and types of cigar operations
     *  
    */

    vector<int> clens = get<0> (cigar);
    vector<char> ctypes = get<1> (cigar);
    int matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
        }
    }

    return matches;
}


int getMatches(vector<int> &clens, vector<char> &ctypes) {
    /*
     * calculates the number of matches from the lengths and types of cigar operations
     * @param clens vector of lengths of cigar operations
     * @param ctypes vector of types of cigar operations
     *  @param matches the number of matches
     *  @param longest_match the length of the longest match
    */
    int matches = 0;
    for (int _=0; _<clens.size(); _++) {
        if (ctypes[_] == '=' || ctypes[_] == 'M') {
            matches += clens[_];
        }
    }

    return matches;
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


void extendTillMatch(string &flank_cigar, string &cigar, int &position, bool end) {
    /*
     * extends a repeat till matches for within the non-overlapping region
     *  @param cigar cigar string for the non-overlapping part of the repeat
     *  @param end extend the end of the cigar; if false extends the start of the cigar
    */
    string length = "";
    if (end) {
        for (int i = 0; i<flank_cigar.length(); i++) {
            if (isdigit(flank_cigar[i])) length += flank_cigar[i];
            else {
                if (flank_cigar[i] == 'M' || flank_cigar[i] == '=') {
                    cigar += length + 'M'; position += stoi(length);
                } return;
            }
        }
    }

    else {
        if (flank_cigar[flank_cigar.length()-1] == 'M' || flank_cigar[flank_cigar.length()-1] == '=') {
            for (int i = flank_cigar.length()-2; i >= 0; i--) {
                if (isdigit(flank_cigar[i])) length += flank_cigar[i];
                else { cigar = length + 'M' + cigar; position -= stoi(length); return; }
            }
        }
    }

    return;
}


void extendTillMatch(tuple<vector<int>, vector<char>> &flank_cigar, string &cigar, int &position, bool end) {
    /*
     * extends a repeat till matches for within the non-overlapping region
     *  @param cigar cigar string for the non-overlapping part of the repeat
     *  @param end the end position of the repeat
    */
    string length = "";
    vector<int> clens = get<0> (flank_cigar);
    vector<char> ctypes = get<1> (flank_cigar);
    int n = clens.size();

    if (end) {
        if (ctypes[0] == 'M' || ctypes[0] == '=') {
            cigar = cigar + to_string(clens[0]) + 'M';
            position = position + clens[0]; return;
        }
    }

    else {
        if (ctypes[ctypes.size()-1] == 'M' || ctypes[ctypes.size()-1] == '=') {
            cigar = to_string(clens[clens.size()-1]) + 'M' + cigar; position -= clens[clens.size()-1]; return;
        }
    }

    return;
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


void separateCigarsOverlappingLoci(int upstart, int upend, tuple<vector<int>, vector<char>> &up_cigarvalues,
                                   int dnstart, int dnend, tuple<vector<int>, vector<char>> &dn_cigarvalues,
                                   tuple<vector<int>, vector<char>> &up_olseg_cigarvalues, tuple<vector<int>, vector<char>> &dn_olseg_cigarvalues,
                                   tuple<vector<int>, vector<char>> &up_nolseg_cigarvalues, tuple<vector<int>, vector<char>> &dn_nolseg_cigarvalues) {
    /*
     * separates the overlapping cigars of two loci into overlapping and non-overlapping cigars
     *  @param upstart start position of the upstream locus
     *  @param upend end position of the upstream locus
     *  @param upcigar CIGAR of the upstream locus
     *  @param dnstart start position of the downstream locus
     *  @param dnend end position of the downstream locus
     *  @param dncigar CIGAR of the downstream locus
     *  @param up_olseg_cigarvalues the lengths of continuous cigar operations and continuous cigar operations of the upstream locus
     *  @param dn_olseg_cigarvalues the lengths of continuous cigar operations and continuous cigar operations of the downstream locus
     *  @param up_nolseg_cigarvalues the lengths of non-overlapping cigar operations and non-overlapping cigar operations of the upstream locus
     *  @param dn_nolseg_cigarvalues the lengths of non-overlapping cigar operations and non-overlapping cigar operations of the downstream locus
    */

    int uppos = upstart, dnpos = dnstart;
    vector<int>  up_clens  = get<0> (up_cigarvalues);
    vector<char> up_ctypes = get<1> (up_cigarvalues);
    vector<int>  dn_clens  = get<0> (dn_cigarvalues);
    vector<char> dn_ctypes = get<1> (dn_cigarvalues);

    int i, clen; char ctype;
    int rpos = upstart, end = dnstart;
    for (i=0; i<up_clens.size(); i++) {
        clen = up_clens[i]; ctype = up_ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }
        get<0>(up_nolseg_cigarvalues).push_back(clen);
        get<1>(up_nolseg_cigarvalues).push_back(ctype);
        if (rpos == end) { 
            if (i+1 < up_ctypes.size() && up_ctypes[i+1] == 'D') {
                i += 1;
                get<0>(up_nolseg_cigarvalues).push_back(up_clens[i]);
                get<1>(up_nolseg_cigarvalues).push_back(up_ctypes[i]);
            }
            break;
        }
        if (rpos > end) break;
    }
    if (rpos > end) {
        get<0>(up_nolseg_cigarvalues)[i] -= (rpos - end);
        up_clens[i] = rpos - end; i = i - 1;
    }
    up_olseg_cigarvalues = { vector<int>(up_clens.begin() + i + 1, up_clens.end()), vector<char>(up_ctypes.begin() + i + 1, up_ctypes.end()) };

    rpos = dnstart; end = upend;
    for (i=0; i<dn_clens.size(); i++) {
        clen = dn_clens[i]; ctype = dn_ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }
        get<0>(dn_olseg_cigarvalues).push_back(clen);
        get<1>(dn_olseg_cigarvalues).push_back(ctype);
        if (rpos >= end) break;
    }
    if (rpos > end) {
        get<0>(dn_olseg_cigarvalues)[i] -= (rpos - end);
        dn_clens[i] = rpos - end; i = i - 1;
    }
    dn_nolseg_cigarvalues = { vector<int>(dn_clens.begin() + i + 1, dn_clens.end()), vector<char>(dn_ctypes.begin() + i + 1, dn_ctypes.end()) };
}

