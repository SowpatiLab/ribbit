#include "output_utils.h"

using namespace std;
using namespace boost;


void printIsolatedRepeats(vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, ostream *out) {
    /*
     *  prints the repeats to output file merging the overlapping ones and reporting the nested ones in a separate column
     *  @param repeat_loci vector of tuples containing the repeat loci information
     *  @param out the output file stream
     */

    cleanCigar(get<5> (repeat_loci[0]));
    *out << get<0>(repeat_loci[0]) << "\t" << get<1>(repeat_loci[0]) << "\t" << get<2>(repeat_loci[0]) << "\t" 
         << get<3>(repeat_loci[0]) << "\t" << get<4>(repeat_loci[0]) << "\t" << get<6>(repeat_loci[0]) << "\t"
         << get<7>(repeat_loci[0]) << "\t" << get<8>(repeat_loci[0]) << "\tI";
    if (CIGAROUTPUT) { *out << ":" << get<5>(repeat_loci[0]); }
    *out << "\n";
}


void printMergedRepeats(vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, ostream* out) {
    /*
     *  prints the repeats to output file merging the overlapping ones and reporting the nested ones in a separate column
     *  @param repeat_loci vector of tuples containing the repeat loci information
     *  @param out the output file stream
     */

    cleanCigar(get<5> (repeat_loci[0]));
    *out << get<0>(repeat_loci[0]) << "\t" << get<1>(repeat_loci[0]) << "\t" << get<2>(repeat_loci[0]) << "\t" 
         << get<3>(repeat_loci[0]) << "\t" << get<4>(repeat_loci[0]) << "\t" << get<6>(repeat_loci[0]) << "\t" 
         << get<7>(repeat_loci[0]) << "\t" << get<8>(repeat_loci[0]) << "\t";

    *out << "M:";
    if (CIGAROUTPUT) { *out << get<5>(repeat_loci[0]) << ":"; }
    *out << repeat_loci.size()-1 << ":";
    for (int j=1; j<repeat_loci.size(); j++) {
        cleanCigar(get<5> (repeat_loci[j]));
        *out << get<1> (repeat_loci[j]) << "-" << get<2> (repeat_loci[j]) << "-" << get<6> (repeat_loci[j]) << "-" << get<4> (repeat_loci[j]);
        if (j != repeat_loci.size() - 1) { *out << ","; }
    }
    *out << ":";
    for (int j=1; j<repeat_loci.size(); j++) {
        *out << get<3> (repeat_loci[j]);
        if (j != repeat_loci.size() - 1) { *out << ","; }
    }
    if (CIGAROUTPUT) {
        *out << ":";
        for (int j=1; j<repeat_loci.size(); j++) {
            *out << get<5> (repeat_loci[j]);
            if (j != repeat_loci.size() - 1) { *out << ","; }
        }
    }
    *out << "\n";
}


bool checkNonSupportCoverage(vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {


    sort(repeat_loci.begin(), repeat_loci.end(), [](const tuple<string,int,int,string,double,string,int,int,int> &a, const tuple<string,int,int,string,double,string,int,int,int> &b) {
        if (get<1>(a) == get<1>(b)) {
            if (get<2>(a) == get<2>(b)) {
                if (get<6>(a) == get<6>(b)) {
                    // higher purity first
                    return get<4>(a) > get<4>(b);
                }
                // lower motif length first
                return get<6>(a) < get<6>(b);
            }
            // higher end position first
            return get<2>(a) > get<2>(b);
        }

        // lower start position first
        return get<1>(a) < get<1>(b);
    });

    int repeat_start = get<1>(repeat_loci[0]), repeat_end = get<2>(repeat_loci[0]);
    int repeat_length = get<7>(repeat_loci[0]);
    int motif_length = get<6>(repeat_loci[0]);
    int threshold;
    if      (repeat_length < 100) threshold = repeat_length;
    else if (repeat_length < 1000) threshold = repeat_length - 20;
    else                         threshold = repeat_length - 100;

    int nstart, nend, nmlen;
    unordered_map<int, int> coverage_map;
    unordered_map<int, int> nmlen_end;
    for (int i=1; i<repeat_loci.size(); i++) {
        nstart = get<1>(repeat_loci[i]);
        nend   = get<2>(repeat_loci[i]);
        nmlen  = get<6>(repeat_loci[i]);
        if (nmlen == motif_length) continue;
        // accumulate non-support coverage per motif length, ensuring we don't double-count overlaps
        auto it = nmlen_end.find(nmlen);
        if (it == nmlen_end.end()) {
            nmlen_end[nmlen] = nend;
            coverage_map[nmlen] = nend - nstart;
        }
        else {
            if (nstart >= it->second) {
                coverage_map[nmlen] += nend - nstart;
                nmlen_end[nmlen] = nend;
            }
            else if (nend > it->second) {
                coverage_map[nmlen] += nend - it->second;
                nmlen_end[nmlen] = nend;
            }
        }
    }

    for (auto const& [nmlen, cov] : coverage_map) {
        if (cov >= threshold) { return false; }
    }
    return true;
}


void mergeRepeats(vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, 
                  ostream *out, int &end_index) {
    /*
     *  merges the repeats in the repeat loci and reports the nested ones in a separate column
     *  @param repeat_loci vector of tuples containing the repeat loci information
     *  @param end_index the index in the repeat loci vector till which the repeats should be merged 
     */

    vector<tuple<string, int, int, string, double, string, int, int, int>> overlapping_repeats;

    *out << std::setprecision(2) << fixed;
    int didx, repeat_start, repeat_end, motif_length;
    string motif;
    double purity;
    int merge_start = -1, merge_end = -1;
    for (int i=0; i < repeat_loci.size(); i++) {
        repeat_start = get<1>(repeat_loci[i]);
        repeat_end = get<2>(repeat_loci[i]);
        motif = get<3>(repeat_loci[i]);
        motif_length = get<6>(repeat_loci[i]);

        if (i == 0) {
            overlapping_repeats.push_back(repeat_loci[i]);
            merge_start = repeat_start; merge_end = repeat_end;
            continue;
        }

        if (repeat_end <= get<2>(overlapping_repeats[0]) ) {
            // add repeats to the set of overlapping repeats
            overlapping_repeats.push_back(repeat_loci[i]);
        }

        else {

            if (overlapping_repeats.size() == 1) printIsolatedRepeats(overlapping_repeats, out);
            else if (overlapping_repeats.size() > 1) {
                if (checkNonSupportCoverage(overlapping_repeats)) printMergedRepeats(overlapping_repeats, out);
                else {
                    overlapping_repeats.erase(overlapping_repeats.begin());
                    int idx = overlapping_repeats.size() - 1;
                    mergeRepeats(overlapping_repeats, out, idx);
                }
            }

            overlapping_repeats.clear();
            if (i <= end_index) {
                overlapping_repeats.push_back(repeat_loci[i]);
                merge_start = repeat_start; merge_end = repeat_end;
            }
            else { end_index = i-1; break; }
        }

    }

    if (overlapping_repeats.size() == 1) printIsolatedRepeats(overlapping_repeats, out);
    else if (overlapping_repeats.size() > 1) printMergedRepeats(overlapping_repeats, out);
}


void printRepeatsToOutput(ostream *out, vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                          int end_index) {
    /*
     *  print repeats from the recorded repeat loci and clears the repeat loci
     *  @param out outstream of the output file
     *  @param repeat_loci vector of recorded repeat loci sorted based on the start
     *  @param end_index index to which point repeats should be printed to the output
     */

    mergeRepeats(repeat_loci, out, end_index);
    repeat_loci.erase(repeat_loci.begin(), repeat_loci.begin() + end_index + 1);
}


int absolute(int x) {
    /*
     *  returns the absolute value of an integer
     *  @param x integer
     *  @return int absolute value of the integer
     */
    if (x < 0) return -x;
    return x;
}


bool checkCyclicalVariation(string query, string ref) {
    /*
     *  checks if a query motif is cyclical variation of reference motif
     *  @param query sequence of the query motif 
     *  @param ref sequence of reference motif
     *  @return bool if the query motif is a cyclical variation of the reference motif
     */
    if (query.length() != ref.length()) return false;

    string cycle;
    for (int _=0; _<query.length(); _++) {
        cycle =  query.substr(_, query.length() - _) + query.substr(0,_);
        if (ref == cycle) return true;
    }

    return false;
}


int leastDistance(string reference, string query) {
    /*
     *  calculates the least distance between two strings
     *  @param reference string to be compared with
     *  @param query string to be compared
     *  @return int the least distance between the two strings
     */
    StripedSmithWaterman::Aligner   aligner;
    StripedSmithWaterman::Filter    filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(query.c_str(), reference.c_str(), reference.length(), filter, &alignment, 15);
    string cigar = alignment.cigar_string;
    return alignment.mismatches + ((reference.length()/2) - getAlignmentLength(cigar));
}


tuple<vector<int>, vector<char>> extractNonOverlapCigar(int a_end, int b_start, vector<int> dn_clens, vector<char> dn_ctypes) {
    /*
     *  extracts the non-overlapping cigar from the downstream locus of two overlapping loci
     *  @param a_end     end position of the upstream locus
     *  @param b_start   start position of the downstream locus
     *  @param dn_clens   the lengths of the cigar operations for the downstream locus
     *  @param dn_ctypes  the consecutive cigar operations of the downstream locus
     *  @return tuple<vector<int>, vector<char>> the lengths of continuous cigar operations and continuous cigar operations
     */
    vector<int> nol_clens; vector<char> nol_ctypes;
    int rpos = b_start, i = 0;
    int clen; char ctype;
    for (i=0; i<dn_clens.size(); i++) {  // pass through the locus
        clen = dn_clens[i]; ctype = dn_ctypes[i];
        if (ctype == '=' || ctype == 'M' || ctype == 'X' || ctype == 'I') {
            rpos += clen;
        }

        if (rpos == a_end) {
            // if reached the end of the upstream overlapping locus; record the remaining cigar
            nol_clens = vector<int>(dn_clens.begin() + i + 1, dn_clens.end());
            nol_ctypes = vector<char>(dn_ctypes.begin() + i + 1, dn_ctypes.end());
            return {nol_clens, nol_ctypes};
        }
        else if (rpos > a_end) {
            // if reached beyond the end of the upstream overlapping locus
            if (i < dn_clens.size() - 1) {
                nol_clens = {rpos-a_end};
                for (int _=i+1; _<dn_clens.size(); _++) { nol_clens.push_back(dn_clens[_]); }
                nol_ctypes = vector<char>(dn_ctypes.begin() + i, dn_ctypes.end());
                return {nol_clens, nol_ctypes};
            }
            else if (i == dn_clens.size() - 1){
                nol_clens.push_back(rpos-a_end);
                nol_ctypes = vector<char>(dn_ctypes.begin() + i, dn_ctypes.end());
                return {nol_clens, nol_ctypes};
            }
        }
    }

    if (rpos >= a_end) {
        // if the end is reached only at the last cigar operation
        if (rpos == a_end) return {nol_clens, nol_ctypes};
        else if (rpos > a_end) {
            nol_clens.push_back(rpos-a_end);
            nol_ctypes.push_back(dn_ctypes[i]);
            return {nol_clens, nol_ctypes};
        }
    }

    return {nol_clens, nol_ctypes};
}


void alignSequenceWithPerfectRepeat(string &sequence, string &motif, int &repeat_start, int &repeat_end,
                                    double &purity, string &new_cigar) {
    /*
     *  aligns the sequence with a perfect repeat of the motif
     *  @param sequence the sequence to be aligned
     *  @param motif the motif to be aligned with
     *  @param repeat_start start position of the repeat
     *  @param repeat_end end position of the repeat
     *  @param purity purity of the alignment
     *  @param new_cigar cigar string of the alignment
     */

    string perfect_repeat = "";
    int ppr_length = 0;
    int units = sequence.length()/motif.length();
    for (int _=0; _ <= units + 2; _++) {
        perfect_repeat += motif;
        ppr_length += motif.length();
    }
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(sequence.c_str(), perfect_repeat.c_str(), ppr_length, filter, &alignment, 15);
    string cigar = alignment.cigar_string;

    tuple<vector<int>, vector<char>> csplit = cigarSplit(cigar);
    vector<int>  clens = get<0> (csplit);
    vector<char> ctypes = get<1> (csplit);

    // initialise the repeat coordinates to the seed coordinates
    int alignment_length = 0;

    char ctype; int clength, cidx = 0;
    int qpos = 0;
    int matches = 0;

    bool mismatch_continue = false;         // to check if there are contiguous mismatches
    new_cigar = "";
    for ( ; cidx < clens.size(); cidx++) {
        clength = clens[cidx]; ctype = ctypes[cidx];

        switch (ctype) {
            case 'S':
                // soft clip: edit the repeat start and end
                if (cidx == 0) { repeat_start += clength; }
                else { repeat_end -= clength; }
                qpos += clength;
                break;

            case 'X':
                qpos += clength; alignment_length += clength;
                new_cigar += to_string(clength) + ctype;
                break;
            case 'I':
                qpos += clength; alignment_length += clength;
                new_cigar += to_string(clength) + ctype;
                break;
            case 'D':
                alignment_length += clength;
                new_cigar += to_string(clength) + ctype;
                break;
            case '=': case 'M':
                qpos += clength; alignment_length += clength;
                matches += clength;

                new_cigar += to_string(clength) + ctype;
                break;
            default: break;
        }
    }

    purity = (double)matches / (double)alignment_length;
}


bool onlyMatches(vector<char> &ctypes) {
    /*
     *  checks if the cigar operations are only matches
     *  @param ctypes vector of cigar operations
     *  @return bool if the cigar operations are only matches
     */
    if (ctypes.size() == 0 || ctypes.size() > 1) return false;
    if (ctypes[0] == '=' || ctypes[0] == 'M') return true;
    return true;
}


int getBoundaryOverlappingLoci(int upstart, int upend, string &upmotif, tuple<vector<int>, vector<char>> &up_cigarvalues, string &upcigar,
                               int dnstart, int dnend, string &dnmotif, tuple<vector<int>, vector<char>> &dn_cigarvalues, string &dncigar) {
    /*
     *  gets the boundary of the overlapping loci
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


void compareOverlappingLoci(int &upstart, int &upend, string &upcigar, string &upmotif, double &uppurity, bool &up_drop, bool &up_update,
                            int &dnstart, int &dnend, string &dncigar, string &dnmotif, double &dnpurity, bool &dn_drop, bool &dn_update) {
    /*
     *  compares the overlapping cigar of two overlapping loci and merges them
     *  @param upstart start position of the upstream locus
     *  @param upend end position of the upstream locus
     *  @param upcigar CIGAR of the upstream locus
     *  @param dnstart start position of the downstream locus
     *  @param dnend end position of the downstream locus
     *  @param dncigar CIGAR of the downstream locus
     */

    string full_sequence = "", new_cigar = "";

    tuple<vector<int>, vector<char>> up_cigarvalues = cigarSplit(upcigar);
    tuple<vector<int>, vector<char>> dn_cigarvalues = cigarSplit(dncigar);
    tuple<vector<int>, vector<char>> up_olseg_cigarvalues;
    tuple<vector<int>, vector<char>> up_nolseg_cigarvalues;
    tuple<vector<int>, vector<char>> dn_olseg_cigarvalues;
    tuple<vector<int>, vector<char>> dn_nolseg_cigarvalues;
    separateCigarsOverlappingLoci(upstart, upend, up_cigarvalues, dnstart, dnend, dn_cigarvalues,
                                  up_olseg_cigarvalues, dn_olseg_cigarvalues, up_nolseg_cigarvalues,
                                  dn_nolseg_cigarvalues);

    vector<int>  up_olclens  = get<0> (up_olseg_cigarvalues);
    vector<char> up_olctypes = get<1> (up_olseg_cigarvalues);
    vector<int>  dn_olclens  = get<0> (dn_olseg_cigarvalues);
    vector<char> dn_olctypes = get<1> (dn_olseg_cigarvalues);

    assert(getRepeatLength(up_olseg_cigarvalues) == getRepeatLength(dn_olseg_cigarvalues));

    if (onlyMatches(up_olctypes) && onlyMatches(dn_olctypes)) {
        // keep both the loci as they are
        // Nothing changes
    }

    else if (onlyMatches(up_olctypes) && !onlyMatches(dn_olctypes)) {
        // if the upstream locus has only matches and downstream locus has mismatches or indels
        // then
        dnstart = upend;
        dncigar = buildCigar(dn_nolseg_cigarvalues);
        extendTillMatch(dn_olseg_cigarvalues, dncigar, dnstart, false);
        int dn_nonoverlap_length = dnend - dnstart;
        if (dn_nonoverlap_length > MINIMUM_LENGTH[dnmotif.length()] && dn_nonoverlap_length > 2*dnmotif.length()) {
            // separate both the loci
            dnpurity = ((double) getMatches(dncigar)) / ((double) getAlignmentLength(dncigar));
            dn_update = true; return;
        }
        else {
            // merge the loci
            // get the genome sequence within the coordinates and align with the upmotif repeat
            dn_drop = true;
            full_sequence = SEQUENCE.substr(upstart, dnend-upstart);
            upend = dnend;
            alignSequenceWithPerfectRepeat(full_sequence, upmotif, upstart, upend, uppurity, new_cigar);
            upcigar = new_cigar;
            up_update = true; return;
        }
    }

    else if (!onlyMatches(up_olctypes) && onlyMatches(dn_olctypes)) {
        // if the downstream locus has only matches and upstream locus has mismatches or indels
        // then
        upend   = dnstart;
        upcigar = buildCigar(up_nolseg_cigarvalues);
        extendTillMatch(up_olseg_cigarvalues, upcigar, upend, true);
        int up_nonoverlap_length = dnstart - upstart;
        if (up_nonoverlap_length > MINIMUM_LENGTH[upmotif.size()] && up_nonoverlap_length > 2*upmotif.size()) {
            // separate both the loci
            uppurity = ((double) getMatches(upcigar)) / ((double) getAlignmentLength(upcigar));
            up_update = true; return;
        }
        else {
            // merge the loci
            up_drop = true;
            full_sequence = SEQUENCE.substr(upstart, dnend-upstart);
            dnstart = upstart;
            alignSequenceWithPerfectRepeat(full_sequence, dnmotif, dnstart, dnend, dnpurity, new_cigar);
            dncigar = new_cigar;
            dn_update = true; return;
        }
    }

    else {
        // if both the loci have mismatches or indels calculate an optimal boundary
        int boundary = getBoundaryOverlappingLoci(upstart, upend, upmotif, up_cigarvalues, upcigar,
                                                  dnstart, dnend, dnmotif, dn_cigarvalues, dncigar);

        uppurity = ((double) (getMatches(upcigar))) / ((double) (getAlignmentLength(upcigar)));
        dnpurity = ((double) (getMatches(dncigar))) / ((double) (getAlignmentLength(dncigar)));

        upend = boundary; dnstart = boundary;
        up_update = true; dn_update = true;
        if (((upend - upstart) < MINIMUM_LENGTH[upmotif.size()]) || ((upend - upstart) < 2*upmotif.size())) {
            up_drop = true;
        }
        if (((dnend - dnstart) < MINIMUM_LENGTH[dnmotif.size()]) || ((dnend - dnstart) < 2*dnmotif.size())) {
            dn_drop = true;
        }
        return;
    }

}


// Define a comparison function for tuples (e.g., comparing the first element)
bool compareRepeatLoci(const tuple<string, int, int, string, double, string, int, int, int> &a,
                       const tuple<string, int, int, string, double, string, int, int, int> &b) {
    /*
     *  compares repeat loci based on the start positions; if starts are same, returns the one with larger end position first
     *  @param a tuple of the first repeat location
     *  @param b tuple of the second repeat location
     *  @returns bool bool value indicating if the first repeat is before second repeat
     */
    if (get<1> (a) == get<1> (b)) {
        // if the start positions are same, compare based on the end positions
        return get<2> (a) > get<2> (b);
    }
    return get<1> (a) < get<1> (b); // Compare based on the first element (int)
}


tuple<string, double> mergeRepeatsIdenticalMotif(int upend, int dnstart, string upcigar, string dncigar) {
    /*
     *  merges the two overlapping loci with identical motifs
     *  @param upend end position of the upstream locus
     *  @param dnstart start position of the downstream locus
     *  @param upcigar CIGAR of the upstream locus
     *  @param dncigar CIGAR of the downstream locus
     *  @return tuple<string, double> merged CIGAR and purity of the merged locus
     */

    tuple<vector<int>, vector<char>> up_cigarvalues = cigarSplit(upcigar);
    vector<int>  up_clens  = get<0> (up_cigarvalues);
    vector<char> up_ctypes = get<1> (up_cigarvalues);
    tuple<vector<int>, vector<char>> dn_cigarvalues = cigarSplit(dncigar);
    vector<int>  dn_clens  = get<0> (dn_cigarvalues);
    vector<char> dn_ctypes = get<1> (dn_cigarvalues);

    tuple<vector<int>, vector<char>> nol_cigarvalues = extractNonOverlapCigar(upend, dnstart, dn_clens, dn_ctypes);
    vector<int>  nol_clens   = get<0> (nol_cigarvalues);
    vector<char> nol_ctypes  = get<1> (nol_cigarvalues);

    vector<int> merged_clens;
    for (int _=0; _< up_clens.size(); _++)     { merged_clens.push_back(up_clens[_]); }
    for (int _=0; _< nol_clens.size(); _++) { merged_clens.push_back(nol_clens[_]); }

    vector<char> merged_ctypes;
    for (int _=0; _< up_ctypes.size(); _++)     { merged_ctypes.push_back(up_ctypes[_]); }
    for (int _=0; _< nol_ctypes.size(); _++) { merged_ctypes.push_back(nol_ctypes[_]); }

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


bool handleOverlapRelation(string &sequence_id, int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                           int motif_length, int repeat_length, int repeat_units, ostream *out,
                           vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, vector<int> &remove_loci,
                           int &recursion_level, vector<tuple<string, int, int, string, double, string, int, int, int>> &new_repeat_loci) {
    /*
     *  handles nested-parent relationship of the current repeat with the recorded repeats
     *  @param repeat_start start position of the current repeat
     *  @param repeat_end end position of the current repeat
     *  @param motif motif of the current repeat
     *  @param motif_length length of the motif of the current repeat
     *  @param purity purity of the current repeat
     *  @param cigar_string CIGAR of the current repeat
     *  @param repeat_loci vector of recorded repeat loci sorted based on the start positions
     */

    int last_start, last_end, last_length, last_mlen, last_units;
    int segment_length; double segment_purity;
    bool up_trim, down_trim;
    string last_seqid,last_motif, last_cigar; double last_purity;
    int i = 0;

    for (i=repeat_loci.size()-1; i >= 0; i--) {
        if (i >= repeat_loci.size()) { i = repeat_loci.size()-1; }
        // start comparing repeats from the end
        last_start  = get<1> (repeat_loci[i]);
        last_end    = get<2> (repeat_loci[i]);

        if ((repeat_end < last_start) || (repeat_start > last_end)) {
            // continue if repeat doesn't overlap with the repeat
            continue;
        }

        last_seqid  = get<0> (repeat_loci[i]);
        last_motif  = get<3> (repeat_loci[i]);
        last_mlen   = get<6> (repeat_loci[i]);
        last_purity = get<4> (repeat_loci[i]);
        last_cigar  = get<5> (repeat_loci[i]);
        last_length = get<7> (repeat_loci[i]);

        // identical coordinates
        if (repeat_start == last_start && repeat_end == last_end) { continue; }

        // nested
        else if (last_start < repeat_start && repeat_end < last_end) { continue; }

        // parent
        else if (repeat_start < last_start && last_end < repeat_end) { continue; }

        else {
            // Only considering the repeats with overlap relationships
            if (motif_length > SMALL_MLEN_LIMIT && last_mlen > SMALL_MLEN_LIMIT) {
                // if both the motifs are larger than the small motif length limit
                if (motif_length != last_mlen) {

                    // the repeats are not of the same motif
                    // we calculate an optimal boundary point and either merge or separate the repeats
                    bool fail = false, last_fail = false;       // if current repeat or last repeat should be dropped
                    bool last_update = false, update = false;   // if current repeat or last repeat should be updated
                    // if the repeats are overlapping
                    if (last_start < repeat_start && repeat_start < last_end) {     // last-repeat is upstream of current repeat
                        compareOverlappingLoci(last_start, last_end, last_cigar, last_motif, last_purity, last_fail, last_update,
                                               repeat_start, repeat_end, cigar_string, motif, purity, fail, update);
                    }
    
                    else if (repeat_start < last_start && last_start < repeat_end) {  // last-repeat is downstream of current repeat
                        compareOverlappingLoci(repeat_start, repeat_end, cigar_string, motif, purity, fail, update,
                                               last_start, last_end, last_cigar, last_motif, last_purity, last_fail, last_update);
                    }
    
                    if (last_update && update) {
                        // if both needs to be updated
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        if (!last_fail) {
                            last_length = last_end - last_start;
                            last_units = last_length / last_mlen;
                            recursion_level += 1;
                            new_repeat_loci.push_back(make_tuple(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                                                 last_length, last_units));
                            addLocusToOutput(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                             last_length, last_units, out, repeat_loci, recursion_level, new_repeat_loci);
                        }
                        if (!fail) {
                            repeat_length = repeat_end - repeat_start;
                            repeat_units = repeat_length / motif_length;
                            recursion_level += 1;
                            new_repeat_loci.push_back(make_tuple(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length,
                                                                 repeat_length, repeat_units));
                            addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length, repeat_length,
                                             repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
                        }
                        return false;
                    }
    
                    else if (update & !(last_update)) {
                        // if only current repeat is updated
                        if (last_fail) { remove_loci.push_back(i); }
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        repeat_length = repeat_end - repeat_start;
                        repeat_units = repeat_length / motif_length;
                        recursion_level += 1;
                        new_repeat_loci.push_back(make_tuple(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length,
                                                             repeat_length, repeat_units));
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length, repeat_length,
                                         repeat_units, out, repeat_loci, recursion_level, new_repeat_loci);
                        return false;
                    }
    
                    else if (last_update & !(update)) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        remove_loci.clear();
                        last_length = last_end - last_start;
                        last_units = last_length / last_mlen;
                        recursion_level += 1;
                        new_repeat_loci.push_back(make_tuple(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                                             last_length, last_units));
                        addLocusToOutput(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                         last_length, last_units, out, repeat_loci, recursion_level, new_repeat_loci);
                        if (fail) return false;
                    }
                }

                else {
                    if (repeat_start == last_start || repeat_end == last_end) { continue; }
                    int upstart, upend, dnstart, dnend, uplength, dnlength;
                    string upcigar, dncigar, upmotif, dnmotif;
                    double uppurity, dnpurity;
                    tuple<vector<int>, vector<char>> up_cigarvalues, dn_cigarvalues;

                    int merge_start = 0, merge_end = 0, merge_length, merge_units;
                    string merge_cigar = ""; double merge_purity = 0.0;
                    string merge_motif = "";

                    if (last_start < repeat_start && repeat_start < last_end) {
                        upstart = last_start; upend = last_end; upcigar = last_cigar; upmotif = last_motif; uppurity = last_purity;
                        dnstart = repeat_start; dnend = repeat_end; dncigar = cigar_string; dnmotif = motif; dnpurity = purity;
                    }
                    else if (repeat_start < last_start && last_start < repeat_end) {
                        dnstart = last_start; dnend = last_end; dncigar = last_cigar; dnmotif = last_motif; dnpurity = last_purity;
                        upstart = repeat_start; upend = repeat_end; upcigar = cigar_string; upmotif = motif; uppurity = purity;
                    }
                    else if (last_end == repeat_start || last_start == repeat_end) {
                        if (last_start < repeat_start) {
                            merge_start = last_start; merge_end = repeat_end;
                            merge_cigar = last_cigar + cigar_string;
                        }
                        else if (repeat_start < last_start) {
                            merge_start = repeat_start; merge_end = last_end;
                            merge_cigar = cigar_string + last_cigar;
                        }
                        cleanCigar(merge_cigar);
                        merge_length = merge_end - merge_start;
                        merge_units  = merge_length / motif_length;
                        merge_motif = ((repeat_end-repeat_start)*purity >= (last_end-last_start)*last_purity) ? motif : last_motif;
                        merge_purity = ((double) (getMatches(merge_cigar))) / ((double) (getAlignmentLength(merge_cigar)));
                        recursion_level += 1;
                        new_repeat_loci.push_back(make_tuple(sequence_id, merge_start, merge_end, merge_motif, merge_purity, merge_cigar, motif_length,
                                                             merge_length, merge_units));
                        addLocusToOutput(sequence_id, merge_start, merge_end, merge_motif, merge_purity, merge_cigar,
                                         motif_length, merge_length, merge_units, out, repeat_loci, recursion_level, new_repeat_loci);
                        continue;
                    }

                    up_cigarvalues = cigarSplit(upcigar);
                    dn_cigarvalues = cigarSplit(dncigar);

                    int boundary = getBoundaryOverlappingLoci(upstart, upend, upmotif, up_cigarvalues, upcigar,
                                                              dnstart, dnend, dnmotif, dn_cigarvalues, dncigar);

                    uplength = boundary - upstart; dnlength = dnend - boundary;
                    uppurity = ((double) (getMatches(upcigar))) / ((double) (getAlignmentLength(upcigar)));
                    dnpurity = ((double) (getMatches(dncigar))) / ((double) (getAlignmentLength(dncigar)));

                    if (uplength * uppurity > dnlength * dnpurity) { merge_motif = upmotif; }
                    else { merge_motif = dnmotif; }
                    
                    
                    merge_cigar = upcigar + dncigar; // simply concatenate the cigars
                    cleanCigar(merge_cigar);
                    merge_purity = ((double) (getMatches(merge_cigar))) / ((double) (getAlignmentLength(merge_cigar)));
                    
                    merge_start = upstart; merge_end = dnend;
                    merge_length = merge_end - merge_start;
                    merge_units  = merge_length / motif_length;
                    assert(merge_length == getRepeatLength(merge_cigar));
                    recursion_level += 1;
                    new_repeat_loci.push_back(make_tuple(sequence_id, merge_start, merge_end, merge_motif, merge_purity, merge_cigar, motif_length,
                                                         merge_length, merge_units));
                    addLocusToOutput(sequence_id, merge_start, merge_end, merge_motif, merge_purity, merge_cigar,
                                     motif_length, merge_length, merge_units, out, repeat_loci, recursion_level, new_repeat_loci);
                }
            }

            else if (motif_length <= SMALL_MLEN_LIMIT && (motif == last_motif || checkCyclicalVariation(motif, last_motif))) {

                // if the repeats are of the same motif
                if (repeat_start == last_start || repeat_end == last_end) { continue; }

                tuple <string, double> merge_values;
                int merge_start = repeat_start, merge_end = repeat_end, merge_length, merge_units;
                string merge_cigar = "";
                double merge_purity = 0.0;

                // if the repeats are having the same motif then we merge them
                if (repeat_start == last_end || (last_start < repeat_start && repeat_start < last_end)) {
                    merge_values = mergeRepeatsIdenticalMotif(last_end, repeat_start, last_cigar, cigar_string);
                    merge_start = last_start;
                }

                else if (last_start == repeat_end || (repeat_start < last_start && last_start < repeat_end)) {
                    merge_values = mergeRepeatsIdenticalMotif(repeat_end, last_start, cigar_string, last_cigar);
                    merge_end = last_end;
                }

                merge_length = merge_end - merge_start;
                merge_units  = merge_length / motif_length;
                merge_cigar = get<0> (merge_values); merge_purity = get<1> (merge_values);
                recursion_level += 1;
                new_repeat_loci.push_back(make_tuple(sequence_id, merge_start, merge_end, motif, merge_purity, merge_cigar, motif_length,
                                                     merge_length, merge_units));
                addLocusToOutput(sequence_id, merge_start, merge_end, motif, merge_purity, merge_cigar, motif_length,
                                 merge_length, merge_units, out, repeat_loci, recursion_level, new_repeat_loci);
            }
        }
    }

    return true;
}


bool handleNestedParentRelation(int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                                int motif_length, int repeat_length, int repeat_units,
                                vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci, vector<int> &remove_loci) {
    /*
     *  handles nested-parent relationship of the current repeat with the recorded repeats
     *  @param repeat_start start position of the current repeat
     *  @param repeat_end end position of the current repeat
     *  @param motif motif of the current repeat
     *  @param motif_length length of the motif of the current repeat
     *  @param purity purity of the current repeat
     *  @param cigar_string CIGAR of the current repeat
     *  @param repeat_loci vector of recorded repeat loci sorted based on the start positions
     */

    int last_start, last_end, last_length, last_mlen, last_units;
    int segment_length; double segment_purity;
    bool up_trim, down_trim;
    bool up_equal, down_equal;
    string last_seqid,last_motif, last_cigar; double last_purity;
    int i = 0;

    for (i=repeat_loci.size()-1; i >= 0; i--) {
        if (i >= repeat_loci.size()) { i = repeat_loci.size()-1; }
        // start comparing repeats from the end
        last_start  = get<1> (repeat_loci[i]);
        last_end    = get<2> (repeat_loci[i]);

        if ((repeat_end < last_start) || (repeat_start > last_end)) {
            // continue if repeat doesn't overlap with the repeat
            continue;
        }

        last_seqid  = get<0> (repeat_loci[i]);
        last_motif  = get<3> (repeat_loci[i]);
        last_mlen   = get<6> (repeat_loci[i]);
        last_purity = get<4> (repeat_loci[i]);
        last_cigar  = get<5> (repeat_loci[i]);
        last_length = get<7> (repeat_loci[i]);

        // identical coordinates
        if (repeat_start == last_start && repeat_end == last_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) { return false; }
            else if (purity <= last_purity) return false;
            else if (purity > last_purity)  remove_loci.push_back(i);
        }

        // nested
        else if (last_start <= repeat_start && repeat_end <= last_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) {
                if (purity <= last_purity) return false;
            }

            else if (motif_length < last_mlen && motif_length > SMALL_MLEN_LIMIT && purity*repeat_length > last_purity*last_length) {
                // if the current motif is smaller than the last motif and the current repeat has more number of matching bases
                // than the last repeat, then we keep the current repeat
                remove_loci.push_back(i);

            }
            
            else {
                tuple<vector<int>, vector<char>> segment_cigarvalues = extractRegionCigar(last_cigar, repeat_start-last_start, repeat_end-last_start);
                double segment_purity = ((double) getMatches(segment_cigarvalues)) / ((double) getAlignmentLength(segment_cigarvalues));
                if (purity < segment_purity) return false;
            }
        }

        // parent
        else if (repeat_start <= last_start && last_end <= repeat_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) {
                if (last_purity <= purity) remove_loci.push_back(i);
            }

            else if (last_mlen < motif_length && last_mlen > SMALL_MLEN_LIMIT && last_purity*last_length > purity*repeat_length) {
                // if the last motif is smaller than the current motif and the last repeat has more number of matching bases
                // than the current repeat, then we keep the last repeat
                return false;
            }

            else {
                tuple<vector<int>, vector<char>> segment_cigarvalues = extractRegionCigar(cigar_string, last_start-repeat_start, last_end-repeat_start);
                segment_purity = ((double) getMatches(segment_cigarvalues)) / ((double) getAlignmentLength(segment_cigarvalues));
                if (last_purity < segment_purity) remove_loci.push_back(i);
            }
        }

        else {
            // handled in overlap function
        }

    }

    return true;
}


bool checkDuplicateEntry(string sequence_id, int repeat_start, int repeat_end, string motif,
                         vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
    /*
     *  checks if the current repeat locus is a duplicate entry in the recorded repeats
     *  @param sequence_id the name of the sequence the repeat is found in
     *  @param repeat_start the start position of the repeat
     *  @param repeat_end the end position of the repeat
     *  @param motif the sequence of the repeat motif
     *  @param repeat_loci the vector of tuples of all repeat loci
     *  @returns bool bool value indicating if the current repeat is a duplicate entry
     */

    for (int i=repeat_loci.size()-1; i >= 0; i--) {
        if (i >= repeat_loci.size()) { i = repeat_loci.size()-1; }

        // if the start position of the recorded repeat is less than current repeat, break
        if (get<1> (repeat_loci[i]) < repeat_start) { break; }
    
        // start comparing repeats from the end
        if (get<0> (repeat_loci[i]) != sequence_id) { break; } // if sequence ids are different, break

        if (get<1> (repeat_loci[i]) == repeat_start && get<2> (repeat_loci[i]) == repeat_end &&
            (get<3> (repeat_loci[i]) == motif || checkCyclicalVariation(get<3> (repeat_loci[i]), motif))) {
            return false;
        }
    }
    return true;
}


void addLocusToOutput(string sequence_id, int repeat_start, int repeat_end, string motif, double purity, string cigar_string,
                      int motif_length, int repeat_length, int repeat_units, ostream *out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci,
                      int &recursion_level, vector<tuple<string, int, int, string, double, string, int, int, int>> &new_repeat_loci) {
    /*
     * adds the repeat locus to set of all the repeats
     *  @param sequence_id the name of the sequence the repeat is found in
     *  @param repeat_start the start position of the repeat
     *  @param repeat_end the end position of the repeat
     *  @param motif the sequence of the repeat motif
     *  @param purity the purity value of the repeat
     *  @param cigar_string CIGAR of the repeat sequence alignment with perfec repeat
     *  @param motif_length the length of the repeating motif
     *  @param repeat_length the total length of th repeat locus
     *  @param repeat_units number of units of the motif in the repeat
     *  @param out the output file stream
     *  @param repeat_loci the vector of tuples of all repeat loci
     *  @returns void
     */

    int recursion_limit = 100;

    string last_seqid = "";
    if (repeat_loci.size() > 0) { last_seqid = get<0> (repeat_loci[repeat_loci.size()-1]); }
    else if (repeat_loci.size() == 0) {
        tuple<string, int, int, string, double, string, int, int, int> repeat_locus;
        repeat_locus = { sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                         motif_length, repeat_length, repeat_units };
        repeat_loci.push_back(repeat_locus);
        return;
    }

    // if the last sequence id is different
    if (last_seqid != "" && sequence_id != last_seqid) {
        // if the output stream is not null, print the repeats to the output
        // output stream is null in the case of the python module
        if (out) { printRepeatsToOutput(out, repeat_loci, repeat_loci.size()-1); }
        tuple<string, int, int, string, double, string, int, int, int> repeat_locus;
        repeat_locus = { sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                         motif_length, repeat_length, repeat_units };
        auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
        // Insert the new element at the found position
        repeat_loci.insert(it, repeat_locus);
        return;
    }

    if (recursion_level > recursion_limit) {
        set<string> added_repeats;
        for (int j = 0; j < new_repeat_loci.size(); j++) {
            string repeat_key = get<0>(new_repeat_loci[j]) + "_" + to_string(get<1>(new_repeat_loci[j])) + "_" +
                                to_string(get<2>(new_repeat_loci[j])) + "_" + get<3>(new_repeat_loci[j]);
            if (added_repeats.find(repeat_key) != added_repeats.end()) { continue; }
            added_repeats.insert(repeat_key);
            auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), new_repeat_loci[j], compareRepeatLoci);
            // Insert the new element at the found position
            repeat_loci.insert(it, new_repeat_loci[j]);
        }
        new_repeat_loci.clear();
    }

    bool addRepeat = checkDuplicateEntry(sequence_id, repeat_start, repeat_end, motif, repeat_loci);
    if (!addRepeat) return;

    vector<int> remove_loci;
    addRepeat = handleNestedParentRelation(repeat_start, repeat_end, motif, purity, cigar_string,
                                           motif_length, repeat_length, repeat_units, repeat_loci, remove_loci);
    if (!addRepeat) return;
    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1); }
    remove_loci.clear();

    addRepeat = handleOverlapRelation(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                      motif_length, repeat_length, repeat_units, out, repeat_loci, remove_loci,
                                      recursion_level, new_repeat_loci);
    if (!addRepeat) return;
    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1); }
    remove_loci.clear();

    addRepeat = handleNestedParentRelation(repeat_start, repeat_end, motif, purity, cigar_string,
                                           motif_length, repeat_length, repeat_units, repeat_loci, remove_loci);
    if (!addRepeat) return;
    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1); }

    int end_breakpoint = 10000; // 10kb
    for (int j=repeat_loci.size()-1; j>=0; j--) {
        // if the repeat is 10kb away from the last repeat
        // print the repeats to the output
        if (repeat_start - get<2> (repeat_loci[j]) > end_breakpoint) {
            // if the output stream is not null, print the repeats to the output
            // output stream is null in the case of the python module
            if (out) { printRepeatsToOutput(out, repeat_loci, j); }
            break;
        }
    }

    tuple<string, int, int, string, double, string, int, int, int> repeat_locus = { sequence_id, repeat_start, repeat_end,
                                                                                    motif, purity, cigar_string,
                                                                                    motif_length, repeat_length, repeat_units };
    auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
    // Insert the new element at the found position
    repeat_loci.insert(it, repeat_locus);
}
