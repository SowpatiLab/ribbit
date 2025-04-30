#include "output_utils.h"

using namespace std;
using namespace boost;


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
            << get<3> (repeat_loci[i]) << "\t" << get<4> (repeat_loci[i]) << "\t+\t" << get<6> (repeat_loci[i]) << "\t" 
            << get<7> (repeat_loci[i]) << "\t" << get<8> (repeat_loci[i]);
        if (CIGAROUTPUT) { out << "\t" << get<5> repeat_loci[i]; }
        out << "\n";
    }

    repeat_loci.erase(repeat_loci.begin(), repeat_loci.begin() + end_index + 1);
}


int absolute(int x) {
    /*
     * returns the absolute value of an integer
     * @param x integer
     * @return int absolute value of the integer
    */
    if (x < 0) return -x;
    return x;
}


bool checkCyclicalVariation(string query, string ref) {
    /*
     * checks if a query motif is cyclical variation of reference motif
     * @param query sequence of the query motif
     * @param ref sequence of reference motif
     * @return bool if the query motif is a cyclical variation of the reference motif
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
     * calculates the least distance between two strings
     * @param reference string to be compared with
     * @param query string to be compared
     * @return int the least distance between the two strings
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
     * extracts the non-overlapping cigar from the downstream locus of two overlapping loci
     * @param a_end     end position of the upstream locus
     * @param b_start   start position of the downstream locus
     * @param dn_clens   the lengths of the cigar operations for the downstream locus
     * @param dn_ctypes  the consecutive cigar operations of the downstream locus
     * @return tuple<vector<int>, vector<char>> the lengths of continuous cigar operations and continuous cigar operations
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


void alignSequenceWithPerfectRepeat(string &sequence, string &motif, int &repeat_start, int &repeat_end,
                                    double &purity, string &new_cigar) {
    /*
     * aligns the sequence with a perfect repeat of the motif
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
     * checks if the cigar operations are only matches
     * @param ctypes vector of cigar operations
     * @return bool if the cigar operations are only matches
    */
    if (ctypes.size() == 0 || ctypes.size() > 1) return false;
    if (ctypes[0] == '=' || ctypes[0] == 'M') return true;
    return true;
}


void compareOverlappingLoci(int &upstart, int &upend, string &upcigar, string &upmotif, double &uppurity, bool &up_drop, bool &up_update,
                            int &dnstart, int &dnend, string &dncigar, string &dnmotif, double &dnpurity, bool &dn_drop, bool &dn_update) {
    /*
     * compares the overlapping cigar of two overlapping loci and merges them
     * @param upstart start position of the upstream locus
     * @param upend end position of the upstream locus
     * @param upcigar CIGAR of the upstream locus
     * @param dnstart start position of the downstream locus
     * @param dnend end position of the downstream locus
     * @param dncigar CIGAR of the downstream locus
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


tuple<string, double> mergeRepeatsIdenticalMotif(int upend, int dnstart, string upcigar, string dncigar) {
    /*
     * merges the two overlapping loci with identical motifs
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


// Define a comparison function for tuples (e.g., comparing the first element)
bool compareRepeatLoci(const tuple<string, int, int, string, double, string, int, int, int> &a,
                       const tuple<string, int, int, string, double, string, int, int, int> &b) {
    /*
     * compares repeat loci based on the start positions
     *  @param a tuple of the first repeat location
     *  @param b tuple of the second repeat location
     *  @returns bool bool value indicating if the first repeat is before second repeat
    */
    return get<1> (a) < get<1> (b); // Compare based on the first element (int)
}


void addLocusToOutput(string &sequence_id, int repeat_start, int repeat_end, string motif, double purity, string &cigar_string,
                      int motif_length, int repeat_length, int repeat_units, ostream &out,
                      vector<tuple<string, int, int, string, double, string, int, int, int>> &repeat_loci) {
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

    string last_seqid = "";
    if (repeat_loci.size() > 0) { last_seqid = get<0> (repeat_loci[repeat_loci.size()-1]); }

    if (last_seqid != "" && sequence_id != last_seqid) {
        // if the last sequence id is different
        printRepeatsToOutput(out, repeat_loci, repeat_loci.size()-1);
        tuple<string, int, int, string, double, string, int, int, int> repeat_locus;
        repeat_locus = { sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                         motif_length, repeat_length, repeat_units };
        auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
        // Insert the new element at the found position
        repeat_loci.insert(it, repeat_locus);
        return;
    }

    int last_start, last_end, last_length, last_mlen, last_units;
    int overlap_length; string overlap_seq;
    string last_motif, last_cigar; double last_purity;
    vector<int> remove_loci;
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

        last_motif  = get<3> (repeat_loci[i]);
        last_mlen   = get<6> (repeat_loci[i]);
        last_purity = get<4> (repeat_loci[i]);
        last_cigar  = get<5> (repeat_loci[i]);
        last_length = get<7> (repeat_loci[i]);

        // identical coordinates
        if (repeat_start == last_start && repeat_end == last_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) { return; }
            else if (purity <= last_purity) return;
            else if (purity > last_purity)  remove_loci.push_back(i);
        }

        // nested
        else if (last_start <= repeat_start && repeat_end <= last_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) {
                return;
            }

            else {
                tuple<vector<int>, vector<char>> segment_cigarvalues = extractRegionCigar(last_cigar, repeat_start-last_start, repeat_end-last_start);
                double segment_purity = ((double) getMatches(segment_cigarvalues)) / ((double) getAlignmentLength(segment_cigarvalues));
                if (purity < segment_purity) return;

                if (last_mlen >= 2*motif_length) {
                    double mcomp_purity = 0.0; int mcomp_start = 0, mcomp_end = 2*last_mlen; string mcomp_cigar = "";
                    string full_sequence = last_motif + last_motif;
                    alignSequenceWithPerfectRepeat(full_sequence, motif, mcomp_start, mcomp_end, mcomp_purity, mcomp_cigar);
                    if (mcomp_purity > 0.9 && mcomp_end-mcomp_start > last_mlen) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        repeat_start = last_start;
                        repeat_end = last_end;
                        string full_sequence = SEQUENCE.substr(repeat_start, repeat_end-repeat_start);
                        alignSequenceWithPerfectRepeat(full_sequence, motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                         motif_length, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                }

                if (segment_purity > PURITY_THRESHOLD) {
                    // keep both repeats as it is
                }

                else {
                    int upflank_length = repeat_start - last_start;
                    int downflank_length = last_end - repeat_end;
                    if (upflank_length < 2*last_mlen && downflank_length < 2*last_mlen) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        repeat_start = last_start;
                        repeat_end = last_end;
                        string full_sequence = SEQUENCE.substr(last_start, last_end-last_start);
                        alignSequenceWithPerfectRepeat(full_sequence, motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                         motif_length, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                    else if (upflank_length < 2*last_mlen) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        tuple<vector<int>, vector<char>> dnflank_cigarvalues = extractUpCigar(last_cigar, last_end, repeat_end);
                        string dnflank_cigar = buildCigar(dnflank_cigarvalues);
                        double dnflank_purity = ((double) getMatches(dnflank_cigar)) / ((double) getAlignmentLength(dnflank_cigar));
                        addLocusToOutput(sequence_id, repeat_end, last_end, last_motif, dnflank_purity, dnflank_cigar,
                                         last_mlen, last_end-repeat_end, (last_end-repeat_end)/last_mlen, out, repeat_loci);
                        repeat_start = last_start;
                        string full_sequence = SEQUENCE.substr(last_start, repeat_end-last_start);
                        alignSequenceWithPerfectRepeat(full_sequence, motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                         motif_length, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                    else if (downflank_length < 2*last_mlen) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        tuple<vector<int>, vector<char>> upflank_cigarvalues = extractDownCigar(last_cigar, last_start, repeat_start);
                        string upflank_cigar = buildCigar(upflank_cigarvalues);
                        double upflank_purity = ((double) getMatches(upflank_cigar)) / ((double) getAlignmentLength(upflank_cigar));
                        addLocusToOutput(sequence_id, last_start, repeat_start, last_motif, upflank_purity, upflank_cigar,
                                         last_mlen, repeat_start-last_start, (repeat_start-last_start)/last_mlen, out, repeat_loci);
                        repeat_end = last_end;
                        string full_sequence = SEQUENCE.substr(repeat_start, last_end-repeat_start);
                        alignSequenceWithPerfectRepeat(full_sequence, motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                         motif_length, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                }
            }
            // if (purity <= last_purity || motif.length() > last_mlen || motif == last_motif || checkCyclicalVariation(motif, last_motif)
            //     || (motif.length() > SMALL_MLEN_LIMIT && last_mlen > SMALL_MLEN_LIMIT && absolute(motif.length() - last_mlen) <= 2)) {
            //     return;
            // }
        }

        // parent
        else if (repeat_start <= last_start && last_end <= repeat_end) {
            if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) {
                remove_loci.push_back(i);
            }
            else {
                tuple<vector<int>, vector<char>> segment_cigarvalues = extractRegionCigar(cigar_string, last_start-repeat_start, last_end-repeat_start);
                double segment_purity = ((double) getMatches(segment_cigarvalues)) / ((double) getAlignmentLength(segment_cigarvalues));

                if (last_purity < segment_purity) remove_loci.push_back(i);

                else {
                    // if the unique region of the parent is more than twice of the motif length, the parent repeat is retained
                    int upflank_length = last_start - repeat_start;
                    int downflank_length = repeat_end - last_end;
                    if (upflank_length < 2*motif_length && downflank_length < 2*motif_length) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        string full_sequence = SEQUENCE.substr(repeat_start, repeat_length);
                        alignSequenceWithPerfectRepeat(full_sequence, last_motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / last_mlen;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, last_motif, purity, cigar_string,
                                         last_mlen, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                    else if (upflank_length < 2*motif_length) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        tuple<vector<int>, vector<char>> dnflank_cigarvalues = extractUpCigar(cigar_string, repeat_end, last_end);
                        string dnflank_cigar = buildCigar(dnflank_cigarvalues);
                        double dnflank_purity = ((double) getMatches(dnflank_cigar)) / ((double) getAlignmentLength(dnflank_cigar));
                        repeat_length = repeat_end - last_end; repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, last_end, repeat_end, motif, dnflank_purity, dnflank_cigar,
                                         motif_length, repeat_length, repeat_units, out, repeat_loci);
                        repeat_end = last_end;
                        string full_sequence = SEQUENCE.substr(repeat_start, last_end-repeat_start);
                        alignSequenceWithPerfectRepeat(full_sequence, last_motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / last_mlen;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, last_motif, purity, cigar_string,
                                         last_mlen, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                    else if (downflank_length < 2*motif_length) {
                        remove_loci.push_back(i);
                        for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                        tuple<vector<int>, vector<char>> upflank_cigarvalues = extractDownCigar(cigar_string, repeat_start, last_start);
                        string upflank_cigar = buildCigar(upflank_cigarvalues);
                        double upflank_purity = ((double) getMatches(upflank_cigar)) / ((double) getAlignmentLength(upflank_cigar));
                        addLocusToOutput(sequence_id, repeat_start, last_start, motif, upflank_purity, upflank_cigar,
                                         motif_length, last_start-repeat_start, (last_start-repeat_start)/motif_length, out, repeat_loci);
                        repeat_start = last_start;
                        string full_sequence = SEQUENCE.substr(repeat_start, repeat_end-last_start);
                        alignSequenceWithPerfectRepeat(full_sequence, last_motif, repeat_start, repeat_end, purity, cigar_string);
                        repeat_length = repeat_end - repeat_start; repeat_units = repeat_length / last_mlen;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, last_motif, purity, cigar_string,
                                         last_mlen, repeat_length, repeat_units, out, repeat_loci);
                        return;
                    }
                }
            }
            // if (last_purity <= purity || last_mlen > motif.length() || motif == last_motif || checkCyclicalVariation(motif, last_motif)
            //     || (motif.length() > SMALL_MLEN_LIMIT && last_mlen > SMALL_MLEN_LIMIT && absolute(motif.length() - last_mlen) <= 2)) {
            //     remove_loci.push_back(i);
            // }
        }

        // if not the above options the repeats should be overlapping
        else if (motif == last_motif || checkCyclicalVariation(motif, last_motif)) {
            // if the repeats are of the same motif
            tuple <string, double> merge_values;

            // if the repeats are having the same motif then we merge them
            if (repeat_start == last_end || (last_start < repeat_start && repeat_start < last_end)) {
                merge_values = mergeRepeatsIdenticalMotif(last_end, repeat_start, last_cigar, cigar_string);
                repeat_start = last_start;
            }
            
            else if (last_start == repeat_end || (repeat_start < last_start && last_start < repeat_end)) {
                merge_values = mergeRepeatsIdenticalMotif(repeat_end, last_start, cigar_string, last_cigar);
                repeat_end = last_end;
            }

            repeat_length = repeat_end - repeat_start;
            cigar_string = get<0> (merge_values); purity = get<1> (merge_values);
            remove_loci.push_back(i);
            for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
            assert(repeat_length == getRepeatLength(cigar_string));
            addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                motif_length, repeat_length, repeat_units, out, repeat_loci); return;
        }

        else {

            int levenshtein_distance = 0, distance_threshold = 0;
            if ((motif.length() > 2*SMALL_MLEN_LIMIT && last_mlen > 2*SMALL_MLEN_LIMIT && absolute(motif.length() - last_mlen) <= 2)) {
                if (motif.length() > last_mlen) {
                    string reference = motif + motif;
                    levenshtein_distance = leastDistance(reference, last_motif);
                    distance_threshold = motif.length()/ 10;
                }
                else {
                    string reference = last_motif + last_motif;
                    levenshtein_distance = leastDistance(reference, motif);
                    distance_threshold = last_mlen/ 10;
                }
            }

            if (((motif.length() > 2*SMALL_MLEN_LIMIT && last_mlen > 2*SMALL_MLEN_LIMIT && absolute(motif.length() - last_mlen) <= 2))
                && levenshtein_distance <= (distance_threshold)) {

                string full_sequence;
                if (last_start < repeat_start) {     // last-repeat is upstream of current repeat
                    full_sequence = SEQUENCE.substr(last_start, repeat_end-last_start);
                    repeat_start = last_start;
                }
    
                else if (repeat_start < last_start) {  // last-repeat is downstream of current repeat
                    full_sequence = SEQUENCE.substr(repeat_start, last_end-repeat_start);
                    repeat_end = last_end;
                }

                double cmotif_purity = 0.0; int cmotif_start = repeat_start, cmotif_end = repeat_end; string cmotif_cigar = "";
                double pmotif_purity = 0.0; int pmotif_start = repeat_start, pmotif_end = repeat_end; string pmotif_cigar = "";
                alignSequenceWithPerfectRepeat(full_sequence, motif, cmotif_start, cmotif_end, cmotif_purity, cmotif_cigar);
                alignSequenceWithPerfectRepeat(full_sequence, last_motif, pmotif_start, pmotif_end, pmotif_purity, pmotif_cigar);
                if (cmotif_purity > pmotif_purity) {
                    repeat_start = cmotif_start; repeat_end = cmotif_end;
                    cigar_string = cmotif_cigar; purity = cmotif_purity;
                }
                else {
                    repeat_start = pmotif_start; repeat_end = pmotif_end;
                    motif = last_motif; motif_length = last_mlen;
                    cigar_string = pmotif_cigar; purity = pmotif_purity;
                }
                remove_loci.push_back(i);
                for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                repeat_length = repeat_end - repeat_start;
                repeat_units = repeat_length / motif_length;
                assert(repeat_length == getRepeatLength(cigar_string));
                addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string,
                                 motif_length, repeat_length, repeat_units, out, repeat_loci);
                return;
            } 

            else {

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
                        addLocusToOutput(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                         last_length, last_units, out, repeat_loci);
                    }
                    if (!fail) {
                        repeat_length = repeat_end - repeat_start;
                        repeat_units = repeat_length / motif_length;
                        addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length, repeat_length,
                                         repeat_units, out, repeat_loci);
                    }
                    return;
                }
                
                else if (update & !(last_update)) {
                    // if only current repeat is updated
                    if (last_fail) { remove_loci.push_back(i); }
                    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                    repeat_length = repeat_end - repeat_start;
                    repeat_units = repeat_length / motif_length;
                    addLocusToOutput(sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length, repeat_length,
                                     repeat_units, out, repeat_loci);
                    return;
                }
    
                else if (last_update & !(update)) {
                    remove_loci.push_back(i);
                    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
                    remove_loci.clear();
                    last_length = last_end - last_start;
                    last_units = last_length / last_mlen;
                    addLocusToOutput(sequence_id, last_start, last_end, last_motif, last_purity, last_cigar, last_mlen,
                                     last_length, last_units, out, repeat_loci);
                    if (fail) return;
                }
            }
        }
    }

    for (int j: remove_loci) { repeat_loci.erase(repeat_loci.begin() + j, repeat_loci.begin() + j + 1);}
    for (int j=repeat_loci.size()-1; j>=0; j--) {
        if (repeat_start - get<2> (repeat_loci[j]) > 50000) {
            // if the repeat is 50kb away from the last repeat
            // print the repeats to the output
            printRepeatsToOutput(out, repeat_loci, j);
            break;
        }
    }

    tuple<string, int, int, string, double, string, int, int, int> repeat_locus = {sequence_id, repeat_start, repeat_end, motif, purity, cigar_string, motif_length, repeat_length, repeat_units};
    auto it = lower_bound(repeat_loci.begin(), repeat_loci.end(), repeat_locus, compareRepeatLoci);
    // Insert the new element at the found position
    repeat_loci.insert(it, repeat_locus);
}
