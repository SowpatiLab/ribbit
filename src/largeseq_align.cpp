#include "largeseq_align.h"

using namespace std;

const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;
const int GAP_OPEN = -2;
const int GAP_EXTEND = -1;

enum State { ST_M = 0, ST_I = 1, ST_D = 2 };

// Utility to get match/mismatch score
inline int sub_score(char a, char b) {
    return (a == b) ? MATCH_SCORE : MISMATCH_SCORE;
}

// Needleman-Wunsch with affine gap penalties and full traceback -> CIGAR
string nw_affine_cigar(const string &s1, const string &s2, int &distance) {
    const int n = (int)s1.size();
    const int m = (int)s2.size();
    // Use a comfortably small NEG_INF
    const int NEG_INF = numeric_limits<int>::min() / 4;

    distance = 0;

    // DP matrices: M, I (insertion in s1 => gap in s1 => consumes s2), D (deletion in s1 => gap in s2)
    vector<vector<int>> M(n+1, vector<int>(m+1, NEG_INF));
    vector<vector<int>> I(n+1, vector<int>(m+1, NEG_INF));
    vector<vector<int>> D(n+1, vector<int>(m+1, NEG_INF));

    // Backpointers: for each cell and state, store which previous state we came from
    // We'll store as char values 'M','I','D'
    vector<vector<std::array<char,3>>> from(n+1, vector<std::array<char,3>>(m+1));

    // Initialization
    M[0][0] = 0;
    I[0][0] = D[0][0] = NEG_INF;
    from[0][0] = {0,0,0};

    for (int i = 1; i <= n; ++i) {
        // At (i,0): only deletions (gaps in s2) possible
        D[i][0] = GAP_OPEN + (i-1) * GAP_EXTEND;
        M[i][0] = NEG_INF;
        I[i][0] = NEG_INF;
        from[i][0][ST_D] = (i==1 ? 'M' : 'D'); // first D comes from M (open), then from D (extend)
    }
    for (int j = 1; j <= m; ++j) {
        // At (0,j): only insertions (gaps in s1) possible
        I[0][j] = GAP_OPEN + (j-1) * GAP_EXTEND;
        M[0][j] = NEG_INF;
        D[0][j] = NEG_INF;
        from[0][j][ST_I] = (j==1 ? 'M' : 'I'); // first I from M, then I from I
    }

    // Fill DP
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            // M(i,j) from max of M/I/D at (i-1,j-1) plus match/mismatch
            int m_from_M = M[i-1][j-1];
            int m_from_I = I[i-1][j-1];
            int m_from_D = D[i-1][j-1];
            int best_m_prev = m_from_M;
            char best_m_from = 'M';
            if (m_from_I > best_m_prev) { best_m_prev = m_from_I; best_m_from = 'I'; }
            if (m_from_D > best_m_prev) { best_m_prev = m_from_D; best_m_from = 'D'; }
            M[i][j] = best_m_prev + sub_score(s1[i-1], s2[j-1]);
            from[i][j][ST_M] = best_m_from;

            // I(i,j): insertion in s1 -> comes from (i, j-1)
            // either extend previous I, or open from M or D
            int i_extend = I[i][j-1] + GAP_EXTEND;
            int i_open_from_M = M[i][j-1] + GAP_OPEN;
            int i_open_from_D = D[i][j-1] + GAP_OPEN;
            int best_i_prev = i_extend;
            char best_i_from = 'I';
            if (i_open_from_M > best_i_prev) { best_i_prev = i_open_from_M; best_i_from = 'M'; }
            if (i_open_from_D > best_i_prev) { best_i_prev = i_open_from_D; best_i_from = 'D'; }
            I[i][j] = best_i_prev;
            from[i][j][ST_I] = best_i_from;

            // D(i,j): deletion in s1 -> comes from (i-1, j)
            int d_extend = D[i-1][j] + GAP_EXTEND;
            int d_open_from_M = M[i-1][j] + GAP_OPEN;
            int d_open_from_I = I[i-1][j] + GAP_OPEN;
            int best_d_prev = d_extend;
            char best_d_from = 'D';
            if (d_open_from_M > best_d_prev) { best_d_prev = d_open_from_M; best_d_from = 'M'; }
            if (d_open_from_I > best_d_prev) { best_d_prev = d_open_from_I; best_d_from = 'I'; }
            D[i][j] = best_d_prev;
            from[i][j][ST_D] = best_d_from;
        }
    }

    // Final state: choose best among M[n][m], I[n][m], D[n][m]
    int final_score = M[n][m];
    char state = 'M';
    if (I[n][m] > final_score) { final_score = I[n][m]; state = 'I'; }
    if (D[n][m] > final_score) { final_score = D[n][m]; state = 'D'; }

    // Traceback building using vector of (op, count) so we can reverse order safely
    vector<pair<char,int>> runs;
    int i = n, j = m;
    while (i > 0 || j > 0) {
        char op;
        if (state == 'M') {
            // M consumes both
            op = 'M'; // use 'M' for alignment (match or mismatch)
            // determine previous state (from[i][j][ST_M])
            if (s1[i-1] != s2[j-1]) {
                op = 'X'; // mismatch
                distance += 1;
            }
            char prev = from[i][j][ST_M];
            i--; j--;
            state = prev;
        } else if (state == 'I') {
            // I consumes j only (insertion to s1 -> characters in s2)
            op = 'I';
            distance += 1; // count insertion as mismatch
            char prev = from[i][j][ST_I];
            j--;
            state = prev;
        } else { // 'D'
            op = 'D';
            distance += 1; // count deletion as mismatch
            char prev = from[i][j][ST_D];
            i--;
            state = prev;
        }

        if (!runs.empty() && runs.back().first == op) {
            runs.back().second += 1;
        } else {
            runs.emplace_back(op, 1);
        }
    }

    // runs currently in reverse (from end->start). Reverse to produce left->right CIGAR
    reverse(runs.begin(), runs.end());

    // Format as CIGAR string like "5M2I3M1D"
    ostringstream out_cigar;
    for (auto &p : runs) {
        out_cigar << p.second << p.first;
    }

    return out_cigar.str();
}


string alignLargeSequence(string &sequence, string &motif, int motif_length, StripedSmithWaterman::Aligner &aligner,
                          StripedSmithWaterman::Filter &filter, StripedSmithWaterman::Alignment &alignment) {
    /*
     *  aligns a large sequence to a motif sequence using Smith-Waterman algorithm
     *  @param sequence the large sequence to be aligned
     *  @param motif_sequence the motif sequence to be aligned to
     *  @param aligner the Smith-Waterman aligner
     *  @param filter the Smith-Waterman filter
     *  @param alignment the Smith-Waterman alignment
     *  @return string the CIGAR string of the alignment
     */

    int alength = 2000;
    int ppr_length = alength + 2 * motif_length + (int)((1 - PURITY_THRESHOLD) * alength);
    string ppr_sequence = "";
    while(ppr_sequence.length() < ppr_length) { ppr_sequence += motif; }
    int astart = 0;
    string asequence = "", cigar = "", clipped_cigar = "";
    int last_motif_end = 0, new_motif_start = 0;
    int ref_gap = 0, itr = 0;
    int distance = 0;

    while(astart < sequence.length()) {
        alength = (alength < sequence.length() - astart) ? alength : (sequence.length() - astart);
        asequence = sequence.substr(astart, alength);

        aligner.Align(asequence.c_str(), ppr_sequence.c_str(), ppr_length, filter, &alignment, 15);

        if (itr > 0) {
            if (alignment.query_begin > 0 && alignment.ref_begin == 0) {
                if (alignment.query_begin < motif_length/2) {
                    cigar += std::to_string(alignment.query_begin) + "I";
                }
                else {
                    distance = 0;
                    // checking the inserted sequence against the whole motif
                    clipped_cigar = nw_affine_cigar(ppr_sequence.substr(0, motif_length), asequence.substr(0, alignment.query_begin), distance);
                    if (distance < alignment.query_begin) { cigar += clipped_cigar; }
                    else { cigar += std::to_string(alignment.query_begin) + "I"; }
                }
            }
            else if (alignment.ref_begin > 0 && alignment.query_begin == 0) {
                ref_gap = alignment.ref_begin % motif_length;
                cigar += std::to_string(ref_gap) + "D";
            }
            else if (alignment.ref_begin > 0 && alignment.query_begin > 0) {
                distance = 0;
                clipped_cigar = nw_affine_cigar(ppr_sequence.substr(0, alignment.ref_begin), asequence.substr(0, alignment.query_begin), distance);
                cigar += clipped_cigar;
            }
            last_motif_end  = alignment.ref_end % motif_length;
            ppr_sequence = motif.substr(last_motif_end) + ppr_sequence;
            motif = ppr_sequence.substr(0, motif_length);
        }

        if (itr == 0) { cigar += trimSoftClipsinCigar(alignment.cigar_string, "end"); }
        else if (astart + asequence.length() >= sequence.length()) { cigar += trimSoftClipsinCigar(alignment.cigar_string, "start"); }
        else { cigar += trimSoftClipsinCigar(alignment.cigar_string, "both"); }

        if (astart + asequence.length() >= sequence.length()) {
            // last iteration
            break;
        }

        astart += (alignment.query_end + 1);
        itr += 1;
    }

    return cigar;
}
