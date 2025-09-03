#include "binomial_thresholds.h"


long double combination(int n, int r) {
    /*
     * Returns the binomial coefficient C(n, r) = n! / (r! * (n - r)!)
     * @param n Total number of trials.
     * @param r Number of successes.
     * @return long double Binomial coefficient C(n, r).
    */

    if (r > n) return 0;
    if (r == 0 || r == n) return 1;

    long double result = 1;
    for (int i = 1; i <= r; ++i) {
        result *= (n - (r - i));
        result /= i;
    }
    return result;
}

long double probWithRunApprox(int n, int x, int r, long double p) {
    /*
     * Returns the probability of having at least one run of length r given x successes in n trials.
     * @param n total number of trials
     * @param x number of successes
     * @param r length of the run
     * @param p probability of success in each trial
     * @return long double approximate probability of having at least one run of length r
    */

    long double binom_probability = combination(n, x) * pow(p, x) * pow(1 - p, n - x);

    // Approximate conditional probability of having ≥1 run of length r given x successes
    long double run_probability = 1.0 - expl(-(n - r + 1) * pow((long double)x / n, r));

    // Final approximation
    return binom_probability * run_probability;
}


int minimumNumberOfSuccesses(int n, int r, long double p) {
    /*
     * Returns the minimum number of successes required to exceed a given probability threshold.
     * @param n total number of trials
     * @param r continuous run of successes
     * @param p probability of success in each trial
     * @return minimum number of successes required
    */

    /*
     * This function calculates the threshold number of 1s in an anchor seed
     * n i.e., the total number of trials is the length of the anchor seed
     * r i.e., the continuous run of successes is the required minimum continuous 
     *         matches in the anchor seed which is also the length of minimum number of 1s
     *         when considering anchors from the neighboring shifts
     * p i.e., the probability of success in a trial is the purity threshold of the repeat
    */

    // Check if n is in THRESHOLD_BITS and return its value if found
    if (THRESHOLD_BITS.find(n) != THRESHOLD_BITS.end()) {
        // the threshold number of 1s for an anchor seed length is stored as an unordered map
        // this threshold only depends on the length of the anchor seed and agnostic to motif length
        return THRESHOLD_BITS.at(n);
    }

    // For the total number of trials n, we calculated the probability for x number of successes with
    // x ranging from 0 to n with at least one run of continuous successes of length r
    vector<long double> probabilities;
    for (int x = n; x >= 0; x--) {
        probabilities.push_back(probWithRunApprox(n, x, r, p));
    }

    // Threshold number of successes is defined as the value x where cumulative probability from x to n
    // values is >= 0.98
    // This is analougous to 98% of the repeat sequences of purity p will have an anchor seed of length n
    // with at least x number of 1s
    for (int x = 0; x < probabilities.size(); x++) {
        long double sumProb = 0.0;
        for (int j = 0; j <= x; j++) { sumProb += probabilities[j]; }
        if (sumProb >= 0.98) { THRESHOLD_BITS[n] = n-x; return n-x; }
    }
}


// Function to calculate cumulative probability
long double cumulativeBinomialProbability(int n, int k, long double p) {
    /*
     * Returns the cumulative probability of getting at least k successes in n trials.
     * @param n total number of trials
     * @param k minimum number of successes
     * @param p probability of success in each trial
     * @return long double cumulative probability of getting at least k successes
    */
    long double cumulativeProb = 0.0;
    for (int i = k; i <= n; ++i) {
        cumulativeProb += combination(n, i) * pow(p, i) * pow(1 - p, n - i);
    }
    return cumulativeProb;
}


unordered_map<int, pair<int, int>> getWindowThresholds(int minimum_mlen, int maximum_mlen) {
    /*
     * Returns a map of motif lengths to their respective window thresholds.
     * @param minimum_mlen: Minimum motif length.
     * @param maximum_mlen: Maximum motif length.
     * @return unordered_map<int, int>: Map of motif lengths to window thresholds.
    */
    unordered_map<int, pair<int, int>> thresholds;
    int window_length;
    for (int m = minimum_mlen; m <= maximum_mlen; m++) {
        if (m <= 16) { window_length = 8; }
        else { window_length = m/2; }
        for (int i = window_length; i >= 0; i--) {
            if (cumulativeBinomialProbability(window_length, i, 0.85) >= 0.9) {
                // If m is divisible by i, set the threshold
                thresholds[m] = std::make_pair(window_length, i);
                break;
            }
        }
    }

    return thresholds;
}
