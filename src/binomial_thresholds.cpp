#include "binomial_thresholds.h"
#include <boost/math/distributions/binomial.hpp>


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

    // cout << n << "\t" << x << "\t" << nCr(n,x) << "\t" << combination(n, x) << "\n";
    // long double binom_probability = combination(n, x) * pow(p, x) * pow(1 - p, n - x);
    boost::math::binomial_distribution<long double> dist(n, p);
    long double binom_probability = boost::math::pdf(dist, x);

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

        
    // For the total number of trials n, we calculated the probability for x number of successes with
    // x ranging from 0 to n with at least one run of continuous successes of length r
    long double sumProb = 0.0;
    for (int x = 0; x <= n; x++) {
        sumProb += probWithRunApprox(n, x, r, p);
        // if (sumProb >= 0.1) { threshold_bits[n] = x; return x; }
        if (sumProb >= 0.02) { return x; }
    }

    // Threshold number of successes is defined as the value x where cumulative probability from x to n
    // values is >= 0.98
    // This is analougous to 98% of the repeat sequences of purity p will have an anchor seed of length n
    // with at least x number of 1s
}


long double probabilityOfSuccesses(int n, int r, long double p, int successes) {
    /*
     *  Returns the probability of number of matches found in an anchor seed
     *  @param n total number of trials (length of the anchor seed)
     *  @param p probability of success in each trial (purity threshold of the repeat)
     *  @param r continuous run of successes (length of the longest continuous matches)
     *  @param successes the total number of success (total number of matches in the anchor seed)
     *  @return minimum number of successes required
     */

        
    // For the total number of trials n, we calculated the probability for x number of successes with
    // x ranging from 0 to n with at least one run of continuous successes of length r
    return probWithRunApprox(n, successes, r, p);
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


void calculateWindowThresholds() {
    /*
     *  calculates the window lengths and thresholds for all motif lengths
     *  @param none
     *  @return void updates the global WINDOW_LENGTHS and WINDOW_THRESHOLDS variables
    */

    int window_length;
    float fraction_motif = 0.7;
    int motif_cutoff = 8 / fraction_motif;
    double cumm_prob_threshold = 0.9, cumm_prob = 0.0;
    for (int mlen = MINIMUM_MLEN; mlen <= MAXIMUM_MLEN; mlen++) {
        window_length = (mlen <= motif_cutoff) ? 8 : mlen*fraction_motif;
        if (window_length < 8) window_length = 8;
        if (window_length > motif_cutoff) { cumm_prob_threshold = 0.9; }
        else { cumm_prob_threshold = 0.9; }
        
        WINDOW_LENGTHS[mlen - MINIMUM_MLEN] = window_length;
        for (int i = window_length; i >= 0; i--) {
            cumm_prob = cumulativeBinomialProbability(window_length, i, PURITY_THRESHOLD);
            // cumm_prob = cumulativeBinomialProbability(window_length, i, 0.85);
            if (cumm_prob >= cumm_prob_threshold) {
                // If m is divisible by i, set the threshold
                WINDOW_THRESHOLDS[mlen - MINIMUM_MLEN] = i;
                break;
            }
        }
    }
}
