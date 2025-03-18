#include <boost/dynamic_bitset.hpp>

#include "global_variables.h"

using namespace std;
using namespace boost::multiprecision;

uint32_t calculateRepeatClass(boost::dynamic_bitset<> &window, int &motif_length);

int calculateAtomicity(uint256_t &motif, int &motif_length);

int calculateAtomicity(uint32_t &motif, int &motif_length);

int calculateAtomicityLongMotif(uint256_t &motif, int &motif_length);

string calculateMotif(uint256_t motif_unit, int motif_length);
