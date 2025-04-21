#include <iostream>
#include <fstream>
#include <stdio.h>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <string>
#include <vector>

#include "global_variables.h"
#include "fasta_utils.h"

using namespace std; 

void parseVCF(string &input_file, int &window_length, int &window_bitcount_threshold,
              int &anchor_length, int &cones_threshold, string &out_file);