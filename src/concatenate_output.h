#include "global_variables.h"
#include "output_utils.h"

using namespace std;

void concatenateOutputs(string out_file, vector<string>seq_names, int THREADS);

void concatenateThreadOutputs(const vector<string> &temp_files, ofstream &out);
