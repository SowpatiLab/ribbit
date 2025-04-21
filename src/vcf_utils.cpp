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
              int &anchor_length, int &cones_threshold, string &out_file) {

    // assigns the output to either a file or standard output
    streambuf * buf; ofstream outstream;

    // if the output file is not given by default: input file + ".ribbit"
    if (out_file == "") { out_file = input_file + ".ribbit"; }
    outstream.open(out_file);
    buf = outstream.rdbuf();    // output file buffer is created
    ostream out(buf);

    // char *vcf_file = input_file.c_str();
    // Open the VCF file
    htsFile *fp = bcf_open(input_file.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Cannot open VCF file: " << input_file << std::endl;
        return;
    }

    // Read the header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: Cannot read VCF header" << std::endl;
        bcf_close(fp);
        return;
    }

    // Read records
    bcf1_t *rec = bcf_init();
    if (!rec) {
        std::cerr << "Error: Cannot initialize VCF record" << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        return;
    }
    
    int nrecords = 0;
    // bcf_unpack(rec, BCF_UN_FLT);
    while (bcf_read(fp, hdr, rec) >= 0) {

        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int64_t position = rec->pos;
        
        bcf_unpack(rec, BCF_UN_ALL); // bcf_read does not population d in bcf_record; this function populates d in bcf_record

        string allele = "", record_id = "";
        
        string motifs = "";
        char* af = nullptr;
        int naf = 0;
        if (bcf_get_info_string(hdr, rec, "MOTIFS", &af, &naf) > 0) {
            for (int i = 0; i < naf - 1; ++i) motifs += af[i];
            free(af); // Free the memory allocated by `bcf_get_info_float`
        }
        else {
            // std::cerr << "INFO: AF not found.\n";
        }
        
        if (nrecords > 0) out << "\n";
        
        for (int i=0; i<rec->n_allele; i++) {
            record_id = chrom;
            
            if (i==0) { record_id += ":" + to_string(position) + ":Ref"; }
            else { record_id += ":" + to_string(position) + ":Alt-" + to_string(i); }
            allele = rec->d.allele[i];
            
            out << ">" << record_id << ";ALLELE_LENGTH=" << allele.length() << ";MOTIFS=" << motifs << "\n";
            out << allele << "\n";
            SEQUENCE_ID = record_id; SEQUENCE = allele;
            processSequence(record_id, allele, window_length, window_bitcount_threshold, anchor_length, cones_threshold, out);
        }

        nrecords += 1;
    }

    // Clean up
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    outstream.close();

    return;
}