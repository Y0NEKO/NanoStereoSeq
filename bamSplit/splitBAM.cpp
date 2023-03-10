#include <iostream>
#include <fstream>
#include <vector>
#include "htslib/sam.h"

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " input.bam ref.txt output_prefix\n";
        return 1;
    }

    std::string input_filename(argv[1]);
    std::string ref_filename(argv[2]);
    std::string output_prefix(argv[3]);

    // Open the input BAM file using htslib
    samFile *input_fp = sam_open(input_filename.c_str(), "r");
    if (!input_fp) {
        std::cerr << "Failed to open input BAM file: " << input_filename << "\n";
        return 1;
    }

    // Open the output BAM files using htslib
    std::vector<samFile*> output_fps;
    std::ifstream ref_file(ref_filename);
    std::string line;
    while (std::getline(ref_file, line)) {
        std::string output_filename = output_prefix + "." + line + ".bam";
        samFile *output_fp = sam_open(output_filename.c_str(), "wb");
        if (!output_fp) {
            std::cerr << "Failed to open output BAM file: " << output_filename << "\n";
            return 1;
        }
        output_fps.push_back(output_fp);
    }

    // Read the input BAM file and write the records to the appropriate output BAM files
    bam1_t *record = bam_init1();
    while (sam_read1(input_fp, bam_hdr_read(input_fp->header), record) >= 0) {
        std::string read_id = bam_get_qname(record);
        size_t dot_pos = read_id.find('.');
        if (dot_pos != std::string::npos) {
            read_id = read_id.substr(0, dot_pos);
        }
        int output_index = -1;
        ref_file.clear();
        ref_file.seekg(0);
        while (std::getline(ref_file, line)) {
            if (line == read_id) {
                output_index = ref_file.tellg();
                break;
            }
        }
        if (output_index < 0 || output_index >= output_fps.size()) {
            std::cerr << "Failed to find output BAM file for read ID: " << read_id << "\n";
            return 1;
        }
        sam_write1(output_fps[output_index], bam_hdr_read(input_fp->header), record);
    }

    // Clean up
    bam_destroy1(record);
    for (samFile *output_fp : output_fps) {
        sam_close(output_fp);
    }
    sam_close(input_fp);

    return 0;
}
