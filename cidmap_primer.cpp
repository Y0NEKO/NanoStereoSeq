#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <ranges>
#include <htslib/sam.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/max_error.hpp>

using namespace std;
using namespace seqan3;

struct ReadInfo {
    string id;
    string seq;
};
std::string reverseComplement(std::string mStr) {
    std::string str(mStr.length(), 0);
    for (int c = 0; c < mStr.length(); c++) {
        char base = mStr[c];
        switch (base) {
        case 'A':
        case 'a':
            str[mStr.length() - c - 1] = 'T';
            break;
        case 'T':
        case 't':
            str[mStr.length() - c - 1] = 'A';
            break;
        case 'C':
        case 'c':
            str[mStr.length() - c - 1] = 'G';
            break;
        case 'G':
        case 'g':
            str[mStr.length() - c - 1] = 'C';
            break;
        default:
            str[mStr.length() - c - 1] = 'N';
        }
    }
    return str;
}

void readBam(const string& bamFile, vector<ReadInfo>& reads) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    bam1_t* record = bam_init1();

    while (sam_read1(fp, header, record) >= 0) {
        if (record->core.flag & BAM_FUNMAP) {
            continue;
        }

        ReadInfo read;
        read.id = bam_get_qname(record);

        uint8_t* seq = bam_get_seq(record);
        std::string seq_str;

        for (int i = 0; i < record->core.l_qseq; ++i) {
            seq_str.push_back(seq_nt16_str[bam_seqi(seq, i)]);
        }
        read.seq = seq_str;
        reads.emplace_back(read);
    }

    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);
}

void findShortSequence(const vector<ReadInfo>& reads, const vector<string>& targets, const string& outputFileName, uint32_t thread, double errorate) {
    ofstream outputFile(outputFileName);
    // Build FM-index of reads
    cout << "read convert..." << endl;
    std::vector<vector<dna5>> seqs;
    std::vector<string> seqstrs, seqstrs_rc;
    for (auto& read : reads)
    {
        int sublength = 150;
        if (read.seq.size() < 150) {
            sublength = read.seq.size();
        }
        string seq = read.seq.substr(0, sublength);
        string seq_rc = reverseComplement(read.seq).substr(0, sublength);
        seqstrs.push_back(seq);
        seqstrs_rc.push_back(seq_rc);
        
    }

    // Build FM index
    cout << "FM index build..." << endl;
    fm_index fm{ seqstrs };
    fm_index fm_rc{ seqstrs_rc };

    // Search for target sequence matches
    seqan3::configuration const cfg1 = seqan3::search_cfg::parallel{ thread } | seqan3::search_cfg::hit_all{} |
        seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_rate{errorate} };
    //seqan3::configuration const cfg1 = seqan3::search_cfg::hit_all_best{} | seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_count{2} };

    cout << "search seq..." << endl;
    auto searchResults = search(targets, fm, cfg1);
    auto searchResults_rc = search(targets, fm_rc, cfg1);
    //auto searchResults = search(targets_dna5, fm, cfg1);
    //auto searchResults_rc = search(targets_dna5, fm, cfg1);

    // Write results to output file
    for (const auto& searchResult : searchResults) {
        auto queryId{ searchResult.query_id() };
        auto readId{ searchResult.reference_id() };
        auto matchPos{ searchResult.reference_begin_position() };
        outputFile << queryId << '\t' << reads[readId].id << '\t' << matchPos << '\t' << "+" <<"\t" << reads[readId].seq <<endl;
    }

    for (const auto& searchResult : searchResults_rc) {
        auto queryId{ searchResult.query_id() };
        auto readId{ searchResult.reference_id() };
        auto matchPos{ searchResult.reference_begin_position() };
        outputFile << queryId << '\t' << reads[readId].id << '\t' << matchPos << '\t' << "-" << "\t" << reads[readId].seq << endl;
    }

    outputFile.close();
}


int main(int argc, char* argv[]) {
    string bamFileA = argv[1];
    string outprex = argv[2];
    uint32_t thread = atoi(argv[3]);
    double errorate = atof(argv[4]);
    vector<string> targets{ "ATGGCGACCTTATCAG", "CTGCTGACGTACTGAGAGGC", "TTGTCTTCCTAAGACCG"};
    vector<ReadInfo> reads;
    cout << "load bam file" << endl;
    readBam(bamFileA, reads);
    cout << "Finding seq" << endl;
    findShortSequence(reads, targets, outprex, thread, errorate);

    return 0;
}