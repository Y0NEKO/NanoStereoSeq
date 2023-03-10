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
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>

using namespace std;
using namespace seqan3;
vector<string> anchor{ "ATGGCGACCTTATCAG", "TTGTCTTCCTAAGACCG" };
vector<int> mywidth{ 20, 30 , 25};


struct ReadInfo {
    string id;
    string seq;
};

struct SearchInfo {
    string id;
    string seq;
    vector <uint32_t> cidPrimerPos;
    vector <uint32_t> captureOligoPos;
    int isPlus;
    string finalCID;
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

std::string CidExtract(std::string qseq, double errorate) {
    string cid = "";
    
    int sublength = 200;
    if (qseq.size() < 200) {
        sublength = qseq.size();
    }
    string seq = qseq.substr(0, sublength);
    string seq_rc = reverseComplement(qseq).substr(0, sublength);
    vector <string> seqvec = { seq, seq_rc };
    
    seqan3::configuration const cfg1 = seqan3::search_cfg::hit_single_best{} | seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_rate{errorate} };
    
    //fm search
    bi_fm_index fm{ seqvec };
    auto searchResults = search(anchor, fm, cfg1);

    uint32_t cidPrimerPos = -1, captureOligoPos = -1, cidPrimerPosp = -1, captureOligoPosp = -1, cidPrimerPosm = -1, captureOligoPosm = -1, start = -1, width = 0;
    string trueQseq = seq;
    int hitcount = 0;
    int isp = 0, ism = 0;

    for (const auto& trueResult : searchResults) {
        hitcount++;
        auto queryId{ trueResult.query_id() };
        auto readId{ trueResult.reference_id() };
        auto matchPos{ trueResult.reference_begin_position() };
        cout << queryId << "\t" << matchPos << endl;
        if (readId == 0) {
            isp = 1;
            if (queryId == 0) {
                cidPrimerPosp = matchPos;
            }
            else if (queryId == 1) {
                captureOligoPosp = matchPos;
            }
        }
        else {
            ism = 1;
            if (queryId == 0) {
                cidPrimerPosm = matchPos;
            }
            else if (queryId == 1) {
                captureOligoPosm = matchPos;
            }
        }   
    }

    if (hitcount == 0 || (isp == 1 && ism == 1)) {
        return(cid);
    }

    if (isp == 1) {
        cidPrimerPos = cidPrimerPosp;
        captureOligoPos = captureOligoPosp;
    }
    else if (ism == 1) {
        cidPrimerPos = cidPrimerPosm;
        captureOligoPos = captureOligoPosm;
        trueQseq = seq_rc;
    }

    if (cidPrimerPos != -1 && captureOligoPos != -1 && captureOligoPos > cidPrimerPos) {
        start = cidPrimerPos + anchor[0].size();
        width = captureOligoPos - start;
        if (width > 35) {
            width = 25;
        }
    }
    else if (cidPrimerPos != -1) {
        start = cidPrimerPos + anchor[0].size();
        width = 25;
    }
    else if (captureOligoPos != -1) {
        start = captureOligoPos - anchor[1].size();
        width = 25;
    }
    if (start != -1 && start < 180) {
        cout << start << "\t" << width << endl;
        cid = trueQseq.substr(start, width);
    }
    return(cid);
}

void findCIDSequence(const vector<ReadInfo>& reads, const string& outputFileName, uint32_t thread, double errorate {
    ofstream outputFile(outputFileName);
    // Build FM-index of reads
    //seqan3::configuration const cfg1 = seqan3::search_cfg::hit_all_best{} | seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_count{2} };
    //trans query to dna5
    vector <string> cidvector;
    for (auto& read : reads)
    {
        cout << read.id << ":" << endl;
        string cid = CidExtract(read.seq, errorate);
        if (cid != "") {
            outputFile << read.id << '\t' <<  cid << "\t"  << read.seq << endl;
        }
    }
    
    outputFile.close();
}


int main(int argc, char* argv[]) {
    string bamFileA = argv[1];
    string outprex = argv[2];
    uint32_t thread = atoi(argv[3]);
    double errorate = atof(argv[4]);
    //part of cidprimer, capture oligo
    
    vector<ReadInfo> reads;
    readBam(bamFileA, reads);
    findCIDSequence(reads,  outprex, thread, errorcount);

    return 0;
}