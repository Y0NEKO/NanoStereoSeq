#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <future>
#include <thread>
#include <map>
#include <sstream>
#include <inttypes.h>


using namespace std;
using namespace seqan3;
vector<string> anchor{ "ATGGCGACCTTATCAG", "TTGTCTTCCTAAGACCG" };
auto left_dna5 = anchor[0] | seqan3::views::char_to<seqan3::dna5>;
auto right_dna5 = anchor[1] | seqan3::views::char_to<seqan3::dna5>;
auto method = seqan3::align_cfg::method_local{};
seqan3::align_cfg::scoring_scheme scheme{ seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}} };
seqan3::align_cfg::gap_cost_affine gap_costs{ seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-1} };
auto aligncfg = method | scheme | gap_costs;

vector<int> mywidth{ 20, 30 , 25};
struct ReadInfo {
    string id;
    string seq;
    string gene;
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

void loadGeneRanges(string gtfFile, map <string, map<string, pair<uint64_t, uint64_t>>> & geneRanges) {
    ifstream inFile(gtfFile);
    if (!inFile) {
        cerr << "Error: failed to open " << gtfFile << endl;
        exit(EXIT_FAILURE);
    }
 
    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string seqname, genename;
        uint64_t start, end;
        getline(ss, seqname, '\t');
        getline(ss, genename, '\t');
        string field;
        getline(ss, field, '\t');
        start = stoull(field);
        getline(ss, field, '\t');
        end = stoull(field);
        geneRanges[seqname][genename] = make_pair(start, end);
    }
}
std::string CidExtractPairwiseAlign(std::string qseq, int errorcount, string& alignstat) {
    string cid = "";

    int sublength = 150;
    if (qseq.size() < 150) {
        sublength = qseq.size();
    }
    string seq = qseq.substr(0, sublength);
    string seq_rc = reverseComplement(qseq).substr(0, sublength);
    //cout << seq << endl;
    //cout << seq_rc << endl;

    auto seq_dna5 = seq | seqan3::views::char_to<seqan3::dna5>;
    auto seq_rc_dna5 = seq_rc | seqan3::views::char_to<seqan3::dna5>;

    auto left_results = seqan3::align_pairwise(std::tie(left_dna5, seq_dna5), aligncfg);
    auto right_results = seqan3::align_pairwise(std::tie(right_dna5, seq_dna5), aligncfg);
    auto left_rc_results = seqan3::align_pairwise(std::tie(left_dna5, seq_rc_dna5), aligncfg);
    auto right_rc_results = seqan3::align_pairwise(std::tie(right_dna5, seq_rc_dna5), aligncfg);
    auto& left_res = *left_results.begin(); auto& right_res = *right_results.begin(); auto& left_rc_res = *left_rc_results.begin(); auto& right_rc_res = *right_rc_results.begin();
    uint32_t leftPosEP = left_res.sequence2_end_position(); uint32_t rightPosSP = right_res.sequence2_begin_position();
    uint32_t leftPosEM = left_rc_res.sequence2_end_position(); uint32_t rightPosSM = right_rc_res.sequence2_begin_position();
    uint32_t leftPosE, rightPosS;
    string tureseq;
    int width, hit = 0;

    int minscore1 = anchor[0].size() - errorcount;
    int minscore2 = anchor[1].size() - errorcount;
    //cout << "left score: " << left_res.score() << "right score: " << right_res.score() << "minscore " << minscore << endl;
    //cout << leftPosEP << "\t" << rightPosSP << "\t" << leftPosEM << "\t" << rightPosSM << endl;
    int leftscore, rightscore;
    if ((left_res.score() >= minscore1 || right_res.score() >= minscore2) && (left_rc_res.score() < minscore1 && right_rc_res.score() < minscore2)) {
        hit++;
        leftPosE = leftPosEP;
        rightPosS = rightPosSP;
        tureseq = seq;
        leftscore = left_res.score();
        rightscore = right_res.score();
        alignstat = "+\t";
    }
    else if ((left_res.score() < minscore1 && right_res.score() < minscore2) && (left_rc_res.score() >= minscore1 || right_rc_res.score() >= minscore2)) {
        hit++;
        leftPosE = leftPosEM;
        rightPosS = rightPosSM;
        tureseq = seq_rc;
        leftscore = left_rc_res.score();
        rightscore = right_rc_res.score();
        alignstat = "-\t";
    }
    if (hit > 0) {
        if (leftscore >= minscore1 && rightscore >= minscore2) {
            width = rightPosS - leftPosE;
            if (width >= mywidth[0] && width <= mywidth[1]) {
                cid = tureseq.substr(leftPosE, width);
                alignstat = alignstat + to_string(leftPosE) + "\t" + "2_anchor";
            }
        }
        else if (rightscore >= minscore1 && rightPosS >= mywidth[2]) {
            cid = tureseq.substr(rightPosS - mywidth[2], mywidth[2]);
            alignstat = alignstat + to_string(rightPosS - mywidth[2]) + "\t" + "right_anchor";
        }
        else if (leftscore >= minscore2) {
            cid = tureseq.substr(leftPosE, mywidth[2]);
            alignstat = alignstat + to_string(leftPosE) + "\t" + "left_anchor";
        }
    }

    return(cid);
}
void handleBamChunk(const string& bamFile, uint64_t start, uint64_t end, map <string, map<string, pair<uint64_t, uint64_t>>>& geneRanges, int errorcount, ofstream &outputFile, mutex& mtx) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    bam1_t* record = bam_init1();
   
    uint64_t count = 0;
    while (sam_read1(fp, header, record) >= 0) {
        if (record->core.flag & BAM_FUNMAP) {
            continue;
        }

        if (++count < start) {
            continue;
        }

        if (count > end) {
            break;
        }
        if (count % 10000 == 0) {
            cout << count << endl;
        }
        ReadInfo read;
        read.id = bam_get_qname(record);

        uint8_t* seq = bam_get_seq(record);
        std::string seq_str;

        for (int i = 0; i < record->core.l_qseq; ++i) {
            seq_str.push_back(seq_nt16_str[bam_seqi(seq, i)]);
        }
        read.seq = seq_str;
        
        int tid = record->core.tid;
        const char* chrName = header->target_name[tid];
        string chrNameStr(chrName);
        uint64_t pos = record->core.pos;
        uint64_t end = bam_endpos(record);
        //cout << chrNameStr << "\t" << pos << "\t" << end << endl;
        read.gene = "none";
        for (const auto& gene : geneRanges[chrNameStr]) {
            if (gene.first.empty()) {
                continue;
            }
            if (pos <= gene.second.second  && end >= gene.second.first) {
                //cout << gene.first << "\t" << gene.second.first << "\t" << gene.second.second << endl;
                if (read.gene == "none") {
                    read.gene = gene.first;
                }
                else {
                    read.gene = read.gene + "," + gene.first;
                }
                
            }
        }

        //extract cid
        string alignstat;
        string cid = CidExtractPairwiseAlign(read.seq, errorcount, alignstat);
        mtx.lock();
        if (cid != "") {
            // Lock the output file stream to prevent multiple threads writing to it simultaneously
            outputFile << read.id << '\t' << alignstat << "\t" << cid << "\t" << read.gene << "\t" << pos << "\t" << end << endl;
        }
        mtx.unlock();
    }

    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);
}

void handleBam(const string& bamFile, const string& outputFileName, map <string, map<string, pair<uint64_t, uint64_t>>>& geneRanges, int errorcount, int numThreads, uint64_t recordnum) {
    //cout << "bam open " << endl;
    ofstream outputFile(outputFileName);
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    //cout << "bam open index " << endl;
    hts_idx_t* idx = sam_index_load(fp, bamFile.c_str());
    uint64_t mapped, unmapped;
    //cout << "bam get stat " << endl;
    if (hts_idx_get_stat(idx, 0, &mapped, &unmapped) < 0) {
        cerr << "Error: failed to get index statistics" << endl;
    }
    uint64_t n_no_coor = hts_idx_get_n_no_coor(idx);
    cout << "The BAM file has " << n_no_coor << " records with no coordinates, "
        << mapped << " records with coordinates, and " << unmapped << " unmapped records." << endl;

    uint64_t n_records;
    if (recordnum != 0) {
        n_records = recordnum;
    }
    else {
        n_records = mapped;
    }
    cout << n_records << endl;
    sam_close(fp);
    bam_hdr_destroy(header);

    vector<thread> threads;
    mutex mtx;

    uint64_t chunk_size = n_records / numThreads;
    uint64_t remainder = n_records % numThreads;

    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(handleBamChunk, bamFile, start, end, ref(geneRanges), errorcount, ref(outputFile), ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
}

void findCIDSequence(const vector<ReadInfo>& reads, const string& outputFileName, int errorcount) {
    ofstream outputFile(outputFileName);
    vector <string> cidvector;
    for (auto& read : reads)
    {
        cout << read.id << ":" << endl;
        string alignstat;
        string cid = CidExtractPairwiseAlign(read.seq, errorcount, alignstat);
        if (cid != "") {
            outputFile << read.id << '\t' << alignstat << "\t" << cid << "\t"  << read.seq << endl;
        }
    }
    
    outputFile.close();
}

void findCIDSequenceThread(const vector<ReadInfo>& reads, const string& outputFileName, int numThreads, int errorcount) {
    ofstream outputFile(outputFileName);
    mutex outputFileMutex;

    // Define a lambda function that will be executed by each thread
    auto worker = [&](const int threadId) {
        // Calculate the range of reads for this thread
        const uint64_t chunkSize = reads.size() / numThreads;
        const uint64_t startIdx = threadId * chunkSize;
        const uint64_t endIdx = (threadId == numThreads - 1) ? reads.size() : (threadId + 1) * chunkSize;

        // Process reads in the range
        for (uint64_t i = startIdx; i < endIdx; i++) {
            cout << startIdx << "-" << endIdx << endl;
            const auto& read = reads[i];
            //cout << read.id;
            string alignstat;
            string cid = CidExtractPairwiseAlign(read.seq, errorcount, alignstat);
            if (i % 10000 == 0) {
                cout << "processed reads:" << i << endl;
            }
            if (cid != "") {
                //cout << "\t"<<"hit" << endl;
                // Lock the output file stream to prevent multiple threads writing to it simultaneously
                lock_guard<mutex> lock(outputFileMutex);
                outputFile << read.id << '\t' << alignstat << "\t" << cid << "\t" << read.gene  << endl;
            }
            else {
                //cout << "\t" << "nohit" << endl;
            }
        }
    };

    // Start multiple threads
    vector<thread> threads;
    for (int i = 0; i < numThreads; i++) {
        threads.emplace_back(worker, i);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    outputFile.close();
}

int main(int argc, char* argv[]) {
    cout << "usage:" << endl;
    string bamFileA = argv[1];
    string geneAn = argv[2];
    string outFile = argv[3];
    int thread = atoi(argv[4]);
    int errorcount = atoi(argv[5]);
    uint64_t recordnum = 0;
    if (argc == 7) {
        string  recordnumst = argv[6];
        recordnum = stoull(recordnumst);
    }
    //double errorate = atof(argv[4]);
    //part of cidprimer, capture oligo
    cout << "load bam: " << endl;
    //vector<ReadInfo> reads;
    map <string, map<string, pair<uint64_t, uint64_t>>> geneRanges;
    loadGeneRanges(geneAn, geneRanges);
    handleBam(bamFileA, outFile, geneRanges, errorcount, thread, recordnum);
    //findCIDSequenceThread(reads,  outprex, thread, errorcount);

    return 0;
}