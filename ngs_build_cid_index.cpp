#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <future>
#include <thread>
#include <map>
#include <sstream>


using namespace std;
struct ReadInfo {
    string id;
    string seq;
    string pos;
    string gene;
};

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

void readBamChunk(const string& bamFile, uint64_t start, uint64_t end, map <string, map<string, pair<uint64_t, uint64_t>>>& geneRanges, vector<ReadInfo>& reads, mutex& mtx, ofstream & outputFile) {
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
        std::string cidseq, cidpos = "";
        cidseq = read.id.substr(37, 25);
        for (int k = 0; read.id[37 + 28 + k] != '|'; k++) {
            cidpos = cidpos + read.id[37 + 28 + k];
        }
        read.seq = cidseq;
        read.pos = cidpos;
        read.id = read.id.substr(0, 29);

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
        // ¼ÓËø
        mtx.lock();
        outputFile << read.id << '\t' << read.seq << "\t" << read.pos << "\t" << read.gene << endl;
        reads.emplace_back(read);
        mtx.unlock();
    }

    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);
}

void readBam(const string& bamFile, vector<ReadInfo>& reads, map <string, map<string, pair<uint64_t, uint64_t>>>& geneRanges,  const string outputFileName, int numThreads, uint64_t recordnum) {
    //cout << "bam open " << endl;
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

    ofstream outputFile(outputFileName);
    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(readBamChunk, bamFile, start, end, ref(geneRanges),ref(reads), ref(mtx), ref(outputFile));
    }
    
    for (auto& t : threads) {
        t.join();
    }
}

int main(int argc, char* argv[]) {
    cout << "usage:" << endl;
    string bamFileA = argv[1];
    string geneAn = argv[2];
    string outfile = argv[3];
    int thread = atoi(argv[4]);
    uint64_t recordnum = 0;
    if (argc == 6) {
        string  recordnumst = argv[5];
        recordnum = stoull(recordnumst);
    }
    //double errorate = atof(argv[4]);
    //part of cidprimer, capture oligo
    cout << "load bam: " << endl;
    vector<ReadInfo> reads;
    map <string, map<string, pair<uint64_t, uint64_t>>> geneRanges;
    loadGeneRanges(geneAn, geneRanges);
    readBam(bamFileA, reads, geneRanges, outfile, thread, recordnum);

    return 0;
}