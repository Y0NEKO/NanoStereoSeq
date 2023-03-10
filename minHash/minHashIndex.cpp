#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include "MurmurHash3.h"
#include<time.h>
#include<map>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <htslib/sam.h>
#include <set>
#include <thread>
using namespace std;

size_t kmerLen = 5;
size_t bucketNum = 4;

void loadGeneRanges(string gtfFile, vector<pair <string, string>>& geneRanges) {
    ifstream inFile(gtfFile);
    if (!inFile) {
        cerr << "Error: failed to open " << gtfFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string seqname, genename;
        string start, end;
        getline(ss, seqname, '\t');
        getline(ss, genename, '\t');
        getline(ss, start, '\t');
        getline(ss, end, '\t');
        string rg = seqname + ":" + start + "-" + end;
        geneRanges.push_back(make_pair(genename, rg));
    }
}

uint64_t murmurHash3(string seq, uint64_t seed) {
    uint64_t out;
    MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
    return out;
}

//get min hash of different k-mer

uint64_t min_hash(string seq, size_t len, uint64_t seed) {
    vector<uint64_t> hashvec;
    for (size_t i = 0; i <= (seq.size() - len); i++)
    {
        string seqsub = seq.substr(i, len);
        uint64_t hasht = murmurHash3(seqsub, seed);
        hashvec.push_back(hasht);
    }
    auto min_it = min_element(hashvec.begin(), hashvec.end());
    uint64_t min_val = static_cast<uint64_t>(*min_it);
    return(min_val);
}

void kmer_hash(string seq, size_t len, uint64_t seed, vector<uint64_t>& hashvec) {
    for (size_t i = 0; i <= (seq.size() - len); i++)
    {
        string seqsub = seq.substr(i, len);
        uint64_t hasht = murmurHash3(seqsub, seed);
        hashvec.push_back(hasht);
    }
}

/*
char rand_base() {
    static const char bases[] = { 'A', 'C', 'G', 'T' };
    return bases[rand() % 4];
}

//generate random reads
void generate_reads(string& read, int index, const int n, const int k, vector<string>& reads) {
    if (index == n) {
        reads.push_back(read);
        return;
    }
    for (int i = 0; i < k; i++) {
        read[index] = rand_base();
        generate_reads(read, index + 1, n, k, reads);
    }
}
*/


// build buckets for sequence vector
void bucket_build(map<uint64_t, map<uint64_t, vector<uint64_t>>>& buckets, map<uint64_t, string> &seqHash, map<uint64_t, string>& posHash, 
    map<uint64_t, pair<uint64_t, uint64_t>> &rangeHash,
    map<uint64_t, vector<uint64_t>>& kmerHash, 
    string seq, string pos, pair<uint64_t, uint64_t> &range) {

    uint64_t hashTotal = murmurHash3(seq, 0);
    vector<uint64_t> hashVec;
    kmer_hash(seq, 3, 0, hashVec);
    seqHash[hashTotal] = seq;
    posHash[hashTotal] = pos;
    kmerHash[hashTotal] = hashVec;
    rangeHash[hashTotal] = range;

    for (size_t seed = 0; seed < bucketNum; seed++)
    {
        uint64_t hashMin = min_hash(seq, kmerLen, seed);
        buckets[seed][hashMin].push_back(hashTotal);
    }

}

void map_stat(map<uint64_t, vector<uint64_t>> &myHash) {
    int numKeys = myHash.size();
    std::vector<int> value_sizes;
    for (auto const& [key, value] : myHash) {
        value_sizes.push_back(value.size());
    }
    std::sort(value_sizes.begin(), value_sizes.end());
    int min_size = value_sizes.front();
    int median_size = value_sizes[value_sizes.size() / 2];
    int max_size = value_sizes.back();
    std::cout << "keys: " << numKeys <<  "\tMin: " << min_size << "\tMedian: " << median_size << "\tMax: " << max_size << std::endl;
}

void gene_index_precess(vector<pair <string, string>>& geneRanges, string bamFile, string outFold, uint64_t gstart, uint64_t gend) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    bam1_t* record = bam_init1();
    hts_idx_t* idx = sam_index_load(fp, bamFile.c_str());

    for (int gi = gstart; gi <= gend; gi++) {
        //cout << "extract gene" << endl;
        string geneid = geneRanges[gi].first;
        string rg = geneRanges[gi].second;
        const char* c_rg = rg.c_str();
        hts_itr_t* iter = sam_itr_querys(idx, header, c_rg);

        // extract cid seq and pos
        set<string> uniqueSet;
        // extract read align pos
        map<uint64_t, map<uint64_t, vector<uint64_t>>> buckets;
        map<uint64_t, string> seqHash;
        map<uint64_t, string> posHash;
        map<uint64_t, pair<uint64_t, uint64_t>> rangeHash;
        map<uint64_t, vector<uint64_t>> kmerHash;
        //cout << "process reads" << endl;
        while (sam_itr_next(fp, iter, record) >= 0) {
            
            std::string cidseq, cidpos;
            string rqns = bam_get_qname(record);
            cidseq = rqns.substr(37, 25);
            for (int k = 0; rqns[37 + 28 + k] != '|'; k++) {
                cidpos = cidpos + rqns[37 + 28 + k];
            }
            //cout << rqns << endl;
            if (uniqueSet.count(cidseq) == 0) {
                uniqueSet.insert(cidseq);
                //cout << "find range" << endl;
                uint64_t rpos = record->core.pos;
                uint64_t rend = bam_endpos(record);
                pair<uint64_t, uint64_t> cidrange = make_pair(rpos, rend);
                //cout << "build bucket" << endl;
                bucket_build(buckets, seqHash, posHash, rangeHash, kmerHash, cidseq, cidpos, cidrange);
            }
        }
        hts_itr_destroy(iter);

        if (uniqueSet.size() == 0) {
            continue;
        }
        std::cout << geneid << " unique_CIDseq_size\t" << uniqueSet.size() << std::endl;

        //Build minhash index
        cout << "buckets size: " << buckets.size() << endl;
        for (size_t bi = 0; bi < buckets.size(); bi++)
        {
            cout << "bucket" << bi << "\t";
            map_stat(buckets[bi]);
        }

        //output index
        string outFm = outFold + "/" + geneid + ".bin";
        std::ofstream os{ outFm, std::ios::binary };
        cereal::BinaryOutputArchive oarchive{ os };
        oarchive(buckets, seqHash, posHash, kmerHash, rangeHash);
        os.close();

        //mtx.lock();
        //posVec[geneid] = cidPosls;
       // mtx.unlock();
        buckets.clear();
        seqHash.clear();
        posHash.clear();
        kmerHash.clear();
        rangeHash.clear();

    }

    hts_idx_destroy(idx);
    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);


}

void gene_index_process_thread(vector<pair <string, string>>& geneRanges, string bamFile, int numThreads, string outFold) {
    uint64_t genecount = geneRanges.size();
    uint64_t chunk_size = genecount / numThreads;
    uint64_t remainder = genecount % numThreads;
    vector<thread> threads;
    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(gene_index_precess, ref(geneRanges), bamFile, outFold, start, end - 1);
    }

    for (auto& t : threads) {
        t.join();
    }


}

int main(int argc, char* argv[]) {
    string bamFile = argv[1];
    string geneFile = argv[2];
    string outFold = argv[3];
    int numThreads = atoi(argv[4]);
    
    if (!filesystem::exists(outFold)) // if dir exist
    {
        if (filesystem::create_directory(outFold)) // create dir
        {
            cout << "Create fold" << std::endl;
        }
        else
        {
            cerr << "Create fold fail£¡" << std::endl;
        }
    }

    vector<pair <string, string>> geneRanges;
    cout << "loading readlist" << endl;
    loadGeneRanges(geneFile, geneRanges);
    cout << "building index" << endl;
    gene_index_process_thread(geneRanges, bamFile, numThreads, outFold);
    cout << "All index has been built successfully !" << endl;
    return 0;
}
