#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <future>
#include <thread>
#include <map>
#include <sstream>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/utility.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include<time.h>
#include "MurmurHash3.h"
#include <unordered_set>
using namespace std;
using namespace seqan3;

size_t bucketNum = 4;
size_t kmerLen = 5;
auto alnconfig = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

uint64_t murmurHash3(string seq, uint64_t seed) {
    uint64_t out;
    MurmurHash3_x64_128(seq.data(), seq.size(), seed, &out);
    return out;
}

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

void loadNanoReads(string readFile, map<string, vector<string>> & readsIdls, map<string, vector<string>> & readsSeqls,
    map<string, vector<pair<uint64_t, uint64_t>>>& readsRangels,
    vector<string> &genels) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " <<readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string readid, strand, mpos, mstat, seq, genes, startst, endst;
        getline(ss, readid, '\t');
        getline(ss, strand, '\t');
        getline(ss, mpos, '\t');
        getline(ss, mstat, '\t');
        getline(ss, seq, '\t');
        getline(ss, genes, '\t');
        getline(ss, startst, '\t');
        getline(ss, endst, '\t');
        uint64_t start, end;
        start = stoull(startst);
        end = stoull(endst);

        if (genes != "none") {
            istringstream genesS(genes);
            string gene;
            while (getline(genesS, gene, ',')) {
                //cout << gene << endl;
                readsIdls[gene].push_back(readid);
                readsSeqls[gene].push_back(seq);
                readsRangels[gene].push_back(make_pair(start, end));
                genels.push_back(gene);
            }
        }

        
    }
    cout << genels.size() << endl;
    auto new_end = std::unique(genels.begin(), genels.end());
    genels.erase(new_end, genels.end());
    cout << genels.size() << endl;

}

int jaccardCount(vector<uint64_t> &v1, vector<uint64_t> &v2) {
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    vector<uint64_t> result;
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(result));
    return result.size();
}

int pairwiseAlign(string seq1, string seq2) {
    auto seq1_dna5 = seq1 | seqan3::views::char_to<seqan3::dna5>;
    auto seq2_dna5 = seq2 | seqan3::views::char_to<seqan3::dna5>;
    auto results = seqan3::align_pairwise(std::tie(seq1_dna5, seq2_dna5), alnconfig);
    auto& result = *results.begin();
    return result.score();
}

void matchCidChunk(map<string, vector<string>>& readsIdls, map<string, vector<string>>& readsSeqls, map<string, vector<pair<uint64_t, uint64_t>>>& readsRangels,
    vector<string> & genels, string indexFold, ofstream &outputFile, int errorcount, int start, int end, mutex& mtx) {
    cout << "start" << endl;

    for (int i = start; i <= end; i++)
    {
        //cout << genels.size() << endl;
        string gene = genels[i];
        //cout << "load fm index" << endl;
        string indexFile = indexFold + "/" + gene + ".bin";
        //bi_fm_index<char, seqan3::text_layout::collection> bifm;
        map<uint64_t, map<uint64_t, vector<uint64_t>>> buckets;
        map<uint64_t, string> seqHash;
        map<uint64_t, string> posHash;
        map<uint64_t, vector<uint64_t>> kmerHash;
        map<uint64_t, pair<uint64_t, uint64_t>> rangeHash;
        std::ifstream is(indexFile, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(buckets, seqHash, posHash, kmerHash, rangeHash);

        //cout << "search cid" << endl;
        for (uint64_t k = 0; k < readsSeqls[gene].size(); k++)
        {
            double startTime = clock();
            string seqi = readsSeqls[gene][k];
            string readid = readsIdls[gene][k] + "\t" + readsSeqls[gene][k];
            pair<uint64_t, uint64_t> readRange = readsRangels[gene][k];

            string resPos, resSeq;
            pair<uint64_t, uint64_t> resRange;

            if (seqi.size() < 20 | seqi.size() > 30) {
                continue;
            }
            uint64_t hashTotal = murmurHash3(seqi, 0);

            //if exact match;
            auto it = posHash.find(hashTotal);
            if (it != posHash.end()) {
                cout << "exactMatch\t" << readid << endl;
                resPos = it->second;
                resSeq = seqHash[hashTotal];
                resRange = rangeHash[hashTotal];
                double endTime = clock();
                double duration = (endTime - startTime) / CLOCKS_PER_SEC;
                mtx.lock();
                outputFile << readid << "\t" << resPos << "\t" << resSeq << "\t" << 0 << "\t"  << gene << "\t" << readRange.first << "\t" << readRange.second<< "\t" << resRange.first << "\t" << resRange.second << "\t" << duration << endl;
                mtx.unlock();
            }
            else {
                //build hash vec for compare
                
                vector<uint64_t> hashVec;
                kmer_hash(seqi, 3, 0, hashVec);

                int hit = 0;
                //check each bucket collection
                unordered_set<uint64_t> appeared;
                for (size_t seed = 0; seed < bucketNum; seed++)
                {
                    //cout << seed << ":\t";
                    uint64_t hashMin = min_hash(seqi, kmerLen, seed);
                    auto bktit = buckets[seed].find(hashMin);

                    //if found bucket in this seed
                    if (bktit != buckets[seed].end()) {
                        hit = 1;
                        //cout << "minHash found\t";
                        vector<uint64_t> targetBkt = buckets[seed][hashMin];

                        //iterate through bucket
                        int maxscore = -100;
                        for (size_t tgti = 0; tgti < targetBkt.size(); tgti++)
                        {
                            uint64_t hashTotali = targetBkt[tgti];

                            if(appeared.find(hashTotali) == appeared.end()) {
                                appeared.insert(hashTotali);

                                //compare hashVec
                                vector<uint64_t> hashVeci = kmerHash[hashTotali];
                                int matchCount = jaccardCount(hashVec, hashVeci);

                                //if Vec is similar, process pairwise alignment
                                if (matchCount >= 7) {
                                    string refseqi = seqHash[hashTotali];
                                    int score = pairwiseAlign(seqi, refseqi);
                                    // cout << "match score: " << score << endl;
                                    if (score > maxscore && score >= -errorcount) {
                                        maxscore = score;
                                        resSeq = refseqi;
                                        resPos = posHash[hashTotali];
                                        resRange = rangeHash[hashTotali];
                                    }
                                }
                            }
                            
                        }

                        if (maxscore != -100) {
                            double endTime = clock();
                            double duration = (endTime - startTime) / CLOCKS_PER_SEC;
                            mtx.lock();
                            outputFile << readid << "\t" << resPos << "\t" << resSeq << "\t" << 0 << "\t" << gene << "\t" << readRange.first << "\t" << readRange.second << "\t" << resRange.first << "\t" << resRange.second << "\t" << duration << endl;
                            mtx.unlock();
                            break;
                        }
                        
                    }
                }

                if (hit == 0) {
                    cout << "noHash\t" << readid << endl;
                }
                else {
                    cout << "hitHash\t" << readid << endl;
                }
            }
            
        }

        
    }

}
void matchCidThread(map<string, vector<string>>& readsIdls, map<string, vector<string>>& readsSeqls, map<string, vector<pair<uint64_t, uint64_t>>>& readsRangels,
    vector<string>& genels, string indexFold, string outputFileName, unsigned char errorcount, int numThreads) {
    vector<thread> threads;
    mutex mtx;
    int n_records = genels.size();
    int chunk_size = n_records / numThreads;
    int remainder = n_records % numThreads;
    
    ofstream outputFile(outputFileName);
    outputFile << "readid\t" << "cidPos\t" << "cidSeq\t" << "editDi\t" << "gene\t" << "readStart\t" << "readEnd\t" << "refStart\t" << "refEnd\t" << "time" << endl;
    int start = 0, end = 0;
    for (int i = 0; i < numThreads; i++) {
        start = end;
        if (start != 0) {
            start++;
        }
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(matchCidChunk, ref(readsIdls),  ref(readsSeqls), ref(readsRangels), ref(genels), indexFold, ref(outputFile), errorcount, start, end, ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }
}

int main(int argc, char* argv[]) {
    string readFile = argv[1];
    string indexFold = argv[2];
    string outFile = argv[3];
    int thread = atoi(argv[4]);
    int errorcount = atoi(argv[5]);

    map<string, vector<string>> readsIdls;
    map<string, vector<string>> readsSeqls;
    map<string, vector<pair<uint64_t, uint64_t>>> readsRangels;
    vector<string> genels;
    loadNanoReads(readFile, readsIdls, readsSeqls, readsRangels, genels);
    matchCidThread(readsIdls, readsSeqls, readsRangels, genels, indexFold, outFile, errorcount, thread);

    return 0;
}
