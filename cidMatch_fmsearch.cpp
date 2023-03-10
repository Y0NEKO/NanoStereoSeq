#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <ranges>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/all.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <future>
#include <thread>
#include <map>
#include <sstream>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include<time.h>
using namespace std;
using namespace seqan3;

void loadNanoReads(string readFile, map<string, vector<string>> & readsIdls, map<string, vector<string>> & readsSeqls, vector<string> &genels) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " <<readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(inFile, line)) {
        istringstream ss(line);
        string readid, strand, mpos, mstat, seq, genes;
        getline(ss, readid, '\t');
        getline(ss, strand, '\t');
        getline(ss, mpos, '\t');
        getline(ss, mstat, '\t');
        getline(ss, seq, '\t');
        getline(ss, genes, '\t');

        if (genes != "none") {
            istringstream genesS(genes);
            string gene;
            while (getline(genesS, gene, ',')) {
                //cout << gene << endl;
                readsIdls[gene].push_back(readid);
                readsSeqls[gene].push_back(seq);
                genels.push_back(gene);
            }
        }
        
        
    }
    cout << genels.size() << endl;
    auto new_end = std::unique(genels.begin(), genels.end());
    genels.erase(new_end, genels.end());
    cout << genels.size() << endl;

}

void loadPosIndex(map<string, vector<string>> &posVec, string indexFold) {
    string posIndexFile = indexFold + "/" + "cidPos.bin";
    std::ifstream ifs(posIndexFile, std::ios::binary);
    cereal::BinaryInputArchive ar(ifs);
    ar(posVec);

}

void matchCidChunk(map<string, vector<string>>& readsIdls, map<string, vector<string>>& readsSeqls, map<string, vector<string>>& posVec,
    vector<string> & genels, string indexFold, ofstream &outputFile, unsigned char errorcount, int start, int end, mutex& mtx) {
    cout << "start" << endl;
    seqan3::configuration cfg1 = seqan3::search_cfg::max_error_total{ seqan3::search_cfg::error_count{errorcount} }
    | seqan3::search_cfg::hit{ seqan3::search_cfg::hit_all_best{} };

    for (int i = start; i <= end; i++)
    {
        //cout << genels.size() << endl;
        string gene = genels[i];
        //cout << "load fm index" << endl;
        string fmFile = indexFold + "/" + gene + ".bin";
        //bi_fm_index<char, seqan3::text_layout::collection> bifm;
        bi_fm_index<dna5, seqan3::text_layout::single> bifm;
        {
            std::ifstream is{ fmFile, std::ios::binary };
            cereal::BinaryInputArchive iarchive{ is };
            iarchive(bifm);
        }
        //cout << "search cid" << endl;
        for (uint64_t k = 0; k < readsSeqls[gene].size(); k++)
        {
            double startTime = clock();
            string seqi = readsSeqls[gene][k];
            if (seqi.size() < 20 | seqi.size() > 30) {
                continue;
            }
            auto seqdna5 = seqi | seqan3::views::char_to<seqan3::dna5>;
            auto searchResults = search(seqdna5, bifm, cfg1);

            string cidpos = "none";
            vector<string> cidposvec;
            for (const auto& searchResult : searchResults)
            {
                //for collect fm
                //auto rid = searchResult.reference_id();
                //string cidpos = posVec[gene][rid];
                //for single fm
                auto rpos = searchResult.reference_begin_position();
                uint64_t rid = rpos / 30;
                uint64_t res = rpos % 30;
                //cout << rid << "\t" << res << endl;
                if (res > 25) {
                    rid = rid + 1;
                }

                if (cidpos == "none") {
                    cidpos = posVec[gene][rid];
                    cidposvec.push_back(cidpos);
                }
                else{
                    string cidposi = posVec[gene][rid];
                    auto it = std::find(cidposvec.begin(), cidposvec.end(), cidposi);
                    if (it == cidposvec.end())
                    {
                        cidposvec.push_back(cidposi);
                        cidpos = cidpos + "," + cidposi;
                    }
                    
                }


            }
            double endTime = clock();
            double duration = (endTime - startTime) / CLOCKS_PER_SEC;
            // std::cout << "Build index " << " time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
            string readid = readsIdls[gene][k] + "\t" + readsSeqls[gene][k];
            mtx.lock();
            outputFile << readid << "\t" << cidpos << "\t" << cidposvec.size() << "\t" << gene  << "\t"<< duration << endl;
            mtx.unlock();
        }

        
    }

}
void matchCidThread(map<string, vector<string>>& readsIdls, map<string, vector<string>>& readsSeqls, map<string, vector<string>>& posVec,
    vector<string>& genels, string indexFold, string outputFileName, unsigned char errorcount, int numThreads) {
    vector<thread> threads;
    mutex mtx;
    int n_records = genels.size();
    int chunk_size = n_records / numThreads;
    int remainder = n_records % numThreads;
    
    ofstream outputFile(outputFileName);
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
        threads.emplace_back(matchCidChunk, ref(readsIdls),  ref(readsSeqls), ref(posVec), ref(genels), indexFold, ref(outputFile), errorcount, start, end, ref(mtx));
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
    unsigned char errorcount = atoi(argv[5]);

    map<string, vector<string>> posVec;
    map<string, vector<string>> readsIdls;
    map<string, vector<string>> readsSeqls;
    vector<string> genels;
    loadPosIndex(posVec, indexFold);
    loadNanoReads(readFile, readsIdls, readsSeqls, genels);
    matchCidThread(readsIdls, readsSeqls, posVec, genels, indexFold, outFile, errorcount, thread);


    return 0;
}
