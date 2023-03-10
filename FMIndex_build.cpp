#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <map>
#include <unordered_set>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>
#include <htslib/sam.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <filesystem>
using namespace std;
using namespace seqan3;

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
        geneRanges.push_back(make_pair(genename,rg));
    }
}

void GenerateFasta(vector<pair <string, string>>& geneRanges, int gstart, int gend, string bamFile, string outFold,
    map<string, vector<string>>& posVec , mutex & mtx) {
    samFile* fp = sam_open(bamFile.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    bam1_t* record = bam_init1();
    hts_idx_t* idx = sam_index_load(fp, bamFile.c_str());

    for (int gi = gstart; gi < gend; gi++) {
        string geneid = geneRanges[gi].first;
        string rg = geneRanges[gi].second;
        const char* c_rg = rg.c_str();

        /*=============================================Make CID list================================================ */
        hts_itr_t* iter = sam_itr_querys(idx, header, c_rg);
        //cout << "extract cidseq" << endl;
       // std::vector<std::vector<seqan3::dna4>>  cidSeqDnals;
        
        std::vector<pair<string, string>>  cidInfols;
        //std::cout << "load cid" << thread << "\n";
        while (sam_itr_next(fp, iter, record) >= 0) {
            std::string cidseq, cidpos;
            string rqns = bam_get_qname(record);
            cidseq = rqns.substr(37, 25);
            for (int k = 0; rqns[37 + 28 + k] != '|'; k++) {
                cidpos = cidpos + rqns[37 + 28 + k];
            }
            //auto cidseq_dna5 = cidseq | seqan3::views::char_to<seqan3::dna4>;
            //cidSeqDnals.push_back(cidseq_dna5);
            cidInfols.push_back(make_pair(cidseq, cidpos));
        }
        hts_itr_destroy(iter);

        if (cidInfols.size() == 0) {
            continue;
        }
        //build index

        std::cout << geneid << " CID raw size" << cidInfols.size() << std::endl;

        auto new_end = std::unique(cidInfols.begin(), cidInfols.end());
        cidInfols.erase(new_end, cidInfols.end());
        std::cout << geneid << " CID unique size" << cidInfols.size() << std::endl;
        std::vector<string>  cidSeqls(cidInfols.size());
        std::vector<string>  cidPosls(cidInfols.size());

        transform(cidInfols.begin(), cidInfols.end(), cidSeqls.begin(),
            [](const pair<string, string>& p) { return p.first; });

        transform(cidInfols.begin(), cidInfols.end(), cidPosls.begin(),
            [](const pair<string, string>& p) { return p.second; });

        //using collection bi fm index
        //bi_fm_index  bifm{ cidSeqls };

        string cidSeq = cidSeqls[0];
        for (size_t i = 1; i < cidSeqls.size(); i++)
        {
            cidSeq = cidSeq + "NNNNN" + cidSeqls[i];
        }
        //using single bi fm index
        auto cidSeq_dna5 = cidSeq | seqan3::views::char_to<seqan3::dna5>;
        bi_fm_index  bifm{ cidSeq_dna5 };
        
        string outFm = outFold + "/" + geneid + ".bin";
        std::ofstream os{ outFm, std::ios::binary };
        cereal::BinaryOutputArchive oarchive{ os };
        oarchive(bifm);

        mtx.lock();
        posVec[geneid] = cidPosls;
        mtx.unlock();
        cidSeqls.clear();
        cidPosls.clear();
        cidInfols.clear();
        std::vector<pair<string, string>>().swap(cidInfols);
        std::vector <string>().swap(cidSeqls);
        std::vector <std::string>().swap(cidPosls);

    }

    hts_idx_destroy(idx);
    sam_close(fp);
    bam_destroy1(record);
    bam_hdr_destroy(header);

}

void GenerateFastaThread(vector<pair <string, string>>& geneRanges, string bamFile, int numThreads, string outFold) {
    uint64_t genecount = geneRanges.size();
    uint64_t chunk_size = genecount / numThreads;
    uint64_t remainder = genecount % numThreads;
    vector<thread> threads;
    mutex mtx;
    uint64_t start = 0, end = 0;

    map<string, vector<string>> posVec;
    map<string, vector<string>> seqVec;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        if (start != 0) {
            start++;
        }
        end = start + chunk_size;
        cout << start << "-" << end << endl;
        threads.emplace_back(GenerateFasta, ref(geneRanges), start, end, bamFile, outFold, ref(posVec),  ref(mtx));
    }

    for (auto& t : threads) {
        t.join();
    }

    // 序列化向量到文件中
    string posout = outFold + "cidPos.bin";

    std::ofstream os2{ posout, std::ios::binary };
    cereal::BinaryOutputArchive oarchive2{ os2 };
    oarchive2(posVec);
    os2.close();
}

/*
void loadReadsList(string readFile, map <string, vector<string>>& geneCidls, map <string, vector<string>>& geneCidPosls, unordered_set<std::string>& gene_set) {
    ifstream inFile(readFile);
    if (!inFile) {
        cerr << "Error: failed to open " << readFile << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    uint64_t count = 0;
    while (getline(inFile, line)) {
        count++;
        if (count % 1000000 == 0) {
            cout << count << endl;
        }
        istringstream ss(line);
        string readid, cid, cidpos, genes, gene;
        getline(ss, readid, '\t');
        getline(ss, cid, '\t');
        getline(ss, cidpos, '\t');
        getline(ss, genes, '\t');
        if (genes != "none") {
            istringstream geneS(genes);
            while (getline(geneS, gene, ',')) {
                geneCidls[gene].push_back(cid);
                geneCidPosls[gene].push_back(cidpos);
                if (gene_set.count(gene) == 0)
                {
                    gene_set.insert(gene);
                }
            }
        }

    }

}

void BuildFmIndex(map <string, vector<string>> &geneCidls, map <string, vector<string>>& geneCidPosls, vector<bi_fm_index> &fm_vector, uint64_t start, uint64_t end, mutex& mtx, vector<string> &gene_vector) {
    for (uint64_t i = start; i < end; i++)
    {
        string gene = gene_vector[i];
        auto new_end = std::unique(geneCidPosls[gene].begin(), geneCidPosls[gene].end());
        geneCidls[gene].erase(new_end, geneCidls[gene].end());
        geneCidPosls[gene].erase(new_end, geneCidPosls[gene].end());
        bi_fm_index fm{ geneCidls[gene] };
        mtx.lock();
        fm_vector.push_back(fm);
        mtx.unlock();
    }
}
void BuildFmIndexThread(int numThreads, map <string, vector<string>>& geneCidls, map <string, vector<string>>& geneCidPosls, vector<string>& gene_vector, string outPrex) {
    vector<bi_fm_index> fm_vector;
    uint64_t genecount = gene_vector.size();
    uint64_t chunk_size = genecount / numThreads;
    uint64_t remainder = genecount % numThreads;
    vector<thread> threads;
    mutex mtx;
    uint64_t start = 0, end = 0;

    std::unordered_map<std::string, sdsl::csa_wt<>> fm_index_map;
    vector<bi_fm_index> fm_vector;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(BuildFmIndex, ref(geneCidls),ref(geneCidPosls), ref(fm_vector), start, end, ref(mtx), ref(gene_vector));
    }

    for (auto& t : threads) {
        t.join();
    }

    // 序列化向量到文件中
    std::ofstream os{ "index.file", std::ios::binary };
    cereal::BinaryOutputArchive oarchive{ os };
    oarchive(index);

    string outSeq = outPrex + "_seq.bin";
    ofstream os{ outSeq, ios::binary };
    serialize(os, fm_vector);
    os.close();

    string outPos = outPrex + "_pos.bin";
    ofstream ofs(outPos, ios::binary);
    cereal::BinaryOutputArchive ar(ofs);
    ar(geneCidPosls);
    ofs.close();

    // 从文件中反序列化向量
    /*std::vector<bi_fm_index<seqan3::dna5_vector>> fm_index_vector_from_file;
    std::ifstream is{ "fm_index_vector.bin", std::ios::binary };
    seqan3::deserialize(is, fm_index_vector_from_file);
    is.close();

    ifstream ifs("my_map.bin", ios::binary);
    cereal::BinaryInputArchive iar(ifs);
    unordered_map<string, int> new_map;
    iar(new_map);

}

*/

int main(int argc, char* argv[]) {
    cout << "usage:" << endl;
    string readFile = argv[1];
    string geneFile = argv[2];
    string outFile = argv[3];
    int thread = atoi(argv[4]);

    if (!filesystem::exists(outFile)) // 判断目录是否存在
    {
        if (filesystem::create_directory(outFile)) // 创建目录
        {
            std::cout << "Create fold" << std::endl;
        }
        else
        {
            std::cerr << "Create fold fail！" << std::endl;
        }
    }

    vector<pair <string, string>> geneRanges;
    cout << "load readlist: " << endl;
    loadGeneRanges(geneFile, geneRanges);
    GenerateFastaThread(geneRanges, readFile, thread, outFile);
    //loadReadsList(readFile, geneCidls, geneCidPosls, gene_set);
    //vector<string> gene_vector(gene_set.begin(), gene_set.end());
    //BuildFmIndexThread(thread, geneCidls, geneCidPosls, gene_vector, outFile);
    return 0;
}