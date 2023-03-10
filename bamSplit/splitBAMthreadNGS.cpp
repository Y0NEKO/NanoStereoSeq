#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <sstream>
#include <set>
#include <omp.h>
using namespace std;

// Function to split BAM file
void splitBAM(const string& inputBAM,   map<string, string>& readIdToOutput, map<string, BGZF*> &outputFiles,  mutex& mtx, uint64_t start, uint64_t end) {
    samFile* fp = sam_open(inputBAM.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(fp);
    cout << "assigning" << endl;
    bam1_t* record = bam_init1();
    /*
    #pragma omp parallel for
    for (uint64_t i = 0; i < recordnum; i++) {
        bam1_t* aln = bam_init1();
#pragma omp critical
        {
            auto readres = sam_read1(in, header, aln);
        }
        if (readres >= 0) {
            // 解析BAM记录
            readnm++;
            if (readnm % 10000 == 0) {
                cout << readnm << endl;
            }
            string readId = bam_get_qname(aln);

            // for stereo seq
            string cidpos = "";
            for (int k = 0; readId[37 + 28 + k] != '|'; k++) {
                cidpos = cidpos + readId[37 + 28 + k];
            }
            string outputBAM;
#pragma omp critical // 临界区，保证reads的线程安全
            if (readIdToOutput.count(cidpos) > 0) {
                outputBAM = readIdToOutput[cidpos];
                //cout << "hit!" << " cidpos:" << cidpos << ", output: " << outputBAM << endl;
                auto writeres = bam_write1(outputFiles[outputBAM], aln);
            }
        }
        bam_destroy1(aln);
    }
    */
    
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
        string readId = bam_get_qname(record);

        // for stereo seq
        string cidpos = "";
        for (int k = 0; readId[37 + 28 + k] != '|'; k++) {
            cidpos = cidpos + readId[37 + 28 + k];
        }
        string outputBAM;
        mtx.lock();
        if (readIdToOutput.count(cidpos) > 0) {
            outputBAM = readIdToOutput[cidpos];
            //cout << "hit!" << " cidpos:" << cidpos << ", output: " << outputBAM << endl;
            auto writeres = bam_write1(outputFiles[outputBAM], record);
        }
        mtx.unlock();

    }
    
    bam_destroy1(record);
    bam_hdr_destroy(header);
    sam_close(fp);
    
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <input.bam> <ref.txt> <output_prefix> <thread_num> <record_num>" << endl;
        return 1;
    }

    // Input parameters
    string inputBAM = argv[1];
    string outputPrefix = argv[3];
    string refFile = argv[2];
    int numThreads = atoi(argv[4]);
    string recordst = argv[5];
    uint64_t recordnum = stoull(recordst);
    //load ref
    ifstream myref(refFile);
    string line;
    map<string, string> readIdToOutput;
    vector <string> areals;

    cout << "loading ref" << endl;
    while (getline(myref, line)) {
        vector<string> cols;
        stringstream iss(line);
        string col;
        while (getline(iss, col, '\t')) {
            cols.push_back(col);
        }
        if (cols.size() != 2) {
            cerr << "Invalid reference file format" << endl;
            exit(1);
        }
        readIdToOutput[cols[0]] = outputPrefix + "." + cols[1] + ".bam";
        areals.push_back(readIdToOutput[cols[0]]);
    }
    myref.close();

    //open output bam file
    cout << "open bams" << endl;
    set <string> stmp(areals.begin(), areals.end());
    areals.assign(stmp.begin(), stmp.end());

    map<string, BGZF*> outputFiles;
    BGZF* fp;
    samFile* in = sam_open(inputBAM.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(in);

    //build split file
    for (size_t i = 0; i < areals.size(); i++)
    {
        cout << areals[i] << endl;
        outputFiles[areals[i]] = bgzf_open(areals[i].c_str(), "wb");
        int restmp = bam_hdr_write(outputFiles[areals[i]], header);
    }

    // Create threads
    //splitBAM(inputBAM,  readIdToOutput, outputFiles,  recordnum, numThreads);

    mutex mtx;
    uint64_t chunk_size = recordnum / numThreads;
    uint64_t remainder = recordnum % numThreads;    
    cout << "create thread" << endl;
    vector<thread> threads;

    uint64_t start = 0, end = 0;
    for (uint64_t i = 0; i < numThreads; i++) {
        start = end;
        end = start + chunk_size;
        if (i == numThreads - 1) {
            end += remainder;
        }
        cout << start << "-" << end << endl;
        threads.emplace_back(splitBAM, std::ref(inputBAM), std::ref(readIdToOutput), std::ref(outputFiles), std::ref(mtx), start, end);
    }

    for (auto& t : threads) {
        t.join();
    }

    for (auto& kv : outputFiles) {
        bgzf_close(kv.second);
    }

    return 0;
}
