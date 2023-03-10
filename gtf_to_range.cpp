#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <htslib/sam.h>

using namespace std;

// GTF��ÿ�м�¼��Ϣ
struct GtfRecord {
    string seqname;  // Ⱦɫ������
    string source;   // ע������Դ
    string feature;  // ע�����ͣ�����gene��transcript�ȣ�
    uint64_t start;       // ��ʼλ��
    uint64_t end;         // ��ֹλ��
    string score;    // ��ֵ
    string strand;     // Ⱦɫ����������
    string frame;      // ������
    map<string, string> attributes;  // ����
};

// ��GTF�ļ��ж�ȡע����Ϣ
void readGtf(const string& gtfFile, vector<GtfRecord>& records) {
    ifstream inFile(gtfFile);
    if (!inFile) {
        cerr << "Error: failed to open " << gtfFile << endl;
        exit(EXIT_FAILURE);
    }
    
    string line;
    while (getline(inFile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        //cout << line << endl;
        istringstream ss(line);
        GtfRecord record;
        getline(ss, record.seqname, '\t');
        //cout << record.seqname << endl;
        getline(ss, record.source, '\t');
        getline(ss, record.feature, '\t');
        string field;
        getline(ss, field, '\t');
        record.start = stoull(field);
        getline(ss, field, '\t');
        record.end = stoull(field);
        getline(ss, record.score, '\t');
        getline(ss, record.strand, '\t');
        getline(ss, record.frame, '\t');
        string attrStr;
        getline(ss, attrStr);
        //cout << record.seqname << "\t" <<record.strand << "\t" << record.frame << endl;
        attrStr = "\t" + attrStr;
        istringstream attrSS(attrStr);
        string attr;
        

        while (getline(attrSS, attr, ';')) {
            //cout << attr << endl;
            size_t pos = attr.find_first_of(" \t");
            if (pos != string::npos) {
                attr.erase(pos, 1);
            }
            pos = attr.find_last_not_of(" \t");
            if (pos != string::npos) {
                attr.erase(pos + 1);
            }
            pos = attr.find_first_of("\"\';");
            if (pos != string::npos) {
                attr.erase(pos, 1);
            }
            pos = attr.find_last_not_of("\"\';");
            if (pos != string::npos) {
                attr.erase(pos + 1);
            }
            if (!attr.empty()) {
                size_t sepPos = attr.find_first_of(" \t");
                if (sepPos != string::npos) {
                    string key = attr.substr(0, sepPos);
                    string value = attr.substr(sepPos + 1);
                    //cout << key << "\t" << value << endl;
                    record.attributes[key] = value;
                }
            }
        }

        records.emplace_back(record);
    }
    inFile.close();
}

// ����GTF��Ϣ�����������Ƶ�λ�÷�Χ��ӳ��
map <string, map<string, pair<uint64_t, uint64_t>>> buildGeneRanges(const vector<GtfRecord>& records, string outfile) {
    ofstream outputFile(outfile);
    //map<string, pair<uint64_t, uint64_t>> geneRanges;
    map <string, map<string, pair<uint64_t, uint64_t>>> seqGeneRanges;
    for (const auto& record : records) {
        if (record.feature == "gene") {
            string geneName = record.attributes.at("gene_id");
            string seqname = record.seqname;
            if (seqGeneRanges[seqname].count(geneName) == 0) {
                seqGeneRanges[seqname][geneName] = make_pair(record.start, record.end);
            }
            else {
                auto& range = seqGeneRanges[seqname][geneName];
                range.first = min(range.first, record.start);
                range.second = max(range.second, record.end);
            }
        }
    }
    for (const auto& seqgenes : seqGeneRanges) {
        string seqname = seqgenes.first;
        auto genes = seqgenes.second;
        for (const auto& gene : genes)
        {
            if (gene.first.empty()) {
                continue;
            }
            outputFile << seqname << "\t" << gene.first << "\t" << gene.second.first << "\t" << gene.second.second << endl;
        }
        
    }

    return seqGeneRanges;
}

int main(int argc, char* argv[]) {
    string gtf = argv[1];
    string outfile = argv[2];

    vector<GtfRecord> records;
    cout << "load gtf " << endl;
    readGtf(gtf, records);
    cout << "build range" << endl;
    auto geneRanges = buildGeneRanges(records, outfile);
    return 0;

}