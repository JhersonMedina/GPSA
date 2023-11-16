#include <fstream>

vector<string> getSequences(string fileName) {
    ifstream file(fileName);
    vector<string> sequences;
    string line;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') sequences.push_back("");
        else sequences.back() += line;
    }

    file.close();
    return sequences;
}
