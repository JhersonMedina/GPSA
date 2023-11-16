#include <iostream>
#include "Utilities/Aligner.h"
#include "Entities/AlignmentResults.h"
#include "Utilities/Parser.h"

using namespace std;

int main(int argv, char** argc) {
    if (argv != 5) {
        cout << "Enter a single FASTA file with both sequences" << endl;
        exit(0);
    }

    string fileName = argc[1];
    vector<string> sequences = getSequences(fileName);
    double matchCost = stod(argc[2]), missMatchCost = stod(argc[3]), indelCost = stod(argc[4]);

    auto alingmentResults = getBestAlingment({sequences[0], sequences[1], matchCost, missMatchCost, indelCost});

    cout << "Score: " << alingmentResults.score << endl;
    cout << "Alignments: " << endl << alingmentResults.firstSequence << endl << alingmentResults.secondSequence << endl;
    return 0;
}

