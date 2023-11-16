#include "AlignmentInformation.h"
#include "AlignmentResults.h"
#include <vector>

#pragma once

#define PAIR make_pair

using namespace std;

const double WORST_SCORE = 1e18;

double matchCost;
double missMatchCost;
double indelCost;

string firstSequence;
string secondSequence;

vector<vector<double>> bestAlingment;
vector<vector<bool>> memoized;

double findBestAlignmentScore(int, int);
pair<string, string> buildBestAlignment(int, int);

double getMatchCost(char firstNucleotide, char secondNucleotide) {
    return firstNucleotide == secondNucleotide ? matchCost : missMatchCost;
}

AlignmentResults getBestAlingment(AlignmentInformation alingmentInformation) {
    firstSequence = alingmentInformation.firstSequence;
    secondSequence = alingmentInformation.secondSequence;


    matchCost = alingmentInformation.matchCost;
    missMatchCost = alingmentInformation.missMatchCost;

    indelCost = alingmentInformation.indelCost;

    bestAlingment.assign(firstSequence.size()+1, vector<double> (secondSequence.size()+1, 0));
    memoized.assign(firstSequence.size()+1, vector<bool> (secondSequence.size()+1, 0));

    pair<string, string> alignedSequences = buildBestAlignment(0, 0);
    return { alignedSequences.first, alignedSequences.second, findBestAlignmentScore(0, 0) };
}

pair<string, string> buildBestAlignment(int firstSequencePosition, int secondSequencePosition) {
    if (firstSequencePosition >= (int)firstSequence.size() && secondSequencePosition >= (int)secondSequence.size()) {
        return PAIR("", "");
    }

    double bestScore = findBestAlignmentScore(firstSequencePosition, secondSequencePosition);

    if (firstSequencePosition < (int)firstSequence.size() &&
        indelCost + findBestAlignmentScore(firstSequencePosition+1, secondSequencePosition) == bestScore) {
            pair<string, string> alingment = buildBestAlignment(firstSequencePosition+1, secondSequencePosition);
            return PAIR("_" + alingment.first, alingment.second);
    }

    if (secondSequencePosition < (int)secondSequence.size() &&
        indelCost + findBestAlignmentScore(firstSequencePosition, secondSequencePosition+1) == bestScore) {
            pair<string, string> alingment = buildBestAlignment(firstSequencePosition, secondSequencePosition+1);
            return PAIR(alingment.first, "_" + alingment.second);
    }

    pair<string, string> alingment = buildBestAlignment(firstSequencePosition+1, secondSequencePosition+1);
    return PAIR(firstSequence[firstSequencePosition] + alingment.first, secondSequence[secondSequencePosition] + alingment.second);

}

double findBestAlignmentScore (int firstSequencePosition, int secondSequencePosition) {
    if (firstSequencePosition >= (int)firstSequence.size() && secondSequencePosition >= (int)secondSequence.size()) {
        return 0;
    }

    if (memoized[firstSequencePosition][secondSequencePosition]) {
        return bestAlingment[firstSequencePosition][secondSequencePosition];
    }

    bestAlingment[firstSequencePosition][secondSequencePosition] = -WORST_SCORE;

    if (firstSequencePosition < (int)firstSequence.size()) {
        bestAlingment[firstSequencePosition][secondSequencePosition] =
            max(bestAlingment[firstSequencePosition][secondSequencePosition],
                indelCost + findBestAlignmentScore(firstSequencePosition+1, secondSequencePosition));
    }

    if (secondSequencePosition < (int)secondSequence.size()) {
        bestAlingment[firstSequencePosition][secondSequencePosition] =
            max(bestAlingment[firstSequencePosition][secondSequencePosition],
                indelCost + findBestAlignmentScore(firstSequencePosition, secondSequencePosition+1));
    }

    if (firstSequencePosition < (int)firstSequence.size() && secondSequencePosition < (int)secondSequence.size()) {
        bestAlingment[firstSequencePosition][secondSequencePosition] =
            max(bestAlingment[firstSequencePosition][secondSequencePosition],
                getMatchCost(firstSequence[firstSequencePosition], secondSequence[secondSequencePosition]) +
                    findBestAlignmentScore(firstSequencePosition+1, secondSequencePosition+1));
    }

    memoized[firstSequencePosition][secondSequencePosition] = true;
    return bestAlingment[firstSequencePosition][secondSequencePosition];
}
