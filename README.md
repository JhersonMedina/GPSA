# Global Pairwise Sequence Aligner

This is the C++ implementation of de well-known global sequence aligner algorithm using dynamic programing.

# Usage

To use you must have C++ installed on your computer.
Simply compile with 

```
g++ main -o main
```
This will compile an create the execution file. Then align two sequences run.

```
./main [sequences] [matchCost] [missMatchCost] [indelCost]
```

where:
- ```Sequences```: is fasta file containing exactly two sequences. If there are more it will only consider the first two.
- ```matchCost```: is the associated cost for match.
- ```misMatchCost```: is the associated cost for a mismatch.
- ```indelCost```: is the associated cost for an insertion or deletion (indel).

The output looks like this:
```
Score: -1
Alignments:
GCATGCG
GATTACA
```

Where score is the overall score of the alignment and the alignments are the resulting alignment for each sequence.

# Further work and next steps

Right now I'm working on making this a CLI application and in future itereations I expect to inlcude PAM and BLOSUM matrixes as part of the cost.

Beyond that, extend this to multiple sequence alignments although that must require different approaches.
