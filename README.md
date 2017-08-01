# The Motif Dataset

## Purpose
This dataset is meant to work as a starting point for comparison of motif finding algorithms. As the length of sequence can affect the performance of the algorithms, this dataset provides the start and stop points of the motifs such that the user may pull the sequence from NCBI-Entrez +/- some indent to vary sequence length.

## Schema
- **pfm.csv:**
    - contains: Point Frequency Matrices for motifs
    - columns: motif-id, organism, pfm
    - rows: ?
- **motif-data.csv:**
    - contains: location data of motifs
    - columns: motif, organism, chromosome, start, stop, query
    - rows: ?

## How To Start
TODO
