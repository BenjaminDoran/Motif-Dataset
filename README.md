# The Motif Extracter

## Details

This repo works as a glue script between JASPER db and NCBI. JASPER holds motifs with loci in the UCSC genomes. This script takes JASPER's loci for the motif instances from `.bed` files and extracts from NCBI genome sequences containing the motif instance. 

Due to the reasons for creating this script (See [Unimodality Based Clustering Applied to Motif Discovery](https://github.com/BenjaminDoran/motif-paper)) the motifs are randomly located along the returned sequences. See the `mstart` and `mstop` columns in `.csv` file for position labels.

### Example use:

```
python3.6 extract_sequences.py --email <your email> --bedfile http://jaspar.genereg.net/download/bed_files/MA0852.2.bed --outfile MA0852.2SEQS1k
```
The email arguement is required for interfacing with the NCBI API, so please by honest.

Using the perameters above, the script will output `MA0852.2SEQS1k.csv` and `MA0852.2SEQS1k.fasta` files. Both files contain the sequnces, but the `.csv` file also contains extra annotation for each sequence.

This script currently works for human motifs from the hg38 genome. Extending the script is possibe by adding more conversion rules to the `id_swticher.py` file. Pull requests to that end are encouraged and welcomed.

