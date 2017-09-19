"""
Author: Benjamin Doran
Name: ncbi_import.py
Date: July 2017
Purpose: take organism + chromosome + posistion info ->
         return dataset of sequences
Output: .csv file (sequences.csv) columns: label, sequence
Data Source: NCBI Entrez -- https://www.ncbi.nlm.nih.gov/
"""
from os import path, makedirs
from urllib.error import HTTPError
from time import sleep
from random import randint
from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna
import pandas as pd

class NCBIimporter:
    """ TODO """
    def __init__(self, user_email: str, in_file: str, indent: int=0, rand_lim: int=0):
        """ Requires user_email and file of motif locations """
        # Setup
        Entrez.email = user_email
        self.locdat = pd.read_csv(in_file, sep=",")
        # sort locdat by chromosome
        self.chromosomes = self.locdat[['organism', 'chromosome']] \
                           .groupby(['organism', 'chromosome']) \
                           .count()

        self.seqdat = pd.DataFrame()
        self.indent = indent
        self.rand_lim = rand_lim


        # make sequences directory
        if not path.isdir('sequences'):
            makedirs('sequences')

    id_switch = {
        # conversion table for organism + chromosome to RefSeq or GI id
        # it does not seem to feasible to do this with
        # biopython searching programmatically
        # thale cress (plant)
        'arabidopsis thaliana':
            lambda chrom:
            "CP00268{:01d}".format(int(chrom)+3) if
            chrom.isnumeric() and (0 < int(chrom) < 6) else
            None,
        # c. elegans (worm)
        'caenorhabditis elegans':
            lambda chrom:
            "NC_003279" if chrom == "I" else
            "NC_003280" if chrom == "II" else
            "NC_003281" if chrom == "III" else
            "NC_003282" if chrom == "IV" else
            "NC_003283" if chrom == "V" else
            "NC_003284" if chrom == "X" else
            None,
        # fruit fly
        'drosophila melanogaster':
            lambda chrom:
            "NC_004354" if chrom == "X" else
            "NT_033778" if chrom == "2R" else
            "NT_033777" if chrom == "3R" else
            "NT_033779" if chrom == "2L" else
            "NT_037436" if chrom == "3L" else
            "NC_004353" if chrom == "4" else
            None,
        # e. coli (bacteria) - no chromosomes
        'escherichia coli': lambda chrom=None: "NC_000913",
        # human
        'homo sapiens':
            lambda chrom:
            "NC_0000{:02d}".format(int(chrom)) if
            chrom.isnumeric() and (0 < int(chrom) < 23) else
            "NC_000023" if chrom == "X" else
            "NC_000024" if chrom == "Y" else
            None,
        # house mouse
        'mus musculus':
            lambda chrom:
            "NC_0000{:02d}".format(int(chrom)+66) if
            chrom.isnumeric() and (0 < int(chrom) < 20) else
            "NC_000086" if chrom == "X" else
            "NC_000087" if chrom == "Y" else
            None
    }

    def _org_to_id(self, org: str, chrom: str) -> str:
        """ get RefSeq or GI id from organism and chromosome """
        return self.id_switch[org](chrom)

    def _search_chromosome(self, 
                           chrom: SeqIO.SeqRecord,
                           start: pd.Series,
                           stop: pd.Series, 
                           strand: pd.Series) -> tuple:
        """ returns (motif, mstart, mstop) """
        # get randomized indents within a set range
        rstart = start.apply(lambda d: d - (self.indent - randint(0, self.rand_lim)))
        rstop = stop.apply(lambda d: d + (self.indent - randint(0, self.rand_lim)))

        # length of motif
        len_motifs = (stop - start) + 1 # plus 1 because 0 vs. 1 indexing

        # select motif +/- indents
        motifs = pd.concat([rstart, rstop, strand], axis=1)\
                   .apply(lambda r: str(chrom[r["start"]-1:r["stop"]].seq)
                                    if r["strand"] == "+" else
                                    str(chrom[r["start"]-1:r["stop"]].seq.complement()),
                                    axis=1)
        
        # return motif, start index from selected sequence, and
        # stop index from selected sequence
        return (motifs, start - rstart, start - rstart + len_motifs)

    def _fetch_chromosome(self, iden: str):
        """ fetch fasta sequence from NCBI given id """
        try:
            with Entrez.efetch(db='nucleotide', id=iden, rettype='fasta') as handle:
                rec = SeqIO.read(handle, 'fasta')
                rec.Alphabet = generic_dna
                return rec
        except HTTPError:
            return None

    def save_sequences(self) -> int:
        """save chromosome sequences as fasta files in sequences folder
        returns: number of sequences successfully saved
        """
        saved = 0
        # for each chromosome:
        for org, chrom in self.chromosomes.index:

            # get id for looking in NCBI
            refseq_id = self._org_to_id(org, chrom.strip('chr'))

            # fetch sequence from NCBI
            big_sequence = self._fetch_chromosome(refseq_id)
            if not big_sequence:
                print('error: in {} {} refseq_id {}'.format(org, chrom, refseq_id))

            else:
                # save sequence
                outfilename = "sequences/" + org.replace(" ", "-") +\
                              chrom.strip('chr') + ".fasta"
                SeqIO.write(big_sequence, outfilename, 'fasta')

                # update status
                saved += 1
                print('saved', big_sequence.description)

                # ensure requests do not go to NCBI too fast
                sleep(1)
        return saved

    def load_sequences(self) -> int:
        """ TODO """
        loaded = 0
        dirname = "sequences/"
        for org, chrom in self.chromosomes.index:

            chrom_file = dirname + org.replace(" ", "-") +\
                         chrom.strip("chr") + ".fasta"
            chrom_record = SeqIO.read(chrom_file, 'fasta')

            # get rows for organism and chromosome
            startstops = self.locdat.loc[(self.locdat['organism'] == org) &
                                         (self.locdat['chromosome'] == chrom)]
            # retrive motif + indent
            motifs, mstarts, mstops = self._search_chromosome(chrom_record,
                                                              startstops["start"],
                                                              startstops["stop"],
                                                              startstops["strand"])
            rows = pd.concat([startstops, motifs, mstarts, mstops], axis=1)
            rows.columns = ["motif-id", "organism", "chromosome", "start",
                           "stop", "strand", "seq", "mstart", "mstop"]
            self.seqdat = self.seqdat.append(rows, ignore_index=True)
            loaded += 1

        return loaded


if __name__ == "__main__":
    pass
