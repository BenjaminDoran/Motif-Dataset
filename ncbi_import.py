"""
Author: Benjamin Doran
Name: ncbi_import.py
Date: July 2017
Purpose: take organism + chromosome + posistion info ->
         return dataset of sequences
Output: .csv file (sequences.csv) columns: label, sequence
Data Source: NCBI Entrez -- https://www.ncbi.nlm.nih.gov/
"""

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import pandas as pd

class NCBIimporter:
    """ TODO """

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

    def __init__(self, user_email: str, in_file: str):
        """ Requires user_email and file of motif locations """
        # Setup
        Entrez.email = user_email
        self.locdat = pd.read_csv(in_file, sep=",")
        self.seqdat = pd.DataFrame()

    def org_to_id(self, org: str, chrom: str) -> str:
        """ get RefSeq or GI id from organism and chromosome """
        return self.id_switch[org](chrom)

    def _search_chromosome(self) -> bool:
        """ TODO """
        return False

    def _load_chromosome(self, iden: str):
        """ fetch fasta sequence from NCBI given id """
        with Entrez.efetch(db='nucleotide', id=iden, rettype="fasta") as handle:
            return Seq(SeqIO.read(handle, 'fasta'))

    def make_sequences(self) -> bool:
        """ TODO """
        # sort locdat by chromosome
        # for each chromosome:
        #   org_to_id
        #   _load_chromosome
        #   _search_chromosome
        #   _add_label_seq
        # export_sequences
        return False

    def export_sequences(self) -> bool:
        """ TODO """
        return False


if __name__ == "__main__":
    print("TODO")
