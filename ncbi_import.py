"""
Author: Benjamin Doran
Name: ncbi_import.py
Date: July 2017
Purpose: take organism + chromosome + posistion info -> return dataset of sequences
Output: .csv file (sequences.csv)
Data Source: NCBI Entrez -- https://www.ncbi.nlm.nih.gov/
"""

from Bio import Entrez, SeqIO
import pandas as pd

""" handle = Entrez.esearch(db='nucleotide', term='{}[Orgn] AND chromosome {}') -> id 
    record = Entrez.read(handle)
    record["IdList"]
"""
""" Entrez.efetch(db='nucleotide',
                  id='returned id',
                  rettype="gb",
                  retmode='text')
"""

class NCBI_Importer:
    # conversion table for organism + chromosome to Accesion number
    # at this point it does not seem to feasible to do this with biopython
    self.org_to_id = {
        'arabidopsis thaliana': lambda n: \
                                    "CP00268{:01d}".format(int(n)+3) if \
                                        (n.isnumeric()) & (0 < int(n) < 6) else \
                                    None,
        'caenorhabditis elegans': lambda n: \
                                    "NC_003279" if n == "I"   else \
                                    "NC_003280" if n == "II"  else \
                                    "NC_003281" if n == "III" else \
                                    "NC_003282" if n == "IV"  else \
                                    "NC_003283" if n == "V"   else \
                                    "NC_003284" if n == "X"   else \
                                    None,
        'drosophila melanogaster': lambda n: \
                                    "NC_004354" if n == "X"  else \
                                    "NT_033778" if n == "2R" else \
                                    "NT_033777" if n == "3R" else \
                                    "NT_033779" if n == "2L" else \
                                    "NT_037436" if n == "3L" else \
                                    "NC_004353" if n == "4"  else \
                                    None,
        'escherichia coli': lambda n = None: "NC_000913",
        'homo sapiens': lambda n: \
                        "NC_0000{:02d}".format(int(n)) if \
                            (n.isnumeric()) & (0 < int(n) < 23) else \
                        "NC_000023" if n == "X" else \
                        "NC_000024" if n == "Y" else \
                        None,
        'mus musculus': lambda n: \
                        "NC_0000{:02d}".format(int(n)+66) if \
                            (n.isnumeric()) & (0 < int(n) < 20) else \
                        "NC_000086" if n == "X" else \
                        "NC_000087" if n == "Y" else \
                        None
    }

    def __init__(self, user_email:str, in_file:str):
        Entrez.email = user_email
        self.df = pd.read_csv(in_file, sep=",")

    def search_chromosomes(self) -> bool:
        """ TODO """
        return False

    def load_chromosomes(self) -> bool:
        """ TODO """
        return False

    def make_sequences(self) -> bool:
        """ TODO """
        return False

    def export_sequences(self) -> bool:
        """ TODO """
        return False


if __name__ == "__main__":
    print("TODO")



