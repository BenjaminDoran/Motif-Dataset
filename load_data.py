"""
Author: Benjamin Doran
Name: load_data.py
Date: July 2017
Purpose: make a usable motif finding dataset
Output:
Data Source: "motif-data.csv", NCBI-Entrez
"""
import sys, os

from jasper_import import JasperImporter
from ncbi_import import NCBIimporter

JASPER_URL = "http://jaspar.genereg.net/download/bed_files.tar.gz"
REGULON_URL = "http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt"

def import_jasper_loc_data(url, outfile, createnew=False) -> bool:
    """ fetch motif locations from jasper database """
    jas_imp = JasperImporter(url, outfile)

    if not os.path.exists("./tmp/bed_files/"):
        jas_imp.download_jasper_bed_files()
    jas_imp.load_jasper_bed_files()
    jas_imp.export()


def ncbi_save_sequences(email, infile):
    """ fetch NCBI chromosome sequences and save in sequences folder """
    ncbi_imp = NCBIimporter(email, infile)
    return ncbi_imp.save_sequences()


def load(email, locfile):
    """ Loads Jasper Bed files, NCBI human genome, and extracts sequences. """
    # fetch jasper motif location data
    print("loading jasper bed files")
    import_jasper_loc_data(JASPER_URL, locfile, createnew=True)
    print("loaded jasper location data")

    # # From location data fetch and save chromosome sequences
    # print("saving ncbi sequences for all organisms")
    # ncbi_save_sequences(email, locfile)
    # print("saved ncbi sequences for all organisms")

    print("all tasks completed successfully!")

if __name__ == "__main__":
    # get user email
    email = ""
    while not email:
        email = input("Please enter your email address: ")

    # confirm location file name
    locfile = input("\nPulling from 'motif-data.csv'\n" + 
                    "If this is correct press ENTER (else type filename and press ENTER): ")
    if not locfile:
        locfile = "motif-data.csv"
    load(email, locfile)
