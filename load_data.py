"""
Author: Benjamin Doran
Name: load_data.py
Date: July 2017
Purpose: make a usable motif finding dataset
Output:
Data Source: "motif-data.csv", NCBI-Entrez <>
"""
import sys
import progressbar
from jasper_import import JasperImporter
from regulondb_import import RegulonDBImporter
from ncbi_import import NCBIimporter

JASPER_URL = "http://jaspar2016.genereg.net/html/DOWNLOAD/"
REGULON_URL = "http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt"

def import_jasper_loc_data(url, outfile, createnew=False) -> bool:
    """ fetch motif locations from jasper database """
    jas_imp = JasperImporter(url, outfile)

    if not jas_imp.load_jasper_bed_files():
        print("error loading jasper bedfiles")
        return False

    if not jas_imp.export():
        print("unable to export location data to csv")
        return False
    return True

def import_regulon_loc_data(url, outfile, createnew=False):
    reg_imp = RegulonDBImporter(url, outfile)

    if not reg_imp.load_regulon_data():
        print("unable to load Regulon data")
        return False

    if not reg_imp.export():
        print('unable to export regulon location data to csv')
        return False

    return True


def ncbi_save_sequences(email, infile):
    """ fetch NCBI chromosome sequences and save in sequences folder """
    ncbi_imp = NCBIimporter(email, infile)
    return ncbi_imp.save_sequences()


def load():
    """ TODO """
    bar = progressbar.ProgressBar(maxval=20,widgets=[progressbar.Bar('=', '[', ']'), ' '])

    # get user email
    email = ""
    while not email:
        email = input("Please enter your email address: ")

    # confirm location file name
    locfile = input("pulling from 'motif-data.csv'\n" + 
                    "if this is not correct please type filename and press enter.\n"+
                    "else press enter:\n")
    if not locfile:
        locfile = "motif-data.csv"

    # fetch jasper motif location data
    print("loading jasper bed files")
    bar.start()
    if not import_jasper_loc_data(JASPER_URL, locfile, createnew=True):
        sys.exit(1)
    bar.finish()
    print("loaded jasper location data")

    print("loading regulonDB location data")
    bar.start()
    if not import_regulon_loc_data(REGULON_URL, locfile):
        sys.exit(2)
    bar.finish()
    print("loaded regulonDB location data")

    # From location data fetch and save chromosome sequences
    print("saving ncbi sequences for all organisms")
    bar.start()
    if ncbi_save_sequences(email, locfile) < 1:
        sys.exit(3)
    bar.finish()
    print("saved ncbi sequences for all organisms")

    print("all tasks completed successfully!")

if __name__ == "__main__":
    load()
