"""
Author: Benjamin Doran
Name: load_data.py
Date: July 2017
Purpose: make a usable motif finding dataset
Output:
Data Source: "motif-data.csv", NCBI-Entrez <>
"""
import sys
from jasper_import import JasperImporter
from ncbi_import import NCBIimporter

def import_jasper_loc_data(url, outfile) -> bool:
    """ fetch motif locations from jasper database """
    jas_imp = JasperImporter(url, outfile)

    worked = jas_imp.load_jasper_bed_files()
    if not worked:
        print("error loading jasper bedfiles")
        return False

    worked = jas_imp.export()
    if not worked:
        print("unable to export location data to csv")
        return False
    return True

def ncbi_save_sequences(email, infile):
    """ fetch NCBI chromosome sequences and save in sequences folder """
    ncbi_imp = NCBIimporter(email, infile)
    return ncbi_imp.save_sequences()



def main():
    """ TODO """
    # get user email
    email = ""
    while not email:
        email = input("Please enter your email address: ")

    # confirm location file name
    locfile = input("pulling from 'motif-data.csv' if this is not correct \
                    please type the file name and press enter,\
                    else just press enter:\n")
    if not locfile:
        locfile = "motif-data.csv"

    # hard code jaspers url
    jasper_url = "http://jaspar.genereg.net/html/DOWNLOAD/"

    # fetch jasper motif location data
    if not import_jasper_loc_data(jasper_url, locfile):
        sys.exit(1)

    # From location data fetch and save chromosome sequences
    if ncbi_save_sequences(email, locfile) < 1:
        sys.exit(2)

    # from saved sequences make final data file
    # columns:
    # motif-id: unique id of sequence
    # sequence: actual sequence
    # chromosome-index: in form of <start>:<stop><strand>
    # mstart: index in sequence where motif starts
    # mstop: index in sequence where motif stops

if __name__ == "__main__":
    main()
