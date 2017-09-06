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

def import_jasper_loc_data() -> bool:
    """ fetch motif locations from jasper database """
    jas_imp = JasperImporter("http://jaspar.genereg.net/html/DOWNLOAD/",
                             "motif-data.csv")

    worked = jas_imp.load_jasper_bed_files()
    if not worked:
        print("error loading jasper bedfiles")
        return False

    worked = jas_imp.export()
    if not worked:
        print("unable to export location data to csv")
        return False
    return True


def main():
    """ TODO """
    if not import_jasper_loc_data():
        sys.exit(1)

if __name__ == "__main__":
    main()
