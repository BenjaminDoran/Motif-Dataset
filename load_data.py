"""
Author: Benjamin Doran
Name: load_data.py
Date: July 2017
Purpose: make a usable motif finding dataset
Output: 1 .csv file
Data Source: "motif-data.csv", NCBI-Entrez <>
"""
from jasper_import import JasperImporter
from Bio import Entrez

# main
def make_dataset(newfile:str, fromfile:str = "motif-data.csv",
                 indent:int = 50, **kwargs):
    """ Makes new train/test set from NCBI-Entrez
    Outputs a single csv file"""
    print("TODO")

if __name__ == "__main__":
    # make_dataset()
    print("TODO")
