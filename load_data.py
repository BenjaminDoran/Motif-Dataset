"""
Author: Benjamin Doran
Date: July 2017
Purpose: make a usable motif finding dataset
Output: 1 .csv file
"""

# imports
import requests
from bs4 import BeautifulSoup
import pandas as pd
import re

# functions
def load_jasper_bed_files(str:url) -> bool:
    """ TODO """
    return False

def load_jasper_pfm_files(str:url) -> bool:
    """ TODO """
    return False


# main
def main():
    """ creates a csv file with motif id, locus, and point freqency matrix
    as pulled from the JASPER database (http://jaspar.genereg.net)
    Outputs a single csv file"""
    print("TODO")

if __name__ == "__main__":
    main()
