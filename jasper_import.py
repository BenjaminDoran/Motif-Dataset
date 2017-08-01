"""
Author: Benjamin Doran
Name: jasper_import.py
Date: July 2017
Purpose: pull bed files and pfm files from Jasper database
Output: 2 .csv files (pfm.csv, motif-data.csv)
Data Source: Jasper database, http://jaspar.genereg.net
"""
import pandas as pd
import os.path
import re
import requests
from bs4 import BeautifulSoup

class JasperImporter:
    # startup
    def __init__(self, url = "http://jaspar.genereg.net/html/DOWNLOAD/", bedfilename = "motif-data.csv"):
        """ Setup Importer """
        self.url = url
        self.bedfilename = bedfilename

        if os.path.isfile(bedfilename):
            self.load_csv()
        else:
            self.df = None

    # functions
    def load_csv(self) -> bool:
        """ reads in data.csv """
        try:
            self.df = pd.read_csv(self.bedfilename, sep = ",")
            return True
        except:
            return False

    def load_jasper_bed_files(self) -> bool:
        """loads jasper bed files which contain jasper motif id
        and locus of motif, writes information to .csv file
        """
        if self.df != None:
            print("bed files already loaded from csv file")
            return False
        else:
            bed_url = self.url + "bed_files/"
            self.df = pd.DataFrame([])
            final_cols = ["motif", "organism",
                          "chromosome", "start",
                           "stop", "query"]

            # enter bed directory
            response = requests.get(bed_url)
            soup = BeautifulSoup(response.content, "html.parser")
            # all links to files
            bed_motifs = soup.find_all("a", text = re.compile("MA."))
            # for link to file
            for a in bed_motifs:
                locus_file = pd.read_csv(bed_url + a.text, sep = "\t", header = None)
                try: # files have either 4 or 6 columns
                    locus_file.columns = ["chromosome", "start", "stop", "query", "frame", "strand"]
                except ValueError:
                    locus_file.columns = ["chromosome", "start", "stop", "query"]
                # parse query for organism name (and potentially strand)
                locus_file["organism"] = locus_file.loc[:,"query"].str.split(r"_", expand = True).loc[:,0]
                # motif id is part of file name
                locus_file["motif"] = a.text[:-4]
                self.df = self.df.append(pd.DataFrame(locus_file.loc[:,final_cols]), ignore_index=True)
            self.df.columns = final_cols
        return True

    def load_jasper_pfm_files(self) -> bool:
        """ TODO """
        return False

    def export(self) -> bool:
        """ TODO """
        self.df.to_csv(self.bedfilename, sep = ",", header = True, index = False)
