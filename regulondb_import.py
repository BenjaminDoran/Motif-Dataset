"""
Author: Benjamin Doran
Name: regulondb_import.py
Date: Aug 2017
Purpose: Pull data from RegulonDB database
Output: 2 csv files
Data Source: RegulonDB <http://regulondb.ccg.unam.mx/>
"""
import pandas as pd
import os
import re

DATAFILE = "./motif-data.csv"

BASE_URL = "http://regulondb.ccg.unam.mx/menu/download/datasets/files/BindingSiteSet.txt"

IN_COLUMNS = ["motif-id", "tf-name", "site-id", "start",
 		      "stop", "strand", "tf-gen-id", "tf-unit",
		      "gene-ex", "promoter", "center", "sequence",
		      "evidence", "evi-level"]

OUT_COLUMNS = ["motif-id", "organism", "chromosome",
			   "start", "stop", "strand"]

class RegulonDBImporter:
    def __init__(self, url=BASE_URL, outfile=DATAFILE, createnew=False):
        """ Setup Importer """
        self.url = url
        self.outfile = outfile

        if createnew:
            self.datf = pd.DataFrame([], columns=OUT_COLUMNS)
        elif os.path.isfile(outfile) and not createnew:
            self.datf = pd.read_csv(self.outfile, sep=",")
        else:
            print('No data to load, acting in createnew mode')
            self.datf = pd.DataFrame([], columns=OUT_COLUMNS)  

    def load_regulon_data(self) -> bool:
        """ read data from ReglonDB """
        # read data from ReglonDB
        dat_in = pd.read_csv(self.url,
                             header=None,
                             names=IN_COLUMNS,
                             comment="#",
                             sep="\t").dropna()

        # set columns not contained in infile
        dat_in['organism'] = "escherichia coli"
        dat_in['chromosome'] = "None"

        # modify strand to match jasper identifier
        dat_in["strand"].replace("forward", '+', inplace=True)
        dat_in["strand"].replace("reverse", '-', inplace=True)

        # append to our data
        self.datf = self.datf.append(dat_in[OUT_COLUMNS],
        				             ignore_index=True)

        if self.datf.shape[0] < 3021:
            return False
        return True

    def export(self) -> bool:
        """ export data to csv file """
        try:
            self.datf.to_csv(self.outfile, sep=",", header=True, index=False)
            return True
        except:
            return False