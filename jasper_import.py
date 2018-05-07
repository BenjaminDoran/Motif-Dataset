"""
Author: Benjamin Doran
Name: jasper_import.py
Date: July 2017
Purpose: pull bed files and pfm files from Jasper database
Output: 2 .csv files (pfm.csv, motif-data.csv)
Data Source: Jasper database, http://jaspar.genereg.net
"""

import os.path
from urllib.request import urlretrieve
import tarfile
import re

import pandas as pd

def bed_files(members):
    for tarinfo in members:
        fn = os.path.splitext(tarinfo.name)
        if fn[1] == ".bed" and not "._" in fn[0]:
            yield tarinfo

def switch_org_name(orgs: pd.Series) -> pd.Series:
    """ Converts UCSC style org ref names to Latin Names """
    name_table = {'tair10': 'arabidopsis thaliana',
                  'at': 'arabidopsis thaliana',
                  'hg19': 'homo sapiens',
                  'hg38': 'homo sapiens',
                  'ce': 'caenorhabditis elegans',
                  'mm9': 'mus musculus',
                  'dm': 'drosophila melanogaster'}

    return orgs.apply(lambda x: name_table[x])


class JasperImporter:
    """ Saves Bed files from JasperDB to local tmp folder and
        exports human bed data to csv file.
    """
    # startup
    def __init__(self, url, outfile, createnew=True):
        """ Setup Importer """
        self.url = url
        self.outfile = outfile
        self.final_cols = ["motif-id", "organism", "genome", "chromosome", "start", "stop", "strand"]

        if createnew:
            self.datf = pd.DataFrame([], columns=["motif-id", "organism",
                                                  "chromosome", "start",
                                                  "stop", "strand"])
        elif os.path.isfile(outfile) and not createnew:
            self.load_csv()
        else:
            print('No data to load, acting in createnew mode')
            self.datf = pd.DataFrame([])
 
 
    def _expand_query(self, query: pd.Series) -> pd.DataFrame:
        """ Parses query formatted for UCSC Biologic database """
        # expand query
        parser = re.compile(r"^([a-z]{2,4}\d{,2})(?:\_)(chr\d{,2}[A-Z]{,3})(?:[^\:]*):(\d+)-(\d+)(?:\()([\-\+])(?:\))")
        results = query.apply(parser.findall).apply(lambda r: list(sum(r, ()))).apply(pd.Series)
        return results

    # functions
    def load_csv(self) -> bool:
        """ reads in data.csv """
        self.datf = pd.read_csv(self.outfile, sep=",")
        return True

    def download_jasper_bed_files(self):
        """ download and extract files """
        filename, _ = urlretrieve(self.url)
        tar = tarfile.open(filename)
        tar.extractall("./tmp/", members=bed_files(tar))

    def load_jasper_bed_files(self):
        """ loads jasper bed files which contain
            jasper motif-id and query of motif
        """
        directory = "./tmp/bed_files/"

        for filename in sorted(os.listdir(directory)):

            # ignore files that do not follow format
            if filename in ("MA0548.1.bed", "MA0555.1.bed", "MA0563.1.bed",
                            "MA1417.1.bed", "MA1416.1.bed", "MA0036.3.bed", 
                            "MA1198.1.bed", "MA1414.1.bed", "MA1415.1.bed"):
                continue            

            print(f"Loading {filename}")

            # read csv
            locus_file = pd.read_csv(os.path.join(directory, filename),
                                     sep="\t",
                                     header=None)

            # 3: query -> organism, chormosome, start, stop, strand
            result_df = self._expand_query(locus_file.iloc[:,3])

            # name columns
            result_df.columns = self.final_cols[2:]

            # motif id is part of file name
            result_df["motif-id"] = filename[:-4]

            # switch org names to Latin
            result_df["organism"] = switch_org_name(result_df["genome"])

            # reorder columns
            result_df = result_df[self.final_cols]

            # append new columns
            self.datf = self.datf.append(result_df, ignore_index=True)

    def export(self) -> bool:
        """ save data to csv """
        try:
            self.datf.to_csv(self.outfile, sep=",", header=True, index=False)
            return True
        except:
            return False
