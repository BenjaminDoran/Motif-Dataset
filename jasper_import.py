"""
Author: Benjamin Doran
Name: jasper_import.py
Date: July 2017
Purpose: pull bed files and pfm files from Jasper database
Output: 2 .csv files (pfm.csv, motif-data.csv)
Data Source: Jasper database, http://jaspar.genereg.net
"""

import os.path
import re
import requests
import pandas as pd
from bs4 import BeautifulSoup

BASE_URL = "http://jaspar.genereg.net/html/DOWNLOAD/"
OUTFILE = "./motif-data.csv"


def switch_org_name(orgs: pd.Series) -> pd.Series:
    """ Converts UCSC style org ref names to Latin Names """
    name_table = {'tair10': 'arabidopsis thaliana',
                  'at': 'arabidopsis thaliana',
                  'hg19': 'homo sapiens',
                  'ce': 'caenorhabditis elegans',
                  'mm9': 'mus musculus',
                  'dm': 'drosophila melanogaster'}

    return orgs.apply(lambda x: name_table[x])


def expand_query(query: pd.Series) -> pd.DataFrame:
    """ Parses query formatted for UCSC Biologic database """
    # expand query
    parse_query = re.compile(r"[a-zA-Z0-9]+|\([\+\-]\)")
    results = query.apply(parse_query.findall).apply(pd.Series)

    # get strand without "()"
    restrip = re.compile(r"[\+\-]")
    results[4] = results[4].apply(lambda y: restrip.search(y).group(0))
    return results


class JasperImporter:
    """ TODO """
    # startup
    def __init__(self, url=BASE_URL, outfile=OUTFILE, createnew=True):
        """ Setup Importer """
        self.url = url
        self.outfile = outfile

        if createnew:
            self.datf = pd.DataFrame([], columns=["motif-id", "organism",
                                                  "chromosome", "start",
                                                  "stop", "strand"])
        elif os.path.isfile(outfile) and not createnew:
            self.load_csv()
        else:
            print('No data to load, acting in createnew mode')
            self.datf = pd.DataFrame([])

    # functions
    def load_csv(self) -> bool:
        """ reads in data.csv """
        self.datf = pd.read_csv(self.outfile, sep=",")
        return True

    def load_jasper_bed_files(self) -> bool:
        """
        loads jasper bed files which contain
        jasper motif-id and query of motif
        """
        bed_url = self.url + "bed_files/"
        final_cols = ["motif-id", "organism", "chromosome",
                      "start", "stop", "strand"]

        # enter bed directory
        response = requests.get(bed_url)
        soup = BeautifulSoup(response.content, "html.parser")
        # all links to files
        bed_motifs = soup.find_all("a", text=re.compile("MA."))

        # for link to file
        for link in bed_motifs:
            # ignore files that do not follow format (-400 entries)
            if link.text in ("MA0548.1.bed", "MA0555.1.bed", "MA0563.1.bed"):
                continue

            # read csv
            locus_file = pd.read_csv(bed_url+link.text,
                                     sep="\t",
                                     header=None)

            # query is column 3 (0 indexed)
            locus_file.rename(columns={3: "query"}, inplace=True)

            # query -> organism, chormosome, start, stop, strand
            result_df = expand_query(locus_file["query"])

            # name columns
            result_df.columns = final_cols[1:]

            # motif id is part of file name
            result_df["motif-id"] = link.text[:-4]

            # reorder columns
            result_df = result_df[final_cols]

            # switch org names to Latin
            result_df["organism"] = switch_org_name(result_df["organism"])

            # append new columns
            self.datf = self.datf.append(result_df, ignore_index=True)

        return True

    def export(self) -> bool:
        """ TODO """
        try:
            self.datf.to_csv(self.outfile, sep=",", header=True, index=False)
            return True
        except:
            return False
