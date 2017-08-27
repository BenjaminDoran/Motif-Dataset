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
    def __init__(self, url="http://jaspar.genereg.net/html/DOWNLOAD/",
                 bedfilename="motif-data.csv", createnew=True):
        """ Setup Importer """
        self.url = url
        self.bedfilename = bedfilename

        if createnew:
            self.df = pd.DataFrame([], columns=["motif-id", "organism",
                                                "chromosome", "start",
                                                "stop", "strand"])
        elif os.path.isfile(bedfilename) and not createnew:
            self.load_csv()
        else:
            print('No data to load, acting in createnew mode')
            self.df = pd.DataFrame([])

    # functions
    def load_csv(self) -> bool:
        """ reads in data.csv """
        self.df = pd.read_csv(self.bedfilename, sep=",")
        return True

    def load_jasper_bed_files(self) -> bool:
        """ 
        loads jasper bed files which contain
        jasper motif-id and query of motif
        """
        bed_url = self.url + "bed_files/"
        final_cols = ["motif-id","organism","chromosome","start","stop","strand"]

        # enter bed directory
        response = requests.get(bed_url)
        soup = BeautifulSoup(response.content, "html.parser")
        # all links to files
        bed_motifs = soup.find_all("a", text=re.compile("MA."))
        # for link to file
        for a in bed_motifs:
            # ignore files that do not follow format (-400 entries)
            if a.text not in ("MA0548.1.bed","MA0555.1.bed","MA0563.1.bed"):
                # read csv
                locus_file = pd.read_csv(bed_url+a.text, sep="\t", header=None)

                # query is column 3 (0 indexed)
                locus_file.rename(columns={3: "query"}, inplace=True)

                # query -> organism, chormosome, start, stop, strand
                result_df = self.expand_query(locus_file["query"])

                # name columns
                result_df.columns = final_cols[1:]
                
                # motif id is part of file name
                result_df["motif-id"] = a.text[:-4]

                # reorder columns
                result_df = result_df[final_cols]

                # append new columns
                self.df = self.df.append(result_df, ignore_index=True)
        return True

    def load_jasper_pfm_files(self) -> bool:
        """ TODO """
        return False

    def expand_query(self, query:pd.Series) -> pd.DataFrame:
        """ Parses query formatted for UCSC Biologic database """
        # expand query
        parse_query = re.compile("[a-zA-Z0-9]+|\([\+\-]\)")
        results = query.apply(parse_query.findall).apply(pd.Series)

        # get strand without "()"
        restrip = re.compile("[\+\-]")
        results[4] = results[4].apply(lambda y: restrip.search(y).group(0))
        
        return results

    def switch_org_name(self, orgs:pd.Series) -> pd.Series:
        """ Converts UCSC style org ref names to Latin Names """
        name_table = {'tair10': 'arabidopsis thaliana', 
                      'at': 'arabidopsis thaliana', 
                      'hg19': 'homo sapian', 
                      'ce': 'caenorhabditis elegans ', 
                      'mm9': 'mus muscosis', 
                      'dm': 'drosophila melanogaster'}

        return orgs.apply(lambda x: name_table[x])

    def export(self) -> bool:
        """ TODO """
        self.df.to_csv(self.bedfilename, sep=",", header=True, index=False)
