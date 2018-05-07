# standard
import re, random, argparse
from  os import path, makedirs
from time import sleep

# HTTP
from urllib.error import HTTPError

# dataframes
import pandas as pd
import numpy as np

# bio
from Bio import Entrez, SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# custom
from id_switcher import id_to_refseq, switch_org_name

# Global settings
Entrez.email = None
OUTFILE = './MA0852.2SEQS1k'
BEDFILE = "http://jaspar.genereg.net/download/bed_files/MA0852.2.bed" # change bedfile if needed
CHROM_DIR = "./Genome/" # store chromosomes in Genome folder

VERBOSE = True

# change length of returned sequence
LENGTH = 1000

def expand_query(query: pd.Series) -> pd.DataFrame:
    """ Parses query formatted for UCSC Biologic database """
    # expand query
    parser = re.compile(r"^([a-z]{2,4}\d{,2})(?:\_)(chr\d{,2}[A-Z]{,3})(?:[^\:]*):(\d+)-(\d+)(?:\()([\-\+])(?:\))")
    results = query.apply(parser.findall).apply(lambda r: list(sum(r, ()))).apply(pd.Series)
    return results

def load_jasper_bed_files(filename):
    """ loads jasper bed files which contain
        jasper motif-id and query of motif
    """
    final_cols = ["motif-id", "organism", "genome", "chromosome", "start", "stop", "strand"]
    # read csv
    locus_file = pd.read_csv(filename,
                             sep="\t",
                             header=None)

    # 3: query -> organism, chormosome, start, stop, strand
    result_df = expand_query(locus_file.iloc[:,3])

    # name columns
    result_df.columns = final_cols[2:]
    
    # convert start and stop
    result_df['start'] = result_df['start'].apply(int)
    result_df['stop'] = result_df['stop'].apply(int)

    # motif id is part of file name
    result_df["motif-id"] = filename[-12:-4]

    # switch org names to Latin
    result_df["organism"] = switch_org_name(result_df['genome'])

    # reorder columns
    return result_df[final_cols]

def fetch_chromosome(iden):
    """ fetch fasta sequence from NCBI given id """
    try:
        with Entrez.efetch(db='nucleotide', id=iden, rettype='fasta') as handle:
            rec = SeqIO.read(handle, 'fasta')
            rec.Alphabet = generic_dna
            return rec
    except HTTPError:
        return None

def download_needed_chromosomes(datf) -> int:
    """save chromosome sequences as fasta files in sequences folder
    returns: number of sequences successfully saved
    """
    saved = 0
    # for each chromosome:
    for gen, chrom in datf[['genome', 'chromosome']] \
                           .groupby(['genome', 'chromosome']).count().index:

        # get id for looking in NCBI
        refseq_id = id_to_refseq(gen, chrom.strip('chr'))

        outfilename = CHROM_DIR + gen + "_" + chrom.strip('chr') + ".fasta"
            
        if not path.isdir(CHROM_DIR):
            makedirs(CHROM_DIR)

        if path.isfile(outfilename):
            if VERBOSE:
                print(f"{outfilename} already exists.")
            continue

        # fetch sequence from NCBI
        big_sequence = fetch_chromosome(refseq_id)

        if not big_sequence:
            print('error: in {} {} refseq_id {}'.format(gen, chrom, refseq_id))

        else:
            # save sequence
            SeqIO.write(big_sequence, outfilename, 'fasta')

            # update status
            saved += 1
            print('saved', big_sequence.description)

            # ensure requests do not go to NCBI too fast
            sleep(1)
    return saved

def choose_row(r, chrom):
    if r["strand"] == "+":
        return str(chrom[r["rstart"]-1:r["rstop"]].seq)
    else: 
        return str(chrom[r["rstart"]-1:r["rstop"]].seq.reverse_complement())

def search_chromosome(chrom: SeqIO.SeqRecord, start: pd.Series,
                      stop: pd.Series, strand: pd.Series, length: int=LENGTH) -> tuple:
    """ returns (motif, mstart, mstop) """
    
    # length of motif
    len_motifs = (stop - start) + 1 # plus 1 because 0 vs. 1 indexing
    
    rstart = start - len_motifs.apply(lambda d: np.random.randint(0, length - d))
    rstop = rstart + length
    # get randomized indents within a set range

    # select motif +/- indents
    motifs = pd.concat([rstart, rstop, strand], keys=["rstart", "rstop", 'strand'], axis=1)
    motifs = motifs.apply(lambda r: choose_row(r, chrom), axis=1)

    # return motif, start index from selected sequence, and
    # stop index from selected sequence
    return (motifs, start - rstart, start - rstart + len_motifs)

def load_sequences(datf, length):
    """ load chromosomes and return dataframe of sequences with motifs """
    dirname = CHROM_DIR
    seqdat = pd.DataFrame()
    for gen, chrom in datf[['genome', 'chromosome']] \
                           .groupby(['genome', 'chromosome']).count().index:

        chrom_file = dirname + gen + "_" + chrom.strip("chr") + ".fasta"
        chrom_record = SeqIO.read(chrom_file, 'fasta')

        # get rows for organism and chromosome
        startstops = datf.loc[(datf['genome'] == gen) & (datf['chromosome'] == chrom)]
        # retrive motif + indent
        motifs, mstarts, mstops = search_chromosome(chrom_record,
                                                    startstops["start"],
                                                    startstops["stop"],
                                                    startstops["strand"],
                                                    length)
        rows = pd.concat([startstops, motifs, mstarts, mstops], axis=1)
        rows.columns = ["motif-id", "organism", "genome", "chromosome", "start",
                       "stop", "strand", "seq", "mstart", "mstop"]
        seqdat = seqdat.append(rows, ignore_index=True)

    return seqdat

def convert_to_fasta(outfile):
    anno_motifs = pd.read_csv(f"{outfile}.csv") # annotated motif instances
    # transform to seqRecords
    seqs = [SeqRecord(Seq(s), id=f"Seq{i}", description=f"Sequence from {outfile}") 
                                for i, s in zip(range(anno_motifs.shape[0]), anno_motifs['seq'])]
    SeqIO.write(seqs, f"{outfile}.fasta", 'fasta')

def load(bedfile, outfile, length):
    bed_df = load_jasper_bed_files(bedfile)
    download_needed_chromosomes(bed_df)
    seq_df = load_sequences(bed_df, length)
    seq_df.to_csv(f"{outfile}.csv", sep=",", header=True, index=False)
    convert_to_fasta(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=
        """Extracts Motif Sequences from NCBI using .bed files downloaded from JASPER db, 
        with the motif at a random position within the sequence. 

        Input a link to the bed file from JASPER and this script will attempt to download the NCBI sequences.
        """
    )
    parser.add_argument('--email', '-e', required=True, 
                        help="email address to use with NCBI")
    parser.add_argument('--outfile', '-o', default=OUTFILE, 
                        help="where to output final csv")
    parser.add_argument('--bedfile', '-i', default=BEDFILE, 
                        help="input .bed file")
    parser.add_argument('--chrom_dir', default=CHROM_DIR, 
                        help="folder to place downloaded chromosomes from NCBI")
    parser.add_argument('--verbose', '-v', type=bool, default=VERBOSE, 
                        help="whether to output all messages")
    parser.add_argument('--length', '-l', type=int, default=LENGTH, 
                        help="changes length of returned sequence")

    # parse args
    args = parser.parse_args()
    varz = vars(args)

    VERBOSE = varz['verbose']
    CHROM_DIR = varz['chrom_dir']
    Entrez.email = varz['email']

    load(varz['bedfile'], varz['outfile'], varz['length'])



