""" switches from UCSC genome to Accesion number in NCBI"""
id_switch = {
    'hg38':
        lambda chrom:
        f"NC_0000{(int(chrom)):02d}" if
        chrom.isnumeric() and (0 < int(chrom) < 23) else
        "NC_000023" if chrom == "X" else
        "NC_000024" if chrom == "Y" else
        None,
}


def id_to_refseq(gen, chrom):
    return id_switch[gen](chrom)

""" Converts UCSC style org ref names to Latin Names """
name_table = {'tair10': 'arabidopsis thaliana',
              'at': 'arabidopsis thaliana',
              'hg19': 'homo sapiens',
              'hg38': 'homo sapiens',
              'ce': 'caenorhabditis elegans',
              'mm9': 'mus musculus',
              'dm': 'drosophila melanogaster'}

def switch_org_name(orgs):
    return orgs.apply(lambda x: name_table[x])