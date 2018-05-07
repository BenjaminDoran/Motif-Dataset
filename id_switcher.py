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