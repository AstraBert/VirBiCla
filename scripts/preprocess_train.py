import gzip
from Bio import SeqIO
from utils import *
import pandas as pd


def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if (
            infile.endswith(".fasta.gz")
            or infile.endswith(".fna.gz")
            or infile.endswith(".fsa.gz")
            or infile.endswith(".fa.gz")
        ):
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
        else:
            y.close()
            raise ValueError("File is the wrong format")
    # Read file directly as fasta if it is a not zipped fasta: handle also more uncommon extensions :-)
    elif (
        infile.endswith(".fasta")
        or infile.endswith(".fna")
        or infile.endswith(".fsa")
        or infile.endswith(".fa")
    ):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")


def df_from_listofdicts(listofdifcs: list):
    refdict = listofdifcs[0]
    keys = list(refdict.keys())
    newdict = {
        key: [listofdifcs[i][key] for i in range(len(listofdifcs))] for key in keys
    }
    df = pd.DataFrame.from_dict(newdict)
    return df


MAPPING_DOMAINS = {
    "train/16S_ribosomal_RNA": "NV",
    "train/28S_fungal_sequences": "NV",
    "train/SSU_eukaryote_rRNA": "NV",
    "train/ref_viruses_rep_genomes": "V",
}

if __name__ == "__main__":
    csvpath = "train/viral-vs-nonviral_train.csv"
    hugelist = []
    for fsa in list(MAPPING_DOMAINS.keys()):
        print(f"Loading data from {fsa}...")
        fastafile = f"{fsa}.fsa.gz"
        fasta = load_data(fastafile)
        print("Done")
        print("Writing csv...")
        c = 0
        for i in list(fasta.keys()):
            c += 1
            domain = MAPPING_DOMAINS[fsa]
            smalldict = process_dna(domain, fasta[i])
            hugelist.append(smalldict)
            if c % 500 == 0:
                print(f"Processed {c} reads")
        print("Done")
    df = df_from_listofdicts(hugelist)
    df.to_csv(csvpath, index=False)
