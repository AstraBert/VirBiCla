import gzip
from Bio import SeqIO
import pandas as pd
import os


def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if (
            infile.endswith(".fasta.gz")
            or infile.endswith(".fna.gz")
            or infile.endswith(".fas.gz")
            or infile.endswith(".fa.gz")
            or infile.endswith(".fsa.gz")    
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
        or infile.endswith(".fas")
        or infile.endswith(".fa")
        or infile.endswith(".fsa")
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


if __name__ == "__main__":
    files = [f for f in os.listdir("./") if f.endswith(".fasta")]
    df = {"SEQUENCE": [], "CLASSIFICATION": []}
    for f in files:
        seqsdict = load_data(f)
        for i in list(seqsdict.keys()):
            df["SEQUENCE"].append(i)
        if f.startswith("ERR"):
            for i in range(len(list(seqsdict.keys()))):
                df["CLASSIFICATION"].append("Non-viral")
        else:
            for i in range(len(list(seqsdict.keys()))):
                df["CLASSIFICATION"].append("Viral")
    df = pd.DataFrame.from_dict(df)
    df.to_csv("preclassification_resized.csv", index=False)            