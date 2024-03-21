from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument(
    "-i",
    "--input_file",
    help="Path to the input file where all the raw reads are stored (must be fasta)",
    required=True,
)

argparse.add_argument(
    "-o", "--output_file", help="Path to the output csv file", required=True
)

argparse.add_argument(
    "-p","--plot", help="Plot final data to see a summary report of viral/non-viral sequences", required=False, default=False, action="store_true"
)

args = argparse.parse_args()


inf = args.input_file
out = args.output_file
plotthings = args.plot

from utils import *
from model import classifier
import gzip
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

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


def df_from_listofdicts(listofdifcs: list):
    refdict = listofdifcs[0]
    keys = list(refdict.keys())
    newdict = {
        key: [listofdifcs[i][key] for i in range(len(listofdifcs))] for key in keys
    }
    df = pd.DataFrame.from_dict(newdict)
    return df


if __name__=="__main__":
    csvpath = out
    hugelist = []
    print(f"Loading data from {inf}...")
    fasta = load_data(inf)
    print("Done")
    print("Writing csv...")
    c = 0
    for i in list(fasta.keys()):
        c += 1
        domain = "ND"
        smalldict = process_dna(domain, fasta[i])
        hugelist.append(smalldict)
        if c % 500 == 0:
            print(f"Processed {c} reads")
    print("Done")
    df = df_from_listofdicts(hugelist)
    X_pred = df.iloc[:, 1:]
    y_pred = classifier.predict(X_pred)
    newdict = {"SEQUENCE": [i.replace(",","-") for i in list(fasta.keys())], "PREDICTED_VALUE": ["Viral" if v == "V" else "Non-viral" for v in y_pred]}
    newdf = pd.DataFrame.from_dict(newdict)
    newdf.to_csv(out, index=False)
    if plotthings:
        dataf = pd.read_csv(out)
        names = ["Viral", "Non-Viral"]
        vals = [list(dataf["PREDICTED_VALUE"]).count("Viral"), list(dataf["PREDICTED_VALUE"]).count("Non-viral")]
        fig, ax = plt.subplots(figsize=(10,5))
        ax.bar(names, vals, color="blue")
        ax.set_title("VIRAL AND NON-VIRAL SEQ COUNT")
        ax.set_ylabel("Count")
        ax.set_xlabel("Predicted value")
        fig.savefig(out.replace(".csv",".png"))
