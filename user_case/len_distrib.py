from preclassify import load_data
import matplotlib.pyplot as plt
import os


plt.style.use("seaborn-v0_8-paper")

if __name__ == "__main__":
    files = [f for f in os.listdir("./") if f.endswith(".fasta.gz")]
    df = {"ERR2576718": []}
    for f in files:
        seqsdict = load_data(f)
        for i in list(seqsdict.keys()):
            df[f.split(".")[0]].append(len(seqsdict[i]))
    fig, ax = plt.subplots(figsize=(10,5))
    ax.hist(df["ERR2576718"], bins=30, color="blue")
    ax.set_title("LEN DISTRIBUTION ERR2576718 - Mycetome")
    # ax[0].set_yscale("log")
    fig.tight_layout()       
    fig.savefig("lendistrib_mycetome_resized.png")
    plt.show()