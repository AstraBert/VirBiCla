from preclassify import load_data
import os


if __name__ == "__main__":
    files = [f for f in os.listdir("./") if f.endswith(".fas")]
    df = {"SRR20766474": {}}
    for f in files:
        seqsdict = load_data(f)
        for i in list(seqsdict.keys()):
            if f.split(".")[0]=="SRR20766474" and len(seqsdict[i])>1500:
                df[f.split(".")[0]].update({i: seqsdict[i]})
            else:
                continue   
    for key in list(df.keys()):
        outfasta = open(f"{key}_resized.fasta", "w")
        for k in list(df[key].keys()):
            outfasta.write(">"+k+"\n")
            outfasta.write(df[key][k]+"\n")
        outfasta.close()
    
