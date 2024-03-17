from preprocess_train import *
import random as r
import time

MAPPING_DOMAINS_TRUE = {"train/18S_fungal_sequences": "NV", "train/LSU_prokaryote_rRNA": "NV", "test2/viral1072": "V"}

if __name__ == "__main__":
    csvpath = "test2/viral-vs-nonviral_test2.csv"
    hugelist = []
    for fsa in list(MAPPING_DOMAINS_TRUE.keys()):
        print(f"Loading data from {fsa}...")
        fastafile = f"{fsa}.fsa.gz"
        fasta = load_data(fastafile)
        keys = list(fasta.keys())
        if MAPPING_DOMAINS_TRUE[fsa]=="NV":
            r.shuffle(keys)
            keys = keys[:600]
        print("Done")
        print("Writing csv...")
        c = 0
        for i in keys:
            c += 1
            domain=MAPPING_DOMAINS_TRUE[fsa]
            smalldict = process_dna(domain, fasta[i])
            hugelist.append(smalldict)
            if c % 500 == 0:
                print(f"Processed {c} reads")
        print("Done")
    df = df_from_listofdicts(hugelist)
    df.to_csv(csvpath, index=False)