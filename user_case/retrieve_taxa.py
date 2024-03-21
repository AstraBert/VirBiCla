from Bio import Entrez
import pandas as pd
import time

Entrez.email = "your_email@example.com"


def fetch_nucleotide_record(accession):
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = handle.read()
    handle.close()
    return record

def extract_taxonomic_info(record):
    splitrecord = record.split("\n")
    for line in splitrecord:
        if "ORGANISM" in line:
            j = splitrecord.index(line)
            taxonomic_info = [splitrecord[j]]
            j+=1
            while "REFERENCE" not in splitrecord[j]:
                taxonomic_info.append(splitrecord[j])
                j+=1
                if "REFERENCE" in splitrecord[j]:
                    break
            string = " - ".join(taxonomic_info) 
            domains = ["Archaea", "Eukaryota", "Bacteria", "Viruses"]
            for domain in domains:
                if domain in string:
                    return domain
                else:
                    continue
            return string


csv = pd.read_csv("data.csv")
accession_nums = list(csv["ID"])
df = {"ID": [], "CLASSIFICATION": []}
k = 0
for num in accession_nums:
    nucleotide_record = fetch_nucleotide_record(num)
    taxonomic_info = extract_taxonomic_info(nucleotide_record)
    df["ID"].append(num)
    if taxonomic_info in ["Archaea", "Eukaryota", "Bacteria"]:
        df["CLASSIFICATION"].append("Non-viral")
    else:
        df["CLASSIFICATION"].append("Viral")
    k+=1
    print(f"Processed {k} reads")
df = pd.DataFrame.from_dict(df)
df.to_csv("classification.csv", index=False)