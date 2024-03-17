import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("test2/viral-vs-nonviral_test2.csv")
domains = list(df["Domain"])
labels = list(set(domains))
data = [domains.count(label) for label in labels]
fig, ax = plt.subplots(figsize=(10,5))
ax.barh(labels, data, color="green")
fig.savefig("plots/test2/Domain.png")
plt.show()
for key in list(df.keys())[1:]:
    values = list(df[key])
    vir = [values[i] for i in range(len(values)) if domains[i]=="V"]
    nonvir = [values[i] for i in range(len(values)) if domains[i]=="NV"]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    ax[0].hist(vir,bins=30,color="blue")
    ax[0].set_title("Viral " + key)
    ax[1].hist(nonvir,bins=30,color="orange")
    ax[1].set_title("Non-viral " + key)
    fig.savefig(f"plots/test2/{key}.png")
    plt.show()