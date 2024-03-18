import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("train/viral-vs-nonviral_train.csv")
domains = list(df["Domain"])
labels = list(set(domains))
data = [domains.count(label) for label in labels]
fig, ax = plt.subplots(figsize=(10, 5))
ax.barh(labels, data, color="green")
fig.savefig("plots/train/Domain.png")
plt.show()
for key in list(df.keys())[1:]:
    values = list(df[key])
    viral = [values[i] for i in range(len(values)) if domains[i] == "V"]
    nonviral = [values[i] for i in range(len(values)) if domains[i] == "NV"]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    ax[0].hist(viral, bins=30, color="blue")
    ax[0].set_title("Viral " + key)
    ax[1].hist(nonviral, bins=30, color="orange")
    ax[1].set_title("Non-viral " + key)
    fig.savefig(f"plots/train/{key}.png")
    plt.show()
