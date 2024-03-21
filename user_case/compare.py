import pandas as pd

csv1 = "metagenome_predicted_VirBiCla.csv"
csv2 = "preclassification.csv"

df1 = pd.read_csv(csv1)
df2 = pd.read_csv(csv2)

preds = list(df1["PREDICTED_VALUE"])
preclass = list(df2["CLASSIFICATION"])

true_positive = 0
true_negative = 0
false_positive = 0
false_negative = 0

for i in range(len(preds)):
    if preds[i]=="Non-viral" and preds[i]==preclass[i]:
        true_negative+=1
    if preds[i]=="Viral" and preds[i]==preclass[i]:
        true_positive+=1 
    if preds[i]=="Non-viral" and preds[i]!=preclass[i]:
        false_negative+=1
    if preds[i]=="Viral" and preds[i]!=preclass[i]:
        false_positive+=1 
accuracy = (true_positive + true_negative)/(true_negative+true_positive+false_negative+false_positive)
print(f"Confusion matrix\nV\tTrue\tFalse\n\t{true_positive}\t{false_positive}\nNV\tTrue\tFalse\n\t{true_negative}\t{false_negative}")
print(f"Accuracy: {accuracy}")
