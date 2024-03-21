from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier, VotingClassifier, BaggingClassifier, HistGradientBoostingClassifier, ExtraTreesClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
import pandas as pd
from datasets import load_dataset

df = load_dataset("as-cle-bert/VirBiCla-training")
df_train = df["train"]
df_train = pd.DataFrame.from_dict(df_train)
# Features
X_train = df_train.iloc[:, 1:]

# Labels
y_train_rev = df_train["Domain"]

print("Training model...")
clf1 = LogisticRegression(multi_class='auto', max_iter=1000, random_state=1)
clf2 = RandomForestClassifier(n_estimators=50, random_state=1)
clf3 = GaussianNB()
clf4 = DecisionTreeClassifier()
clf5 = KNeighborsClassifier(n_neighbors=5)
clf6 = BaggingClassifier()
clf7 = HistGradientBoostingClassifier()
clf8 = ExtraTreesClassifier()
clf9 = SGDClassifier(loss="log_loss")
classifier = VotingClassifier([('lr', clf1), ('rf', clf2), ('gnb', clf3), ('dt', clf4), ('knn', clf5), ('bc', clf6), ('hgb', clf7), ('etc', clf8), ('sgd', clf9)], voting='soft')
classifier = classifier.fit(X_train,y_train_rev)
print("Done")
