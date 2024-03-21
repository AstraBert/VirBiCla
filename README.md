# VirBiCla

A novel ML-based binary classifier to tell viral and non-viral DNA/RNA sequences apart.

## Table of Contents

- [Model Description](#model-description)
  - [Model Components](#model-components)
  - [Ensemble Voting](#ensemble-voting)
  - [Model Training](#model-training)
  - [Usage](#usage)
- [Training](#training)
  - [Data collection](#data-collection)
  - [Data processing](#data-processing)
  - [Training specs](#training-specs)
- [Testing](#testing)
  - [Data collection](#data-collection-1)
  - [Data processing](#data-processing-1)
  - [Tests](#tests)
- [Use advices and limitations](#use-advices-and-limitations)
- [Customization](#customization)
- [License and right of usage](#license-and-right-of-usage)
- [References](#references)

## Model Description

VirBiCla is an ensemble voting classifier model created using scikit-learn in python.

### Model Components
- **Logistic Regression (LR):** A linear model that predicts the probability of a categorical dependent variable.
- **Random Forest (RF):** An ensemble learning method that constructs a multitude of decision trees during training and outputs the mode of the classes.
- **Gaussian Naive Bayes (GNB):** A simple probabilistic classifier based on applying Bayes' theorem with strong (naive) independence assumptions.
- **Decision Tree (DT):** A non-parametric supervised learning method used for classification and regression.
- **K-Nearest Neighbors (KNN):** A non-parametric method used for classification and regression, where the input consists of the k closest training examples in the feature space.
- **Bagging Classifier (BC):** An ensemble meta-estimator that fits base classifiers each on random subsets of the original dataset.
- **Histogram-Based Gradient Boosting Classifier (HGB):** A gradient boosting model that uses histograms to speed up training.
- **Extra Trees Classifier (ETC):** An ensemble learning method that builds multiple decision trees and merges their predictions.
- **Stochastic Gradient Descent (SGD):** An optimization algorithm that updates the parameters iteratively to minimize the loss function.

### Ensemble Voting
- **Voting Method:** Soft voting, where the predicted probabilities for each class are averaged and the class with the highest probability is chosen.
- **Classifiers:** A combination of LR, RF, GNB, DT, KNN, BC, HGB, ETC, and SGD classifiers are combined into a voting ensemble.

### Model Training
- **Training Data:** The model is trained using the datasets you can see [here](./data/train/)
- **Fit Method:** The classifier is fit to the training data using the `fit` method.

### Usage
Once trained, this ensemble classifier can be used to make predictions on new data by calling the `predict` method.

You can use the model by running the [prediction script](./scripts/predict.py) as follows:

```bash
# make sure to include -p if you wish to get out a nice visualization of the viral/non-viral classification stats
python3 predict.py -i input.fasta -o output.csv -p
```

## Training
### Data collection
Training data were collected using [this script](./shell/retrieve_lots_of_blastdb.sh), which is based on the instructions provided by NCBI on [how to retrieve BLAST databases](https://www.ncbi.nlm.nih.gov/books/NBK569850/).

Training datasets are:

- 16S_ribosomal_RNA: 16S ribosomal RNA (Bacteria and Archaea type strains) - 22239 seqs
- 28S_fungal_sequences: 28S ribosomal RNA sequences (LSU) from Fungi type and reference material - 10292 seqs
- SSU_eukaryote_rRNA: Small subunit ribosomal nucleic acid for Eukaryotes - 8784 seqs
- ref_viruses_rep_genomes: Refseq viruses representative genomes - 18688 seqs

Descriptions are taken from a dedicated BLAST+ tools [repository](https://github.com/ncbi/blast_plus_docs?tab=readme-ov-file#blast-databases), while the number of sequences has been retrieved from the metadata files avauilable on [BLAST ftp website](https://ftp.ncbi.nlm.nih.gov/blast/db/).

### Data processing
Data were processed using [this script](./scripts/preprocess_train.py) and plots of the resulting training dataset were made with the aid of [this script](./scripts/plot_train.py).

Here's a step-by-step breakdown of the processing procedure:

1. **Nucleotide Proportions and GC Content Calculation:**
    - The input DNA sequence is analyzed to compute the frequency of each nucleotide (A, T, C, G) in the sequence.
    - Additionally, the GC content of the sequence is calculated, representing the proportion of Guanine (G) and Cytosine (C) bases compared to the total number of bases in the sequence.

2. **Effective Number of Codons (ENC) Computation:**
    - The program calculates the Effective Number of Codons (ENC), which provides insights into codon usage bias within the DNA sequence.
    - ENC is a measure of how efficiently codons are used to encode amino acids, with lower values indicating higher codon bias.

3. **Homopolymeric Region Detection:**
    - Homopolymeric regions, characterized by consecutive repeats of a single nucleotide (e.g., AAAA), are identified within the DNA sequence.
    - The program computes the percentage of homopolymeric loci for each nucleotide (A, T, C, G) present in the sequence.

4. **Information Entropy Estimation:**
    - Information entropy of the DNA sequence is estimated, providing a measure of uncertainty or randomness in the sequence.
    - Higher entropy values suggest greater sequence diversity, while lower values indicate more uniform or repetitive sequences.

5. **Gene Density Calculation:**
    - The program estimates gene density based on coding statistics of reading frames derived from the DNA sequence.
    - It computes the ratio of coding regions (sequences encoding proteins) to the total length of the translated sequences, providing insights into the density of potential genes within the sequence.

6. **Integration of Results:**
    - Finally, the extracted features and metrics, including nucleotide proportions, GC content, ENC, homopolymeric regions, information entropy, and gene density, are integrated into a structured dictionary.
    - The dictionary encapsulates all relevant information about the processed DNA sequence, facilitating further analysis and interpretation, and it is finally integrated with the viral/non-viral information ("V" or "NV") associated with each file.

7. **Results storage**:
    - The resulting dictionaries for each DNA sequences are stored into a list, converted into a DataFrame and, eventually, merged into a [CSV file](./data/train/viral-vs-nonviral_train.csv)

All the functions employed to process data are defined [here](./scripts/utils.py).

All the plots are available [here](./data/plots/train/).

### Training specs
The binary classificator is trained to predict the "Domain" field (in the CSV) starting from all the other parameters, classifying a sequence as either viral ("V") or non-viral ("NV"). 

Training takes about 40s on Google Colab free python engine (12GB RAM).

## Testing
### Data collection
Three test sets were prepared ([test](./data/test/), [test1](./data/test1/) and [test2](./data/test2/)), and they are composed by:

- 401 sequences (ranging from 1000 to 8000 bp) from recently submitted viral specimens, collected from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), 500 out of 4092 Araneae SSU rRNA sequences collected from [SILVA](https://www.arb-silva.de/) and 500 out of 1982 Tumebacillus SSU rRNA sequences also collected from SILVA.    
- 2000 DNA sequences (ranging from 1000 to 10000 bp) from bacteriophages, collected from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), 1000 out of 4092 Araneae SSU rRNA sequences collected from [SILVA](https://www.arb-silva.de/) and 1000 out of 1982 Tumebacillus SSU rRNA sequences also collected from SILVA.
- 1072 DNA sequences (ranging from 1000 to 15000 bp) from human viruses, collected from [NCBI virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), 600 sequences from [LSU_prpkaryote_rRNA](./data/train/LSU_prokaryote_rRNA.fsa.gz) and 600 sequences from [28S_fungal_sequences](./data/train/28S_fungal_sequences.fsa.gz).

### Data processing
Processing followed the same steps descripted for training, with [this script](./scripts/preprocess_test.py). Plots were also generated using a similar [script](./scripts/plot_train.py)

All the plots are available [here](./data/plots/)

### Tests
You can find the tests in the [dedicated notebook](./scripts/VirBiCla.ipynb). 

- `test` dataset had the following output:

    + **Confusion matrix:**
    
    |  |  |
    |---|---|
    | True postive: 986 |  False positive: 14 |
    | False negative: 14 | True negative: 387 |

    + **Accuracy:** 0.9800142755174875

- `test1` dataset had the following output:

    + **Confusion matrix:**
    
    |  |  |
    |---|---|
    | True postive: 1980 |  False positive: 20 |
    | False negative: 32 | True negative: 1968 |

    + **Accuracy:** 0.987

- `test2` dataset had the following output:

    + **Confusion matrix:**
    
    |  |  |
    |---|---|
    | True postive: 1181 |  False positive: 19 |
    | False negative: 1 | True negative: 1071 |

    + **Accuracy:** 0.9911971830985915

Overall, the models seems to perform well in telling viral and non-viral DNA sequences apart.

## Use: advices and limitations
<!-- User case in based on BioProject [PRJEB52499](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB52499), from paper [*Characterizations of the Gut Bacteriome, Mycobiome, and Virome in Patients with Osteoarthritis*](https://journals.asm.org/doi/10.1128/spectrum.01711-22). -->

User case for VirBiCla is [here](./user_case/).

It is based on data downloaded from NCBI SRA, whose accession numbers are retained in their names, thus all the information about sequencing and related research are easily accesible.

The three samples refer to _Aedes aegypti_'s virome ([SRR25317595](./user_case/SRR25317595_resized.fasta)), _Tuta absoluta_'s microbiome ([SRR20766474](./user_case/SRR20766474.fas)) and rRNA genes region from several fungi samples ([ERR2576718](./user_case/ERR2576718.fasta)).

Size-selected virome (>1000 bp) was first combined with microbiome without any size restraint, and the resulting ca. 17000 reads were predicted by VirBiCla, resulting in [50% accuracy](./user_case/microbiome_notselected/). From the [confusion matrix](./user_case/microbiome_notselected/test.stats) you can easily see that the model predicted correctly almost all the actual viral sequences, whereas it misclassified the bacterial ones.

When size-selected virome was combined with size-selected microbiome ([>1000 bp](./user_case/microbiome_selected/gt1000/) and [>1500 bp](./user_case/microbiome_selected/gt1500/)), progressively greater accuracy (53 and 60%) and better confusion matrices were observed (you can find them in test.stats files in each folder).

When size-selected virome was combined with unprocessed mycetome (which had most of the sequences way longer than 2000 bp), [results](./user_case/microbiome_selected/gt2000/) were definitely better. VirBiCla yielded an almost perfect performance, with 98% accuracy and [few misclassfied sequences](./user_case/microbiome_selected/gt2000/test.stats).

As a general advice, we suggest that the user employes the model "as is" only on long-read sequencing metagenomics samples: as you can see, VirBiCla is really good at differentiating between non-viral and viral sequences when put into a context of long reads (>2000-3000 bp). 

## Customization

If you have different data and you wish to custom VirBiCla for your needs, you need to follow these steps:

1. Modify [preprocess_train.py](./scripts/preprocess_test.py), adding your own fasta files to the preprocessing pipeline. For example, let's say that we have two files called "virome.fa" and "microbiome.fa": we should then make these modifications:
```python
MAPPING_DOMAINS = {
    "virome": "V",
    "microbiome": "NV",
}

if __name__ == "__main__":
    csvpath = "viral-vs-nonviral_train.csv"
    hugelist = []
    for fsa in list(MAPPING_DOMAINS.keys()):
        print(f"Loading data from {fsa}...")
        fastafile = f"{fsa}.fa"
```
2. Modify [model.py](./scripts/model.py) with your data. For example, let's say that we have saved our training data in a CSV file called "viral-vs-nonviral_train.csv": we should then make these modifications:
```python
df = pd.read_csv("viral-vs-nonviral_train.csv")
X_train = df_train.iloc[:, 1:]
y_train_rev = df_train["Domain"]
```

## License and right of usage

The program is provided under MIT license (more at [LICENSE](./LICENSE)).

When using VirBiCla, please consider citing the author ([Astra Bertelli](https://astrabert.vercel.app)) and this repository.

## References
- [scikit-learn](https://scikit-learn.org/stable/)
- [codon-bias](https://codon-bias.readthedocs.io/en/latest/)
- [matplotlib](https://matplotlib.org/)
- [biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://scipy.org/)
- NCBI, SILVA and BLAST+-tools docs repository