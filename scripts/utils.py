from codonbias.scores import EffectiveNumberOfCodons
import re
import collections
from scipy.stats import entropy
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import warnings


def calculate_codonbias(sequence: str):
    """
    Calculate the Effective Number of Codons (ENC) for a given DNA sequence.

    Parameters:
    - sequence (str): A DNA sequence

    Returns:
    - dict: A dictionary containing the Effective Number of Codons (ENC)
    """
    ENCcalc = EffectiveNumberOfCodons(bg_correction=True)
    ENC = ENCcalc.get_score(sequence)
    return {"Effective_Number_of_Codons": ENC}


def find_homopolymers(s: str):
    """Find homopolymeric regions and return the percentage of homopolymeric loci for each nucleotide."""
    a = re.compile(r"A{4,}")
    amatches = a.finditer(s)
    t = re.compile(r"T{4,}")
    tmatches = t.finditer(s)
    c = re.compile(r"C{4,}")
    cmatches = c.finditer(s)
    g = re.compile(r"G{4,}")
    gmatches = g.finditer(s)
    return {
        "Homopolymers_A": sum([match.end() - match.start() for match in amatches])
        / len(s),
        "Homopolymers_T": sum([match.end() - match.start() for match in tmatches])
        / len(s),
        "Homopolymers_C": sum([match.end() - match.start() for match in cmatches])
        / len(s),
        "Homopolymers_G": sum([match.end() - match.start() for match in gmatches])
        / len(s),
    }


def nt_proportion(s: str):
    """Calculate the nucleotide proportions and GC content of a DNA sequence."""
    gc = gc_fraction(s)
    return {
        "A_freq": s.count("A") / len(s),
        "T_freq": s.count("T") / len(s),
        "C_freq": s.count("C") / len(s),
        "G_freq": s.count("G") / len(s),
        "GC_content": gc,
    }


def estimate_information_entropy(dna_sequence):
    """Estimate the information entropy of a DNA sequence."""
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    dist = [x / sum(bases.values()) for x in bases.values()]
    entropy_value = entropy(dist, base=2)
    return {"Entropy": entropy_value}


def gene_density(s: str):
    """Estimate gene density based on coding statistics of reading frames."""
    reading_frames = [
        rf0 := s,
        rf1 := s[1:],
        rf2 := s[2:],
        revrf0 := str(Seq(s).reverse_complement()),
        revrf1 := str(Seq(s).reverse_complement())[1:],
        revrf2 := str(Seq(s).reverse_complement())[2:],
    ]
    translations = [str(Seq(rf).translate(stop_symbol="-")) for rf in reading_frames]
    pattern = r"M.{18,}?-"
    codingstats = [
        sum([len(match) for match in re.findall(pattern, translation)])
        / len(translation)
        for translation in translations
    ]
    return {"Gene_density": max(codingstats)}


def process_dna(domain: str, dna_str: str):
    """Process DNA sequence and extract various features."""
    warnings.filterwarnings("ignore")
    basedict = {"Domain": domain}
    listofdicts = [
        dict0 := nt_proportion(dna_str),
        dict1 := calculate_codonbias(dna_str),
        dict2 := find_homopolymers(dna_str),
        dict3 := estimate_information_entropy(dna_str),
        dict4 := gene_density(dna_str),
    ]
    for d in listofdicts:
        basedict.update(d)
    return basedict
