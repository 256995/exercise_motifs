import itertools

import numpy as np
from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
import random

def count_matrix(motifs):
    """
    counts the bases on each position in the motifs and returns as a dictionary (list of counts for each base)
    """
    motif_length = len(motifs[0])
    count_mat = {
        'A': [0 for i in range(motif_length)],
        'C': [0 for i in range(motif_length)],
        'G': [0 for i in range(motif_length)],
        'T': [0 for i in range(motif_length)]
    }
    for motif in motifs:
        for idx, base in enumerate(motif):
            count_mat[base][idx] += 1

    return count_mat


def score(motifs):
    """
    For every column of the count matrix take the largest count; the score is the sum of these maxima
    """
    count_mat = count_matrix(motifs)
    list_max = []
    for idx in range(len(motifs[0])):
        i_max = max(count_mat['A'][idx], count_mat['C'][idx], count_mat['G'][idx], count_mat['T'][idx])
        list_max.append(i_max)
    return sum(list_max)


def consensus(motifs):
    """
    Returns the most frequent base of every position in motif
    """
    count_mat = count_matrix(motifs)
    cons = ['A' for i in range(len(motifs[0]))] # 'A' at all positions
    for base in ['C', 'G', 'T']:
        for idx in range(len(motifs[0])):
            current_base = cons[idx]
            if count_mat[base][idx] > count_mat[current_base][idx]:
                cons[idx] = base
    return cons


def hamming_distance(a, b):
    """
    returns the number of positions at which the two strings differ
    """
    hamming_dist = 0
    for idx in range(len(a)):
        if a[idx] != b[idx]:
            hamming_dist += 1
    return hamming_dist


def total_distance(pattern, sequences):
    """
    For every sequence, find the smallest Hamming distance between pattern and any of its l-mers (windows); then add these minima up over all sequences.
    """
    list_min = []
    for seq in sequences:
        seq_dist = []
        for shift in range(len(seq) - len(pattern) + 1):
            dist = hamming_distance(pattern, seq[shift:shift+len(pattern)])
            seq_dist.append(dist)
        list_min.append(min(seq_dist))
    return sum(list_min)


class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.pseudocount = pseudocount
        self.motifs = motifs
        self.l = len(motifs[0])
        self.ppm = count_matrix(motifs)
        for base in ['A', 'C', 'G', 'T']:
            ppm_array = np.array(self.ppm[base])
            ppm_array = (ppm_array + self.pseudocount) / (len(motifs) + 4 * self.pseudocount)
            self.ppm[base] = ppm_array.tolist()

    def lmer_probability(self, lmer):
        probability = 1
        for idx, base in enumerate(lmer):
            prob = self.ppm[base][idx]
            probability = probability * prob
        return probability

    def most_probable_lmer(self, sequence):
        lmer = []
        max_prob = 0
        for shift in range(len(sequence) - self.l + 1):
            prob = self.lmer_probability(sequence[shift:shift+self.l])
            if prob > max_prob:
                max_prob = prob
                lmer = sequence[shift:shift+self.l]
        return lmer

    def consensus(self):
        count_mat = count_matrix(self.motifs)
        cons = ['A' for i in range(self.l)]  # 'A' at all positions
        for base in ['C', 'G', 'T']:
            for idx in range(self.l):
                current_base = cons[idx]
                if count_mat[base][idx] > count_mat[current_base][idx]:
                    cons[idx] = base
        cons = ''.join(cons) # convert form list to string
        return cons


class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed) # the finder's own random generator
        self.windows = [] # a list with one list of l-mers per sequence: every sequence cut into all its l-mers
        for seq in self.sequences:
            w = [] # list of windows for seq
            for shift in range(len(seq) - self.l + 1):
                w.append(seq[shift:shift+self.l])
            self.windows.append(w)

    def total_distance(self, pattern):
        list_min = [] # list of minimal distance of the pattern from each sequence
        for seq in self.windows:
            seq_dist = []
            for w in seq:
                dist = hamming_distance(pattern, w)
                seq_dist.append(dist)
            list_min.append(min(seq_dist))
        return sum(list_min)

    def median_string(self):
        """
        Returns the pattern with the smallest total distance out of all 4^l patterns of length l, together with that distance
        """
        all_possible_str = []
        for chars in itertools.product("ACGT", repeat=self.l):
            s = ''.join(chars)
            all_possible_str.append(s)
        distances = []
        for str in all_possible_str:
            distances.append(self.total_distance(str))
        idx = distances.index(min(distances))
        return all_possible_str[idx], distances[idx]

    def randomized_search(self):
        motif_matrix = []
        for seq in self.windows:
            motif_matrix.append(self.rng.choice(seq))
        motif_score = score(motif_matrix)
        while True:
            rand_profile = MotifProfile(motif_matrix, pseudocount=1)
            new_matrix = []
            for sequence in self.sequences:
                motif_matrix.append(rand_profile.most_probable_lmer(sequence))
            if score(new_matrix) <= motif_score:
                break
            motif_matrix = new_matrix
            motif_score = score(motif_matrix)
        return motif_matrix, motif_score

    def best_of(self, runs):
        best_matrix = []
        best_score = 0
        for i in range(runs):
            i_matrix, i_score = self.randomized_search()
            if i_score > best_score:
                best_matrix = i_matrix
                best_score = i_score
        return best_matrix, best_score


if __name__ == "__main__":
    ## TASK 1
    lecture_dna = [
        "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
        "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
        "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
        "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
        "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
        "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
    ]
    print(count_matrix(lecture_dna))
    print(score(["ACGT", "ATGT", "CCGA"]))
    print(consensus(["ACGT", "ATGT", "CCGA"]))
    print(hamming_distance("ACGT", "ACCA"))
    print(total_distance("AC", ["GACT", "TTAG"]))

    ## TASK 2
    print('---- MotifProfile ----')
    profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
    print(profile.l)  # 7
    print(profile.ppm["A"])  # [0.5, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5]
    print(round(profile.lmer_probability("ATGCGTA"), 4))
    print(MotifProfile(["GTAC", "TTAA"]).most_probable_lmer("ACTGGATGACCC")) # "TGAC"
    print(profile.consensus()) # "ATGCGTA"

    # check with biopython
    bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
    bio.pseudocounts = 1
    print(bio.consensus)  # ATGCGTA
    print(bio.pwm["A"])  # the same numbers as your profile.ppm["A"]

    ## TASK 3
    sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]
    finder = MotifFinder(lecture_dna, 6)
    median_str = finder.median_string()
    print(median_str)

    # does not work
    print(finder.best_of(100))
