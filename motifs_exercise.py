import numpy as np

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

if __name__ == "__main__":
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