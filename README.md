# Motif Search

## What are we looking for?

A **motif** is a short sequence pattern that occurs, with small variations, in many
sequences because it has a function: for example the binding site of a transcription
factor, which is found in the promoters of the genes it regulates. The copies are never
quite identical (a transcription factor tolerates mismatches), and we do not know where
they are. The task of **motif discovery** is:

> Given *t* DNA sequences, find one *l*-mer in each sequence so that the chosen
> *l*-mers are as similar to each other as possible.

In the lecture you saw several algorithms for this problem. In this exercise you will
implement two of them:

- the **median string** search, which is exact: it tries all 4^l possible patterns;
- **randomized motif search**, which is a heuristic: it starts from random positions,
  improves them step by step, and is restarted many times.

Structure of the session:
1. **Motif matrix, score and distance** (theory) → **Task 1** 
2. **Profiles** (theory) → **Task 2**
3. **Randomized motif search and the `random` module** (theory) → **Task 3**: a longer
   task you work through independently, combining Tasks 1 and 2.

Write your own code from scratch in a single file, `motifs_exercise.py`. Fork this
repository to your own GitHub account, clone your fork, and commit and push your work
to it as you go (e.g. one commit per finished task).

Reading: Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning
Approach*, 2nd ed., Vol. I, Chapter 2 ("Which DNA Patterns Play the Role of Molecular
Clocks?").

---

## Task 1: Warm-up functions

The six sequences from the lecture 2 (slide 19; *t* = 6, *n* = 57, *l* = 6):

```python
lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]
```

Write the following five functions. Later tasks build on them, so reuse them wherever
a function below says so. Wherever a function takes `motifs`, it is a list of equally
long strings (one *l*-mer per sequence); the bases are always `"A"`, `"C"`, `"G"`,
`"T"`.

### The motif matrix: the building blocks of profiles and randomized search

#### `count_matrix(motifs)`
- **Input:** `motifs`, e.g. `["ACGT", "ATGT", "CCGA"]`
- **Output:** a dict with the four bases as keys; each value is a list with one number
  per column, saying how many of the motifs have that base in that column:
  ```python
  count_matrix(["ACGT", "ATGT", "CCGA"])
  # {"A": [2, 0, 0, 1], "C": [1, 2, 0, 0], "G": [0, 0, 3, 0], "T": [0, 1, 0, 2]}
  ```
- **Hint:** start with a list of zeros of length *l* for every base, then go through
  every motif and every position in it and add 1 to the right base and column. The
  four numbers of each column always add up to the number of motifs.

#### `score(motifs)`
- **Input:** `motifs`
- **Output:** `int`. For every column of the count matrix take the largest count; the
  score is the sum of these maxima. The more similar the motifs, the higher the score;
  the maximum is *t* · *l*.
- **Example:** `score(["ACGT", "ATGT", "CCGA"])` → `2 + 2 + 3 + 2` = `9`

#### `consensus(motifs)`
- **Input:** `motifs`
- **Output:** `str` of length *l*, the most frequent base of every position in motif
- **Example:** `consensus(["ACGT", "ATGT", "CCGA"])` → `"ACGT"`
- **Ties:** if several bases share the highest count at a position, take the first of
  them in the order A, C, G, T, e.g. `consensus(["AC", "GT"])` → `"AC"`.
- **Hint:** Go through the bases in the order `"ACGT"` and replace the current best 
  only when a count is **strictly** higher (gives you the tie rule for free).


### Distances: the building blocks of the median string

#### `hamming_distance(a, b)`
- **Input:** two strings of the same length
- **Output:** `int`, the number of positions at which the two strings differ
- **Example:** `hamming_distance("ACGT", "ACCA")` → `2`

#### `total_distance(pattern, sequences)`
- **Input:** a `pattern` of length *l* and a list of `sequences` (longer strings)
- **Output:** `int`. For every sequence, find the smallest Hamming distance between
  `pattern` and any of its *l*-mers (windows); then add these minima up over all
  sequences.
- **Example:** `total_distance("AC", ["GACT", "TTAG"])` → `1`:
  - `"GACT"`: windows `GA`, `AC`, `CT`, distances 2, 0, 2 → minimum `0`
  - `"TTAG"`: windows `TT`, `TA`, `AG`, distances 2, 2, 1 → minimum `1`

---

Check your functions against the lecture 2:

```python
red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

print(score(red))                                # 26
print(consensus(best), score(best))              # AGATAG 34
print(hamming_distance("TAAGTT", "TGAATT"))      # 2
print(total_distance("TGCGTT", lecture_dna))     # 13
```

Lecture slide 24 gets 14 for `TGCGTT`, because there the distance is measured to the
red l-mers of slide 19. `total_distance` takes the **closest** l-mer in every
sequence, and in sequence 1 a closer one exists.

---

## Task 2: `MotifProfile` class

Implement a class `MotifProfile`. It turns a list of motifs into a profile (the
probability of each base in each column) and then scores *l*-mers with it.

#### Constructor `__init__(self, motifs, pseudocount=1)`
- **Input:** `motifs` (a list of equally long strings, as in Task 1) and the
  `pseudocount` added to every count
- **Stores:** `self.l`, the length of the motifs, and `self.ppm`, the position
  probability matrix: a `dict` of lists like the result of `count_matrix` (the four
  bases as keys, one list of *l* numbers per base), but holding probabilities instead
  of counts: (count + pseudocount) / (*t* + 4 · pseudocount), where *t* is the number
  of motifs
- **Example:**
  ```python
  profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
  profile.l           # 7
  profile.ppm["A"]    # [0.5, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5]
  ```
- **Hint:** reuse `count_matrix` from Task 1. The four probabilities of every column
  add up to 1.

#### `lmer_probability(self, lmer)`
- **Input:** an *l*-mer of length `self.l`
- **Output:** `float`, the product of the probabilities of its bases, column by column
- **Example:** `profile.lmer_probability("ATGCGTA")` → `0.0122` (rounded),
  i.e. 4/8 · 4/8 · 4/8 · 5/8 · 4/8 · 5/8 · 4/8

#### `most_probable_lmer(self, sequence)`
- **Input:** a sequence at least `self.l` long
- **Output:** `str`, the *l*-mer (window) of `sequence` with the highest probability
- **Example:** `MotifProfile(["GTAC", "TTAA"]).most_probable_lmer("ACTGGATGACCC")`
  → `"TGAC"` (probability 1/108 ≈ 0.0093)
- **Ties:** if several windows share the highest probability, take the leftmost one,
  e.g. `MotifProfile(["GTAC", "TTAA"]).most_probable_lmer("GGATGA")` → `"GGAT"`
  (`GGAT` and `ATGA` both have 1/216).
- **Hint:** go through the windows as in `total_distance` and call `lmer_probability`
  for each. Replacing the best window only when a probability is **strictly** higher
  gives you the tie rule for free.

#### `consensus(self)`
- **Output:** `str` of length `self.l`, the most probable base of every column
- **Example:** `profile.consensus()` → `"ATGCGTA"`
- **Ties:** the same rule as for `consensus` in Task 1 (the first base in the order
  A, C, G, T).

Check your class against the lecture (slides 13 and 30):

```python
profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.consensus())                       # ATGCGTA
print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
print(round(two.lmer_probability("TGAC"), 4))         # 0.0093
```

**Check against a library.** Biopython has a module for motifs, `Bio.motifs`
([tutorial](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html),
[full reference](https://biopython.org/docs/latest/api/Bio.motifs.html)). Build the
same profile with it and check that every probability agrees with your `MotifProfile`:

```python
from Bio import motifs
from Bio.Seq import Seq

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]
```

From now on you may use `Bio.motifs` for what it is good at (file formats, scanning,
thresholds, comparing motifs); the search algorithms in Task 3 stay yours.

---

## The `random` module

```python
import random

rng = random.Random(1)                 # a random number generator with seed 1
i = rng.randint(0, 50)                 # random integer, 0 <= i <= 50 (both ends included!)
lmer = rng.choice(["ACG", "CGT", "GTA"])   # one random item of a list
```

`random.Random(seed)` creates a generator **object** with its own sequence of random
numbers. The same seed gives the same numbers every time you run your script, and two
generators never influence each other. Without a seed (`random.Random()`), every run
gives different numbers.

Two different programs will still not get the same results from the same seed,
because they ask the generator for a different number of random values. That is why
the Expected result in Task 3 gives ranges for the random parts.

---

## Task 3: `MotifFinder` class

This is the main task. Work through it independently, combining Tasks 1 and 2.

`planted_motif.fasta` contains 10 random sequences of 60 bp. Into each one, the same
7-bp motif was planted, with one random mutation in every copy. Your job is to find it.

Read the sequences with [Biopython](https://biopython.org/):

```python
from Bio import SeqIO

sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]
```

Implement a class `MotifFinder`. It holds a set of sequences and searches them for a
motif of length *l*, with the median string and with randomized motif search.

#### Constructor `__init__(self, sequences, l, seed=None)`
- **Input:** a list of `sequences`, the motif length `l` and an optional `seed`
- **Stores:**
  - `self.sequences` and `self.l`;
  - `self.rng = random.Random(seed)`, the finder's own random number generator. Use
    it for every random choice in the class (never the global `random.randint`), so
    that a finder with a given seed always gives the same results;
  - `self.windows`, a list with one list of *l*-mers per sequence: every sequence cut
    into all its *l*-mers **once**. All methods below work with these windows instead
    of slicing the sequences again and again.

#### `total_distance(self, pattern)`
- **Output:** `int`, like `total_distance` from Task 1, but using `self.windows`

#### `median_string(self)`
- **Output:** `(pattern, distance)`, the pattern with the smallest total distance out
  of all 4^l patterns of length *l*, together with that distance
- **Hint:** `itertools.product("ACGT", repeat=l)` produces all 4^l combinations of
  letters as tuples; `"".join(...)` turns a tuple into a string.

#### `randomized_search(self)`
- **Output:** `(motifs, score)`, the result of **one run** of randomized motif search:
  1. pick a random *l*-mer in every sequence (`self.rng.choice` over its windows) →
     this is the current motif matrix;
  2. build a `MotifProfile` from the current motifs (pseudocount 1);
  3. in every sequence, take the profile-most-probable *l*-mer (`most_probable_lmer`)
     → new motifs;
  4. if the new motifs have a **higher** `score` (Task 1), keep them and go back to
     step 2; otherwise stop and return the current motifs and their score.

#### `best_of(self, runs)`
- **Output:** `(motifs, score)`, the best result of `runs` calls of
  `randomized_search`

### What to do with it

1. **Check on the lecture data.** `MotifFinder(lecture_dna, 6)` must give the median
   string `AGATAG` with total distance 2 (the answer from the lecture), and
   `best_of(100)` should find motifs with the same consensus and score 34.
2. **Find the planted motif.** Run `median_string()` and `best_of(100)` on
   `planted_motif.fasta` with *l* = 7. Print the *l*-mers that `best_of(100)` picked
   in each sequence and compare their consensus with the median string. Do they
   agree?
3. **How many restarts do you need?** On `planted_motif.fasta` (optionally also on
   the lecture data, see the table below), call `randomized_search()` 200 times. A
   run **succeeds** when the consensus of its motifs (`consensus` from Task 1) equals
   the median string; the fraction of successful runs estimates *p*, the success
   rate of a single run. Then, for R = 5, 10, 20 and 50, repeat `best_of(R)` e.g. 50
   times and measure how often it succeeds in the same sense. Print a small table:
   R, measured success rate, predicted 1 − (1 − *p*)^R. Why does a single run fail
   so often?
4. **From the found motif to a scanner.** The *l*-mers returned by `best_of(100)` in
   step 2 describe the motif. Turn them into a position-specific scoring matrix
   (PSSM: the log-odds weights from the lecture) and scan all sequences with it, on
   **both strands** (imports as in the Task 2 check):
   ```python
   found = motifs.create([Seq(m) for m in best_motifs])  # the motifs from best_of(100)
   found.pseudocounts = 1
   pssm = found.pssm
   threshold = pssm.distribution().threshold_fpr(0.01)
   for sequence in sequences:
       for position, hit_score in pssm.search(Seq(sequence), threshold=threshold):
           ...
   ```
   `threshold_fpr(0.01)` returns the score that a random window exceeds with
   probability 1 % (the false positive rate). `search` reports every window above
   it; a **negative** `position` means a hit on the reverse strand: the window is
   `sequence[len(sequence) + position : len(sequence) + position + l]` and the motif
   matches its reverse complement.
   - How many hits do you get in total? How many of them are the *l*-mers your
     search found, and how many lie on the reverse strand? (A hit on the forward
     strand is a found site if `sequence[position:position + l]` is the *l*-mer that
     `best_of` picked in that sequence.)
   - The sequences contain 10 × 54 windows per strand. How many hits above the 1 %
     threshold would you expect in purely random sequences? What does that say
     about the extra hits?
   - Repeat with `threshold_fpr(0.001)`. What changes, and which threshold would you
     use to report binding sites?

**Expected result:**

|                                                   | lecture data (*l* = 6) | `planted_motif.fasta` (*l* = 7) |
|---------------------------------------------------|------------------------|---------------------------------|
| median string                                     | `AGATAG`               | `GCTAAAG`                       |
| total distance                                    | 2                      | 10                              |
| best score (= *t* · *l* − total distance)         | 34                     | 60                              |
| *p*, one randomized run (random, so only roughly) | about 7–9 %            | about 3–4 %                     |
| `best_of(50)` succeeds (roughly)                  | about 98 %             | about 80–90 %                   |

Scanning `planted_motif.fasta` with the PSSM of the planted sites (step 4):

| threshold              | score     | hits | found sites among them | on the reverse strand | expected by chance |
|------------------------|-----------|------|------------------------|-----------------------|--------------------|
| `threshold_fpr(0.01)`  | about 4.0 | 19   | 10                     | 5                     | about 11           |
| `threshold_fpr(0.001)` | about 6.8 | 10   | 10                     | 0                     | about 1            |

The median string and the scores must match exactly, and so must the scan if
`best_of` found the planted sites. The success rates depend on the random numbers
and on small details of your implementation, so yours will differ a little; the
trend is what matters.
