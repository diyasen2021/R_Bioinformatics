# Lesson 3: Alignment, Multiple Sequence Alignment, and the Foundation of Phylogenetics

## Overview

This lesson covers one of the most fundamental concepts in computational biology: **sequence alignment**. We will answer three big questions:

1. What is phylogenetics and why does it matter?
2. What problem does alignment solve, and why does it exist?
3. How do the major alignment algorithms work, and how do you choose between them?

By the end you will understand not just *how* to run an alignment in R, but *why* alignment is necessary, what the algorithms are actually doing, and how alignment quality directly determines the quality of any downstream analysis — including the phylogenetic trees we built in Lesson 2.

---

## Part 1 — What is Phylogenetics?

### The core question

Phylogenetics is the study of **evolutionary relationships among organisms**. It asks: given a set of species (or genes, or populations), how are they related? Who shares a more recent common ancestor with whom? How long ago did lineages diverge?

These questions matter enormously across biology:

- **Medicine** — tracing the spread and evolution of pathogens (HIV, SARS-CoV-2, influenza)
- **Conservation** — understanding which populations are distinct enough to warrant separate protection
- **Drug discovery** — identifying conserved functional regions across species
- **Forensics** — identifying species from trace biological material
- **Ecology** — understanding how communities assembled over time

### Trees as hypotheses

A phylogenetic tree is a **hypothesis** about evolutionary history — not a fact. It is the best explanation the available data supports, given a particular method and model. Different data, different models, or different methods can produce different trees. This is not a flaw — it is how science works. A good phylogenetics paper always reports which methods were used, which model was selected, and what the support values are at each node.

### What a tree represents

```
        ┌── Human
     ┌──┤
     │  └── Chimpanzee
  ───┤
     │  ┌── Gorilla
     └──┤
        └── Orangutan
```

- **Tips** (leaf nodes) — the sequences or taxa you sampled
- **Internal nodes** — hypothetical common ancestors that no longer exist
- **Branches** — lineages through time; branch *length* represents evolutionary change
- **Root** — the most ancient common ancestor of all taxa in the tree
- **Topology** — the branching pattern (who groups with whom), independent of branch lengths

### Molecular phylogenetics

Before DNA sequencing, phylogenies were built from morphological characters — bone shape, body plan, behaviour. Molecular phylogenetics uses **DNA, RNA, or protein sequences** instead. The logic is elegant:

- Sequences change over time through mutation
- Related species share more recent common ancestors, so they have had less time to accumulate differences
- Therefore, **more similar sequences = more closely related**

But — and this is crucial — you can only compare sequences meaningfully if you know which positions are **homologous**: descended from the same position in the ancestral sequence. This is exactly what alignment provides.

### The molecular clock

Mutations accumulate at roughly predictable rates in certain genes. This means sequence divergence is approximately proportional to time since divergence — a **molecular clock**. This allows phylogenetics to do more than just determine topology; it can estimate *when* lineages diverged, in years, if the clock is calibrated with fossil data. This is the basis of divergence time estimation and the entire field of molecular dating.

---

## Part 2 — Why Alignment Exists: The Problem It Solves

### Evolution changes sequence length

Consider two sequences descended from a common ancestor:

```
Ancestor:   ATGCCTGAACGT
Species A:  ATGCCTGAACGT       (no change)
Species B:  ATGCC--GAACGT      (deletion of TG in Species A... or insertion in B?)
```

Over evolutionary time, mutations do not just substitute one base for another — they also **insert** or **delete** bases (collectively called **indels**). This means that after divergence, two related sequences may be different lengths, and position 7 in one sequence may not correspond to position 7 in the other.

If you simply lined up the sequences by position without accounting for indels:

```
Species A:  A T G C C T G A A C G T
Species B:  A T G C C G A A C G T -
```

You would incorrectly infer that position 6 changed from T to G, position 7 changed from G to A, and so on — when actually a single deletion event is responsible. Your entire downstream analysis would be built on a false picture of what happened.

### Homology: the concept alignment operationalises

**Homology** means "descended from a common ancestral character". In sequence terms, two positions are homologous if they both descended from the same position in the ancestral sequence. Alignment is the computational process of **identifying homologous positions** and arranging them in the same column.

```
Aligned:
Species A:  A T G C C T G A A C G T
Species B:  A T G C C - - G A A C G T
                    ↑ ↑
                    gaps inserted to maintain homology
```

Now we can correctly see that positions 1-5 are identical, an indel occurred at positions 6-7, and positions 8-12 are identical again.

### The gap character

In aligned sequences, `-` represents a **gap** — a position where one sequence has a base and another does not, due to an insertion or deletion event. Gaps are not bases; they represent evolutionary events (indels). In phylogenetic analysis, gaps can be treated as missing data or as a fifth character state, depending on the method.

### Alignment is a biological hypothesis

This is perhaps the most important conceptual point in this lesson: **an alignment is not just a computation — it is a biological claim**. When you align two sequences and place a gap at position 47, you are hypothesising that an insertion or deletion occurred at that position in one lineage. The alignment encodes your beliefs about evolutionary history. A wrong alignment produces a wrong tree, no matter how sophisticated your tree-building method is. This is why biologists inspect their alignments visually, trim poorly aligned regions, and sometimes align manually in difficult cases.

---

## Part 3 — Pairwise Alignment

Pairwise alignment compares **two sequences** at a time. It is the foundation on which multiple sequence alignment is built, and it is also directly useful for tasks like finding where a query sequence matches a reference, or comparing two protein domains.

### The two types of pairwise alignment

#### Global alignment — Needleman-Wunsch (1970)

Global alignment aligns sequences **across their entire length**. It finds the single best way to align sequence A to sequence B from end to end.

Use global alignment when:
- Both sequences are approximately the same length
- You expect the sequences to be homologous across their full length
- Example: comparing the same gene from two different species

```
Sequence A:  ATGCCTGAACGT
Sequence B:  ATGCCTAACGT

Global:
A:  ATGCCT-GAACGT
B:  ATGCCT-AACGT-
        ↑       ↑
       gap      gap
```

#### Local alignment — Smith-Waterman (1981)

Local alignment finds the **best matching subsequence** between two sequences — it does not try to align the full length of both. It identifies the region of highest similarity.

Use local alignment when:
- Sequences are different lengths
- You are looking for a shared domain or motif within larger sequences
- One sequence is a short query (e.g. a primer or probe) and the other is a long reference genome
- Example: BLAST uses a heuristic version of Smith-Waterman

```
Sequence A (long):  GGGGATGCCTGAACGTCCCC
Sequence B (short): ATGCCTGAA

Local:
A:  ATGCCTGAA   (best local match extracted from A)
B:  ATGCCTGAA
    |||||||||
    Score: 9 matches
```

### Scoring: how alignment algorithms make decisions

Both algorithms work by assigning scores to alignment choices. Every position in the alignment must be one of:

- **Match** — same base in both sequences → positive score (reward)
- **Mismatch** — different bases → negative score (penalty)
- **Gap opening** — introducing a new gap → large negative penalty
- **Gap extension** — extending an existing gap → smaller negative penalty (gaps tend to occur in runs)

The algorithm uses **dynamic programming** — building up a matrix of optimal sub-alignment scores — to find the globally or locally optimal alignment given the scoring scheme. The key insight of dynamic programming is that the optimal alignment of the full sequence can be built from the optimal alignments of smaller sub-problems.

### Gap penalties: a biological choice

The gap penalties you choose encode biological assumptions. A high gap opening penalty assumes indels are rare; a lower penalty assumes they are common. For protein-coding sequences, you might penalise gaps that break the reading frame more heavily. There is no universally correct scoring scheme — it depends on the biology of your sequences.

### Pairwise alignment in R

```r
library(Biostrings)

seq_a <- DNAString("ATGCCTGAACGT")
seq_b <- DNAString("ATGCCTAACGT")

# Global alignment (Needleman-Wunsch)
global <- pairwiseAlignment(seq_a, seq_b, type = "global")
global
score(global)        # Alignment score
aligned(global)      # The aligned sequences
nmatch(global)       # Number of matches
nmismatch(global)    # Number of mismatches
nindel(global)       # Number of indels

# Local alignment (Smith-Waterman)
local <- pairwiseAlignment(seq_a, seq_b, type = "local")
local
score(local)

# You can customise the scoring scheme
custom_score <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
global_custom <- pairwiseAlignment(
  seq_a, seq_b,
  type            = "global",
  substitutionMatrix = custom_score,
  gapOpening      = -2,
  gapExtension    = -0.5
)
```

### Visualising a pairwise alignment

```r
# Print a readable alignment view
writePairwiseAlignments(global)
```

This prints the classic alignment view with match bars, showing you exactly where matches, mismatches, and gaps occur.

---

## Part 4 — Multiple Sequence Alignment (MSA)

### Why pairwise is not enough for phylogenetics

Pairwise alignment compares two sequences. Phylogenetics needs to compare many sequences simultaneously — and all comparisons must be **consistent with each other**. If A aligns to B one way, and A aligns to C another way, you cannot simply stack those alignments on top of each other. You need a single alignment where every column represents a homologous position across all sequences simultaneously. This is **Multiple Sequence Alignment (MSA)**.

### The computational challenge

Finding the mathematically optimal MSA for N sequences is an **NP-hard problem** — the computation time grows so fast with the number of sequences that an exact solution becomes impossible for even moderate N. For 10 sequences the exact solution is already impractical; for 100 it is completely out of reach.

This is why all practical MSA tools use **heuristics** — algorithms that find good-enough solutions quickly, sacrificing the guarantee of optimality for computational tractability. Different heuristics make different trade-offs between speed and accuracy.

### Progressive alignment: the dominant strategy

Most MSA algorithms use **progressive alignment**:

1. **Compute pairwise distances** between all sequence pairs
2. **Build a guide tree** — a rough tree based on pairwise similarity (not a phylogenetic result — just a computational scaffold)
3. **Align sequences progressively** — starting with the most similar pair, then adding more sequences in order of decreasing similarity, like building up a multiple alignment from the leaves of the guide tree inward

The weakness of progressive alignment is **error propagation**: mistakes made in early pairwise alignments get locked in and carried forward. More distantly related sequences added later may be forced into a suboptimal alignment because the early decisions cannot be revisited. Iterative refinement (used by MUSCLE and MAFFT) addresses this by going back and realigning after the initial progressive pass.

---

## Part 5 — The Major MSA Algorithms

### ClustalW / ClustalOmega

**ClustalW** (1994) was for decades the most widely used MSA tool in biology, and is still one of the most cited papers in all of science. It uses progressive alignment with a neighbour-joining guide tree.

**ClustalOmega** (2011) is the modern successor — much faster, scales to millions of sequences, and generally more accurate. It uses a different guide tree method (mBed) and HMM-based profile alignment.

**When to use:** ClustalOmega is a solid, well-understood default. Its ubiquity means results are easily compared with published work. Well-supported in R via the `msa` package.

**Weakness:** Like all progressive aligners, early errors can propagate. Can struggle with highly divergent sequences.

### MUSCLE

**MUSCLE** (Multiple Sequence Comparison by Log-Expectation, 2004) improved on ClustalW by adding **iterative refinement** — after the initial progressive alignment, it goes back and realigns subgroups to improve the overall score. Generally faster and more accurate than ClustalW on benchmark datasets.

**When to use:** A good default when you want something better than ClustalW without the complexity of MAFFT. Available in R via the `msa` package.

### MAFFT

**MAFFT** (Multiple Alignment using Fast Fourier Transform, 2002, with many updates since) uses Fourier transform methods to rapidly identify homologous regions. It offers multiple alignment strategies of increasing accuracy and computation time:

- `FFT-NS-2` — fast, good for large datasets (>1000 sequences)
- `G-INS-i` — accurate, for sequences with a single global homologous region (e.g. a single gene)
- `L-INS-i` — very accurate, for sequences with one or a few alignable regions embedded in ulignable flanking sequence
- `E-INS-i` — for sequences with multiple independent alignable regions

MAFFT is currently the most widely recommended MSA tool for phylogenetics. It is extremely well benchmarked and handles both closely and distantly related sequences well.

**When to use:** Large datasets, divergent sequences, or when you need the best accuracy. Typically run externally (command line) and results imported into R, though an R wrapper exists.

### DECIPHER

**DECIPHER** is an R-native aligner built specifically for the Bioconductor ecosystem. It uses a database-driven approach and is particularly good at:

- Working directly with `DNAStringSet` objects (no format conversion)
- Aligning sequences within an R workflow without external tools
- Large ribosomal RNA datasets

**When to use:** When you want to stay entirely within R. Excellent for teaching and for workflows where you do not want to install external command-line tools.

### Summary comparison

| Tool | Speed | Accuracy | Best for | R package |
|---|---|---|---|---|
| ClustalOmega | Fast | Good | General use, well-cited | `msa` |
| MUSCLE | Fast | Very good | General use, better than ClustalW | `msa` |
| MAFFT | Fast–slow (mode-dependent) | Excellent | Large or divergent datasets | External / `msa` |
| DECIPHER | Moderate | Good | In-R workflows, rRNA | `DECIPHER` |

---

## Part 6 — MSA in R: Hands-On

### Using DECIPHER (recommended for in-R workflows)

```r
library(rentrez)
library(Biostrings)
library(DECIPHER)

# Fetch cytochrome b sequences (from Lesson 2)
accessions <- c("NC_012920", "NC_001643", "NC_011120", "NC_001646", "NC_005089")
raw <- entrez_fetch(db = "nucleotide", id = accessions, rettype = "fasta", retmode = "text")
tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)
dna <- readDNAStringSet(tmp)
names(dna) <- c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse")

# Extract cytb region
cytb <- subseq(dna, start = 14747, end = 15887)

# Align with DECIPHER
cytb_aligned <- AlignSeqs(cytb)

# View the alignment interactively
BrowseSeqs(cytb_aligned)

# View as a matrix of characters (useful for inspection)
as.matrix(cytb_aligned)
```

### Using the msa package (ClustalOmega, MUSCLE)

```r
library(msa)

# ClustalOmega
aligned_clustal <- msa(cytb, method = "ClustalOmega")
aligned_clustal

# MUSCLE
aligned_muscle <- msa(cytb, method = "Muscle")

# Print a readable alignment
msaPrettyPrint(aligned_clustal, output = "asis")

# Convert to DNAStringSet for downstream use
aligned_ss <- as(aligned_clustal, "DNAStringSet")
```

### Comparing aligners

For important analyses it is worth running two or three aligners and comparing the results. If alignments differ substantially in a region, that region is **ambiguously aligned** and should be treated with caution or trimmed.

```r
# Align with both
aln_decipher <- AlignSeqs(cytb)
aln_muscle   <- as(msa(cytb, method = "Muscle"), "DNAStringSet")

# Check alignment widths (total columns including gaps)
width(aln_decipher)[1]   # Should all be the same within each alignment
width(aln_muscle)[1]

# Different total widths between aligners suggests different gap placement
# Inspect both with BrowseSeqs
BrowseSeqs(aln_decipher)
BrowseSeqs(aln_muscle)
```

---

## Part 7 — Alignment Quality and Its Impact on Phylogenetics

### Garbage in, garbage out

This cannot be overstated: **a phylogenetic tree is only as good as the alignment it was built from**. An incorrect alignment will produce an incorrect tree, regardless of how sophisticated the tree-building method is. Many published phylogenetic errors trace back to poor alignment rather than poor tree inference.

### What to look for in an alignment

When you inspect an alignment visually (use `BrowseSeqs()` or export to a viewer like Jalview), look for:

- **Conserved columns** — positions that are identical across all sequences. These anchor the alignment and are most trustworthy
- **Variable columns** — positions that differ across sequences. These carry phylogenetic signal
- **Gap-rich regions** — stretches with many gaps suggest a difficult-to-align region. These may introduce noise into phylogenetic analysis
- **Terminal gaps** — gaps at the start or end of sequences, often due to sequences of different lengths. Usually safe to trim
- **Internal gap clusters** — dense gap regions in the middle of an alignment are a warning sign of poor alignment or genuine indel-rich regions

### Trimming: removing poorly aligned regions

It is common practice to **trim** an alignment before phylogenetic analysis — removing columns that are poorly aligned or contain too many gaps. Two popular tools for this are:

- **trimAl** — automated trimming with several heuristic modes; can be called from R
- **Gblocks** — selects conserved, less divergent blocks of the alignment

In R, you can do basic trimming manually:

```r
# Convert alignment to a matrix for inspection
aln_matrix <- as.matrix(cytb_aligned)

# Find columns where more than 50% of sequences have a gap
gap_freq <- colMeans(aln_matrix == "-")
keep_cols <- gap_freq < 0.5

# Trim the alignment
aln_trimmed <- DNAStringSet(apply(aln_matrix[, keep_cols], 1, paste, collapse = ""))
names(aln_trimmed) <- names(cytb_aligned)
width(aln_trimmed)   # Should be shorter than the original
```

### Conserved vs variable regions

Not all positions in an alignment are equally informative for phylogenetics:

- **Invariant positions** — identical across all sequences. Contribute no information about relationships (though they do contribute to branch length estimation)
- **Parsimony-informative positions** — variable, but shared by at least two sequences. These are the positions that actually drive tree topology inference
- **Autapomorphies** — unique to one sequence. Inform branch lengths but not topology

A gene with too few variable positions gives a poorly resolved tree. A gene with too many differences (saturated) loses signal because multiple mutations at the same position obscure history. Cytochrome b, our example gene, sits in a useful middle ground for mammalian phylogenetics.

---

## Part 8 — From Alignment to Phylogenetics: The Conceptual Bridge

Now that you understand alignment, we can revisit the phylogenetics pipeline from Lesson 2 with much deeper understanding.

### Why aligned sequences are required

Every phylogenetic method — distance-based, maximum parsimony, maximum likelihood, Bayesian — needs to compare positions across sequences. That comparison only makes sense if each column of the alignment represents a homologous position. The alignment *is* the data matrix that phylogenetic methods operate on.

### What different tree-building methods do with the alignment

| Method | What it does with the alignment | Strength |
|---|---|---|
| Distance (NJ) | Computes pairwise distances, discards the alignment | Fast, simple |
| Maximum parsimony | Finds the tree requiring fewest mutations | Intuitive, but inconsistent in some conditions |
| Maximum likelihood | Finds the tree that makes the observed alignment most probable under a model | Statistically rigorous |
| Bayesian inference | Estimates a distribution of trees given the alignment and a prior | Full uncertainty quantification |

In Lesson 2 we used neighbour-joining — fast and appropriate for learning. Maximum likelihood (e.g. with `phangorn` in R, or IQ-TREE externally) is the current standard for published work.

### The substitution model connects alignment to tree

When we computed `dist.dna(cytb_bin, model = "TN93")` in Lesson 2, the model (TN93) was our assumption about how mutations occur. Substitution models correct for **multiple hits** — the fact that a position may have mutated more than once since divergence, with later mutations overwriting earlier ones. Without this correction, sequence distances are underestimated, and branch lengths in the tree will be too short.

---

## The Full Alignment-to-Tree Pipeline

```r
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(ape)

# 1. Fetch
accessions <- c("NC_012920", "NC_001643", "NC_011120", "NC_001646", "NC_005089")
raw <- entrez_fetch(db = "nucleotide", id = accessions, rettype = "fasta", retmode = "text")
tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)
dna <- readDNAStringSet(tmp)
names(dna) <- c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse")

# 2. Extract gene region
cytb <- subseq(dna, start = 14747, end = 15887)

# 3. Align — this is the critical step
cytb_aligned <- AlignSeqs(cytb)

# 4. Inspect the alignment — do not skip this
BrowseSeqs(cytb_aligned)

# 5. Trim poorly aligned columns (optional but recommended)
aln_matrix <- as.matrix(cytb_aligned)
gap_freq    <- colMeans(aln_matrix == "-")
cytb_trimmed <- DNAStringSet(
  apply(aln_matrix[, gap_freq < 0.5], 1, paste, collapse = "")
)
names(cytb_trimmed) <- names(cytb_aligned)

# 6. Convert to DNAbin
cytb_bin <- as.DNAbin(cytb_trimmed)

# 7. Distance matrix
dist_matrix <- dist.dna(cytb_bin, model = "TN93")

# 8. Tree
tree <- nj(dist_matrix)
tree_rooted <- root(tree, outgroup = "Mouse", resolve.root = TRUE)

# 9. Plot
plot(tree_rooted, main = "Cytochrome b — Great Apes", cex = 1.2)
add.scalebar()
```

---

## Summary

| Concept | Why it matters |
|---|---|
| Phylogenetics asks about evolutionary relationships | The tree is a hypothesis, not a fact |
| Sequences encode evolutionary history | More similar = more recently diverged |
| Alignment solves the indel problem | Without it, position comparisons are meaningless |
| Alignment is a biological hypothesis | A wrong alignment produces a wrong tree |
| Global vs local alignment | Global for full-length homologs; local for finding shared regions |
| Dynamic programming guarantees optimal pairwise alignment | But is too slow for many sequences — hence heuristic MSA |
| Progressive alignment is the dominant MSA strategy | Build a guide tree, align most similar first |
| DECIPHER for in-R workflows; MAFFT for large/divergent datasets | Choose based on your data and workflow |
| Always inspect and trim your alignment | Poorly aligned regions add noise, not signal |
| Substitution models correct for multiple hits | Raw distances underestimate true evolutionary divergence |

---

## Where to Go Next

- **Lesson 4:** Substitution models in depth — what they are, how to choose one, and model selection with `phangorn::modelTest()`
- **Lesson 5:** Maximum likelihood trees with `phangorn` — moving beyond neighbour-joining
- **Lesson 6:** Bootstrap, jackknife, and Bayesian posterior probabilities — quantifying confidence in your tree
- **Further reading:** The original Needleman-Wunsch (1970) and Smith-Waterman (1981) papers are surprisingly readable. The MUSCLE (Edgar 2004) and MAFFT (Katoh et al. 2002) papers explain their algorithms clearly and are worth reading once you are comfortable with the basics.
