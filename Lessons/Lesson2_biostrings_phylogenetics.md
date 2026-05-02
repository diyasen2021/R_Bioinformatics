# Lesson 2: From Biostrings to a Phylogenetic Tree

## Overview

In Lesson 1 we fetched a single DNA sequence from NCBI and manually parsed the FASTA format. Now we go further — fetching **multiple sequences**, importing them into `Biostrings`, and following a complete pipeline all the way to a **plotted phylogenetic tree**.

This lesson uses **cytochrome b (cytb)**, a mitochondrial gene that is one of the most widely used markers in molecular phylogenetics. It evolves at a useful rate — fast enough to distinguish between related species, slow enough to retain signal across deeper evolutionary time. We will compare cytb sequences from several mammals.

### What you will build

```
rentrez::entrez_fetch()          # Fetch multiple FASTA sequences
       ↓
Biostrings::readDNAStringSet()   # Parse into a DNAStringSet object
       ↓
Biostrings operations            # Explore, characterise, understand the object
       ↓
DECIPHER::AlignSeqs()            # Multiple sequence alignment
       ↓
ape::as.DNAbin()                 # Convert to DNAbin for phylogenetics
       ↓
ape::dist.dna()                  # Build a pairwise distance matrix
       ↓
ape::nj()                        # Construct a neighbour-joining tree
       ↓
ape::plot.phylo()                # Visualise the tree
```

Each arrow is a **format handoff** — a key bioinformatics skill. Understanding why each package uses its own object type is as important as knowing how to use the functions.

---

## Prerequisites

Install the required packages if you haven't already. Note that `Biostrings` and `DECIPHER` are Bioconductor packages and are installed differently from CRAN packages:

```r
# CRAN packages
install.packages("rentrez")
install.packages("ape")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
```

Then load everything:

```r
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(ape)
```

---

## Step 1 — Fetch Multiple Sequences from NCBI

In Lesson 1 we fetched a single sequence by accession number. Now we fetch several at once. We will use cytochrome b sequences from five mammals: human, chimpanzee, gorilla, orangutan, and mouse (the mouse serves as an **outgroup** — more on that later).

```r
# Accession numbers for cytochrome b in five species
accessions <- c(
  "NC_012920",   # Homo sapiens (human)
  "NC_001643",   # Pan troglodytes (chimpanzee)
  "NC_011120",   # Gorilla gorilla (gorilla)
  "NC_001646",   # Pongo pygmaeus (orangutan)
  "NC_005089"    # Mus musculus (mouse — outgroup)
)

# Fetch all sequences in a single call — rentrez handles multiple IDs
raw <- entrez_fetch(
  db      = "nucleotide",
  id      = accessions,
  rettype = "fasta",
  retmode = "text"
)
```

### Why fetch multiple sequences?

Phylogenetics requires **comparison**. A single sequence tells you nothing about evolutionary relationships — you need at least three sequences to form a tree, and in practice five or more gives a more meaningful result.

### Why include a mouse?

The mouse is our **outgroup** — a species we know is more distantly related to the others. Including an outgroup allows us to **root the tree**, meaning we can determine which end represents the ancestral lineage and which direction evolution has proceeded. Without an outgroup you get an unrooted tree, which shows relationships but not direction.

---

## Step 2 — Parse into a DNAStringSet

As in Lesson 1, `entrez_fetch()` returns a plain character string. We write it to a temporary file and read it with `Biostrings`:

```r
tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)

dna <- readDNAStringSet(tmp)
dna
```

### What is a DNAStringSet?

A `DNAStringSet` is Biostrings' core container for a **collection of DNA sequences**. Think of it as a named list of sequences, but with:

- Memory-efficient storage (sequences are stored as bytes, not characters)
- Built-in biological validation (only valid DNA bases allowed)
- A rich set of vectorised operations that work across all sequences at once

When you print `dna` you will see something like:

```
DNAStringSet object of length 5:
    width seq                                          names
[1]  16569 GATCACAGGTCTATCACCCTATT...  NC_012920.1 Homo sapiens...
[2]  16554 GATCACAGGTCTGTCACCCTATT...  NC_001643.1 Pan troglodytes...
[3]  16364 GATCACAAGCCTGTCACCCCATT...  NC_011120.1 Gorilla gorilla...
[4]  16499 GTTAACAGGTCTATCACCCCTTT...  NC_001646.1 Pongo pygmaeus...
[5]  16299 GTTAATGTAGCTTAATTAATAAAG... NC_005089.1 Mus musculus...
```

Note: these are full mitochondrial genome sequences. We will extract just the cytochrome b gene region shortly.

---

## Step 3 — Explore the DNAStringSet Object

Before moving on, spend time getting to know the object. This is good practice and builds intuition.

```r
# How many sequences?
length(dna)

# What are the names (from the FASTA headers)?
names(dna)

# How long is each sequence?
width(dna)

# Access a single sequence by index — returns a DNAString (singular)
dna[[1]]

# Access by name
dna[["NC_012920.1 Homo sapiens mitochondrion, complete genome"]]

# Subset to get a smaller DNAStringSet (keeps the container type)
dna[1:3]
```

### The difference between `[` and `[[`

This is an important R concept that Biostrings follows consistently:

- `dna[1]` returns a **DNAStringSet** of length 1 — the container is preserved
- `dna[[1]]` returns a **DNAString** — a single sequence, the element itself

This mirrors how standard R lists behave. When you want to do operations on one sequence as a sequence (e.g. translate it), use `[[`. When you want a subset of sequences to pass to another function, use `[`.

---

## Step 4 — Clean Up the Sequence Names

The FASTA headers from NCBI are long and unwieldy. Tidy them now — short, clean names make the final tree readable.

```r
# Assign short species names
names(dna) <- c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse")
dna
```

Always do this **before** alignment and tree building. The names travel with the sequences through every subsequent step, and they become your tip labels on the final tree.

---

## Step 5 — Core Biostrings Operations

Before moving to alignment, let's use some of Biostrings' key functions. These work on the `DNAStringSet` as a whole — vectorised across all sequences simultaneously.

### Base composition

```r
# Count each base across all sequences
letterFrequency(dna, letters = c("A", "T", "G", "C"))
```

This returns a matrix — one row per sequence, one column per base. Useful for spotting data quality issues: a sequence with unexpected base frequencies may be misidentified or contaminated.

### GC content

```r
# GC content as a proportion (0 to 1)
letterFrequency(dna, letters = "GC", as.prob = TRUE)
```

Mitochondrial GC content varies across species and is itself phylogenetically informative.

### Reverse complement

```r
# Compute reverse complement of all sequences
reverseComplement(dna)
```

You often need this when sequences are reported on different strands. Biostrings handles the complement (A↔T, G↔C) and the reversal in one step.

### Alphabet frequency check

```r
# Are there any ambiguous bases (N, R, Y etc.)?
alphabetFrequency(dna)
```

Columns beyond A, C, G, T show ambiguous or degenerate bases. High counts of N often indicate low-quality sequencing regions.

---

## Step 6 — Extract the Cytochrome b Region

We fetched full mitochondrial genomes. Cytochrome b occupies a specific region within the mitochondrial genome. For simplicity, we extract a comparable region from each sequence.

In a real analysis you would look up the exact coordinates from the GenBank annotation. For this tutorial we use approximate coordinates that capture the cytb gene across these species:

```r
# Approximate cytb coordinates in mammalian mitochondrial genomes
# These vary slightly per species — in a real study check the annotation
cytb <- subseq(dna, start = 14747, end = 15887)
width(cytb)   # Should be 1141 bp for each sequence
```

### Why extract a gene region?

Aligning complete mitochondrial genomes is computationally expensive and complicated by insertions, deletions, and rearrangements across species. Aligning a single gene is simpler and more appropriate when your question is about that gene's evolutionary history.

---

## Step 7 — Multiple Sequence Alignment with DECIPHER

Now we reach the first major handoff. Before building a phylogenetic tree we need to **align** the sequences — arrange them so that homologous positions (positions derived from the same ancestral position) are in the same column.

```r
# AlignSeqs works directly with DNAStringSet
cytb_aligned <- AlignSeqs(cytb)
cytb_aligned
```

### Why do we need alignment?

An unaligned set of sequences is just a collection of strings. Position 100 in one sequence is not necessarily homologous to position 100 in another — especially if there have been insertions or deletions during evolution. Alignment establishes **positional homology**, which is the foundation all phylogenetic methods require.

### Why DECIPHER?

`DECIPHER` is beginner-friendly, well-documented, and performs very well on biological sequences. It works directly with `DNAStringSet` objects, so there is no format conversion needed at this step. Other aligners you may encounter include `msa` (which wraps ClustalW, MUSCLE, and ClustalOmega) and `muscle` — all good options, but DECIPHER requires the least setup for in-R workflows.

### Inspect the alignment

```r
# View the alignment — columns now represent homologous positions
BrowseSeqs(cytb_aligned)   # Opens an interactive viewer in your browser
```

In the aligned output, `-` characters represent **gaps** — positions where one sequence has a base but another does not (due to an insertion or deletion event during evolution). This is normal and expected.

---

## Step 8 — Convert to DNAbin for Phylogenetics

Here is the second major handoff. `ape` and the broader phylogenetics ecosystem in R use their own sequence object type: `DNAbin`. We need to convert from Biostrings' `DNAStringSet` to `ape`'s `DNAbin`.

```r
# Convert aligned DNAStringSet to DNAbin
cytb_bin <- as.DNAbin(cytb_aligned)
class(cytb_bin)   # "DNAbin"
cytb_bin
```

### Why does this conversion exist?

`DNAStringSet` and `DNAbin` store sequences differently:

- `DNAStringSet` stores sequences as raw bytes optimised for **fast string operations** (pattern matching, translation, etc.)
- `DNAbin` stores sequences as raw bytes optimised for **distance and phylogenetic calculations** using the codes defined by the IUPAC standard as implemented in `ape`

Neither format is universally better — they are optimised for different tasks. Understanding that this conversion is necessary (and why) is a core bioinformatics competency. You will encounter similar handoffs throughout R-based genomics workflows.

---

## Step 9 — Build a Distance Matrix

Phylogenetic tree construction begins with quantifying how different each pair of sequences is. We compute a **pairwise distance matrix** using a substitution model.

```r
# Compute pairwise distances using the Tamura-Nei 93 model
dist_matrix <- dist.dna(cytb_bin, model = "TN93")
dist_matrix
```

### What is a distance matrix?

A distance matrix is a table where each cell contains the estimated evolutionary distance between two sequences. Higher values mean more differences have accumulated — implying either more evolutionary time has passed, or faster evolution, or both.

### What is the TN93 model?

TN93 (Tamura-Nei 1993) is a **nucleotide substitution model** that accounts for:

- Different rates of transitions (A↔G, C↔T) vs transversions (A↔C, A↔T, G↔C, G↔T)
- Different base frequencies

It is a good general-purpose model for mitochondrial sequences. Other common models include `K80` (simpler, assumes equal base frequencies) and `GTR` (most flexible, but requires more data to estimate reliably). In a rigorous study you would use model selection tools (like `modeltest` or IQ-TREE's ModelFinder) to choose the best model for your data.

---

## Step 10 — Build the Phylogenetic Tree

With a distance matrix in hand, we can construct the tree using the **neighbour-joining** algorithm:

```r
# Construct a neighbour-joining tree
tree <- nj(dist_matrix)
tree
```

### What is neighbour-joining?

Neighbour-joining (NJ) is a **distance-based** tree-building algorithm. It works by:

1. Starting with all sequences in a "star" topology (all connected to one central node)
2. Iteratively joining the pair of sequences that are most similar to each other (the "neighbours")
3. Repeating until all sequences are placed on the tree

NJ is fast, widely used, and appropriate for teaching and exploratory analysis. For publishable work you would typically use maximum likelihood (e.g. with `phangorn` or IQ-TREE) or Bayesian inference (e.g. MrBayes), which are more statistically rigorous but more computationally demanding.

---

## Step 11 — Root the Tree

An unrooted tree shows relationships but not direction. We root the tree using our outgroup — the mouse:

```r
# Root the tree on the mouse (outgroup)
tree_rooted <- root(tree, outgroup = "Mouse", resolve.root = TRUE)
```

### Why root on the mouse?

We know from broader biological knowledge that rodents (Mus) diverged from the primate lineage much earlier than any of the primates diverged from each other. By telling `root()` that Mouse is the outgroup, we are saying: "the root of this tree — the point representing the common ancestor of all five species — lies on the branch leading to Mouse." This gives the tree a direction and allows us to read it as a hypothesis about evolutionary history.

---

## Step 12 — Plot the Tree

Now for the payoff — visualising the tree:

```r
# Basic plot
plot(tree_rooted,
     main  = "Cytochrome b Phylogeny — Great Apes + Mouse",
     cex   = 1.2)    # cex controls tip label size

# Add a scale bar showing branch length (substitutions per site)
add.scalebar()
```

### Reading the tree

- **Tip labels** are your species names
- **Branch lengths** represent the amount of evolutionary change (substitutions per site) along each branch — longer branches mean more change
- **Nodes** (where branches meet) represent hypothetical common ancestors
- Species that share a more recent common ancestor cluster together — you should see Human, Chimpanzee, and Gorilla forming one cluster, with Orangutan slightly more distant, and Mouse as the outgroup

### Adding bootstrap support (optional but recommended)

Bootstrap values tell you how well-supported each node is — how often that grouping appeared when the analysis was repeated on resampled data:

```r
# Bootstrap the tree (100 replicates — use 1000 for publication)
boot_tree <- boot.phylo(
  tree_rooted,
  cytb_bin,
  function(x) root(nj(dist.dna(x, model = "TN93")), "Mouse", resolve.root = TRUE),
  B = 100
)

# Plot with bootstrap values at nodes
plot(tree_rooted, main = "Cytochrome b Phylogeny with Bootstrap Support", cex = 1.2)
nodelabels(boot_tree, cex = 0.8, frame = "none", adj = c(1.2, -0.3))
add.scalebar()
```

Bootstrap values range from 0 to 100. Values above 70 are generally considered well-supported; above 95 is strong support.

---

## The Full Pipeline

Here is the complete workflow from fetch to tree in one place:

```r
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(ape)

# 1. Fetch
accessions <- c("NC_012920", "NC_001643", "NC_011120", "NC_001646", "NC_005089")
raw <- entrez_fetch(db = "nucleotide", id = accessions, rettype = "fasta", retmode = "text")

# 2. Parse
tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)
dna <- readDNAStringSet(tmp)

# 3. Name
names(dna) <- c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse")

# 4. Explore
width(dna)
letterFrequency(dna, letters = "GC", as.prob = TRUE)

# 5. Extract cytb region
cytb <- subseq(dna, start = 14747, end = 15887)

# 6. Align
cytb_aligned <- AlignSeqs(cytb)

# 7. Convert to DNAbin
cytb_bin <- as.DNAbin(cytb_aligned)

# 8. Distance matrix
dist_matrix <- dist.dna(cytb_bin, model = "TN93")

# 9. Build tree
tree <- nj(dist_matrix)

# 10. Root on outgroup
tree_rooted <- root(tree, outgroup = "Mouse", resolve.root = TRUE)

# 11. Plot
plot(tree_rooted, main = "Cytochrome b Phylogeny", cex = 1.2)
add.scalebar()
```

---

## Understanding the Format Handoffs

The most important conceptual lesson in this pipeline is that **different packages use different object types**, and converting between them is a core skill:

| Step | Object type | Package | Optimised for |
|---|---|---|---|
| Raw fetch | `character` | base R | Text storage |
| After import | `DNAStringSet` | Biostrings | Fast string operations |
| After alignment | `DNAStringSet` | DECIPHER / Biostrings | Aligned string operations |
| After conversion | `DNAbin` | ape | Distance and phylogenetic calculations |
| After tree building | `phylo` | ape | Tree manipulation and plotting |

Each format exists for a reason. When you see an error like "expected DNAbin, got DNAStringSet", you now know exactly what to do.

---

## Where to Go Next

This lesson used neighbour-joining, which is distance-based and fast. More rigorous approaches you can explore from here:

- **Maximum likelihood trees** — `phangorn::phyml()` or the external tool IQ-TREE, which can be called from R
- **Bayesian inference** — MrBayes or BEAST (time-calibrated trees), both external tools with R interfaces
- **Model selection** — `phangorn::modelTest()` to choose the best substitution model for your data
- **Tree annotation** — `ggtree` (Bioconductor) for publication-quality, ggplot2-style tree visualisation
- **More sequences** — try adding sequences from more species, or try a different gene marker

---

## Summary

| Concept | Why it matters |
|---|---|
| Fetch multiple sequences in one call | Phylogenetics requires comparison across sequences |
| Include an outgroup | Allows you to root the tree and read evolutionary direction |
| Alignment before tree building | Establishes positional homology — required by all phylogenetic methods |
| Format handoffs (DNAStringSet → DNAbin → phylo) | Different packages are optimised for different tasks |
| Substitution models | Account for the fact that observed differences underestimate true evolutionary change |
| Bootstrap support | Quantifies confidence in each node of the tree |
| Neighbour-joining is a starting point | Fast and useful for exploration; maximum likelihood or Bayesian methods for rigorous analysis |
