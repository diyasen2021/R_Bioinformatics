# Reading DNA Sequences from NCBI into R: A Beginner's Tutorial

## Introduction

This tutorial walks you through how to fetch DNA sequence data from NCBI (the National Center for Biotechnology Information) using the `rentrez` package in R, and how to work with that data manually — step by step. 

Understanding the manual approach is just as important as using convenience packages. It helps you:

- Understand **what FASTA format actually is** — not magic, just text with a convention
- Debug problems when packages fail or return unexpected results
- Handle **non-standard or malformed sequences** you may encounter in the wild
- Build intuition for **what higher-level packages like `Biostrings` are doing internally**

Once you've done it manually a few times, using `readDNAStringSet()` stops feeling like magic and starts feeling like a sensible shortcut.

---

## Prerequisites

Install and load the `rentrez` package if you haven't already:

```r
install.packages("rentrez")
library(rentrez)
```

---

## Step 1 — Fetch the Raw Sequence

Use `entrez_fetch()` to retrieve a sequence from NCBI. Here we fetch the BRCA1 gene using its accession number:

```r
raw <- entrez_fetch(db = "nucleotide", id = "NM_007294", rettype = "fasta", retmode = "text")
```

### What is this doing?

| Argument | What it means |
|---|---|
| `db = "nucleotide"` | Which NCBI database to search — here the nucleotide database |
| `id = "NM_007294"` | The accession number of the sequence you want |
| `rettype = "fasta"` | Return the data in FASTA format (plain text sequence) |
| `retmode = "text"` | Return as plain text (not XML or JSON) |

### What does it return?

`entrez_fetch()` returns a **plain character string** — it is *not* a special bioinformatics object. It is just `character(1)`, a single string of text, exactly like the contents of a `.fasta` file you might download manually from NCBI.

---

## Step 2 — Inspect the Raw Object

Before you do anything with data, **always look at it first**. This is good practice for any data science or bioinformatics work.

```r
class(raw)     # What type of object is it?
length(raw)    # How many elements does the vector have?
nchar(raw)     # How many total characters (including newlines)?

# Print just the first 200 characters so you can see the structure
substr(raw, 1, 200)
```

### What you'll see

The output will look something like this:

```
>NM_007294.4 Homo sapiens BRCA1 DNA repair associated (BRCA1), transcript variant 1, mRNA
ATGGATTTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGA
GTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAGATT...
```

This is **FASTA format**. It has two parts:

1. A **header line** starting with `>` containing the identifier and description
2. The **sequence itself**, split across multiple lines (usually 60 or 70 characters wide)

The sequence is split across lines purely for human readability — biologically it is one continuous string. This is a crucial thing to understand.

---

## Step 3 — Split the String into Lines

Right now, the entire sequence (header and all) is packed into a single string with newline characters (`\n`) separating the lines. To work with it, we need to break it apart.

```r
lines <- strsplit(raw, "\n")[[1]]
```

### Why `[[1]]`?

`strsplit()` always returns a **list**, even when you split a single string. The `[[1]]` extracts the character vector from inside that list. This is a common R gotcha — without `[[1]]`, you'd be working with a list instead of a plain vector.

### Inspect the result

```r
head(lines, 5)   # Look at the first few lines — should see the header
tail(lines, 5)   # Look at the last few — often there's a trailing empty line
length(lines)    # How many lines in total?
```

You'll now have a character vector where each element is one line of the FASTA file.

---

## Step 4 — Separate the Header from the Sequence

Now we can use the defining feature of FASTA format — the header always starts with `>` — to split the data into its two logical parts.

```r
# The header is the line(s) starting with ">"
header <- lines[startsWith(lines, ">")]
header

# Everything else is sequence — filter out blank lines too
seq_lines <- lines[!startsWith(lines, ">") & nchar(lines) > 0]
head(seq_lines)
```

### Why filter blank lines?

FASTA files frequently have trailing empty lines at the end. If you don't remove them, they'll introduce empty strings into your sequence when you collapse the lines together in the next step. The condition `nchar(lines) > 0` removes any line with zero characters.

### What you have now

- `header`: a character vector with the identifier line(s) — e.g. `">NM_007294.4 Homo sapiens BRCA1..."`
- `seq_lines`: a character vector where each element is a 60-character chunk of the sequence

---

## Step 5 — Collapse the Sequence into a Single String

This is the core operation. The sequence is currently split across many lines — we need to join them back into one continuous string.

```r
sequence <- paste(seq_lines, collapse = "")
```

### Why `collapse = ""`?

`paste()` normally joins strings with a space separator. Setting `collapse = ""` joins all elements of the vector with no separator between them — exactly what we want for a DNA sequence.

### Verify it looks right

```r
nchar(sequence)           # Total number of base pairs
substr(sequence, 1, 60)   # Inspect the first 60 bases
```

You now have a clean, continuous DNA sequence as a plain character string, with no header, no newlines, no extra whitespace — just bases.

---

## Step 6 — Explore and Analyse the Sequence

Now that you have a clean string, you can start asking biological questions about it.

### Count individual bases

```r
bases <- strsplit(sequence, "")[[1]]
table(bases)
```

`strsplit(sequence, "")` splits the sequence into individual characters (one per base). `[[1]]` again extracts the vector from the list. `table()` then counts occurrences of each unique character — you should see counts for A, T, G, and C.

### Calculate GC content

GC content (the proportion of G and C bases) is a fundamental sequence property used in primer design, species identification, and more.

```r
gc <- sum(bases %in% c("G", "C")) / length(bases) * 100
cat("GC content:", round(gc, 2), "%\n")
```

`%in%` checks each base against the vector `c("G", "C")` and returns TRUE/FALSE. `sum()` counts the TRUEs. Dividing by `length(bases)` and multiplying by 100 gives a percentage.

### Count occurrences of a motif

For example, how many ATG codons (potential start codons) are in the sequence?

```r
matches <- gregexpr("ATG", sequence)[[1]]
length(matches)
```

`gregexpr()` finds all positions where a pattern occurs. `[[1]]` extracts from the list. `length()` tells you how many matches were found.

### Extract a subsequence

Pull out a region of interest — for example, the first 300 bases:

```r
substr(sequence, 1, 300)
```

---

## The Full Pipeline Together

Here is the entire manual workflow in one place:

```r
library(rentrez)

# 1. Fetch
raw <- entrez_fetch(db = "nucleotide", id = "NM_007294", rettype = "fasta", retmode = "text")

# 2. Inspect
class(raw)
substr(raw, 1, 200)

# 3. Split into lines
lines <- strsplit(raw, "\n")[[1]]
head(lines, 5)

# 4. Separate header and sequence
header   <- lines[startsWith(lines, ">")]
seq_lines <- lines[!startsWith(lines, ">") & nchar(lines) > 0]

# 5. Collapse sequence into one string
sequence <- paste(seq_lines, collapse = "")
nchar(sequence)

# 6. Analyse
bases <- strsplit(sequence, "")[[1]]
table(bases)
gc <- sum(bases %in% c("G", "C")) / length(bases) * 100
cat("GC content:", round(gc, 2), "%\n")
```

---

## What is FASTA Format, Really?

FASTA is deliberately simple — it is just a **text convention**, not a binary or proprietary format:

```
>IDENTIFIER optional free-text description
SEQUENCELINE1 (usually 60 or 70 characters wide)
SEQUENCELINE2
...
```

The sequence is split across multiple lines **only for human readability**. Biologically it is one continuous string. This means that **collapsing the lines** with `paste(..., collapse = "")` is the fundamental operation that every FASTA parser — including `Biostrings::readDNAStringSet()` — is doing internally.

---

## When to Go Beyond Manual Parsing

Manual parsing is great for understanding, but for real analysis you'll want dedicated packages. Here is a quick guide:

| Goal | Package | Object type |
|---|---|---|
| Sequence manipulation, translation, reverse complement | `Biostrings` | `DNAStringSet` |
| Phylogenetics, evolutionary analysis | `ape` / `phangorn` | `DNAbin` |
| Sequence alignment | `DECIPHER` or `msa` | Works with `DNAStringSet` |
| Gene annotations, features, CDS | `genbankr` | `GBAccession` |
| Just need a plain string | Manual (this tutorial) | `character` |

### Using Biostrings with your fetched data

If you want the power of `Biostrings` without losing the understanding you've built here, you can combine both approaches — fetch with `rentrez`, write to a temp file, then read with `Biostrings`:

```r
library(rentrez)
library(Biostrings)

raw <- entrez_fetch(db = "nucleotide", id = "NM_007294", rettype = "fasta", retmode = "text")

tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)
dna <- readDNAStringSet(tmp)

# Now you have full Biostrings capabilities
width(dna)                  # Sequence length in bp
reverseComplement(dna)      # Reverse complement
translate(dna)              # Translate to amino acids
letterFrequency(dna, "GC")  # GC content
```

The temp file step is the only extra overhead — and now you understand exactly what `readDNAStringSet()` is parsing, because you've seen the raw text yourself.

---

## Summary of Key Concepts

| Concept | Why it matters |
|---|---|
| `entrez_fetch()` returns `character(1)` | It's just text — don't expect a special object |
| FASTA header starts with `>` | This is the universal convention for identifying metadata vs sequence |
| Sequence is split across lines | Cosmetic only — collapse with `paste(..., collapse = "")` |
| `strsplit()` returns a list | Always use `[[1]]` to extract the character vector |
| Filter blank lines | Trailing empty lines silently corrupt your sequence if left in |
| Manual parsing builds intuition | Packages are doing the same thing, just faster and more robustly |
