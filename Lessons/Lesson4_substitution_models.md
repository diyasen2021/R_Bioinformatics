# Lesson 4: Substitution Models and Maximum Likelihood Phylogenetics

## Overview

In Lesson 2 we built a neighbour-joining tree using the TN93 substitution model without really explaining what that model is or why we chose it. In Lesson 3 we built a deeper understanding of alignment. Now we address the question that serious phylogenetics demands: **how do we model sequence evolution, and how do we find the tree that best explains our data?**

This lesson covers:
- What substitution models are and why they are necessary
- The hierarchy of models from simple to complex
- How to select the best model for your data
- Maximum likelihood — what it is and why it produces better trees than neighbour-joining
- Bootstrap support, properly understood
- A complete ML phylogenetics workflow in R using `phangorn`

---

## Part 1 — Why We Need Substitution Models

### The multiple hits problem

Imagine two sequences that diverged from a common ancestor 50 million years ago. At position 100, the ancestral base was A. Over time, in one lineage it mutated to G. Much later, that G mutated to T. When we observe the two sequences today, position 100 shows A in one and T in the other — we see **one difference**, but **two mutations** actually occurred.

This is called **multiple hits** — the same position has been hit by more than one mutation. The more time has passed, the more likely multiple hits become. If we simply count observed differences and call that the evolutionary distance, we will systematically **underestimate** the true amount of change that has occurred.

Substitution models correct for multiple hits by modelling the evolutionary process mathematically, allowing us to estimate the **true** number of substitutions from the **observed** number of differences.

### Transitions vs transversions

Not all substitutions are equally likely. In DNA, there are two classes:

- **Transitions** — substitutions between chemically similar bases: purines (A ↔ G) or pyrimidines (C ↔ T). These are chemically easier and occur roughly twice as often as transversions.
- **Transversions** — substitutions between a purine and a pyrimidine (A ↔ C, A ↔ T, G ↔ C, G ↔ T). Chemically harder, occur less frequently.

The ratio of transitions to transversions in a dataset (the **ti/tv ratio** or **κ**) varies across genes and organisms. Simple models ignore this distinction; more realistic models account for it.

### Base frequencies

In real genomes, bases are not equally frequent. Human mitochondrial DNA is AT-rich. Some bacterial genomes are extremely GC-rich. A model that assumes equal base frequencies will be wrong for these sequences. More realistic models estimate base frequencies from the data.

---

## Part 2 — The Model Hierarchy

Substitution models form a nested hierarchy — each model is a special case of the next, more complex one. Understanding this hierarchy helps you choose the right model and understand what each one assumes.

### JC69 — Jukes-Cantor (1969)

The simplest possible model. Assumes:
- All four bases are equally frequent (25% each)
- All substitutions occur at the same rate

Just one parameter: the overall substitution rate. Rarely appropriate for real data but useful as a baseline and for teaching.

```
All substitution rates equal: A→C = A→G = A→T = C→G = C→T = G→T
Base frequencies: πA = πC = πG = πT = 0.25
```

### K80 — Kimura 2-parameter (1980)

Adds one distinction over JC69:
- Transitions and transversions occur at different rates
- Still assumes equal base frequencies

Two parameters: α (transition rate) and β (transversion rate). Better than JC69 for most biological data.

### HKY85 — Hasegawa, Kishino, Yano (1985)

Adds unequal base frequencies to K80:
- Different transition and transversion rates
- Base frequencies estimated from the data

Four parameters: κ (ti/tv ratio) plus three base frequencies (the fourth is determined by the others summing to 1). A good general-purpose model.

### TN93 — Tamura-Nei (1993)

What we used in Lesson 2. Extends HKY85 by allowing:
- **Two different transition rates** — one for purines (A↔G) and one for pyrimidines (C↔T)
- Unequal base frequencies

This is often a better fit for mitochondrial data, where the two transition types occur at noticeably different rates.

### GTR — General Time Reversible

The most parameter-rich standard model. Allows:
- Six different substitution rates (one for each pair of bases)
- Unequal base frequencies

Twelve parameters in total (six rates + four base frequencies, minus constraints). This is the most flexible model and often fits data best, but requires sufficient data to estimate all parameters reliably. GTR is the most commonly used model in modern phylogenetics.

### Rate variation across sites — Gamma (Γ) and Invariant sites (I)

An important addition to any of the above models: in real sequences, not all positions evolve at the same rate. Some positions are highly conserved (e.g. active sites of enzymes), others evolve rapidly.

Two common extensions account for this:

- **+G (Gamma distribution)** — models rate variation across sites using a Gamma distribution. Controlled by the shape parameter α. Small α = high variation; large α = rates are more uniform. Almost always improves model fit significantly.
- **+I (Invariant sites)** — allows a proportion of sites to be completely invariable (evolutionarily frozen). Can be combined with +G.

In practice you will often see models written as **GTR+G**, **TN93+G+I**, **HKY+G**, etc. The +G extension is almost always worth including.

### The model hierarchy at a glance

```
JC69  (1 param)
  ↓ + ti/tv distinction
K80   (2 params)
  ↓ + unequal base frequencies
HKY85 (4 params)
  ↓ + separate purine/pyrimidine transition rates
TN93  (6 params)
  ↓ + all rates free
GTR   (8 params)
  ↓ + rate variation across sites
GTR+G (9 params)
  ↓ + invariant sites
GTR+G+I (10 params)
```

---

## Part 3 — Model Selection

### The problem of overfitting

A more complex model will always fit the data at least as well as a simpler one — it has more parameters to play with. But a model with too many parameters may be **overfitting** — fitting the noise in your specific dataset rather than the true underlying signal. An overfitted model performs poorly on new data and can produce unstable or misleading trees.

We need a principled way to choose a model that fits well without being unnecessarily complex.

### Information criteria — AIC and BIC

The two most commonly used model selection criteria are:

**AIC (Akaike Information Criterion)**
```
AIC = 2k - 2 ln(L)
```
Where k is the number of parameters and L is the likelihood. Lower AIC is better. AIC penalises complexity but relatively lightly — it tends to favour slightly more complex models.

**BIC (Bayesian Information Criterion)**
```
BIC = k ln(n) - 2 ln(L)
```
Where n is the number of sites in the alignment. BIC penalises complexity more heavily than AIC — it tends to favour simpler models. For phylogenetics BIC is often preferred.

In practice, the model selected by AIC and BIC is often the same or very similar, and both are vastly better than just picking a model arbitrarily.

### Model selection in R with phangorn

```r
library(phangorn)

# Read the aligned sequences (from Lesson 3 workflow)
# phangorn uses its own object type: phyDat
cytb_phyDat <- as.phyDat(cytb_bin)   # cytb_bin is the DNAbin from Lesson 3

# modelTest tests all standard models and ranks them by AIC and BIC
# This takes a few minutes — it fits every model to your data
mt <- modelTest(cytb_phyDat, model = "all", multicore = FALSE)

# View results sorted by BIC
mt[order(mt$BIC), ]

# The top model by BIC
best_model <- mt$Model[which.min(mt$BIC)]
cat("Best model (BIC):", best_model, "\n")
```

### What modelTest does

`modelTest()` fits every model in its repertoire to your alignment using a starting tree (usually NJ), computes the log-likelihood for each, and returns a table with AIC, BIC, and other statistics. You then select the model with the lowest AIC or BIC and use it for your maximum likelihood tree search.

---

## Part 4 — Maximum Likelihood Phylogenetics

### What is maximum likelihood?

Maximum likelihood (ML) is a statistical framework for estimation. In phylogenetics, the question is:

> Given my alignment and my substitution model, what tree (topology + branch lengths) makes the observed data most probable?

Formally, we want to find the tree T that maximises P(data | T, model) — the probability of observing our alignment given the tree and the model. This probability is called the **likelihood**.

### Why is ML better than neighbour-joining?

Neighbour-joining is a **distance method** — it compresses all the information in your alignment into a single number per sequence pair (the distance), then builds a tree from those distances. Information is lost in this compression.

Maximum likelihood uses the **full alignment** — every column, every base, in every sequence — directly in the calculation. It evaluates how well each possible tree explains every position in the alignment simultaneously, under an explicit model of how sequences evolve.

ML is:
- Statistically consistent — as you add more data, it converges on the true tree
- More accurate than NJ, especially for divergent sequences
- Able to use a biologically realistic substitution model throughout
- The standard for modern phylogenetic publications

The trade-off is computation. Finding the ML tree requires searching through an enormous space of possible tree topologies (for 10 sequences there are over 2 million possible unrooted topologies; for 20 sequences, over 2 × 10²⁰). ML tree searches use heuristics to navigate this space efficiently.

### The ML tree search in phangorn

`phangorn` implements ML phylogenetics entirely within R. The typical workflow is:

1. Start with a reasonable starting tree (usually NJ)
2. Use the best-fit model from `modelTest()`
3. Optimise the tree using nearest-neighbour interchange (NNI) or subtree pruning and regrafting (SPR)
4. Assess support with bootstrap

```r
library(phangorn)

# Step 1 — Starting tree (NJ)
dist_matrix  <- dist.dna(cytb_bin, model = "TN93")
nj_tree      <- nj(dist_matrix)

# Step 2 — Set up the model
# Use the best model identified by modelTest
# Here we use GTR+G as an example — substitute your best_model
pml_fit <- pml(nj_tree, data = cytb_phyDat)

# Step 3 — Optimise under GTR+G
ml_fit <- optim.pml(
  pml_fit,
  model       = "GTR",
  optGamma    = TRUE,    # optimise the Gamma shape parameter
  optInv      = FALSE,   # set TRUE to include invariant sites (+I)
  rearrangement = "stochastic",  # tree topology search strategy
  control     = pml.control(trace = 0)  # suppress verbose output
)

ml_fit                     # summary of the fitted model
logLik(ml_fit)             # log-likelihood of the ML tree
ml_tree <- ml_fit$tree     # extract the tree
```

### Understanding the optimisation output

```r
# Model parameter estimates
ml_fit$Q      # rate matrix (GTR rates)
ml_fit$bf     # estimated base frequencies
ml_fit$shape  # Gamma shape parameter (α)

# Compare log-likelihoods of starting vs optimised tree
logLik(pml_fit)    # NJ starting tree
logLik(ml_fit)     # ML optimised tree — should be higher (less negative)
```

---

## Part 5 — Bootstrap Support

### What bootstrap means in phylogenetics

Bootstrap support quantifies **how well-supported each node in the tree is by the data**. The procedure is:

1. Randomly resample columns from your alignment **with replacement** to create a new alignment of the same size (a "bootstrap replicate")
2. Build a tree from this resampled alignment using the same method and model
3. Repeat 100 to 1000 times
4. For each node in your original tree, count what percentage of bootstrap trees contained that same grouping — this is the **bootstrap support value**

A bootstrap value of 95 means that 95% of resampled datasets produced a tree with that node. Values above 70 are generally considered supported; above 95 is strong support. Values below 50 suggest that node is poorly supported by the data.

### Bootstrap in phangorn

```r
# Bootstrap the ML tree — 100 replicates (use 1000 for publication)
bs <- bootstrap.pml(
  ml_fit,
  bs          = 100,
  optNni      = TRUE,     # optimise topology for each replicate
  control     = pml.control(trace = 0)
)

# Plot with bootstrap values
ml_tree_rooted <- root(ml_tree, outgroup = "Mouse", resolve.root = TRUE)

plotBS(
  ml_tree_rooted,
  bs,
  type = "phylogram",
  bs.col = "blue",
  cex  = 1.0,
  main = "Maximum Likelihood Tree — Cytochrome b\nGTR+G, bootstrap support shown"
)
add.scalebar()
```

### Bootstrap values vs Bayesian posterior probabilities

You may encounter **posterior probability** support values in Bayesian phylogenetics (e.g. from MrBayes or BEAST). These are not the same as bootstrap values:

- Bootstrap values are **frequentist** — they measure reproducibility of the node across resampled datasets
- Posterior probabilities are **Bayesian** — they represent the probability that the node is correct given the data and prior

Posterior probabilities tend to be higher than bootstrap values for the same node. A posterior probability of 0.95 and a bootstrap of 70 might support the same node with similar actual confidence. Neither is universally better — both are valid measures of support under their respective frameworks.

---

## Part 6 — Comparing NJ and ML Trees

It is instructive to compare your NJ and ML trees to see how much the method matters for your dataset:

```r
# Root both trees on the outgroup
nj_tree_rooted <- root(nj_tree, outgroup = "Mouse", resolve.root = TRUE)
ml_tree_rooted <- root(ml_tree, outgroup = "Mouse", resolve.root = TRUE)

# Plot side by side
par(mfrow = c(1, 2))

plot(nj_tree_rooted,
     main = "Neighbour-Joining (TN93)",
     cex  = 1.0)
add.scalebar()

plot(ml_tree_rooted,
     main = "Maximum Likelihood (GTR+G)",
     cex  = 1.0)
add.scalebar()

par(mfrow = c(1, 1))
```

For closely related sequences with a good alignment, NJ and ML often agree. For more divergent sequences, or datasets with rate variation and compositional bias, ML with an appropriate model can give substantially different (and more accurate) results.

---

## Part 7 — The Complete ML Workflow

```r
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(ape)
library(phangorn)

# --- 1. Fetch and parse ---
accessions <- c("NC_012920", "NC_001643", "NC_011120", "NC_001646", "NC_005089")
raw        <- entrez_fetch(db = "nucleotide", id = accessions,
                           rettype = "fasta", retmode = "text")
tmp <- tempfile(fileext = ".fasta")
writeLines(raw, tmp)
dna <- readDNAStringSet(tmp)
names(dna) <- c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse")

# --- 2. Extract and align ---
cytb         <- subseq(dna, start = 14747, end = 15887)
cytb_aligned <- AlignSeqs(cytb)

# --- 3. Convert formats ---
cytb_bin    <- as.DNAbin(cytb_aligned)
cytb_phyDat <- as.phyDat(cytb_bin)

# --- 4. Model selection ---
mt          <- modelTest(cytb_phyDat, model = "all",
                         multicore = FALSE)
best_model  <- mt$Model[which.min(mt$BIC)]
cat("Best model:", best_model, "\n")

# --- 5. Starting tree ---
dist_matrix <- dist.dna(cytb_bin, model = "TN93")
nj_tree     <- nj(dist_matrix)

# --- 6. ML optimisation ---
pml_fit <- pml(nj_tree, data = cytb_phyDat)
ml_fit  <- optim.pml(
  pml_fit,
  model         = "GTR",
  optGamma      = TRUE,
  rearrangement = "stochastic",
  control       = pml.control(trace = 0)
)

# --- 7. Bootstrap ---
bs <- bootstrap.pml(ml_fit, bs = 100, optNni = TRUE,
                    control = pml.control(trace = 0))

# --- 8. Plot ---
ml_tree_rooted <- root(ml_fit$tree, outgroup = "Mouse", resolve.root = TRUE)
plotBS(ml_tree_rooted, bs, type = "phylogram", cex = 1.1,
       main = "ML Cytochrome b Phylogeny — GTR+G")
add.scalebar()
```

---

## Summary

| Concept | Why it matters |
|---|---|
| Multiple hits | Observed differences underestimate true evolutionary change — models correct for this |
| Transitions vs transversions | Not all substitutions are equally likely — realistic models account for this |
| Model hierarchy JC69 → GTR | More complex models fit better but risk overfitting — model selection balances this |
| +G (Gamma rate variation) | Sites evolve at different rates — almost always worth including |
| AIC and BIC | Principled ways to choose the best model without overfitting |
| `phangorn::modelTest()` | Tests all standard models and ranks them — always run this |
| Maximum likelihood | Uses the full alignment under an explicit model — more accurate than NJ |
| Bootstrap support | Measures reproducibility of each node — values above 70 are generally supported |
| NJ as a starting point | Fast and useful for exploration — ML is the standard for rigorous work |

---

## What's Next

Lesson 5 covers **tree visualisation** — moving beyond `plot.phylo()` to publication-quality, annotated trees using `ggtree`. With a properly estimated ML tree in hand, you are ready to make figures worthy of a paper.
