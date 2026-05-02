# Lesson 5: Tree Visualisation and Interpretation with ggtree

## Overview

In Lessons 2 and 4 we used `ape`'s `plot.phylo()` to display our trees — functional, but basic. In this lesson we move to `ggtree`, the gold standard for phylogenetic tree visualisation in R. Built on top of `ggplot2`, `ggtree` gives you precise control over every aspect of your tree's appearance and allows you to annotate it with external data like geography, morphological traits, dates, or statistical values.

By the end of this lesson you will be able to:

- Understand the anatomy of a phylogenetic tree plot
- Use `ggtree` to create clean, professional tree figures
- Display bootstrap support values clearly
- Change tree layouts (rectangular, circular, fan, cladogram)
- Annotate tips with metadata (colour, shape, labels)
- Highlight clades with shading and labels
- Combine a tree with a data table or heatmap
- Save publication-quality figures

---

## Prerequisites

```r
# Install ggtree and related packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
BiocManager::install("treeio")

# Also needed
install.packages("ggplot2")
install.packages("dplyr")
```

```r
library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(dplyr)
```

---

## Part 1 — The Anatomy of a Phylogenetic Tree Plot

Before writing any code it helps to be precise about the visual elements of a tree:

```
                    ┌── Human           ← tip (leaf)
               ┌────┤
               │    └── Chimpanzee
          ─────┤                        ← internal node (common ancestor)
               │    ┌── Gorilla
               └────┤
                    └── Orangutan
     ──────────── Mouse                 ← outgroup tip

     |──────|                           ← scale bar (substitutions per site)
     branch                             ← branch (line connecting nodes)
```

- **Tips** — the sequences or taxa you sampled
- **Internal nodes** — hypothetical common ancestors
- **Root** — the most ancient common ancestor in the tree (the base)
- **Branches** — connections between nodes; their lengths represent evolutionary change
- **Clade** — a node and all of its descendants
- **Bootstrap values** — support values placed at internal nodes

---

## Part 2 — Your First ggtree Plot

`ggtree` works like `ggplot2` — you build a plot by adding layers with `+`. The base layer is always `ggtree(tree)`.

```r
# Use the ML tree from Lesson 4
# ml_tree_rooted is a rooted phylo object

# Simplest possible plot
ggtree(ml_tree_rooted)
```

This gives you a tree with no labels — just the topology. Now add tip labels:

```r
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4)
```

### Adding a scale bar

```r
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  geom_treescale(x = 0, y = 0.5, fontsize = 3)
```

### The ggtree philosophy

`ggtree` is designed so that every visual element is a separate layer — tips, nodes, branches, labels, highlights — and each can be controlled independently. This is the same philosophy as `ggplot2`, so if you are already comfortable with ggplot, ggtree will feel natural very quickly.

---

## Part 3 — Tree Layouts

One of `ggtree`'s strengths is that switching between layouts is a single argument change:

```r
# Rectangular (default) — most common for publications
ggtree(ml_tree_rooted, layout = "rectangular") +
  geom_tiplab(size = 4)

# Circular — good for many taxa
ggtree(ml_tree_rooted, layout = "circular") +
  geom_tiplab(size = 3)

# Fan — similar to circular, tips fan outward
ggtree(ml_tree_rooted, layout = "fan", open.angle = 30) +
  geom_tiplab(size = 3)

# Cladogram — ignores branch lengths, shows topology only
ggtree(ml_tree_rooted, branch.length = "none") +
  geom_tiplab(size = 4)

# Slanted (triangular)
ggtree(ml_tree_rooted, layout = "slanted") +
  geom_tiplab(size = 4)
```

### When to use each layout

| Layout | Best for |
|---|---|
| Rectangular | Standard publications, most datasets |
| Circular / fan | Large trees with many taxa (30+) |
| Cladogram | Showing topology when branch lengths are not the focus |
| Slanted | Personal preference — some find it cleaner than rectangular |

---

## Part 4 — Displaying Bootstrap Support

Bootstrap values are stored at internal nodes. There are several ways to display them clearly.

### As node labels

```r
# First attach bootstrap values to the tree
# bs_values is a numeric vector from bootstrap.pml() in Lesson 4
ml_tree_bs <- ml_tree_rooted
ml_tree_bs$node.label <- bs   # attach bootstrap values

# Plot with node labels
ggtree(ml_tree_bs) +
  geom_tiplab(size = 4) +
  geom_nodelab(
    aes(label = label),
    size  = 3,
    nudge_x = -0.002,   # adjust position
    nudge_y =  0.2
  )
```

### As coloured point symbols at nodes

A cleaner approach for publications — use point size or colour to indicate support level:

```r
ggtree(ml_tree_bs) +
  geom_tiplab(size = 4) +
  geom_nodepoint(
    aes(subset = !isTip, color = as.numeric(label)),
    size = 4
  ) +
  scale_color_gradient(
    low  = "red",
    high = "steelblue",
    name = "Bootstrap\nsupport",
    limits = c(0, 100)
  ) +
  geom_treescale()
```

### As threshold categories

A common publication style: show a symbol only for well-supported nodes:

```r
# Add a column indicating support category
ggtree(ml_tree_bs) +
  geom_tiplab(size = 4) +
  geom_nodepoint(
    aes(subset = (!isTip & as.numeric(label) >= 70)),
    shape = 21,
    size  = 3,
    fill  = "steelblue",
    color = "white"
  ) +
  geom_nodepoint(
    aes(subset = (!isTip & as.numeric(label) < 70)),
    shape = 21,
    size  = 3,
    fill  = "tomato",
    color = "white"
  )
```

---

## Part 5 — Annotating Tips with Metadata

The real power of `ggtree` shows when you attach external data to your tree. Suppose you have metadata about each species — habitat, body mass, conservation status, or any other trait.

### Create a metadata table

```r
# Metadata for our five species
metadata <- data.frame(
  label  = c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse"),
  group  = c("Hominini", "Hominini", "Homininae", "Ponginae", "Outgroup"),
  mass_kg = c(70, 45, 140, 75, 0.02),
  status = c("LC", "EN", "EN", "CR", "LC"),
  stringsAsFactors = FALSE
)
```

### Attach metadata and plot

```r
# %<+% is ggtree's operator for attaching metadata to a tree
p <- ggtree(ml_tree_rooted) %<+% metadata

# Now you can use metadata columns in aes()
p +
  geom_tiplab(size = 4) +
  geom_tippoint(
    aes(color = group, size = log(mass_kg + 1)),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Hominini"  = "#E64B35",
      "Homininae" = "#4DBBD5",
      "Ponginae"  = "#00A087",
      "Outgroup"  = "#3C5488"
    ),
    name = "Clade"
  ) +
  scale_size_continuous(name = "log(Body mass)") +
  theme(legend.position = "right")
```

### Colour tip labels by group

```r
p +
  geom_tiplab(aes(color = group), size = 4) +
  scale_color_manual(
    values = c(
      "Hominini"  = "#E64B35",
      "Homininae" = "#4DBBD5",
      "Ponginae"  = "#00A087",
      "Outgroup"  = "#3C5488"
    )
  ) +
  theme(legend.position = "none")
```

---

## Part 6 — Highlighting Clades

Highlighting a clade with a background rectangle is a very common figure element in phylogenetics papers.

### Method 1 — geom_hilight

```r
# First you need the node number of the clade's ancestor
# Use MRCA() to find the most recent common ancestor of a set of tips
hominini_node <- MRCA(ml_tree_rooted, "Human", "Chimpanzee")
ape_node      <- MRCA(ml_tree_rooted, "Human", "Orangutan")

ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  geom_hilight(
    node  = ape_node,
    fill  = "#AED6F1",
    alpha = 0.4
  ) +
  geom_hilight(
    node  = hominini_node,
    fill  = "#A9DFBF",
    alpha = 0.5
  ) +
  geom_treescale()
```

### Method 2 — geom_cladelabel

Add a label alongside a highlighted clade:

```r
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  geom_cladelabel(
    node   = ape_node,
    label  = "Great Apes",
    color  = "#2471A3",
    offset = 0.01,
    fontsize = 3.5
  ) +
  geom_cladelabel(
    node   = hominini_node,
    label  = "Hominini",
    color  = "#1E8449",
    offset = 0.005,
    fontsize = 3.5
  ) +
  geom_treescale()
```

---

## Part 7 — Combining a Tree with a Data Panel

`ggtree` can display a heatmap, bar chart, or data table alongside the tree — a **faceted tree + data** figure. This is extremely useful for showing trait data, sequence features, or gene presence/absence alongside the tree.

### Tree + bar chart

```r
# Create a simple data frame with a numeric trait
trait_data <- data.frame(
  label    = c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse"),
  gc_content = c(44.4, 43.6, 43.7, 42.8, 36.2)   # approximate cytb GC%
)

p <- ggtree(ml_tree_rooted) %<+% trait_data +
  geom_tiplab(size = 4, offset = 0.001)

# facet_plot adds a panel to the right of the tree
facet_plot(
  p,
  panel = "GC content (%)",
  data  = trait_data,
  geom  = geom_bar,
  mapping = aes(x = gc_content, fill = label),
  stat  = "identity",
  orientation = "y"
) +
  theme_tree2() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")
```

### Tree + heatmap

For multiple traits at once, a heatmap alongside the tree is ideal:

```r
# Multiple trait columns
heatmap_data <- data.frame(
  label       = c("Human", "Chimpanzee", "Gorilla", "Orangutan", "Mouse"),
  GC_content  = c(44.4, 43.6, 43.7, 42.8, 36.2),
  Length_kb   = c(1.14, 1.14, 1.14, 1.14, 1.14),
  AT_skew     = c(0.12, 0.14, 0.13, 0.15, 0.22)
)

# Reshape to long format for ggplot
library(tidyr)
heatmap_long <- pivot_longer(
  heatmap_data,
  cols      = -label,
  names_to  = "trait",
  values_to = "value"
)

p2 <- ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4, offset = 0.001)

gheatmap(
  p2,
  heatmap_data[, -1] |> `rownames<-`(heatmap_data$label),
  offset     = 0.02,
  width      = 0.3,
  colnames   = TRUE,
  font.size  = 3
) +
  scale_fill_viridis_c(name = "Value") +
  theme(legend.position = "right")
```

---

## Part 8 — Comparing Two Trees

When you have trees built from different genes or methods, comparing them visually is valuable. `ggtree` supports this with `cophyloplot` via `ape`, or using `ggtree`'s own layout tools.

### Side-by-side with ape

```r
# Compare NJ vs ML tree
par(mfrow = c(1, 2))
plot(root(nj_tree, "Mouse", resolve.root = TRUE),
     main = "Neighbour-Joining", cex = 1)
plot(ml_tree_rooted, main = "Maximum Likelihood", cex = 1)
par(mfrow = c(1, 1))
```

### Tanglegram with ape

A tanglegram connects matching tips between two trees with lines, making it easy to spot topological differences:

```r
# Root both trees
nj_rooted <- root(nj_tree, "Mouse", resolve.root = TRUE)

# Tanglegram
cophyloplot(
  nj_rooted,
  ml_tree_rooted,
  assoc     = cbind(nj_rooted$tip.label, nj_rooted$tip.label),
  space     = 40,
  length.line = 1,
  gap       = 2,
  col       = "steelblue"
)
```

---

## Part 9 — Styling and Theming

`ggtree` provides several built-in themes and full compatibility with `ggplot2` themes:

```r
# Clean white background — good for publications
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  theme_tree2()   # adds axis with branch length scale

# Completely minimal
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  theme_tree()

# Use ggplot2 themes
ggtree(ml_tree_rooted) +
  geom_tiplab(size = 4) +
  theme_minimal()

# Custom fonts and colours
ggtree(ml_tree_rooted, color = "#2C3E50", size = 0.8) +
  geom_tiplab(size = 4, color = "#2C3E50", fontface = "italic") +
  geom_nodepoint(aes(subset = !isTip), size = 2, color = "#E74C3C") +
  theme_tree2() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
```

---

## Part 10 — Saving Publication-Quality Figures

```r
# Build your final figure
final_tree <- ggtree(ml_tree_bs, layout = "rectangular") %<+% metadata +
  geom_tiplab(aes(color = group), size = 4, fontface = "italic") +
  geom_nodepoint(
    aes(subset = (!isTip & as.numeric(label) >= 70)),
    size = 3, shape = 21, fill = "steelblue", color = "white"
  ) +
  geom_hilight(node = ape_node,     fill = "#AED6F1", alpha = 0.3) +
  geom_hilight(node = hominini_node, fill = "#A9DFBF", alpha = 0.3) +
  geom_cladelabel(node = ape_node,      label = "Great Apes",
                  color = "#2471A3", offset = 0.012, fontsize = 3) +
  geom_cladelabel(node = hominini_node, label = "Hominini",
                  color = "#1E8449", offset = 0.006, fontsize = 3) +
  geom_treescale(x = 0, y = 0.5, fontsize = 3) +
  scale_color_manual(
    values = c("Hominini" = "#E64B35", "Homininae" = "#4DBBD5",
               "Ponginae" = "#00A087", "Outgroup"  = "#3C5488")
  ) +
  theme_tree2() +
  theme(legend.position = "none") +
  ggtitle("Cytochrome b Phylogeny — Great Apes",
          subtitle = "Maximum likelihood, GTR+G, bootstrap ≥70 shown")

# Save as PDF (vector — infinitely scalable, ideal for journals)
ggsave("cytb_phylogeny.pdf",
       plot   = final_tree,
       width  = 7,
       height = 5,
       units  = "in")

# Save as PNG (raster — good for presentations, posters)
ggsave("cytb_phylogeny.png",
       plot   = final_tree,
       width  = 7,
       height = 5,
       units  = "in",
       dpi    = 300)   # 300 dpi is the standard for print publication
```

### File format guidance

| Format | Use case |
|---|---|
| PDF | Journal submission, infinitely scalable vector format |
| SVG | Web figures, can be edited in Inkscape or Illustrator |
| PNG (300 dpi) | Presentations, posters, supplementary figures |
| TIFF (300+ dpi) | Some journals require TIFF specifically |

---

## Part 11 — Reading and Plotting Trees from Files

In real workflows you will often have tree files from external tools (IQ-TREE, RAxML, MrBayes). `treeio` handles reading these formats:

```r
library(treeio)

# Read a Newick tree file
tree_nwk <- read.newick("my_tree.nwk")

# Read a NEXUS tree file (from MrBayes or PAUP)
tree_nex <- read.nexus("my_tree.nex")

# Read an IQ-TREE output file (with bootstrap values)
tree_iqtree <- read.iqtree("my_tree.treefile")

# Read a BEAST output file (with divergence times)
tree_beast <- read.beast("my_tree.beast")

# All of these can be passed directly to ggtree
ggtree(tree_nwk) + geom_tiplab()
ggtree(tree_beast) +
  geom_tiplab() +
  geom_range("height_0.95_HPD", color = "red", size = 2, alpha = 0.5)
```

---

## Part 12 — A Complete Annotated Figure

Here is the complete workflow from Lesson 4's ML tree to a fully annotated, publication-ready figure:

```r
library(ggtree)
library(ggplot2)
library(ape)

# Assuming ml_tree_rooted and bs (bootstrap values) from Lesson 4
# and metadata data frame defined earlier in this lesson

# Attach bootstrap to tree
ml_tree_bs              <- ml_tree_rooted
ml_tree_bs$node.label   <- bs

# Find clade nodes
ape_node      <- MRCA(ml_tree_bs, "Human", "Orangutan")
hominini_node <- MRCA(ml_tree_bs, "Human", "Chimpanzee")

# Build the figure
p_final <- ggtree(ml_tree_bs) %<+% metadata +

  # Clade highlights (add before other layers so they sit behind)
  geom_hilight(node = ape_node,      fill = "#D6EAF8", alpha = 0.6) +
  geom_hilight(node = hominini_node, fill = "#D5F5E3", alpha = 0.6) +

  # Clade labels
  geom_cladelabel(node = ape_node,
                  label = "Great Apes", barsize = 1.2,
                  color = "#1A5276", offset = 0.014, fontsize = 3.2) +
  geom_cladelabel(node = hominini_node,
                  label = "Hominini", barsize = 1.2,
                  color = "#1D6A39", offset = 0.007, fontsize = 3.2) +

  # Tip labels coloured by clade
  geom_tiplab(aes(color = group), size = 4, fontface = "italic",
              offset = 0.001) +

  # Bootstrap points at nodes
  geom_nodepoint(
    aes(subset = (!isTip & as.numeric(label) >= 70)),
    size = 3.5, shape = 21, fill = "#2980B9", color = "white", stroke = 0.5
  ) +

  # Tip points coloured by conservation status
  geom_tippoint(aes(color = group), size = 3) +

  # Scale bar
  geom_treescale(x = 0, y = 4.8, width = 0.01,
                 offset = 0.3, fontsize = 3) +

  # Colour scale
  scale_color_manual(
    values = c("Hominini"  = "#C0392B",
               "Homininae" = "#2471A3",
               "Ponginae"  = "#148F77",
               "Outgroup"  = "#566573"),
    guide = "none"
  ) +

  # Theme and labels
  theme_tree2() +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.margin   = margin(10, 60, 10, 10)
  ) +
  ggtitle(
    "Cytochrome b Phylogeny — Hominidae",
    subtitle = "Maximum likelihood · GTR+G · Filled circles = bootstrap ≥ 70"
  )

p_final

# Save
ggsave("cytb_final_figure.pdf", p_final, width = 8, height = 5)
ggsave("cytb_final_figure.png", p_final, width = 8, height = 5, dpi = 300)
```

---

## Summary

| Concept | Why it matters |
|---|---|
| `ggtree` is built on `ggplot2` | Familiar syntax, composable layers, full theme control |
| `%<+%` operator | Attaches metadata to the tree for use in `aes()` |
| Layout options | Circular for many taxa, rectangular for publications |
| Bootstrap as points | Cleaner than text labels, easy to read at a glance |
| `geom_hilight` and `geom_cladelabel` | Visually group and label clades — essential for complex trees |
| `facet_plot` and `gheatmap` | Combine tree with trait data in one figure |
| `treeio` | Reads trees from external tools (IQ-TREE, BEAST, MrBayes) |
| Save as PDF or 300 dpi PNG | Vector for journals, raster for presentations |

---

## What's Next

With four lessons completed you now have a full phylogenetics toolkit:

- Fetch sequences (Lesson 1)
- Parse and manipulate with Biostrings (Lesson 2)
- Align with DECIPHER and MSA (Lesson 3)
- Model selection and ML tree building (Lesson 4)
- Publication-quality visualisation (Lesson 5)

**Lesson 6** extends this foundation into protein sequences — translating DNA to amino acids, working with `AAStringSet` objects, protein alignment, and why protein sequences are sometimes better than DNA for deep evolutionary comparisons.
