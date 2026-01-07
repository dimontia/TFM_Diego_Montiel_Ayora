# Libraries
library(ape)
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyr)
library(tibble)
library(stringr)
library(scales)
library(ggthemes)
library(paletteer)

# Phylo trees
LIM <- "LIM_final_doms_phylogeny.treefile"
LIM_tree <- read.tree(LIM)
CERS <- "CERS_final_doms_phylogeny.treefile"
CERS_tree <- read.tree(CERS)
CUT <- "CUT_final_doms_phylogeny.treefile"
CUT_tree <- read.tree(CUT)
HNF <- "HNF_final_doms_phylogeny.treefile"
HNF_tree <- read.tree(HNF)
POU <- "POU_final_doms_phylogeny.treefile"
POU_tree <- read.tree(POU)
SINE <- "SINE_final_doms_phylogeny.treefile"
SINE_tree <- read.tree(SINE)
SOX <- "SOX_final_doms_phylogeny.treefile"
SOX_tree <- read.tree(SOX)
HOX <- "HOX_final_doms_phylogeny.treefile"
HOX_tree <- read.tree(HOX)
NKL <- "NKL_final_doms_phylogeny.treefile"
NKL_tree <- read.tree(NKL)
ZF <- "ZF_final_doms_phylogeny.treefile"
ZF_tree <- read.tree(ZF)

# Species reference code csv
reference_code <- read.csv("code.csv", header = TRUE, sep = ";")
reference_code[4:30]<-list(NULL)
reference_code <- reference_code[-c(163:1182), ]

reference_code <- rows_insert(reference_code, tibble(Species = "Microstomoeca roanoka", Taxonomy = "Choanoflagellata", Code = "Mroan"))
reference_code <- rows_insert(reference_code, tibble(Code = "Skvev", Species = "Salpingoeca kvevrii", Taxonomy = "Choanoflagellata"))
reference_code <- rows_insert(reference_code, tibble(Code = "Bmino", Species = "Bicosta minor", Taxonomy = "Choanoflagellata"))
reference_code <- rows_insert(reference_code, tibble(Code = "Bmono", Species = "Barroeca monosierra", Taxonomy = "Choanoflagellata"))

# LIM
LIM_leaf <- LIM_tree$tip.label
LIM_dataset <- data.frame(Code=substr(LIM_leaf, 1, 5), Gene="LIM", Header=LIM_leaf)
LIM_dataset <- LIM_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
LIM_dataset <- LIM_dataset[-c(as.numeric(rownames(LIM_dataset[apply(LIM_dataset, 1, function(x) any(is.na(x))), ]))), ]
# CERS
CERS_leaf <- CERS_tree$tip.label
CERS_dataset <- data.frame(Code=substr(CERS_leaf, 1, 5), Gene="CERS", Header=CERS_leaf)
CERS_dataset <- CERS_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
CERS_dataset <- CERS_dataset[-c(as.numeric(rownames(CERS_dataset[apply(CERS_dataset, 1, function(x) any(is.na(x))), ]))), ]
# CUT
CUT_leaf <- CUT_tree$tip.label
CUT_dataset <- data.frame(Code=substr(CUT_leaf, 1, 5), Gene="CUT", Header=CUT_leaf)
CUT_dataset <- CUT_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
CUT_dataset <- CUT_dataset[-c(as.numeric(rownames(CUT_dataset[apply(CUT_dataset, 1, function(x) any(is.na(x))), ]))), ]
# HNF
HNF_leaf <- HNF_tree$tip.label
HNF_dataset <- data.frame(Code=substr(HNF_leaf, 1, 5), Gene="HNF", Header=HNF_leaf)
HNF_dataset <- HNF_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
HNF_dataset <- HNF_dataset[-c(as.numeric(rownames(HNF_dataset[apply(HNF_dataset, 1, function(x) any(is.na(x))), ]))), ]
# POU
POU_leaf <- POU_tree$tip.label
POU_dataset <- data.frame(Code=substr(POU_leaf, 1, 5), Gene="POU", Header=POU_leaf)
POU_dataset <- POU_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
POU_dataset <- POU_dataset[-c(as.numeric(rownames(POU_dataset[apply(POU_dataset, 1, function(x) any(is.na(x))), ]))), ]
# SINE
SINE_leaf <- SINE_tree$tip.label
SINE_dataset <- data.frame(Code=substr(SINE_leaf, 1, 5), Gene="SINE", Header=SINE_leaf)
SINE_dataset <- SINE_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
SINE_dataset <- SINE_dataset[-c(as.numeric(rownames(SINE_dataset[apply(SINE_dataset, 1, function(x) any(is.na(x))), ]))), ]

# Indicate the number of closer nodes it can search
K <- 4

# SOX
queries_SOX <- c("sox")
SOX_leaf_queries <- grep(paste(queries_SOX, collapse = "|"),
                         SOX_tree$tip.label,
                         value = TRUE,
                         ignore.case = TRUE)
unit_tree <- SOX_tree
unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))
D <- dist.nodes(unit_tree)
ntip <- length(unit_tree$tip.label)
dist_to_mrca_edges <- function(tree, tip1, tip2, Dmat) {
  m <- getMRCA(tree, c(tip1, tip2))
  if (is.na(m)) return(list(d1=NA_real_, d2=NA_real_))
  i1 <- match(tip1, tree$tip.label)
  i2 <- match(tip2, tree$tip.label)
  list(d1 = Dmat[i1, m], d2 = Dmat[i2, m])
}
all_tips <- unit_tree$tip.label
pairs <- expand.grid(query = SOX_leaf_queries,
                     tip   = all_tips,      
                     stringsAsFactors = FALSE)
withinK <- logical(nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  q <- pairs$query[i]
  t <- pairs$tip[i]                 
  dd <- dist_to_mrca_edges(unit_tree, q, t, D)
  withinK[i] <- !any(is.na(unlist(dd))) && max(dd$d1, dd$d2) <= K
}
SOX_tips_withinK <- unique(pairs$tip[withinK])
SOX_dataset <- data.frame(Code=substr(SOX_tips_withinK, 1, 5), Gene="SOX", Header=SOX_tips_withinK)
SOX_dataset <- SOX_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
SOX_dataset <- SOX_dataset[-c(as.numeric(rownames(SOX_dataset[apply(SOX_dataset, 1, function(x) any(is.na(x))), ]))), ]

# HOX
queries_HOX <- c("cdx", "evx", "gbx", "gsx", "hoxa", "mnx", "meox", "pdx", "hox", "labial", "lab", "abdominal", "abd", "bicoid", "bcd")
HOX_leaf_queries <- grep(paste(queries_HOX, collapse = "|"),
                         HOX_tree$tip.label,
                         value = TRUE,
                         ignore.case = TRUE)
unit_tree <- HOX_tree
unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))
D <- dist.nodes(unit_tree)
ntip <- length(unit_tree$tip.label)
dist_to_mrca_edges <- function(tree, tip1, tip2, Dmat) {
  m <- getMRCA(tree, c(tip1, tip2))
  if (is.na(m)) return(list(d1=NA_real_, d2=NA_real_))
  i1 <- match(tip1, tree$tip.label)
  i2 <- match(tip2, tree$tip.label)
  list(d1 = Dmat[i1, m], d2 = Dmat[i2, m])
}
all_tips <- unit_tree$tip.label
pairs <- expand.grid(query = HOX_leaf_queries,
                     tip   = all_tips,
                     stringsAsFactors = FALSE)
withinK <- logical(nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  q <- pairs$query[i]
  t <- pairs$tip[i]
  dd <- dist_to_mrca_edges(unit_tree, q, t, D)
  withinK[i] <- !any(is.na(unlist(dd))) && max(dd$d1, dd$d2) <= K
}
HOX_tips_withinK <- unique(pairs$tip[withinK])
HOX_dataset <- data.frame(Code=substr(HOX_tips_withinK, 1, 5), Gene="HOX", Header=HOX_tips_withinK)
HOX_dataset <- HOX_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
HOX_dataset <- HOX_dataset[-c(as.numeric(rownames(HOX_dataset[apply(HOX_dataset, 1, function(x) any(is.na(x))), ]))), ]

# NKL
queries_NKL <- c("BARHL", "BARX", "BSX", "DBX", "DLX", "EMX", "HHEX", "HLX", "LBX", "NANOG", "NKX", "NOTO", "TLX", "VAX", "VENT", "NK", "distal", "empty_spiracles", "scarecrow", "scro")
NKL_leaf_queries <- grep(paste(queries_NKL, collapse = "|"),
                         NKL_tree$tip.label,
                         value = TRUE,
                         ignore.case = TRUE)
unit_tree <- NKL_tree
unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))
D <- dist.nodes(unit_tree)
ntip <- length(unit_tree$tip.label)
dist_to_mrca_edges <- function(tree, tip1, tip2, Dmat) {
  m <- getMRCA(tree, c(tip1, tip2))
  if (is.na(m)) return(list(d1=NA_real_, d2=NA_real_))
  i1 <- match(tip1, tree$tip.label)
  i2 <- match(tip2, tree$tip.label)
  list(d1 = Dmat[i1, m], d2 = Dmat[i2, m])
}
all_tips <- unit_tree$tip.label
pairs <- expand.grid(query = NKL_leaf_queries,
                     tip   = all_tips,
                     stringsAsFactors = FALSE)
withinK <- logical(nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  q <- pairs$query[i]
  t <- pairs$tip[i]
  dd <- dist_to_mrca_edges(unit_tree, q, t, D)
  withinK[i] <- !any(is.na(unlist(dd))) && max(dd$d1, dd$d2) <= K
}
NKL_tips_withinK <- unique(pairs$tip[withinK])
NKL_dataset <- data.frame(Code=substr(NKL_tips_withinK, 1, 5), Gene="NKL", Header=NKL_tips_withinK)
NKL_dataset <- NKL_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
NKL_dataset <- NKL_dataset[-c(as.numeric(rownames(NKL_dataset[apply(NKL_dataset, 1, function(x) any(is.na(x))), ]))), ]

# ZF
queries_ZF <- c("adnp", "tshz", "zeb", "zfhx", "homez", "teashirt", "zinc_finger")
ZF_leaf_queries <- grep(paste(queries_ZF, collapse = "|"),
                        ZF_tree$tip.label,
                        value = TRUE,
                        ignore.case = TRUE)
unit_tree <- ZF_tree
unit_tree$edge.length <- rep(1, nrow(unit_tree$edge))
D <- dist.nodes(unit_tree)
ntip <- length(unit_tree$tip.label)
dist_to_mrca_edges <- function(tree, tip1, tip2, Dmat) {
  m <- getMRCA(tree, c(tip1, tip2))
  if (is.na(m)) return(list(d1=NA_real_, d2=NA_real_))
  i1 <- match(tip1, tree$tip.label)
  i2 <- match(tip2, tree$tip.label)
  list(d1 = Dmat[i1, m], d2 = Dmat[i2, m])
}
all_tips <- unit_tree$tip.label
pairs <- expand.grid(query = ZF_leaf_queries,
                     tip   = all_tips,
                     stringsAsFactors = FALSE)
withinK <- logical(nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  q <- pairs$query[i]
  t <- pairs$tip[i]
  dd <- dist_to_mrca_edges(unit_tree, q, t, D)
  withinK[i] <- !any(is.na(unlist(dd))) && max(dd$d1, dd$d2) <= K
}
ZF_tips_withinK <- unique(pairs$tip[withinK])
ZF_dataset <- data.frame(Code=substr(ZF_tips_withinK, 1, 5), Gene="ZF", Header=ZF_tips_withinK)
ZF_dataset <- ZF_dataset %>%
  full_join(reference_code, by = c("Code" = "Code"));
ZF_dataset <- ZF_dataset[-c(as.numeric(rownames(ZF_dataset[apply(ZF_dataset, 1, function(x) any(is.na(x))), ]))), ]

# Build full dataset
dataset <- rbind(LIM_dataset, CERS_dataset, POU_dataset, SOX_dataset, NKL_dataset, HOX_dataset, ZF_dataset, SINE_dataset, HNF_dataset, CUT_dataset)

# Taxonomy subsets
Non_Bilateria <- c("Porifera", "Ctenophora", "Cnidaria", "Placozoa")
Bilateria <- c("Arthropoda", "Nematoda", "Priapulida", "Tardigrada", "Mollusca", "Annelida", "Platyhelminthes", "Rotifera", "Nemertea", "Bryozoa", "Phoronida", "Brachiopoda", "Chordata", "Chordata", "Echinodermata", "Hemichordata", "Acoelomorpha", "Xenoturbellida")
Unicellular_Holozoa <- c("Choanoflagellata", "Ichthyosporea", "Filasterea", "Corallochytrea")
Fungi <- c("Ascomycota", "Basidiomycota", "Mucoromycota", "Zoopagomycota", "Chytridiomycota", "Blastocladiomycota")
dataset <- dataset %>%
  mutate(Classification = case_when(
    Taxonomy %in% Non_Bilateria        ~ "Non_Bilateria",
    Taxonomy %in% Bilateria            ~ "Bilateria",
    Taxonomy %in% Unicellular_Holozoa  ~ "Unicellular_Holozoa",
    Taxonomy %in% Fungi                ~ "Fungi",
    TRUE                               ~ NA_character_
  ))
dataset <- dataset %>%
  mutate(
    Classification = factor(
      Classification,
      levels = c("Bilateria", "Non_Bilateria", "Unicellular_Holozoa", "Fungi")
    )
  ) %>%
  arrange(Classification, Taxonomy) %>%   # orders rows by Classification first
  mutate(Taxonomy = factor(Taxonomy, levels = unique(Taxonomy)))

# Duplicated headers
duplications <- dataset[duplicated(dataset$Header) | duplicated(dataset$Header, fromLast = TRUE), ]

# Cleaning LIM 
duplications_LIM_clean <- duplications[
  duplications$Header %in% duplications$Header[duplications$Gene == "LIM"] &
    duplications$Gene != "LIM",
]
# Cleaning CUT
duplications_CUT_clean <- duplications[
  duplications$Header %in% duplications$Header[duplications$Gene == "CUT"] &
    duplications$Gene != "CUT",
]
# Cleaning POU
duplications_POU_clean <- duplications[
  duplications$Header %in% duplications$Header[duplications$Gene == "POU"] &
    duplications$Gene != "POU",
]
# Cleaning SINE
duplications_SINE_clean <- duplications[
  duplications$Header %in% duplications$Header[duplications$Gene == "SINE"] &
    duplications$Gene != "SINE",
]
# Cleaning HOX
duplications_HOX <- duplications[
  grepl(paste(queries_HOX, collapse = "|"),
        duplications$Header,
        ignore.case = TRUE),
]
duplications_HOX_clean <- duplications_HOX[duplications_HOX$Gene!="HOX", ]
# Cleaning NKL
duplications_NKL <- duplications[
  grepl(paste(queries_NKL, collapse = "|"),
        duplications$Header,
        ignore.case = TRUE),
]
duplications_NKL_clean <- duplications_NKL[duplications_NKL$Gene!="NKL", ]
# Cleaning ZF
duplications_ZF <- duplications[
  grepl(paste(queries_ZF, collapse = "|"),
        duplications$Header,
        ignore.case = TRUE),
]
duplications_ZF_clean <- duplications_ZF[duplications_ZF$Gene!="ZF", ]

# Joining duplications
duplications_to_clean <- rbind(duplications_LIM_clean, duplications_CUT_clean, duplications_POU_clean, duplications_SINE_clean, duplications_HOX_clean, duplications_ZF_clean, duplications_NKL_clean)

# Delete duplications
dataset_clean <- anti_join(dataset, duplications_to_clean,
                           by = colnames(dataset))
duplications_new <- dataset_clean[duplicated(dataset_clean$Header) | duplicated(dataset_clean$Header, fromLast = TRUE), ]
dataset_clean <- anti_join(dataset_clean, duplications_new,
                           by = colnames(dataset_clean))

