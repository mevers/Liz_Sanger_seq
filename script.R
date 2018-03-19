library(sangerseqR);
library(Biostrings);
source("ggchrom.R");


# From VWF sequencing data 012018/Notes.txt
keys <- list(
    "1" = "Father",
    "2" = "Mother",
    "3" = "Sister",
    "4" = "Brother")


# Sanger sequencing files
fn.ab1 <- list.files(
    path = "VWF sequencing data 012018",
    pattern = "*.ab1$",
    full.names = TRUE,
    recursive = TRUE);


# Reference sequence files
fn.seq <- sub("ab1$", "txt", fn.ab1);


# Parse fn.ab1 to extract exon ID
ID <- gsub("(^\\d-VMF-|_\\w{3}\\.ab1)", "", basename(fn.ab1));
title <- gsub("\\.ab1", "", basename(fn.ab1));


ab1 <- lapply(fn.ab1, readsangerseq);
seq <- lapply(fn.seq, function(x)
    DNAString(paste0(readLines(x), collapse = "")));


# Sanity check
stopifnot(
    all(sapply(ab1, function(x) x@primarySeq@length) == sapply(seq, length)))


# Plot Sanger sequencing results
ggchrom(seq[[1]], ab1[[1]], title[1]);
