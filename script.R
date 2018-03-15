library(sangerseqR);
library(tidyverse);
library(Biostrings);


# From VWF sequencing data 012018/Notes.txt
keys <- list(
    1 = "Father",
    2 = "Mother",
    3 = "Sister",
    4 = "Brother")


# Sanger sequencing files
fn.ab1 <- list.files(
    path = "VWF sequencing data 012018",
    pattern = "*.ab1$",
    full.names = TRUE,
    recursive = TRUE);


# Reference sequence files
fn.seq <- sub("ab1$", "txt", fn.ab1);


ab1 <- lapply(fn.ab1, readsangerseq);
seq <- lapply(fn.seq, function(x)
    DNAString(paste0(readLines(x), collapse = "")));


#lapply(ab1, width)
