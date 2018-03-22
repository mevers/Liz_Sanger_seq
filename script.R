library(sangerseqR);
library(Biostrings);
source("ggchrom.R");


# From VWF sequencing data 012018/Notes.txt
keys <- list(
    "1" = "Father",
    "2" = "Mother",
    "3" = "Sister",
    "4" = "Brother");


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



# Read sequencing and sequence data
ab1 <- lapply(fn.ab1, readsangerseq);
seq <- lapply(fn.seq, function(x)
    DNAString(paste0(readLines(x), collapse = "")));


# Sanity check
stopifnot(length(ab1) == length(seq));
stopifnot(
    all(sapply(ab1, function(x) x@primarySeq@length) == sapply(seq, length)))


# Combine ab1, seq and title in list
lst <- lapply(1:length(ab1), function(i)
    c(sangerseq = ab1[i], seq = seq[i], title = title[i]))


# Plot Sanger sequencing results
lapply(lst, function(x) ggchrom(x$seq, x$sangerseq, x$title));
