library(sangerseqR);
library(Biostrings);
library(GenomicRanges);
library(tidyverse);
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
id.exon <- gsub("^.+-(\\d{2})\\w_\\w\\d{2}\\.ab1", "E\\1", basename(fn.ab1));
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
#lapply(lst, function(x) ggchrom(x$seq, x$sangerseq, x$title));


# pairWise alignment
lst.perExon <- split(lst, id.exon);
ab1.exon <- lapply(lst.perExon, function(x)
    lapply(x, function(y) primarySeq(y$sangerseq)))


# Read reference sequence and annotation of VWF (NM_000552.4)
fa <- readDNAStringSet("ref/NM_000552.4.fa");
gr <- read.delim("ref/NM_000552.4.gff3", header = F, comment.char = "#") %>%
    filter(V3 == "exon") %>%
    separate(V9, into = "id", sep = ";", extra = "drop") %>%
    mutate(id = sub("ID=id", "E", id)) %>%
    select(V1, V4, V5, id) %>%
    makeGRangesFromDataFrame(
        keep.extra.columns = TRUE,
        seqnames.field = "V1",
        start.field = "V4",
        end.field = "V5");


# Extract exon sequences
ir <- ranges(gr);
names(ir) <- gr$id;
seq.exon <- extractAt(unlist(fa), at = ir);



seq.exon <- seq.exon[names(seq.exon) %in% id.exon]


pairwiseAlignment(seq.exon[1], ab1.exon[[1]][[1]], type = "global-local")



tmp <- lapply(ab1.exon, function(x) DNAStringSet(x))
