library(sangerseqR);
library(Biostrings);
library(GenomicRanges);
library(tidyverse);
library(msa);
source("ggchrom.R");
source("get.sangerseq.R");

# From VWF sequencing data 012018/Notes.txt
keys <- list(
    "1" = "Father",
    "2" = "Mother",
    "3" = "Sister",
    "4" = "Brother");


lst.prob <- get.sangerseq(
    path = "2016 sequencing",
    id.regexp = list("^.+-(\\d{1,2})\\w_\\w\\d{2}\\.ab1", "E\\1"));


lst.rltv <- get.sangerseq(
    path = "VWF sequencing data 012018",
    id.regexp = list("^.+-(\\d{1,2})\\w_\\w\\d{2}\\.ab1", "E\\1"));



lst <- get.sangerseq(
    path = "rawdata",
    id.sample.re = list("^(.{1,3})-.+", "\\1"),
    id.exon.re = list("^.+-(\\d{1,2})\\w_\\w\\d{2}\\.ab1", "E\\1"))


# Convert list to long dataframe
df.seq <- do.call(rbind.data.frame, lapply(lst, function(x)
    cbind(id.sample = x$id.sample, id.exon = x$id.exon, seq = x$seq.pri))) %>%
    mutate_if(is.factor, as.character) %>%
    group_by(id.sample, id.exon) %>%
    arrange(id.sample, id.exon) %>%
    mutate(ntot = n()) %>%
    mutate(id.exon.unique = ifelse(
        ntot > 1,
        sprintf("%s_rep%i", id.exon, 1:n()),
        id.exon)) %>%
    ungroup() %>%
    select(-ntot) %>%
    spread(id.sample, seq)



# Plot Sanger sequencing results
#lapply(lst, function(x) ggchrom(x$seq, x$sangerseq, x$title));


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
seq.ref <- extractAt(unlist(fa), at = ir);
df.ref <- as_tibble(as.character(reverseComplement(seq.ref))) %>%
    rownames_to_column("id.exon") %>%
    rename(ref = value)


# Left-join of Sanger and reference sequences; keep only entries where
# we have sequencing data from all source (1:4, QJL, ref)
df <- left_join(df.seq, df.ref, by = "id.exon") %>%
    filter(complete.cases(.)) %>%
    column_to_rownames("id.exon.unique");


# Multiple sequence alignment
seq <- apply(df[, -1], 1, function(x) DNAStringSet(unlist(x)));
res.msa <- lapply(seq, function(x) msa(x, type = "dna", order = "input"));
print(res.msa[[1]], show = "complete");
