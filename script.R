library(sangerseqR);
library(Biostrings);
library(GenomicRanges);
library(tidyverse);
library(msa);
source("ggchrom.R");
source("get.sangerseq.R");


# Samples
keys <- list(
    "QJL" = "Proband",
    "1" = "Father",
    "2" = "Mother",
    "3" = "Sister",
    "4" = "Brother",
    "ref" = "RefSeq");


# Read Sanger sequencing data into list
lst <- get.sangerseq(
    path = "rawdata",
    id.sample.re = list("^(.{1,3})-.+", "\\1"),
    id.exon.re = list("^.+V[MW]F-(\\d{1,2})(-[12])*\\w_\\w\\d{2}\\.ab1", "E\\1"),
    id.strand.re = list("^.+V[MW]F-\\d{1,2}(-[12])*(\\w)_\\w\\d{2}\\.ab1", "\\2"));



# Plot Sanger sequencing results
for (i in 1:length(lst)) ggchrom(
    lst[[i]]$sangerseq,
    lst[[i]]$title);


## Extract primer information
tab.primers <- do.call(rbind.data.frame, lapply(
    lst[sapply(lst,
        function(x) x$id.exon %in% c("E23", "E27", "E33", "E39"))],
    function(x) cbind(x$id.exon, x$id.strand, x$id.sample))) %>%
    spread(V3, V2);


# Convert list to long dataframe
df.seq <- do.call(rbind.data.frame, lapply(lst, function(x)
    cbind(
        id.sample = x$id.sample,
        id.exon = x$id.exon,
        id.strand = x$id.strand,
        seq = x$seq.pri))) %>%
    mutate_if(is.factor, as.character) %>%
    group_by(id.sample, id.exon) %>%
    arrange(id.sample, id.exon) %>%
    mutate(
        ntot = n(),
        id.exon.unique = ifelse(
            ntot > 1,
            sprintf("%s_rep%i", id.exon, 1:n()),
            id.exon)) %>%
    ungroup() %>%
    mutate(seq = ifelse(
        id.strand == "R",
        as.character(reverseComplement(DNAStringSet(seq))),
        seq)) %>%
    select(-ntot, -id.strand) %>%
    spread(id.sample, seq);


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
df.ref <- as_tibble(as.character((seq.ref))) %>%
    rownames_to_column("id.exon") %>%
    rename(ref = value);


# Left-join of Sanger and reference sequences; keep only entries where
# we have sequencing data from all source (1:4, QJL, ref)
df <- left_join(df.seq, df.ref, by = "id.exon") %>%
#    filter(complete.cases(.)) %>%
    column_to_rownames("id.exon.unique");


# Multiple sequence alignment
seq <- apply(df[, -1], 1, function(x) DNAStringSet(unlist(x[!is.na(x)])));
names(seq) <- rownames(df);
res.msa <- lapply(seq, function(x) {
    names(x) <- unlist(keys[names(x)]);
    msa(x, type = "dna", order = "input")
});


# Print results
for (i in 1:length(res.msa)) msaPrettyPrint(
    res.msa[[i]],
    file = sprintf("plots/%s.pdf", names(res.msa)[i]),
    paperWidth = 11.69, paperHeight = 8.27,
    margins = c(0.1, 0.2),
    askForOverwrite = FALSE,
    shadingMode = "identical",
    showConsensus = "none",
    logoColors = "accessible area",
    verbose = FALSE);



# NM_000552.3 transcript != NM_000552.4
#fa <- readDNAStringSet("ref/NM_000552.3.fa");
#df <- data.frame(
#    seqnames = "NM_000552.3",
#    start = c(640, 1701, 3264, 3860, 3942, 3936, 6018, 7489),
#    end = c(640, 1701, 3264, 3860, 3942, 3936, 6018, 7489),
#    id = c("rs2229444", "rs1800378", "rs749285654", "rs769502210", "rs61749368", "rs61749367", "rs778370191", "rs216867"));
#gr <- df %>%
#    makeGRangesFromDataFrame(keep.extra.columns = TRUE);
#ir <- ranges(gr);
#names(ir) <- gr$id;
#snp.seq <- extractAt(unlist(fa), at = ir);
