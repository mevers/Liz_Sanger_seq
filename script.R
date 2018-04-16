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
        seq.pri = x$seq.het1,
        seq.sec = x$seq.het2))) %>%
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
    mutate(
        seq.pri = ifelse(
            id.strand == "R",
            as.character(reverseComplement(DNAStringSet(seq.pri))),
            seq.pri),
        seq.sec = ifelse(
            id.strand == "R",
            as.character(reverseComplement(DNAStringSet(seq.sec))),
            seq.sec)) %>%
    select(-ntot, -id.strand) %>%
    gather(call, seq, seq.pri:seq.sec) %>%
    unite(id, id.sample, call, sep = "_") %>%
    spread(id, seq)



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
    names(x) <- sapply(strsplit(names(x), "_"), function(x)
        ifelse(
            !is.na(x[2]),
            paste(unlist(keys[x[1]]), gsub("seq\\.", "", x[2])),
            unlist(keys[x[1]])));
    msa(x, type = "dna", order = "input")
});
save(res.msa, file = "res.msa.RData");


# Print results
setwd("plots");
for (i in 1:length(res.msa)) msaPrettyPrint(
    res.msa[[i]],
    file = sprintf("%s.pdf", names(res.msa)[i]),
    paperWidth = 11.69, paperHeight = 8.27,
    margins = c(0.1, 0.2),
    askForOverwrite = FALSE,
    shadingMode = "identical",
    showConsensus = "none",
    logoColors = "accessible area",
    verbose = FALSE);
setwd("..");


# Determine number of mismatches
n.mm <- lapply(res.msa, function(x) {
    ss <- unmasked(x);
    idx.ref <- which(names(ss) == "RefSeq");
    pos.ref <- as.character(ss[idx.ref]) %>% str_locate_all("[^-]") %>% range();
    ss <- subseq(ss, pos.ref[1], pos.ref[2]);
    c(
        `Length ref` = diff(pos.ref) + 1,
        mapply(adist, as.character(ss[-idx.ref]), as.character(ss[idx.ref]), SIMPLIFY = F));
})
df.mm <- data.frame(bind_rows(lapply(n.mm, as.data.frame)), row.names = names(res.msa))


# Plot number of mismatches per exon
exons.ordered <- rownames(df.mm)[order(as.numeric(
    gsub("^E(\\d+).*$", "\\1", rownames(df.mm))))]
samples.ordered <- expand.grid(unlist(keys), c("pri", "sec")) %>%
    arrange(Var1, Var2) %>% mutate(id = paste(Var1, Var2, sep = ".")) %>% pull(id);
df.mm %>%
    rownames_to_column("Exon") %>%
    gather(Sample, Number_of_mm, -Exon, -Length.ref) %>%
    mutate(
        Exon = factor(Exon, levels = exons.ordered),
        Sample = factor(Sample, levels = samples.ordered),
        Number_of_mm_binned = cut(
            Number_of_mm,
            breaks = c(-1:10, max(Number_of_mm, na.rm = TRUE)),
            labels = as.character(c(0:10, ">10")))) %>%
    ggplot(aes(x = Sample, y = Exon, fill = Number_of_mm_binned)) +
    geom_tile() +
    geom_text(
        aes(label = Number_of_mm, colour = ifelse(Number_of_mm < 10, "0", "1")),
        size = 2) +
    scale_fill_brewer(palette = "Reds") +
    scale_colour_manual(values = c("black", "white"), guide = F) +
    labs(fill = "Number of mismatches") +
    theme_bw() +
    theme(legend.position = "bottom")
ggsave(file = "mismatches.pdf", height = 8, width = 12)
