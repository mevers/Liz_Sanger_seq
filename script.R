library(sangerseqR);
library(tidyverse);
library(Biostrings);
library(grid);



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


# Parse fn.ab1 to extract exon ID
ID <- gsub("(^\\d-VMF-|_\\w{3}\\.ab1)", "", basename(fn.ab1));


ab1 <- lapply(fn.ab1, readsangerseq);
seq <- lapply(fn.seq, function(x)
    DNAString(paste0(readLines(x), collapse = "")));


# Sanity check
stopifnot(
    all(sapply(ab1, function(x) x@primarySeq@length) == sapply(seq, length)))


split(lapply(ab1, function(x) x@primarySeq), ID)



x <- unlist(strsplit(as.character(seq[[1]]), ""));

set.seed(2018);
df <- cbind.data.frame(nt = x, pos = 1:length(x), val = sample(1:10, length(x), replace = T));


maxN <- 250;
df %>%
    mutate(bin = cut(pos, breaks = seq(0, n() + maxN, by = maxN))) %>%
    group_by(bin) %>%
    mutate(pos.in.bin = factor(1:n(), levels = as.character(1:maxN))) %>%
    complete(pos.in.bin) %>%
    ungroup() %>%
    mutate(pos = 1:n()) %>%
    mutate(
        nt.A = ifelse(nt == "A", as.character(nt), NA),
        nt.C = ifelse(nt == "C", as.character(nt), NA),
        nt.G = ifelse(nt == "G", as.character(nt), NA),
        nt.T = ifelse(nt == "T", as.character(nt), NA)) %>%
    ggplot(aes(x = pos, val)) +
        geom_line() +
        facet_wrap(~ bin, scales = "free_x", ncol = 1) +
        theme_bw() +
        geom_text(
            aes(label = nt.A, x = pos, y = -Inf),
            size = 1.5,
            vjust = -1,
            colour = "green") +
        geom_text(
            aes(label = nt.C, x = pos, y = -Inf),
            size = 1.5,
            vjust = -1,
            colour = "blue") +
        geom_text(
            aes(label = nt.G, x = pos, y = -Inf),
            size = 1.5,
            vjust = -1,
            colour = "black") +
        geom_text(
            aes(label = nt.T, x = pos, y = -Inf),
            size = 1.5,
            vjust = -1,
            colour = "red") + 
    labs(x = "Position [in nt]", y = "Value") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank());
ggsave("example.png", width = 11.69, height = 8.27);
