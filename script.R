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
x2 <- unlist(strsplit(as.character(primarySeq(ab1[[1]])), ""));


# Add mismatch
x2[10] <- "A";



set.seed(2018);
df <- cbind.data.frame(nt = x, nt2 = x2, pos = 1:length(x), val = sample(1:10, length(x), replace = T));


maxN <- 200;
tmp <- df %>%
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
        nt.T = ifelse(nt == "T", as.character(nt), NA),
        flag = ifelse(nt == nt2, NA, 1));



library(gtable);
theme_set(
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black")))


## Different attempt
# https://stackoverflow.com/questions/17492230/how-to-place-grobs-with-annotation-custom-at-precise-areas-of-the-plot-region/17493256#17493256
# https://stackoverflow.com/questions/22818061/annotating-facet-title-as-strip-over-facet

gg.val <- ggplot(tmp, aes(x = pos, val)) +
        geom_line() +
        facet_wrap(~ bin, scales = "free_x", ncol = 1) +
         geom_rect(aes(
            xmin = flag * (pos - 0.5),
            xmax  = flag * (pos + 0.5),
            ymin = -Inf, ymax = +Inf), alpha = 0.4, fill = "grey") +
        labs(x = "Position [in nt]", y = "Value") +
        scale_x_continuous(breaks = seq(0, length(x), by = 20));



gg.seq <- ggplot(tmp, aes(x = pos, y = 1)) +
    facet_wrap(~ bin, scales = "free_x", ncol = 1) +
    geom_text(
        aes(label = nt2, x = pos, y = 1.2),
        size = 1.5,
        vjust = +2) +
    geom_text(
        aes(label = nt.A, x = pos, y = 1),
        size = 1.5,
        vjust = -1,
        colour = "green") +
    geom_text(
        aes(label = nt.C, x = pos, y = 1),
        size = 1.5,
        vjust = -1,
        colour = "blue") +
    geom_text(
        aes(label = nt.G, x = pos, y = 1),
        size = 1.5,
        vjust = -1,
        colour = "black") +
    geom_text(
        aes(label = nt.T, x = pos, y = 1),
        size = 1.5,
        vjust = -1,
        colour = "red");


g.val <- ggplotGrob(gg.val);
g.seq <- ggplotGrob(gg.seq);
idx <- subset(g.val$layout, grepl("panel", g.val$layout$name), select = t:l);
pos <- cumsum(c(idx$t[1] - 1, diff(idx$t) + 1));
for (i in 1:length(pos)) {
    g.val <- gtable::gtable_add_rows(
        x = g.val,
        heights = unit(1, "strwidth", "AC"),
        pos = pos[i]);
}
g.val <- gtable_add_grob(
    x = g.val,
    grobs = g.seq$grobs[grepl("panel", g.seq$layout$name)],
    t = pos + 1,
    l = 4);
grid.newpage();
grid.draw(g.val);


ggsave("example.png", g.val, width = 11.69, height = 8.27);
