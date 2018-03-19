library(sangerseqR);
library(tidyverse);
library(Biostrings);
library(grid);
library(gtable);




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


set.seed(2018);
len <- sapply(seq, length)[1];
df <- cbind.data.frame(
    nt.ref = unlist(strsplit(as.character(seq[[1]]), "")),
    nt.pri = unlist(strsplit(as.character(primarySeq(ab1[[1]])), "")),
    nt.sec = unlist(strsplit(as.character(secondarySeq(ab1[[1]])), "")),
    pos = 1:len,
    val = sample(1:10, len, replace = T));


# Add mismatch
df$nt.pri[10] <- "A";



# Trace matrix
tm <- traceMatrix(ab1[[1]]);
df.tm <- matrix(tm, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T"))) %>%
    as.data.frame() %>%
    mutate(
        n = 1:n(),
        pos.float = n / max(n) * x,
        pos = round(pos.float)) %>%
    group_by(pos) %>%
    summarise_at(vars(A:T), sum)



maxN <- 200;
tmp <- df %>%
    mutate(bin = cut(pos, breaks = seq(0, n() + maxN, by = maxN))) %>%
    group_by(bin) %>%
    mutate(pos.in.bin = factor(1:n(), levels = as.character(1:maxN))) %>%
    complete(pos.in.bin) %>%
    ungroup() %>%
    mutate(pos = 1:n()) %>%
    mutate(flag = ifelse(nt.ref == nt.pri, NA, 1));


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


gg.val <- tmp %>%
    ggplot(aes(x = pos, val)) +
        geom_line() +
        facet_wrap(~ bin, scales = "free_x", ncol = 1) +
         geom_rect(aes(
            xmin = flag * (pos - 0.5),
            xmax  = flag * (pos + 0.5),
            ymin = -Inf, ymax = +Inf), alpha = 0.4, fill = "grey") +
        labs(x = "Position [in nt]", y = "Value") +
        scale_x_continuous(breaks = seq(0, length(x), by = 20));


# Produce sequence panels
gg.seq <- tmp %>%
    select(bin, pos, nt.ref, nt.pri, nt.sec) %>%
    gather(seq, nt, nt.ref, nt.pri, nt.sec) %>%
    mutate(
        seq = factor(seq, levels = rev(c("nt.ref", "nt.pri", "nt.sec"))),
        col = case_when(
            seq == "nt.ref" & nt == "A" ~ "A",
            seq == "nt.ref" & nt == "C" ~ "C",
            seq == "nt.ref" & nt == "G" ~ "G",
            seq == "nt.ref" & nt == "T" ~ "T",
            TRUE ~ "G")) %>%
    ggplot() +
        facet_wrap(~ bin, scales = "free_x", ncol = 1) +
        geom_text(
            aes(x = pos, y = seq, label = nt, colour = col),
            size = 1.5,
            show.legend = FALSE,
            fontface = "bold") +
        scale_colour_manual(
            values = c("A" = "green", "C" = "blue", "G" = "black", "T" = "red")) +
        scale_y_discrete(
            labels = c("nt.ref" = "Ref seq", "nt.pri" = "Pri seq", "nt.sec" = "Sec seq"));


# Combine grobs from gg.val and gg.seq
g.val <- ggplotGrob(gg.val);
g.seq <- ggplotGrob(gg.seq);
# Insert rows at position pos
idx <- subset(g.val$layout, grepl("panel", g.val$layout$name), select = t:l);
pos <- cumsum(c(idx$t[1] - 1, diff(idx$t) + 1));
for (i in 1:length(pos)) {
    g.val <- gtable::gtable_add_rows(
        x = g.val,
        heights = unit(1, "strwidth", "AC"),
        pos = pos[i]);
}
# Add grobs from g.seq into new g.val rows
g.val <- gtable_add_grob(
    x = g.val,
    grobs = g.seq$grobs[grepl("panel", g.seq$layout$name)],
    t = pos + 1,
    l = 4);
grid.newpage();
grid.draw(g.val);


ggsave("example.png", g.val, width = 11.69, height = 8.27);
