require(ggplot2);
require(sangerseqR);
require(tidyverse);
require(Biostrings);
require(grid);
require(gtable);


# Plot Sanger sequencing results
ggchrom <- function(seq, ab1, title = "", maxN = 200, col.ACGT = c("green", "blue", "black", "red")) {

    # Sanity check
    stopifnot(class(seq) %in% "DNAString");
    stopifnot(class(ab1) %in% "sangerseq");
    stopifnot(length(col.ACGT) == 4);


    # Get length of sequence
    len <- length(seq);


    # Reference, primary and secondary sequence
    df.seq <- cbind.data.frame(
        nt.ref = unlist(strsplit(as.character(seq), "")),
        nt.pri = unlist(strsplit(as.character(primarySeq(ab1)), "")),
        nt.sec = unlist(strsplit(as.character(secondarySeq(ab1)), "")),
        pos = 1:len);


    # Add mismatch
    df.seq$nt.pri[10] <- "A";


    # Trace matrix
    tm <- traceMatrix(ab1);
    df.tm <- matrix(tm, ncol = 4, dimnames = list(NULL, c("A", "C", "G", "T"))) %>%
        as.data.frame() %>%
        mutate(
            n = 1:n(),
            pos.float = n / max(n) * len,
            pos = ceiling(pos.float)) %>%
        group_by(pos) %>%
        summarise_at(vars(A:T), sum);


    # Merge sequence and trace data, and reformat for plotting
    df <- left_join(df.seq, df.tm, by = "pos") %>%
        mutate(bin = cut(pos, breaks = seq(0, n() + maxN, by = maxN))) %>%
        group_by(bin) %>%
        mutate(pos.in.bin = factor(1:n(), levels = as.character(1:maxN))) %>%
        complete(pos.in.bin) %>%
        ungroup() %>%
        mutate(pos = 1:n()) %>%
        mutate(flag = ifelse(nt.ref == nt.pri, NA, 1)) %>%
        gather(nt.tr, val, A:T)


    # Set general ggplot theme
    theme_set(
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour = "black")));


    # Trace plot
    gg.val <- df %>%
        ggplot(aes(x = pos, y = val)) +
            geom_line(aes(colour = nt.tr)) +
#            geom_bar(aes(fill = nt.tr), stat = "identity") +
            facet_wrap(~ bin, scales = "free_x", ncol = 1) +
            geom_rect(aes(
                xmin = flag * (pos - 0.5),
                xmax  = flag * (pos + 0.5),
                ymin = -Inf, ymax = +Inf), alpha = 0.4, fill = "grey") +
            labs(x = "Position [in nt]", y = "Signal", title = title) +
            scale_x_continuous(breaks = seq(0, len, by = 20)) +
#            scale_fill_manual(
#                values = c(
#                    "A" = col.ACGT[1],
#                    "C" = col.ACGT[2],
#                    "G" = col.ACGT[3],
#                    "T" = col.ACGT[4])) +
            scale_colour_manual(
                values = c(
                    "A" = col.ACGT[1],
                    "C" = col.ACGT[2],
                    "G" = col.ACGT[3],
                    "T" = col.ACGT[4])) +
            theme(legend.position = "bottom") +
#            guides(fill = guide_legend(title = "Nucleotide"))
            guides(colour = guide_legend(title = "Nucleotide"))


    # Produce sequence panels
    gg.seq <- df %>%
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
                values = c(
                    "A" = col.ACGT[1],
                    "C" = col.ACGT[2],
                    "G" = col.ACGT[3],
                    "T" = col.ACGT[4])) +
            scale_y_discrete(
                labels = c(
                    "nt.ref" = "Ref seq",
                    "nt.pri" = "Pri seq",
                    "nt.sec" = "Sec seq"));


    # Combine grobs from gg.val and gg.seq
    g.all <- ggplotGrob(gg.val);
    g.seq <- ggplotGrob(gg.seq);
    # Insert rows at position pos
    idx <- subset(g.all$layout, grepl("panel", g.all$layout$name), select = t:l);
    pos <- cumsum(c(idx$t[1] - 1, diff(idx$t) + 1));
    for (i in 1:length(pos)) {
        g.all <- gtable::gtable_add_rows(
            x = g.all,
            heights = unit(1, "strwidth", "AC"),
            pos = pos[i]);
    }
    # Add grobs from g.seq into new g.all rows
    g.all <- gtable_add_grob(
        x = g.all,
        grobs = g.seq$grobs[grepl("panel", g.seq$layout$name)],
        t = pos + 1,
        l = 4);
        grid.newpage();
        grid.draw(g.all);

    # Save as plot
    ggsave("example.png", g.all, width = 10, height = 1.5 * nrow(idx));
}
