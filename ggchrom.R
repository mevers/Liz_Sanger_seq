require(ggplot2);
require(sangerseqR);
require(tidyverse);
require(Biostrings);
require(grid);
require(gtable);


# Plot Sanger sequencing results
ggchrom <- function(
    ab1 = NULL,
    title = "",
    maxN = 200,
    col.ACGT = c("green", "blue", "black", "red"),
    ratio = 0.33) {


    # Sanity check
    stopifnot(class(ab1) %in% "sangerseq");
    stopifnot(length(col.ACGT) == 4);


    # Get length of sequence
    len <- length(primarySeq(makeBaseCalls(ab1, ratio = ratio)));


    # Reference, primary and secondary sequence
    df.seq <- cbind.data.frame(
        nt.het1 = unlist(strsplit(as.character(primarySeq(makeBaseCalls(
            ab1, ratio = ratio))), "")),
        nt.het2 = unlist(strsplit(as.character(secondarySeq(makeBaseCalls(
            ab1, ratio = ratio))), "")),
        pos = 1:len);


    # Ensure that primary and secondary sequences only contain
    # A, C, G, T, N
    df.seq <- df.seq %>%
        mutate_at(vars(contains("nt")), function(x) sub("[^ACGTN]", "N", x));


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
        gather(nt.tr, val, A:T);


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
        ggplot(aes(x = pos, y = log10(val))) +
            geom_rect(
                data = df %>%
                    filter(nt.het1 != nt.het2) %>%
                    group_by(pos, bin) %>%
                    summarise(val = 1),
                aes(xmin = pos - 0.5, xmax = pos + 0.5, ymin = 0, ymax = Inf),
                fill = "black", alpha = 0.2) +
            geom_line(aes(colour = nt.tr)) +
            geom_point(aes(colour = nt.tr), size = 0.5) +
            facet_wrap(~ bin, scales = "free_x", ncol = 1) +
            labs(x = "Position [in nt]", y = "log10 Signal", title = title) +
            scale_x_continuous(breaks = seq(0, len, by = 20)) +
            scale_colour_manual(
                values = c(
                    "A" = col.ACGT[1],
                    "C" = col.ACGT[2],
                    "G" = col.ACGT[3],
                    "T" = col.ACGT[4])) +
            theme(legend.position = "bottom") +
            guides(colour = guide_legend(title = "Nucleotide"));


    # Produce sequence panels
    gg.seq <- df %>%
        select(bin, pos, nt.het1, nt.het2) %>%
        gather(seq, nt, nt.het1, nt.het2) %>%
        mutate(
            seq = factor(seq, levels = rev(c("nt.het1", "nt.het2"))),
            col = case_when(
                seq == "nt.het1" & nt == "A" ~ "A",
                seq == "nt.het1" & nt == "C" ~ "C",
                seq == "nt.het1" & nt == "G" ~ "G",
                seq == "nt.het1" & nt == "T" ~ "T",
                TRUE ~ "G")) %>%
        ggplot() +
            facet_wrap(~ bin, scales = "free_x", ncol = 1) +
            geom_text(
#                aes(x = pos, y = seq, label = nt, colour = col),
                aes(x = pos, y = seq, label = nt),
                size = 1.5,
                show.legend = FALSE,
                fontface = "bold") +
            scale_colour_manual(
                values = c(
                    "A" = col.ACGT[1],
                    "C" = col.ACGT[2],
                    "G" = col.ACGT[3],
                    "T" = col.ACGT[4])) +
            scale_x_continuous(breaks = seq(0, len, by = 20)) +
            scale_y_discrete(
                labels = c(
                    "nt.het1" = "Pri seq",
                    "nt.het2" = "Sec seq"));


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
    cat(sprintf("Saving as %s.pdf\n", title));
    ggsave(
        file = sprintf("plots/%s.pdf", title),
        plot = g.all,
        width = 10,
        height = 1.5 * nrow(idx));
}
