# Necessary R libraries
library(sangerseqR);


# This function reads Sanger sequencing files recursively starting from path,
# and extracts two IDs from the filename based on the three regular expressions
# id.sample.re, id.exon.re, and id.strand.re
get.sangerseq <- function(
    path = ".",
    pattern = "*.ab1",
    id.sample.re = NULL,
    id.exon.re = NULL,
    id.strand.re = NULL) {


    # Sanger sequencing files
    fn.ab1 <- list.files(
        path = path,
        pattern = "*.ab1$",
        full.names = TRUE,
        recursive = TRUE);


    # Parse fn.ab1 to extract exon ID
    id.sample <- gsub(id.sample.re[[1]], id.sample.re[[2]], basename(fn.ab1));
    id.exon <- gsub(id.exon.re[[1]], id.exon.re[[2]], basename(fn.ab1));
    id.strand <- gsub(id.strand.re[[1]], id.strand.re[[2]], basename(fn.ab1));
    title <- gsub("\\.ab1", "", basename(fn.ab1));


    # Read sequencing and sequence data
    ab1 <- lapply(fn.ab1, readsangerseq);


    # Combine ab1, seq and title in list
    lst <- lapply(1:length(ab1), function(i)
        c(
            seq.pri = as.character(primarySeq(ab1[[i]])),
            seq.sec = as.character(secondarySeq(ab1[[i]])),
            sangerseq = ab1[i],
            seq.het1 = as.character(primarySeq(makeBaseCalls(ab1[[i]], ratio = 0.33))),
            seq.het2 = as.character(secondarySeq(makeBaseCalls(ab1[[i]], ratio = 0.33))),
            title = title[i],
            id.sample = id.sample[i],
            id.exon = id.exon[i],
            id.strand = id.strand[i]));


    # Return list
    return(lst);
}
