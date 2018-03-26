# Necessary R libraries
library(sangerseqR);


# This function reads Sanger sequencing files recursively starting from path,
# and extracts two IDs from the filename based on the two regular expressions
# id1.regexp and id2.regexp
get.sangerseq <- function(path = ".", pattern = "*.ab1", id.sample.re = NULL, id.exon.re = NULL) {


    # Sanger sequencing files
    fn.ab1 <- list.files(
        path = path,
        pattern = "*.ab1$",
        full.names = TRUE,
        recursive = TRUE);


    # Parse fn.ab1 to extract exon ID
    id.sample <- gsub(id.sample.re[[1]], id.sample.re[[2]], basename(fn.ab1));
    id.exon <- gsub(id.exon.re[[1]], id.exon.re[[2]], basename(fn.ab1));
    title <- gsub("\\.ab1", "", basename(fn.ab1));


    # Read sequencing and sequence data
    ab1 <- lapply(fn.ab1, readsangerseq);


    # Combine ab1, seq and title in list
    lst <- lapply(1:length(ab1), function(i)
        c(
            seq.pri = as.character(primarySeq(ab1[[i]])),
            seq.sec = as.character(secondarySeq(ab1[[i]])),
            sangerseq = ab1[i],
            title = title[i],
            id.sample = id.sample[i],
            id.exon = id.exon[i]));


    # Return list
    return(lst);
}
