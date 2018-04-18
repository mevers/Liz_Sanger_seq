library(readxl);

# Read mutation file from
# http://www.ragtimedesign.com/vwf/mutation/
df <- read_excel("docs/vwfmutations.xls", skip = 1) %>%
    rename(
        VWD_class = `VWD Classification`,
        exon_number = `Exon Number`);


# Plot number of mutations per exon per type
df %>%
    count(VWD_class, exon_number) %>%
    complete(VWD_class, exon_number) %>%
    ggplot(aes(VWD_class, n, fill = VWD_class)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ exon_number, scales = "free")


# Plot number of mutations per exon (x-axis) per VWD classification
df %>%
    filter(grepl("^\\d+$", exon_number)) %>%
    count(VWD_class, exon_number) %>%
    complete(VWD_class, exon_number) %>%
    mutate(exon_number = factor(
        exon_number,
        levels = as.character(sort(as.numeric(unique(exon_number)))))) %>%
    ggplot(aes(x = exon_number, y = n, fill = VWD_class)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ VWD_class, ncol = 2) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
        labs(
            x = "Exon",
            y = "Number of mutations",
            title = "Number of mutations per exon per VWD classifaction",
            caption = "[Source: http://www.ragtimedesign.com/vwf/mutation]") +
        scale_fill_brewer(palette = "Set2");
ggsave(file = "mutations_per_exon_per_VWD.pdf", height = 16, width = 12);




# Get VWD mutation list from Leiden Open Variation Database (LOVD)
if (!file.exists("df.LOVD.VWF.var.RData")) {
    library(rvest);
    url <- paste0(paste0(
        "https://grenada.lumc.nl/LOVD2/VWF/variants.php",
        "?select_db=VWF",
        "&action=view_all",
        "&limit=1000",
        "&order=Variant%2FDNA%2CASC"),
        c("&page=1", "&page=2"))
    lst.table <- lapply(url, function(x) {
        tables <- html_nodes(read_html(x), "table");
        idx.hdr <- 10;                                                           # Header is table 10
        hdr <- html_table(tables[idx.hdr], fill = TRUE)[[1]][, 1];               # Extract header
        if (identical(hdr[1], hdr[2])) hdr <- hdr[-1];                           # Remove duplicate 1st, 2nd entries
        hdr <- gsub("\\s", "_", hdr);                                            # Replace whitespace with "_"
        hdr <- gsub("#", "No", hdr);                                              # Replace "#" with "No"
        idx.data <- length(tables) - 1;                                          # Data is second last table
        data <- html_table(tables[idx.data], fill = TRUE)[[1]][-1, ];            # Extract data
        data[] <- lapply(data, as.character);                                    # Convert all columns to character
        setNames(data, hdr);
    })
    df.LOVD.VWF.var <- bind_rows(lst.table);
    save(df.LOVD.VWF.var, file = "df.LOVD.VWF.var.RData");
} else {
    load("df.LOVD.VWF.var.RData");
}


# Clean VWF LOVD variant data
df <- df.LOVD.VWF.var %>%
    select(Exon, DNA_change, Protein_change, VWD_type, No_Reported) %>%
    filter(
        str_detect(Exon, "^\\d+$"),
        !str_detect(VWD_type, "-")) %>%
    mutate(
        No_Reported = as.numeric(as.character(No_Reported)),
        VWD_type = gsub("^(type [1-3]\\w*).*$", "\\1", tolower(VWD_type)));


# Plot
df %>%
    mutate_if(is.character, as.factor) %>%
    group_by(VWD_type, Exon) %>%
    summarise(count = sum(No_Reported)) %>%
    complete(Exon) %>%
    mutate(Exon = factor(Exon, levels = as.character(sort(as.numeric(unique(as.character(Exon))))))) %>%
    ggplot(aes(x = Exon, y = count, fill = VWD_type)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ VWD_type, ncol = 2) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "bottom") +
        labs(
            x = "Exon",
            y = "Number of mutations",
            title = "Number of mutations per exon per VWD classifaction",
            caption = "[Source: https://grenada.lumc.nl/LOVD2/VWF/variants.php]") +
        scale_fill_brewer(palette = "Set2");
ggsave(file = "mutations_per_exon_per_VWD_LOVD.pdf", height = 16, width = 12);



df %>%
    complete(Exon, VWD_type) %>%
    count()
