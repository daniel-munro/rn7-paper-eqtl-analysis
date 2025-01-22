library(tidyverse)

columns <- c("Sequence-Name", "Sequence-Role", "Assigned-Molecule",
             "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship",
             "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name")

names_rn6 <- read_tsv("data/GCA_000001895.4_Rnor_6.0_assembly_report.txt",
                      comment = "#", col_names = columns, col_types = "ccccccccic",
                      na = "na") %>%
    select(name_rn6 = "Sequence-Name",
           assigned_rn6 = "Assigned-Molecule",
           chrom_rn6 = "GenBank-Accn")

names_rn7 <- read_tsv("data/GCA_015227675.2_mRatBN7.2_assembly_report.txt",
                      comment = "#", col_names = columns, col_types = "ccccccccic",
                      na = "na") %>%
    select(name_rn7 = "Sequence-Name",
           assigned_rn7 = "Assigned-Molecule",
           chrom_rn7 = "GenBank-Accn")

chrs <- c(1:20, "X", "Y", "MT")

d <- read_tsv("chr_counts.txt", col_types = "ccci",
                  col_names = c("sample", "chrom_rn6", "chrom_rn7", "n")) %>%
    left_join(names_rn6, by = "chrom_rn6") %>%
    left_join(names_rn7, by = "chrom_rn7") %>%
    mutate(name_rn6 = str_replace(name_rn6, "chr", ""),
           name_rn6 = factor(name_rn6, levels = chrs),
           name_rn7 = factor(name_rn7, levels = chrs))

# d <- d_all %>%
#     mutate(
#         tmp = str_replace(name_rn6, "chr", ""),
#         chr_rn6 = if_else(!(tmp %in% chrs) | is.na(tmp), "Other", tmp),
#         chr_rn7 = if_else(!(name_rn7 %in% chrs) | is.na(name_rn7),
#                           "Other",
#                           name_rn7)
#     ) %>%
#     group_by(chr_rn6, chr_rn7) %>%
#     summarise(n = sum(n), .groups = "drop")

d %>%
    ggplot(aes(x = name_rn6, y = name_rn7, fill = n)) +
    facet_wrap(~ sample) +
    geom_tile() +
    # scale_fill_continuous(low = "white", high = "black") +
    scale_fill_viridis_c(trans = "log",
                         breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6),
                         labels = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
    ggtitle("Reads mapped to each chromosome")

# Without reads that didn't move:
d %>%
    filter(name_rn6 != name_rn7) %>%
    ggplot(aes(x = name_rn6, y = name_rn7, fill = n)) +
    facet_wrap(~ sample) +
    geom_tile() +
    # scale_fill_continuous(low = "white", high = "black") +
    scale_fill_viridis_c(trans = "log",
                         breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6),
                         labels = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
    ggtitle("Just relocated reads")

# Averaged:
d %>%
    group_by(name_rn6, name_rn7) %>%
    summarise(n = sum(n) / n_distinct(d$sample), .groups = "drop") %>%
    ggplot(aes(x = name_rn6, y = name_rn7, fill = n)) +
    geom_tile() +
    # scale_fill_continuous(low = "white", high = "black") +
    scale_fill_viridis_c(trans = "log",
                         breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6),
                         labels = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
    ggtitle("Reads mapped to each chromosome (mean across samples)")

d %>%
    filter(name_rn6 != name_rn7) %>%
    group_by(name_rn6, name_rn7) %>%
    summarise(n = sum(n) / n_distinct(d$sample), .groups = "drop") %>%
    ggplot(aes(x = name_rn6, y = name_rn7, fill = n)) +
    geom_tile() +
    # scale_fill_continuous(low = "white", high = "black") +
    scale_fill_viridis_c(trans = "log",
                         breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6),
                         labels = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6)) +
    theme_minimal() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) +
    ggtitle("Just relocated reads (mean across samples)")

# Total reads per chromosome
d %>%
    group_by(name_rn6) %>%
    summarise(n = sum(n)) %>%
    ggplot(aes(x = name_rn6, y = n)) +
    geom_col() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

d %>%
    group_by(name_rn7) %>%
    summarise(n = sum(n)) %>%
    ggplot(aes(x = name_rn7, y = n)) +
    geom_col() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))

# Fraction of reads remapped
with(d, sum(n[name_rn6 != name_rn7]) / sum(n))
# By sample
d %>%
    mutate(remapped = name_rn6 != name_rn7) %>%
    group_by(sample) %>%
    summarise(frac_remapped = sum(n[remapped]) / sum(n),
              .groups = "drop") %>%
    ggplot(aes(x = 1, y = frac_remapped)) +
    geom_boxplot(outlier.shape = NA) +
    ggbeeswarm::geom_beeswarm(size = 1, color = "blue")
