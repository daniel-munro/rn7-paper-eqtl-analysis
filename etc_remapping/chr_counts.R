suppressPackageStartupMessages(library(tidyverse))

# Chromosomes:
keep6 <- c("NC_005100.4", "NC_005101.4", "NC_005102.4", "NC_005103.4", "NC_005104.4",
           "NC_005105.4", "NC_005106.4", "NC_005107.4", "NC_005108.4", "NC_005109.4",
           "NC_005110.4", "NC_005111.4", "NC_005112.4", "NC_005113.4", "NC_005114.4",
           "NC_005115.4", "NC_005116.4", "NC_005117.4", "NC_005118.4", "NC_005119.4",
           "NC_005120.4", "NC_024475.1", "NC_001665.2")
keep7 <- c("NC_051336.1", "NC_051337.1", "NC_051338.1", "NC_051339.1", "NC_051340.1",
           "NC_051341.1", "NC_051342.1", "NC_051343.1", "NC_051344.1", "NC_051345.1",
           "NC_051346.1", "NC_051347.1", "NC_051348.1", "NC_051349.1", "NC_051350.1",
           "NC_051351.1", "NC_051352.1", "NC_051353.1", "NC_051354.1", "NC_051355.1",
           "NC_051356.1", "NC_051357.1", "NC_001665.2")

chrom_counts <- function(sample_id) {
    # Get only main chromosomes, and remove reads with multiple mappings among those.
    d1 <- read_tsv(str_glue("data/read_chrs_rn6/{sample_id}.txt.gz"),
                   col_types = "cc",
                   col_names = c("read", "chrom_rn6")) |>
        filter(chrom_rn6 %in% keep6)
    dup <- unique(d1$read[duplicated(d1$read)])
    d1 <- filter(d1, !(read %in% dup))
    d2 <- read_tsv(str_glue("data/read_chrs_rn7/{sample_id}.txt.gz"),
                   col_types = "cc",
                   col_names = c("read", "chrom_rn7")) |>
        filter(chrom_rn7 %in% keep7)
    dup <- unique(d2$read[duplicated(d2$read)])
    d2 <- filter(d2, !(read %in% dup))
    inner_join(d1, d2, by = "read") |>
        count(chrom_rn6, chrom_rn7)
}

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
outfile <- args[2]

# Since this takes a while, just do one sample per run.
counts <- tibble(sample = sample_id) |>
    summarise(chrom_counts(sample), .by = sample)

write_tsv(counts, outfile, col_names = FALSE) # no colnames to allow concat
