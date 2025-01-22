suppressPackageStartupMessages(library(tidyverse))

## GenBank chromosomes:
# keep6 <- c("CM000072.5", "CM000073.5", "CM000074.5", "CM000075.5", "CM000076.5",
#            "CM000077.5", "CM000078.5", "CM000079.5", "CM000080.5", "CM000081.5",
#            "CM000082.5", "CM000083.5", "CM000084.5", "CM000085.5", "CM000086.5",
#            "CM000087.5", "CM000088.5", "CM000089.5", "CM000090.5", "CM000091.5",
#            "CM000092.5", "CM002824.1", "AY172581.1")
# keep7 <- c("CM026974.1", "CM026975.1", "CM026976.1", "CM026977.1", "CM026978.1",
#            "CM026979.1", "CM026980.1", "CM026981.1", "CM026982.1", "CM026983.1",
#            "CM026984.1", "CM026985.1", "CM026986.1", "CM026987.1", "CM026988.1",
#            "CM026989.1", "CM026990.1", "CM026991.1", "CM026992.1", "CM026993.1",
#            "CM026994.1", "CM026995.1", "AY172581.1")
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
    d1 <- read_tsv(str_glue("data/read_chrs_rn6/{sample_id}.txt.gz"), col_types = "cc",
                   col_names = c("read", "chrom_rn6")) %>%
        filter(chrom_rn6 %in% keep6)
    dup <- unique(d1$read[duplicated(d1$read)])
    d1 <- filter(d1, !(read %in% dup))
    d2 <- read_tsv(str_glue("data/read_chrs_rn7/{sample_id}.txt.gz"), col_types = "cc",
                   col_names = c("read", "chrom_rn7")) %>%
        filter(chrom_rn7 %in% keep7)
    dup <- unique(d2$read[duplicated(d2$read)])
    d2 <- filter(d2, !(read %in% dup))
        # group_by(read) %>%
        # filter(n() == 1) %>%
        # ungroup()
    inner_join(d1, d2, by = "read") %>%
        count(chrom_rn6, chrom_rn7)
}

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
outfile <- args[2]

# samples <- c("00077E67B5_Acbc", "00077E8336_Acbc")
# Since this takes a while, just do one sample per run.

counts <- tibble(sample = sample_id) %>%
    group_by(sample) %>%
    summarise(chrom_counts(sample), .groups = "drop")

write_tsv(counts, outfile, col_names = FALSE) # no colnames to allow concat
