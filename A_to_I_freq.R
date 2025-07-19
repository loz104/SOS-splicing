
library(Rsamtools)
library(GenomicRanges)
library(tidyverse)

fasta <- FaFile("reference.fasta")

ref_seq <- scanFa(fasta)[[1]]
ref_pos <- seq_along(ref_seq)
ref_base <- as.character(ref_seq)
is_A <- ref_base == "A"
a_pos <- ref_pos[is_A]

target_gr <- GRanges("reference", IRanges(start = a_pos, width = 1))

bamfile <- "sample.bam"
pup <- pileup(
  bamfile,
  index = paste0(bamfile, ".bai"),
  scanBamParam = ScanBamParam(which = target_gr),
  pileupParam = PileupParam(
    distinguish_nucleotides = TRUE,
    min_base_quality = 0,
    include_insertions = FALSE,
    include_deletions = FALSE
  )
)

# Calculate A-to-G frequencies
freq_df <- pup %>%
  filter(nucleotide %in% c("A", "G")) %>%
  group_by(pos, nucleotide) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  group_by(pos) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

all_pos <- tibble(pos = ref_pos)
full_df <- all_pos %>%
  crossing(nucleotide = c("A", "G")) %>%
  left_join(freq_df, by = c("pos", "nucleotide")) %>%
  mutate(freq = replace_na(freq, 0))

ggplot(full_df, aes(x = pos, y = freq, fill = nucleotide)) +
  geom_bar(stat = "identity", position = "stack")

