
library(tidyverse)
library(ggtranscript)

# Set exon1 range
end_of_exon1 <- targetlocus

# Read BED files
beds <- list.files("sample.bed")
names(beds) <- basename(beds)

bed_combined <- lapply(beds, function(bed) {
  read_tsv(bed, col_types = cols(.default = col_guess()))
}) %>% bind_rows(.id = "sample_name")

# Process exon1 intervals
bed_alt_exon1 <- bed_combined %>% 
  separate_longer_delim(c(blockStarts, blockSizes), delim = ",") %>%
  mutate(
    blockStarts = as.numeric(blockStarts),
    blockSizes = as.numeric(blockSizes),
    start = blockStarts + 1,
    end = blockStarts + blockSizes
  ) %>%
  select(-blockStarts, -blockSizes) %>%
  filter(end <= end_of_exon1) %>%
  group_by(sample_name, name) %>%
  summarize(
    starts = paste0(start, collapse = ","),
    ends = paste0(end, collapse = ","),
    .groups = "drop"
  )

# Summarize isoform usage
bed_alt_exon1_stats <- bed_alt_exon1 %>% 
  group_by(sample_name) %>% 
  mutate(reads_per_sample = n()) %>% 
  group_by(starts, ends, sample_name, reads_per_sample) %>% 
  summarize(count_per_sample = n(), .groups = "drop") %>%
  mutate(pct_reads_per_sample = count_per_sample / reads_per_sample * 100)

# Filter by minimal percentage
min_pct <- 2

bed_alt_exon1_stats_filt <- bed_alt_exon1_stats %>% 
  filter(pct_reads_per_sample > min_pct) %>%
  arrange(desc(pct_reads_per_sample)) %>%
  group_by(starts, ends) %>% 
  mutate(exon1_isoform_id = paste0("isoform", cur_group_id())) %>%
  mutate(exon1_isoform_id = factor(exon1_isoform_id, levels = rev(unique(exon1_isoform_id)))) %>%
  separate_longer_delim(c(starts, ends), delim = ",") %>%
  mutate(
    start = as.numeric(starts),
    end = as.numeric(ends)
  ) %>%
  ungroup() %>%
  select(-starts, -ends)

# Order isoforms
isoform_order <- bed_alt_exon1_stats_filt %>%
  select(exon1_isoform_id, pct_reads_per_sample) %>%
  distinct() %>%
  arrange(pct_reads_per_sample) %>%
  pull(exon1_isoform_id)

# Prepare data for plotting
for_ggtranscript <- bed_alt_exon1_stats_filt %>%
  select(start, end, exon1_isoform_id) %>%
  distinct() %>%
  mutate(exon1_isoform_id = factor(exon1_isoform_id, levels = isoform_order))

# Transcript structure plots
p <- ggplot(for_ggtranscript, aes(xstart = start, xend = end, y = exon1_isoform_id)) +
  geom_intron(data = to_intron(for_ggtranscript, "exon1_isoform_id")) +
  geom_line(aes(x = start, y = exon1_isoform_id)) +
  theme_minimal()

beginning <- p + coord_cartesian(xlim = targetlocus)
end <- p + coord_cartesian(xlim = targetlocus)

# Isoform usage plot
isoform_usage <- ggplot(
  bed_alt_exon1_stats_filt %>%
    select(exon1_isoform_id, pct_reads_per_sample, sample_name) %>%
    distinct(),
  aes(x = pct_reads_per_sample, y = reorder(exon1_isoform_id, pct_reads_per_sample))
) + 
  geom_col() +
  facet_wrap(~sample_name) +
  theme_minimal()

# Combine plots
combined_plot <- beginning + end + isoform_usage + patchwork::plot_layout(widths = c(1.5, 1.5, 1))
