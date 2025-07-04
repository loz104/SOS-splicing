
library(tidyverse)
library(tools)
library(ggpubr)

rmats_dir <- "path/to/rmats_output/"

rmats_files <- list.files(rmats_dir, pattern = "JC\\.", full.names = TRUE)
rmats_files <- rmats_files[!grepl("raw", rmats_files)]
names(rmats_files) <- file_path_sans_ext(basename(rmats_files))

rmats <- lapply(rmats_files, read_tsv) %>% bind_rows(.id = "filename")

rmats_for_boxplot <- rmats %>%
  mutate(event = gsub("\\.MATS.JC", "", filename)) %>%
  select(event, ID = ID...1, GeneID, IncLevel1, IncLevel2) %>%
  separate_longer_delim(c(IncLevel1, IncLevel2), delim = ",") %>%
  filter(!is.na(IncLevel1), !is.na(IncLevel2), IncLevel1 != "NA", IncLevel2 != "NA") %>%
  mutate(
    IncLevel1 = as.numeric(IncLevel1),
    IncLevel2 = as.numeric(IncLevel2)
  ) %>%
  group_by(event, ID, GeneID) %>%
  summarize(target = mean(IncLevel1), NC = mean(IncLevel2), .groups = "drop") %>%
  group_by(event) %>%
  mutate(event = paste0(event, ", n = ", n())) %>%
  pivot_longer(c(target, NC), names_to = "group", values_to = "PSI") %>%
  mutate(group = factor(group, levels = c("NC", "target")))

p <- ggplot(rmats_for_boxplot, aes(x = group, y = PSI)) +
  geom_boxplot() +
  facet_wrap(~event, nrow = 1) +
  stat_compare_means() +
  ylim(0, 1.1) +
  theme_minimal()

print(p)
