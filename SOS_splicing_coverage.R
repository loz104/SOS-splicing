
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)

indexBam("sample.bam")
als <- list(sample = readGAlignments("sample.bam"))

all_cov <- lapply(names(als), function(sample) {
  coverage <- coverage(als[[sample]])
  cov <- as.vector(coverage[[1]])
  data.frame(pos = seq_along(cov),
             cov = cov,
             sample = sample) %>%
    mutate(covn = cov / max(cov))
}) %>% bind_rows()

p <- ggplot(all_cov, aes(x = pos, y = covn)) +
  geom_area() +
  theme_minimal()

print(p)


