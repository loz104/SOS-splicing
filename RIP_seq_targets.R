
library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(limma)

#-----------------------------------------------
# Function: calculate enrichment between IP and Input

calc_enrich_pc <- function(df, input_col, ip_col, out_prefix, pseudocount = 1) {
  counts_mat <- as.matrix(df[, c(input_col, ip_col)]) + pseudocount
  rownames(counts_mat) <- df$gene
  
  colData <- DataFrame(condition = factor(c("Input", "IP")),
                       row.names = c(input_col, ip_col))
  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                colData = colData, design = ~condition)
  dds <- estimateSizeFactors(dds)
  norm_ct <- counts(dds, normalized = TRUE)
  
  enrichment <- norm_ct[, ip_col] / norm_ct[, input_col]
  
  res <- data.frame(
    gene = rownames(norm_ct),
    enrichment = enrichment,
    row.names = NULL
  )
  names(res)[2] <- paste0(out_prefix, "_enrich")
  res
}

# Calculate enrichment for each pair of IP and Input (nt: notag; F: Flag)

df_enr1 <- calc_enrich_pc(filtered, "ntinput1_count", "ntip1_count", "nt1")
df_enr2 <- calc_enrich_pc(filtered, "ntinput2_count", "ntip2_count", "nt2")
df_enr3 <- calc_enrich_pc(filtered, "Finput1_count", "Fip1_count", "F1")
df_enr4 <- calc_enrich_pc(filtered, "Finput2_count", "Fip2_count", "F2")

enr_all <- filtered %>%
  select(gene) %>%
  left_join(df_enr1, by = "gene") %>%
  left_join(df_enr2, by = "gene") %>%
  left_join(df_enr3, by = "gene") %>%
  left_join(df_enr4, by = "gene")

ranked_final <- enr_all %>%
  mutate(
    nt_mean    = (nt1_enrich + nt2_enrich) / 2,
    F_mean     = (F1_enrich  + F2_enrich)  / 2,
    ratio_F_nt = F_mean / nt_mean,
    R_nt2      = pmax(0, 1 - abs(log2(nt1_enrich) - log2(nt2_enrich))),
    R_F2       = pmax(0, 1 - abs(log2(F1_enrich)  - log2(F2_enrich))),
    reproducibility = R_nt2 * R_F2,
    final_score     = ratio_F_nt * reproducibility
  )

# Statistical significance (P-values) was assessed using the limma

expr <- df %>%
  select(gene, nt1_enrich, nt2_enrich, F1_enrich, F2_enrich) %>%
  column_to_rownames("gene") %>%
  as.matrix() %>%
  `+`(1) %>%
  log2()

group <- factor(c("noTag", "noTag", "tag", "tag"))
design <- model.matrix(~ group)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = 2, number = Inf) %>%
  rownames_to_column("gene")

#-----------------------------------------------
# Visualization of enriched genes (log2FC vs reproducibility)
# Assume 'df' contains log2FC and reproducibility info

df_pos <- subset(df, log2FC > 0)

p <- ggplot(df_pos, aes(x = log2FC, y = reproducibility)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 1.0, linetype = "dashed") +
  theme_bw(base_size = 14)

sig_genes <- subset(df_pos, reproducibility > 0.55 & log2FC > 1)

p <- p +
  geom_point(data = sig_genes) +
  geom_label_repel(data = sig_genes, aes(label = gene), max.overlaps = 10000)

print(p)


