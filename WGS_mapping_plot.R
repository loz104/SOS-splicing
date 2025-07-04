
library(vcfR)
library(ggplot2)

vcf <- read.vcfR("sample.vcf")
vcf_data <- as.data.frame(vcf@fix)

vcf_data$calAF <- sapply(vcf_data$INFO, function(info) {
  dp4_match <- regmatches(info, regexpr("DP4=[0-9,]+", info))
  if (length(dp4_match) > 0) {
    dp4_values <- as.numeric(strsplit(sub("DP4=", "", dp4_match), ",")[[1]])
    if (length(dp4_values) == 4) {
      return((dp4_values[3] + dp4_values[4]) / sum(dp4_values))
    }
  }
  return(NA)
})

#filter by quality 
vcf_data$QUAL <- as.numeric(vcf_data$QUAL)
filtered_vcf_data <- vcf_data[vcf_data$QUAL > 30 & !is.na(vcf_data$calAF), ]


filtered_vcf_data$CHROM <- as.factor(filtered_vcf_data$CHROM)

pdf("SNP_Allele_Frequency.pdf")
ggplot(filtered_vcf_data, aes(x = as.numeric(POS) / 1e6, y = calAF)) +
  facet_wrap(~ CHROM, scales = "free_x") +
  theme_minimal()
dev.off()
