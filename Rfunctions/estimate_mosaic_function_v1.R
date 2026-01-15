estimate_mosaic <- function(df, chr, start, end) {
  #' @description estimate mosaic fraction m in each CN=3 segment from BAF band spacing
  #' @input df: sample genotype dataframe with at least columns SNP, Chr, Position, GType, Log.R.Ration and B.Allele.Freq
  #' @input chr: string, chromosome
  #' @input start: string, genomic position along chromosome that is the segment start
  #' @input end: string, genomic position along chromosome that is the segment end
  #' @return mos_df: dataframe with lower and upper BAF estimates and the mosaic fraction and number of heterozygotes per index
  #' 
  idx <- which(df$Position >= start & df$Position <= end & df$Chr==chr)
  x   <- df$B.Allele.Freq[idx]
  het <- x[x > 0.05 & x < 0.95]
  upper <- median(het[het > 0.55], na.rm = TRUE)
  lower <- median(het[het < 0.45], na.rm = TRUE)
  m_est <- 5.88 * (upper - 0.5) 
  mos_df <- data.frame(baf_lower=lower, baf_upper=upper, cn_fraction_m=round(m_est,3), n_het=length(het))
  return(mos_df)
}
