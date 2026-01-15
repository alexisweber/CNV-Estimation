plot_baf_mosaic <- function(df, seg_df) {
  #' @description plot BAF distributions by BAF score count histogram 
  #' @input df: sample genotype dataframe with at least columns SNP, Chr, Position, GType, Log.R.Ration and B.Allele.Freq
  #' @input seg_df: sample segments dataframe with at least columns ID, chrom, loc.start, loc.end, num.mark, seg.mean, CN where CN is triplicated
  #' @return baf_plot: histogram count of each BAF score, separated by upper and lower clusters
  #' 
  
  if (nrow(seg_df) == 0) {
    return(ggplot() + theme_void() +
             geom_text(aes(0,0,label="No CN=3 region on chr21")) +
             xlab(NULL) + ylab(NULL))
  }
  
  # concatenate BAF from all CN=3 regions
  idx <- unlist(mapply(function(s,e) which(df$Position >= s & df$Position <= e), 
                       seg_df$loc.start, seg_df$loc.end, SIMPLIFY = FALSE))
  baf <- df$B.Allele.Freq[idx]
  het <- baf[baf > 0.05 & baf < 0.95]
  upper <- median(het[het > 0.5], na.rm = TRUE)
  lower <- median(het[het < 0.5], na.rm = TRUE)
  m_est <- 6 * (upper - 0.5)
  
  df <- data.frame(BAF = het)
  
  baf_plot <- ggplot(df, aes(x = BAF)) +
    geom_histogram(bins = 120) +
    geom_vline(xintercept = c(lower, 0.5, upper), linetype = c("dashed","solid","dashed")) +
    annotate("label", x = upper, y = Inf, vjust = 1.5,
             label = sprintf("upper≈%.3f\nm≈%.2f", upper, m_est), size = 3) +
    annotate("label", x = lower, y = Inf, vjust = 3.2,
             label = sprintf("lower≈%.3f", lower), size = 3) +
    labs(title = "BAF Score Distribution",
         subtitle = "Mosaic fraction m ≈ 6 × (upper band − 0.5)",
         x = "BAF (heterozygote cluster)", y = "Count") +
    theme_bw(base_size = 11)
  
  return(baf_plot)
}
