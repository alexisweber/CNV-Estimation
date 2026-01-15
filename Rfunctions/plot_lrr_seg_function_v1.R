plot_lrr_seg <- function(df, seg_df, prop_gain) {
  #' @description plot scatter plot with trend line of log-r-ratio scores across genomic positions
  #' @input df: sample genotype dataframe with at least columns SNP, Chr, Position, GType, Log.R.Ration and B.Allele.Freq
  #' @input seg_df: sample segments dataframe with at least columns ID, chrom, loc.start, loc.end, num.mark, seg.mean, CN where CN is triplicated
  #' @input prop_gain: integer, proportion of total segment that are CN=3
  #' @return p: LRR scatter plot of lrr score along genomic positions per chromosome 
  #' 
  p <- ggplot(df, aes(x = Position, y = Log.R.Ratio)) +
    { if (nrow(seg_df) > 0)
      geom_rect(data = seg_df,
                aes(xmin = loc.start, xmax = loc.end, ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "grey90", alpha = 0.8)
    } +
    geom_point(alpha = 0.35, size = 0.4) +
    #geom_smooth(se = FALSE, span = 0.01, size = 0.5, color = "black") +
    labs(title = "Chr21 Log R Ratio (LRR)",
         subtitle = sprintf("Shaded = CN=3 segments • Proportion chr21 gain = %.2f", prop_gain),
         x = "Genomic position (bp)", y = "LRR") +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_bw(base_size = 11)
  return(p)
}