gain_coverage <- function(gt_df, seg_df, chr){
  #' @description input sample genotype dataframe and lrr segment dataframe and output the proportion of gain per segment
  #' @input gt_df: sample genotype dataframe with at least columns SNP, Chr, Position, GType, Log.R.Ration and B.Allele.Freq
  #' @input seg_df: sample segments dataframe with at least columns ID, chrom, loc.start, loc.end, num.mark, seg.mean, CN
  #' @input chr: string of chromosome studying
  #' @return gain_prop: integer, proportion of chromosome that is triplicated
  
  gt_total_len <- max(gt_df$Position) - min(gt_df$Position)
  gain_df <- seg_df %>% filter(chrom == chr & CN == 3)
  gain_prop = sum(gain_df$loc.end - gain_df$loc.start) / gt_total_len
  return(gain_prop)
}