plot_baf <- function(gt_df) {
  #' @description plot BAF distributions by genomic region by trisomy, facet wrap all chromosomes
  #' @input gt_df:sample genotype dataframe with at least columns SNP, Chr, Position, GType, Log.R.Ration and B.Allele.Freq
  #' @return p: facet wrapped individual chr BAF distribution by genomic region, colored by trisomy
  #' 
  
  
  gt_df <- gt_df %>%
    mutate(mean_lrr = mean(Log.R.Ratio[Chr == "21"], na.rm = TRUE), 
           col = case_when(Chr == "21" & mean_lrr > 0.025 ~ "red", TRUE ~ "darkgrey")) %>%
    dplyr::select(-mean_lrr) %>%
    mutate(Chr = factor(Chr, levels=c(1:22, "X", "Y")))
  
  p <- ggplot(gt_df, aes(x = Position, y = B.Allele.Freq, color = col)) +
    geom_point(size = 0.3) +
    facet_wrap(~ Chr, scales = "free_x", nrow = 1) +
    scale_color_identity() +
    theme_classic() +
    labs(x = "Position on Chromosome", y = "B Allele Frequency") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=20),
          axis.title.x = element_text(size=25),
          axis.title.y = element_text(size=25),
          strip.text = element_text(size = 20, face="bold"))
  
  return(p)
    
    
  
}
  
  #data <- gt3[c("Chr", "Position", paste0(k, ".B.Allele.Freq"))]
  #data$chr_position <- as.factor(paste0("chr", data$Chr))
  #data$sample <- gtDonor3[which(gtDonor3$DonorID == k), "DonorID"]
  #data$status <- gtDonor3[which(gtDonor3$DonorID == k), "Status"]
  #data$color <- ifelse(data$Chr == "21" & data$status == "T21", "red", "darkgrey")
  #data$Chr <- factor(data$Chr, levels=c(1:22, "X", "Y"))
  

