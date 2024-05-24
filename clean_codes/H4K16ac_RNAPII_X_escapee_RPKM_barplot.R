#RPKM of each antibody (H4K16ac and Pol II Ser2p) from SPC and SPD was summarized in the region of interest (+/- 3kb from TSS) of selected genes
#The escapee genes are defined as genes with Pol II Ser2p peaks in SPD, because all the X-linked genes were inactive in SPC
#The control genes are defined as genes with H3K4me3 peaks at promoter but no peaks of Pol II Ser2p
#Chr1 genes are used as comparison to indicate the situation where the chromosomes have no dramatic changes of transcription status

#RPKM summarised into regions of interest of selected genes
#The tools and parameters used are the same as as those in summarize_RPKM_individual_chromosome.smk
SPD_Pol_X <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/intersect_counts/counts/H4K16ac_SPD_RNAPII_chrX_promoter_counts.txt",
                                header = TRUE)

SPD_H3K4me3_spec_X <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/intersect_counts/counts/H4K16ac_SPD_H3K4me3_spec_chrX_promoter_counts.txt",
                                         header = TRUE)

SPD_Pol_chr1 <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/intersect_counts/counts/SPD_chr1_Pol_gene_promoter.txt",
                        header = TRUE)

SPD_H3K4me3_spec_chr1 <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/intersect_counts/counts/SPD_chr1_H3K4me3_spec_gene_promoter.txt",
                                 header = TRUE)


#now we use X-linked genes which escaped from X chromosome inactivation as targets, to explore the H4K16ac intensity changes
SPC_H4K16ac_to_SPD_X_pol <- (SPD_Pol_X$SPC_H4K16ac_2.filtered)/(SPD_Pol_X$SPC_IgG.filtered+1)
SPC_H4K16ac_to_SPD_X_H3K4me3 <- (SPD_H3K4me3_spec_X$SPC_H4K16ac_2.filtered)/(SPD_H3K4me3_spec_X$SPC_IgG.filtered+1)

SPD_H4K16ac_to_SPD_X_pol <- (SPD_Pol_X$SPD_H4K16ac_2.filtered)/(SPD_Pol_X$SPD_IgG.filtered+1)
SPD_H4K16ac_to_SPD_X_H3K4me3 <- (SPD_H3K4me3_spec_X$SPD_H4K16ac_2.filtered)/(SPD_H3K4me3_spec_X$SPD_IgG.filtered+1)

#chr1 genes for comparison
SPC_H4K16ac_to_SPD_chr1_pol <- (SPD_Pol_chr1$SPC_H4K16ac_2.filtered)/(SPD_Pol_chr1$SPC_IgG.filtered+1)
SPC_H4K16ac_to_SPD_chr1_H3K4me3 <- (SPD_H3K4me3_spec_chr1$SPC_H4K16ac_2.filtered)/(SPD_H3K4me3_spec_chr1$SPC_IgG.filtered+1)

SPD_H4K16ac_to_SPD_chr1_pol <- (SPD_Pol_chr1$SPD_H4K16ac_2.filtered)/(SPD_Pol_chr1$SPD_IgG.filtered+1)
SPD_H4K16ac_to_SPD_chr1_H3K4me3 <- (SPD_H3K4me3_spec_chr1$SPD_H4K16ac_2.filtered)/(SPD_H3K4me3_spec_chr1$SPD_IgG.filtered+1)

# create a dataset
#H4K16ac chrX
med_SPC_control <- median(SPC_H4K16ac_to_SPD_X_H3K4me3)
SPC_H4K16ac <- data.frame(
  genes=c(rep("Transcribe",length(SPC_H4K16ac_to_SPD_X_pol)), rep("Control",length(SPC_H4K16ac_to_SPD_X_H3K4me3))),
  value=c(SPC_H4K16ac_to_SPD_X_pol/(med_SPC_control+1), 
          SPC_H4K16ac_to_SPD_X_H3K4me3/(med_SPC_control+1)),
  stage = rep("SPC", length(SPC_H4K16ac_to_SPD_X_pol)+length(SPC_H4K16ac_to_SPD_X_H3K4me3))
)

med_SPD_control <- median(SPD_H4K16ac_to_SPD_X_H3K4me3)
SPD_H4K16ac <- data.frame(
  genes=c(rep("Transcribe",length(SPD_H4K16ac_to_SPD_X_pol)), rep("Control",length(SPD_H4K16ac_to_SPD_X_H3K4me3))),
  value=c(SPD_H4K16ac_to_SPD_X_pol/(med_SPD_control+1), 
          SPD_H4K16ac_to_SPD_X_H3K4me3/(med_SPD_control+1)),
  stage = rep("SPD", length(SPD_H4K16ac_to_SPD_X_pol)+length(SPD_H4K16ac_to_SPD_X_H3K4me3))
)

H4K16ac_dat <- rbind(SPC_H4K16ac, SPD_H4K16ac)

my_sum <- H4K16ac_dat %>%
  group_by(genes,stage) %>%
  summarise( 
    n=n(),
    mean=mean(log2(value+1)),
    median=median(log2(value+1)),
    sd=sd(log2(value+1))
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

my_sum$chr <- factor(my_sum$stage, levels = c("SPC", "SPD"))

ggplot(my_sum, aes(x=genes, y=median, fill=genes, colour = genes)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.92), alpha = 0.6) + # Use dodge to separate genders
  facet_wrap(~stage, scales="free_x") + # Separate plots for each marker
  scale_fill_brewer(palette="Pastel1") + # Use a color palette that is easily distinguishable
  labs(title="Chromosome Comparison by Marker and Cell type",
       x="genes",
       y="Median",
       fill="genes") +
  scale_fill_manual(values=c("Transcribe"="#ffc48a", "Control"="#a369b0")) + # Set custom colors
  scale_colour_manual(values=c("Transcribe"="#ffc48a", "Control"="#a369b0")) + # Set custom colors for the error bars
  geom_errorbar(aes(ymin=median-ic, ymax=median+ic, group=interaction(genes, stage), colour = genes), 
                position=position_dodge(width=0.85), width=0.25, alpha=0.9, size=1.5) +
  theme_minimal() + # Use a minimal theme for a cleaner look
  theme(axis.text.x = element_text(angle=45, hjust=1)) # Rotate x labels for clarity

#H4K16ac chr1
med_SPC_control <- median(SPC_H4K16ac_to_SPD_chr1_H3K4me3)
SPC_H4K16ac <- data.frame(
  genes=c(rep("Transcribe",length(SPC_H4K16ac_to_SPD_chr1_pol)), rep("Control",length(SPC_H4K16ac_to_SPD_chr1_H3K4me3))),
  value=c(SPC_H4K16ac_to_SPD_chr1_pol/(med_SPC_control+1), 
          SPC_H4K16ac_to_SPD_chr1_H3K4me3/(med_SPC_control+1)),
  stage = rep("SPC", length(SPC_H4K16ac_to_SPD_chr1_pol)+length(SPC_H4K16ac_to_SPD_chr1_H3K4me3))
)

med_SPD_control <- median(SPD_H4K16ac_to_SPD_chr1_H3K4me3)
SPD_H4K16ac <- data.frame(
  genes=c(rep("Transcribe",length(SPD_H4K16ac_to_SPD_chr1_pol)), rep("Control",length(SPD_H4K16ac_to_SPD_chr1_H3K4me3))),
  value=c(SPD_H4K16ac_to_SPD_chr1_pol/(med_SPD_control+1), 
          SPD_H4K16ac_to_SPD_chr1_H3K4me3/(med_SPD_control+1)),
  stage = rep("SPD", length(SPD_H4K16ac_to_SPD_chr1_pol)+length(SPD_H4K16ac_to_SPD_chr1_H3K4me3))
)

#The method for cleaning data and plotting is the same as above used for chrX genes
H4K16ac_dat <- rbind(SPC_H4K16ac, SPD_H4K16ac)

#The same analysis was done for Pol II Ser2p on X-linked escapees and chr1 genes



