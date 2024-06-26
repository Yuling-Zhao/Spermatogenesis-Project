#The R script takes summarized RPKM values (from the output of summarize_RPKM_individual_chr.smk) at each chromosome and normalizes the counts by chromosome size (bp). 
#The normalized RPKM at chromosome level is calculated for H4K16ac, H3K9me3 and IgG in spermatogonia (SSC), spermatocytes (SPC) and spermatids (SPD). 
#The normalized RPKM is reffered as chromosome coverage in the script annotation, but is different from the traditional coverage, as it's not the percentage of each chromosome covered by reads. 
#To compare the coverage of autosomes and X chromosome, the normalized RPKM from H4K16ac and H3K9me3 was normalized again by the normalized RPKM of IgG, to get rid of the potential bias from copy number and accessibility of chromosomes.


#compare RPKM coverage of chrX and autosome, between H4K16ac and H3K9me3, from SPG, SPC and SPD####
all_cov <- read.table("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/marker_compare/chr_coverage/coverage/SSC_SPC_SPD_seqnorm_cov_merge.txt",
                      header = TRUE)
rownames(all_cov) <- all_cov$chr
chr <- c(1:19, "X", "Y")
all_cov <- all_cov[chr, ]
k16ac_cols <- c("SSC_IgG", "SPG_H4K16ac_1", "SPG_H4K16ac_2", "SPG_H4K16ac_3", 
                "MI_IgG_1", "MI_IgG_2", "SPC_H4K16ac_1", "SPC_H4K16ac_2",
                "MII_IgG_1", "MII_IgG_2", "SPD_H4K16ac_1", "SPD_H4K16ac_2", "SPD_H4K16ac_3")

k9me3_cols <- c("SSC_IgG", "SSC_H3K9me3_1", "SSC_H3K9me3_2", "SSC_H3K9me3_3", 
                "MI_IgG_1", "MI_IgG_2", "SPC_H3K9me3_1", "SPC_H3K9me3_2", "SPC_H3K9me3_3",
                "MII_IgG_1", "MII_IgG_2", "SPD_H3K9me3_1","SPD_H3K9me3_3")

k16ac_coverage <- all_cov[, k16ac_cols]
k9me3_coverage <- all_cov[, k9me3_cols]

#calculate sum counts for all chromosome
k16ac_cov_sum <- sapply(k16ac_coverage, sum)
k9me3_cov_sum <- sapply(k9me3_coverage, sum)

#calculate each chromosome coverage and normalize by chromosome size
k16ac_cov_norm <- sweep(k16ac_coverage, 2, k16ac_cov_sum, `/`) / all_cov$end * sum(all_cov$end)
k9me3_cov_norm <- sweep(k9me3_coverage, 2, k9me3_cov_sum, `/`) / all_cov$end * sum(all_cov$end)

#calculate the average autosome coverage
k16ac_cov_auto <- sapply(k16ac_cov_norm[1:19, ], mean)
k9me3_cov_auto <- sapply(k9me3_cov_norm[1:19, ], mean)

k16ac_cov_chrX <- k16ac_cov_norm["X", ]
k9me3_cov_chrX <- k9me3_cov_norm["X", ]

#prepare data for plotting
k16ac_dat <- data.frame(
  cov_norm = c(t(k16ac_cov_auto), t(k16ac_cov_chrX)),
  chr = rep(c("autosome", "chrX"), each = length(k16ac_cov_auto)),
  stage = rep(c(rep("SSC", 4), rep("SPC", 4), rep("SPD", 5)), 2),
  marker = rep(c("IgG", "H4K16ac", "H4K16ac", "H4K16ac", "IgG", "IgG", "H4K16ac", "H4K16ac", "IgG", "IgG", "H4K16ac", "H4K16ac", "H4K16ac"), 2)
)

k16ac_dat$stage <- factor(k16ac_dat$stage,
                          levels = c("SSC", "SPC", "SPD"))
  
k9me3_dat <- data.frame(
  cov_norm = c(t(k9me3_cov_auto), t(k9me3_cov_chrX)),
  chr = rep(c("autosome", "chrX"), each = length(k9me3_cov_auto)),
  stage = rep(c(rep("SSC", 4), rep("SPC", 5), rep("SPD", 4)), 2),
  marker = rep(c("IgG", "H3K9me3", "H3K9me3", "H3K9me3", "IgG", "IgG", "H3K9me3", "H3K9me3", "H3K9me3", "IgG", "IgG", "H3K9me3", "H3K9me3"), 2)
)
k9me3_dat$stage <- factor(k9me3_dat$stage,
                          levels = c("SSC", "SPC", "SPD"))

#normalize the coverage by IgG
#I am manually typing the IgG data from the above calculation
IgG_cov <- data.frame(
  cov_norm = c(1.0415272, 1.0673300, 1.0652466, 0.3989144, 0.3536154, 0.2824583),
  chr = c(rep("autosome", 3), rep("chrX", 3)),
  stage = c("SSC", "SPC", "SPD"),
  marker = c(rep("IgG", 6))
)

#the coverage from the calculation is not really how much percentage of each chromosome is covered, so it's unnecessary to be less than 1;
#It reflects the distribution of normalized counts to each chromosome, considering the bias from chromosome size (bigger chromosomes in principle may get more reads by chance);


k16ac_dat_norm <- c("cov_norm", "stage", "chr")
for (i in c("SSC", "SPC", "SPD")) {
  IgG_auto <- filter(IgG_cov, chr == "autosome" & stage %in% i)
  k16ac_auto <- filter(k16ac_dat, chr == "autosome" & stage %in% i & marker == "H4K16ac")
  k16ac_auto_norm <- k16ac_auto$cov_norm / IgG_auto$cov_norm
  
  IgG_chrX <- filter(IgG_cov, chr == "chrX" & stage %in% i)
  k16ac_chrX <- filter(k16ac_dat, chr == "chrX" & stage %in% i & marker == "H4K16ac")
  k16ac_chrX_norm <- k16ac_chrX$cov_norm / IgG_chrX$cov_norm
  
  cov_norm <- c(k16ac_auto_norm, k16ac_chrX_norm)
  stage <- rep(i, length(cov_norm))
  chr <- rep(c("autosome", "chrX"), each = length(k16ac_auto_norm))
  k16ac_stage <- data.frame(cov_norm, stage, chr)
  k16ac_dat_norm <- rbind(k16ac_dat_norm, k16ac_stage)
}
k16ac_dat_norm$cov_norm <- as.numeric(k16ac_dat_norm$cov_norm)
k16ac_dat_norm$stage <- factor(k16ac_dat_norm$stage,
                               levels = c("SSC", "SPC", "SPD"))

k9me3_dat_norm <- c("cov_norm", "stage", "chr")
for (i in c("SSC", "SPC", "SPD")) {
  IgG_auto <- filter(IgG_cov, chr == "autosome" & stage %in% i)
  k9me3_auto <- filter(k9me3_dat, chr == "autosome" & stage %in% i & marker == "H3K9me3")
  k9me3_auto_norm <- k9me3_auto$cov_norm / IgG_auto$cov_norm
  
  IgG_chrX <- filter(IgG_cov, chr == "chrX" & stage %in% i)
  k9me3_chrX <- filter(k9me3_dat, chr == "chrX" & stage %in% i & marker == "H3K9me3")
  k9me3_chrX_norm <- k9me3_chrX$cov_norm / IgG_chrX$cov_norm
  
  cov_norm <- c(k9me3_auto_norm, k9me3_chrX_norm)
  stage <- rep(i, length(cov_norm))
  chr <- rep(c("autosome", "chrX"), each = length(k9me3_auto_norm))
  k9me3_stage <- data.frame(cov_norm, stage, chr)
  k9me3_dat_norm <- rbind(k9me3_dat_norm, k9me3_stage)
}
k9me3_dat_norm$cov_norm <- as.numeric(k9me3_dat_norm$cov_norm)
k9me3_dat_norm$stage <- factor(k9me3_dat_norm$stage,
                               levels = c("SSC", "SPC", "SPD"))


#code example for dotplot using ggplot2
my_colors <- c("autosome" = "grey", "chrX" = "red")

library(ggplot2)
ggplot(k9me3_dat_norm[-1, ], aes(x=stage, y=cov_norm, fill=chr)) +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  labs(title="Coverage comparison between autosome and chrX",
       x="Stage",
       y="Coverage",
       fill="Chromosome")  +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
  ) +
  labs(fill = "Chromosome") +
  scale_fill_manual(values = my_colors)
