#This is Rscript annotating H3K4me3 peaks (at promoter) and Pol II Ser2p peaks (on gene body) to relative genes, based on the distance of peaks to genes;
#peaks are called by SEACR for each stage, using SEACR_nonctr.smk
#Only identical gene names are counted when calculate the ratio of genes on each chromosome;
#When calculate the overlap ratio between H3K4me3 and Pol II Ser2p peaked genes, using either total genes with H3K4me3 or Pol II Ser2p peaks as reference was plotted


#annotate H3K4me3 and Pol II Ser2p peaks to genes####
#for H3K4me3 SPG####
#load peaks
SPG_H3K4me3 <- read.table("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/chromHMM/peak_all/H3K4me3_ChIPSeq_AGSC_inter.stringent.bed")
SPG_H3K4me3 <- SPG_H3K4me3[, 1:3]
names(SPG_H3K4me3) <- c("chr", "start", "end")

#annotate with
library(GenomicFeatures)
txdb <- loadDb("/data/akhtar/group2/zhao/reference_data/GRCm38_ensembl_m9_txdb")

library(GenomicRanges)
SPG_H3K4me3_Grange <- GRanges(seqnames = SPG_H3K4me3$chr,
                             ranges = IRanges(start = SPG_H3K4me3$start, end = SPG_H3K4me3$end))

library(ChIPseeker)
#plot annotated regions
library(ggplot2)
SPG_H3K4me3_anno <- annotatePeak(SPG_H3K4me3_Grange, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = F)
SPG_H3K4me3_anno_dat <- as.data.frame(SPG_H3K4me3_anno@anno)
#get the genes with H3K4me3 at promoter
SPG_H3K4me3_genes <- SPG_H3K4me3_anno_dat[SPG_H3K4me3_anno_dat$annotation == "Promoter", ]
#for H3K4me3 SPC####
#load peaks
MI_H3K4me3 <- read.table("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/chromHMM/peak_all/MI_H3K4me3_nonctr_inter_sub_all.stringent.bed")
MI_H3K4me3 <- MI_H3K4me3[, 1:3]
names(MI_H3K4me3) <- c("chr", "start", "end")


library(GenomicRanges)
MI_H3K4me3_Grange <- GRanges(seqnames = MI_H3K4me3$chr,
                             ranges = IRanges(start = MI_H3K4me3$start, end = MI_H3K4me3$end))

library(ChIPseeker)
#plot annotated regions
library(ggplot2)
MI_H3K4me3_anno <- annotatePeak(MI_H3K4me3_Grange, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = F)
MI_H3K4me3_anno_dat <- as.data.frame(MI_H3K4me3_anno@anno)
#get the genes with H3K4me3 at promoter
SPC_H3K4me3_genes <- MI_H3K4me3_anno_dat[MI_H3K4me3_anno_dat$annotation == "Promoter", ]

#for H3K4me3 SPD####
#load peaks
MII_H3K3me3 <- read.table("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/chromHMM/peak_all/MII_H3K4me3_nonctr_inter_sub_all.stringent.bed")
MII_H3K3me3 <- MII_H3K3me3[, 1:3]
names(MII_H3K3me3) <- c("chr", "start", "end")

MII_H3K3me3_Grange <- GRanges(seqnames = MII_H3K3me3$chr,
                             ranges = IRanges(start = MII_H3K3me3$start, end = MII_H3K3me3$end))

MII_H3K3me3_anno <- annotatePeak(MII_H3K3me3_Grange, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = F)
MII_H3K3me3_anno_dat <- as.data.frame(MII_H3K3me3_anno@anno)
#get the genes with H3K4me3 at promoter
SPD_H3K4me3_genes <- MII_H3K3me3_anno_dat[MII_H3K3me3_anno_dat$annotation == "Promoter", ]

#for Pol II SPC####
#load peaks
SPC_Pol <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/SPC_RNAPIIser2_nonctr_inter_sub.stringent.bed")
SPC_Pol <- SPC_Pol[, 1:3]
names(SPC_Pol) <- c("chr", "start", "end")

SPC_Pol_Grange <- GRanges(seqnames = SPC_Pol$chr,
                             ranges = IRanges(start = SPC_Pol$start, end = SPC_Pol$end))

SPC_Pol_anno <- annotatePeak(SPC_Pol_Grange, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = F)
SPC_Pol_anno_dat <- as.data.frame(SPC_Pol_anno@anno)
#get the genes with Pol II on genes
SPC_Pol_anno@annoStat
#remove distal Intergenic, and consider the rest regions as gene related
SPC_Pol_genes <- SPC_Pol_anno_dat[!SPC_Pol_anno_dat$annotation == "Distal Intergenic", ]

table(SPC_Pol_genes$seqnames)
SPC_K4me3_Pol_olp <- SPC_Pol_genes[SPC_Pol_genes$geneId %in% SPC_H3K4me3_genes$geneId, ]

#for Pol II SPD####
#load peaks
SPD_Pol <- read.table("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/SPD_RNAPIIser2_nonctr_inter_sub.stringent.bed")
SPD_Pol <- SPD_Pol[, 1:3]
names(SPD_Pol) <- c("chr", "start", "end")

SPD_Pol_Grange <- GRanges(seqnames = SPD_Pol$chr,
                          ranges = IRanges(start = SPD_Pol$start, end = SPD_Pol$end))

SPD_Pol_anno <- annotatePeak(SPD_Pol_Grange, TxDb = txdb, tssRegion = c(-1000, 1000), verbose = F)
SPD_Pol_anno_dat <- as.data.frame(SPD_Pol_anno@anno)
#get the genes with Pol II on genes
SPD_Pol_anno@annoStat
#remove distal Intergenic, and consider the rest regions as gene related
SPD_Pol_genes <- SPD_Pol_anno_dat[!SPD_Pol_anno_dat$annotation == "Distal Intergenic", ]

table(SPD_Pol_genes$seqnames)
SPD_K4me3_Pol_olp <- SPD_Pol_genes[SPD_Pol_genes$geneId %in% SPD_H3K4me3_genes$geneId, ]


#the X:A ratio is much lower in Pol II (<0.1%) than in H3K4me3 in both stages
#plot the ratio of each chromosome with genes targeted by H3K4me3 and Pol II####
library(dplyr)
#there are already no unconventional chromosomes
#SPG H3K4me3

#SPC H3K4me3
#calculate percentage
SPG_H3K4me3_gene_dat <- SPG_H3K4me3_genes[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPG_H3K4me3_gene_dat$pro <- SPG_H3K4me3_gene_dat$n/sum(SPG_H3K4me3_gene_dat$n) * 100

#add reference for to color ChrX and autosome (Chr1)
ref_chr <- data.frame(
  seqnames = c(1:19, "X"),
  ref = c("Chromosome 1", rep("Autosome", 18), "Chromosome X")
)

SPG_H3K4me3_gene_dat_ref <- merge(SPG_H3K4me3_gene_dat, ref_chr, by = "seqnames")

SPG_H3K4me3_gene_dat_ref$marker <- rep("H3K4me3", nrow(SPG_H3K4me3_gene_dat_ref))
SPG_H3K4me3_gene_dat_ref$stage <- rep("SPG", nrow(SPG_H3K4me3_gene_dat_ref))


#calculate percentage
SPC_H3K4me3_gene_dat <- SPC_H3K4me3_genes[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPC_H3K4me3_gene_dat$pro <- SPC_H3K4me3_gene_dat$n/sum(SPC_H3K4me3_gene_dat$n) * 100

SPC_H3K4me3_gene_dat_ref <- merge(SPC_H3K4me3_gene_dat, ref_chr, by = "seqnames")

SPC_H3K4me3_gene_dat_ref$marker <- rep("H3K4me3", nrow(SPC_H3K4me3_gene_dat_ref))
SPC_H3K4me3_gene_dat_ref$stage <- rep("SPC", nrow(SPC_H3K4me3_gene_dat_ref))

#SPD H3K4me3
#calculate percentage
SPD_H3K4me3_gene_dat <- SPD_H3K4me3_genes[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPD_H3K4me3_gene_dat$pro <- SPD_H3K4me3_gene_dat$n/sum(SPD_H3K4me3_gene_dat$n) * 100

SPD_H3K4me3_gene_dat_ref <- merge(SPD_H3K4me3_gene_dat, ref_chr, by = "seqnames")

SPD_H3K4me3_gene_dat_ref$marker <- rep("H3K4me3", nrow(SPD_H3K4me3_gene_dat_ref))
SPD_H3K4me3_gene_dat_ref$stage <- rep("SPD", nrow(SPD_H3K4me3_gene_dat_ref))

all_H3K4me3_peak <- rbind(SPG_H3K4me3_gene_dat_ref, SPC_H3K4me3_gene_dat_ref, SPD_H3K4me3_gene_dat_ref)
all_H3K4me3_peak$stage <- factor(all_H3K4me3_peak$stage, 
                                 levels = c("SPG", "SPC", "SPD"))

#SPC Pol II
SPC_Pol_genes_dat <- SPC_Pol_genes[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPC_H3K4me3_Pol_overlap_dat <- SPC_K4me3_Pol_olp[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))


#calculate percentage
SPC_Pol_genes_dat$pro <- SPC_Pol_genes_dat$n/sum(SPC_Pol_genes_dat$n) * 100

#add reference for to color ChrX and autosome (Chr1)
ref_chr <- data.frame(
  seqnames = c(1:19, "X"),
  ref = c("Chromosome 1", rep("Autosome", 18), "Chromosome X")
)

SPC_Pol_genes_dat_ref <- merge(SPC_Pol_genes_dat, ref_chr, by = "seqnames")

SPC_Pol_genes_dat_ref$marker <- rep("H3K4me3", nrow(SPC_Pol_genes_dat_ref))
SPC_Pol_genes_dat_ref$stage <- rep("SPC", nrow(SPC_Pol_genes_dat_ref))

#SPD H3K4me3
#calculate percentage
SPD_Pol_genes_dat <- SPD_Pol_genes[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPD_H3K4me3_Pol_overlap_dat <- SPD_K4me3_Pol_olp[, c("seqnames", "geneId")] %>% 
  group_by(seqnames) %>% 
  summarise(n = n_distinct(geneId))

SPD_Pol_genes_dat$pro <- SPD_Pol_genes_dat$n/sum(SPD_Pol_genes_dat$n) * 100

#add reference for to color ChrX and autosome (Chr1)
ref_chr <- data.frame(
  seqnames = c(1:19, "X"),
  ref = c("Chromosome 1", rep("Autosome", 18), "Chromosome X")
)

SPD_Pol_genes_dat_ref <- merge(SPD_Pol_genes_dat, ref_chr, by = "seqnames")

SPD_Pol_genes_dat_ref$marker <- rep("Pol_II", nrow(SPD_Pol_genes_dat_ref))
SPD_Pol_genes_dat_ref$stage <- rep("SPD", nrow(SPD_Pol_genes_dat_ref))

all_Pol_peak <- rbind(SPC_Pol_genes_dat_ref, SPD_Pol_genes_dat_ref)
#plot one marker with different stages
library(ggplot2)
my_colors <- c("Autosome" = "grey", "Chromosome 1" = "blue", "Chromosome X" = "red")

ggplot(all_H3K4me3_peak, aes(x = factor(stage), y = pro, fill = factor(ref))) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.5) +
  theme_classic() +
  ylim(0, 12) +
  xlab("") +
  ylab("Proportion of peaks on each chromosome") +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
  ) +
  labs(fill = "Chromosome") +
  scale_fill_manual(values = my_colors)

#plot the overlap (between H3K4me3 and Pol II) in the background of H3K4me3 and Pol II all genes####
SPC_H3K4me3_Pol_overlap_dat$pro_K4me3 <- SPC_H3K4me3_Pol_overlap_dat$n/SPC_H3K4me3_gene_dat$n* 100
SPC_H3K4me3_Pol_overlap_dat$pro_Pol <- SPC_H3K4me3_Pol_overlap_dat$n/SPC_Pol_genes_dat$n* 100
SPC_H3K4me3_Pol_overlap_ref <- merge(SPC_H3K4me3_Pol_overlap_dat, ref_chr, by = "seqnames")

SPD_H3K4me3_Pol_overlap_dat$pro_K4me3 <- SPD_H3K4me3_Pol_overlap_dat$n/SPD_H3K4me3_gene_dat$n* 100
SPD_H3K4me3_Pol_overlap_dat$pro_Pol <- SPD_H3K4me3_Pol_overlap_dat$n/SPD_Pol_genes_dat$n* 100
SPD_H3K4me3_Pol_overlap_ref <- merge(SPD_H3K4me3_Pol_overlap_dat, ref_chr, by = "seqnames")

overlap_pro <- rbind(SPC_H3K4me3_Pol_overlap_ref, SPD_H3K4me3_Pol_overlap_ref)
overlap_pro$stage <- rep(c("SPC", "SPD"), each = nrow(SPC_H3K4me3_Pol_overlap_ref))

ggplot(overlap_pro, aes(x = factor(stage), y = pro_K4me3, fill = factor(ref))) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2) +
  theme_classic() +
  xlab("") +
  ylim(1, 75) +
  ylab("Proportion of overlap genes on each chromosome") +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 12)
  ) +
  labs(fill = "seqnames") +
  scale_fill_manual(values = my_colors)
