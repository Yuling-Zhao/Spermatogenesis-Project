#This script is run on bash
#This script shows example codes to get intersect of peaks from replicates, which was done separately on X chromosome and autosomes
#IgG peaks were substract from the peaks of experimental groups
#Total peaks of whole genome was from the merge of peaks from X chromosome and autosomes

module load bedtools2
#intersect
#chrX
bedtools intersect -a /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_1_nonctr_chrX.stringent.bed \
-b /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_2_nonctr_chrX.stringent.bed \
-f 0.1 -r > /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_nonctr_chrX_inter.stringent.bed

#autosome
bedtools intersect -a /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_1_nonctr_auto.stringent.bed \
-b /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_2_nonctr_auto.stringent.bed \
-f 0.1 -r > /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_nonctr_auto_inter.stringent.bed


#subtract IgG
#chrX
bedtools subtract -A -a /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_nonctr_chrX_inter.stringent.bed \
-b /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_IgG_nonctr_chrX.stringent.bed \
> /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_nonctr_chrX_inter_sub.stringent.bed

#autosome
bedtools subtract -A -a /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_nonctr_auto_inter.stringent.bed \
-b /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_IgG_nonctr_auto.stringent.bed \
> /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_nonctr_auto_inter_sub.stringent.bed

#merge chrX and autosome peaks
cat /data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/auto/SPC_RNAPIIser2_nonctr_auto_inter_sub.stringent.bed \
/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/chrX/SPC_RNAPIIser2_nonctr_chrX_inter_sub.stringent.bed \
>/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_3142/SEACR/peak/SPC_RNAPIIser2_nonctr_inter_sub.stringent.bed
