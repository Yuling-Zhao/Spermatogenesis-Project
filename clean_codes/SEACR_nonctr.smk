#This is a snakefile using SEACR to call peaks for CUT&TAG experiments
#The X chromosome and autosomes were called peaks separately, to avoid the bias from copy number difference
#The IgG was not used as control in SEACR but rather also called peaks;
#IgG peaks were subtract from the peaks of experimental groups later by bedtools

#SEACR
#call peaks seperately in autosomes and chrX
import glob
import os
file_paths = glob.glob("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/SEACR/bed3/*.filtered.fragments.bedgraph")
file_names_without_suffix = [os.path.basename(file_path)[:-28] for file_path in file_paths]

wdir = "/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967"

rule all:
    input:
        expand("{wdir}/SEACR/bed3/auto_bed/{file_names}_auto.filtered.fragments.bedgraph", 
               wdir=wdir, file_names = file_names_without_suffix),
        expand("{wdir}/SEACR/bed3/chrX_bed/{file_names}_chrX.filtered.fragments.bedgraph", 
               wdir=wdir, file_names = file_names_without_suffix),
        expand("{wdir}/SEACR/peak/auto/{file_names}_nonctr_auto.stringent.bed", 
               wdir=wdir, file_names = file_names_without_suffix),
        expand("{wdir}/SEACR/peak/chrX/{file_names}_nonctr_chrX.stringent.bed", 
               wdir=wdir, file_names = file_names_without_suffix)

rule sep_auto_x:
    input:
        bedgraph="{wdir}/SEACR/bed3/{file_names}.filtered.fragments.bedgraph"
    output:
        auto_bed="{wdir}/SEACR/bed3/auto_bed/{file_names}_auto.filtered.fragments.bedgraph",
        x_bed="{wdir}/SEACR/bed3/chrX_bed/{file_names}_chrX.filtered.fragments.bedgraph"
    shell:
        """
        awk '$1 >= 1 && $1 <= 19' {input.bedgraph} > {output.auto_bed}
        awk '$1 == "X"' {input.bedgraph} > {output.x_bed}
        """     
        
rule SEACR_auto:
    input:
        bedgraph="{wdir}/SEACR/bed3/auto_bed/{file_names}_auto.filtered.fragments.bedgraph"
    output:
        bed="{wdir}/SEACR/peak/auto/{file_names}_nonctr_auto.stringent.bed"
    params:
        file_names=lambda wildcards: wildcards.file_names,
        outdir="{wdir}/SEACR/peak/auto"
    shell:
        """
        module load R
        bash /data/akhtar/group2/zhao/SEACR/SEACR_1.3.sh \
        {input.bedgraph} \
        0.05 non stringent \
        {params.outdir}/{params.file_names}_nonctr_auto
        """

rule SEACR_x:
    input:
        bedgraph="{wdir}/SEACR/bed3/chrX_bed/{file_names}_chrX.filtered.fragments.bedgraph"
    output:
        bed="{wdir}/SEACR/peak/chrX/{file_names}_nonctr_chrX.stringent.bed"
    params:
        file_names=lambda wildcards: wildcards.file_names,
        outdir="{wdir}/SEACR/peak/chrX"
    shell:
        """
        module load R
        bash /data/akhtar/group2/zhao/SEACR/SEACR_1.3.sh \
        {input.bedgraph} \
        0.05 non stringent \
        {params.outdir}/{params.file_names}_nonctr_chrX
        """
