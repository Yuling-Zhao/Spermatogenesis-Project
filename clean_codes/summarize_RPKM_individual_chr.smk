#snakefile
import glob
import os
file_paths = glob.glob("/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/marker_compare/chr_coverage/seqnorm_bed/*.filtered.bedgraph")
file_names_without_suffix = [os.path.basename(file_path)[:-26] for file_path in file_paths]

wdir = "/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967"

rule all:
    input:
        expand("{wdir}/marker_compare/chr_coverage/coverage/{file_names}_seqnorm_cov.txt", 
               wdir=wdir, file_names = file_names_without_suffix),
        expand("{wdir}/marker_compare/chr_coverage/coverage/SSC_SPC_SPD_seqnorm_cov_merge.txt", wdir=wdir)

rule coverage_cal:
    input:
        bed = "{wdir}/marker_compare/chr_coverage/seqnorm_bed/{file_names}_seqnorm.filtered.bedgraph",
        genome = "/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/marker_compare/chr_coverage/genome_chr.bed"
    output:
        "{wdir}/marker_compare/chr_coverage/coverage/{file_names}_seqnorm_cov.txt"
    params:
        file_names=lambda wildcards: wildcards.file_names
    shell:
        """
        module load bedtools2
        bedtools intersect -a {input.genome} \
        -b {input.bed} \
        -wao | cut -f 1-3,7 | bedtools groupby -g 1-3 -c 4 -o sum > {params.file_names}_temp.txt
        (echo -e "chr\tstart\tend\t{params.file_names}"; cat {params.file_names}_temp.txt) > {output}
        rm {params.file_names}_temp.txt
        """
        
rule coverage_merge:
    input:
        coverage = expand("{wdir}/marker_compare/chr_coverage/coverage/{file_names}_seqnorm_cov.txt", 
                          wdir=wdir, file_names = file_names_without_suffix),
        genome = "/data/akhtar/group6/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/marker_compare/chr_coverage/genome_chr.bed"
    output:
        "{wdir}/marker_compare/chr_coverage/coverage/SSC_SPC_SPD_seqnorm_cov_merge.txt"
    shell:
        """
        (echo -e "chr\tstart\tend"; cat {input.genome}) > tmp_regions.txt
        paste {input.coverage} | awk '{{for(i=4;i<=NF;i+=4) printf "%s ", $i; print "\t"}}' > tmp_counts.txt
        paste tmp_regions.txt tmp_counts.txt > {output}
        rm tmp_regions.txt tmp_counts.txt
        """
