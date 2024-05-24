#This is a snakefile which generates bedgraph for peak calling with SEACR and further used for coverage of individual chromosomes analysis

wdir = "/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967"

import glob
import os
file_paths = glob.glob("/data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/ChIP-seq_ba/filtered_bam/*.bam")
file_names_without_suffix = [os.path.basename(file_path)[:-13] for file_path in file_paths]


rule all:
    input:
        expand("{wdir}/SEACR/bedgraph/{file_names_without_suffix}_seqnorm.filtered.bedgraph", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix),
        expand("{wdir}/SEACR/bamtobed/{file_names_without_suffix}.filtered.sorted.bam", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix),
        expand("{wdir}/SEACR/bamtobed/{file_names_without_suffix}.filtered.sorted.bed", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix),
        expand("{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.clean.bed", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix),
        expand("{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.fragments.bed", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix),
        expand("{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.fragments.bedgraph", 
               wdir=wdir, file_names_without_suffix=file_names_without_suffix)

rule bedgraph:
    input:
        bam="{wdir}/ChIP-seq_ba/filtered_bam/{file_names_without_suffix}.filtered.bam",
    output:
        bedgraph="{wdir}/SEACR/bedgraph/{file_names_without_suffix}_seqnorm.filtered.bedgraph"
    shell:
        """
        bamCoverage -b {input.bam} \
        -o {output.bedgraph} \
        -p 'max' \
        --outFileFormat "bedgraph" \
        --normalizeUsing RPKM
        """

rule bamtobed:
    input:
        bam="{wdir}/ChIP-seq_ba/filtered_bam/{file_names_without_suffix}.filtered.bam"
    output:
        sorted_bam="{wdir}/SEACR/bamtobed/{file_names_without_suffix}.filtered.sorted.bam",
        bed="{wdir}/SEACR/bamtobed/{file_names_without_suffix}.filtered.sorted.bed"
    shell:
        """
        samtools sort -n {input.bam} \
        -o {output.sorted_bam}
        bedtools bamtobed -bedpe \
        -i {output.sorted_bam} \
        > {output.bed}
        """
        
rule bed_prep:
    input:
        bedpe = "{wdir}/SEACR/bamtobed/{file_names_without_suffix}.filtered.sorted.bed",
        genome = "/data/repository/organisms/GRCm38_ensembl/genome_fasta/genome.fa.fai"
    output:
        clean_bed = "{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.clean.bed",
        fragments_bed = "{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.fragments.bed",
        fragments_bedgraph = "{wdir}/SEACR/bed3/{file_names_without_suffix}.filtered.fragments.bedgraph"
    shell:
        """
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {input.bedpe} > {output.clean_bed}
        cut -f 1,2,6 {output.clean_bed} | sort -k1,1 -k2,2n -k3,3n > {output.fragments_bed}
        bedtools genomecov -bg -i {output.fragments_bed} -g {input.genome} > {output.fragments_bedgraph}
        """
