#All of the CUT&TAG data was mapped using snakePipe and options as shown below

#bash
#snakePipes/2.5.1: DNA-mapping
DNA-mapping \
-i /data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/CUTandTag-seq_mouse/FASTQ/ \
-o /data/akhtar/group2/zhao/spermatogenic_CUTnX/SSCs_CUTnTAG_2967/ChIP-seq_ba/ \
-j 16 --DAG --fastqc --properPairs --trim \
--alignerOpts '--local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 60 -X 700' \
mm10
