#!/bin/bash
reads_path=project/reads_files

aln_fa=resource/clade_refseq_flu_h1n1pdm_ha_MW626062.fasta
aln_target=h1n1pdm_ha_MW626062
barcode_file=resource/clade_barcodes_flu_h1n1pdm_ha_MW626062.csv

bam_path=project/bam_files_${aln_target}
vcf_path=project/vcf_files_${aln_target}
stat_path=project/stat_files

demix_result=project/demix_output/demix_result_${aln_target}.tsv
raw_sample=project/demix_output/PP_raw_${aln_target}.tsv
filtered_sample=project/demix_output/PPFF_${aln_target}.tsv
filtered_barcode=project/demix_output/MMFF_${aln_target}.tsv
potential_sites=project/demix_output/potential_sites_${aln_target}.tsv

THREADS=24

bwa index $aln_fa
samtools faidx $aln_fa
mkdir $bam_path
mkdir $vcf_path
mkdir $stat_path

for sample_name in `ls $reads_path`
do
    sample_dir=$reads_path/$sample_name
    echo handling: $sample_dir
    reads1=$sample_dir/${sample_name}_1.fastq.gz
    reads2=$sample_dir/${sample_name}_2.fastq.gz
    bam_fname=$bam_path/${sample_name}_${aln_target}.bam
    vcf_fname=$vcf_path/${sample_name}_${aln_target}.vcf
    depth_file=$stat_path/${sample_name}_${aln_target}.depth.txt
    coverage_file=$stat_path/${sample_name}_${aln_target}.coverage.txt
    bwa mem -t $THREADS $aln_fa $reads1 $reads2 -v 2 | samtools view -bS --threads $THREADS | samtools sort --threads $THREADS -o $bam_fname
    samtools index $bam_fname
    samtools depth -a $bam_fname > $depth_file
    samtools coverage $bam_fname > $coverage_file
    freebayes --haplotype-length 0 --pooled-continuous -m 20 -p 1 --strict-vcf -q 13 --min-coverage 10 -F 0.4 -f $aln_fa $bam_fname | \
    # bcftools view --include 'QUAL>=10 && (FMT/AO)/(FMT/DP)>=0.01' | \
    bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT | \
    bcftools norm -f $aln_fa -s -m +any -Ov > $vcf_fname
done

python /home/gzy/NextVpower/NextVpower.py -i ${vcf_path} -b ${barcode_file} -o ${demix_result} \
    --vcsample ${raw_sample} --fbarcode ${filtered_barcode} --fsample ${filtered_sample} --potentials ${potential_sites}
