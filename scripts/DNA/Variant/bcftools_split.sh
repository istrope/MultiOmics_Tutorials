#!/bin/bash
#SBATCH --job-name=filter_bcf    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load bcftools

#THIS SCRIPT TAKES FILES WHERE ALL SAMPLE'S ARE COMBINED INTO ONE FILE AND SPLITS THEM INTO INDIVIDUAL SAMPLES

for file in *.vcf; do
    merged_vcf=$file
    # Output directory for the split VCF files
    sample_dir=$(basename $file .vcf)
    output_dir=${sample_dir}/bcfcalls
    mkdir -p $sample_dir
    mkdir -p $output_dir
    # List all samples from the VCF file
    samples=$(bcftools query -l $merged_vcf)
    echo "Extracting Sample: $sample_dir"
    # Loop through each sample and extract it into a separate VCF file
    for sample in $samples; do
        echo "Extracting cell: $sample"
        out=$(basename $sample)
        bcftools view --threads 4 -s $sample -Oz -o "$output_dir/${out}.vcf.gz" $merged_vcf
        bcftools index --threads 4 "$output_dir/${out}.vcf.gz"
    done
    echo "All samples have been extracted."
done

