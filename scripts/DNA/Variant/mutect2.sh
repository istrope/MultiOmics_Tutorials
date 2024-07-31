#!/bin/bash
#SBATCH --job-name=SC136KCCIS8_mutect2    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail
#SBATCH --ntasks=2                    # Run on a single CPU
#SBATCH --mem=35gb                     # Job memory request
#SBATCH --time=150:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
module load java
module load samtools
module load gatk

for file in SC136_KCIS8_Post/bam/*.bam
do
name=${file##*/}
unfiltered=${name/_final.bam/_unfiltered.vcf}
filtered=${name/_final.bam/_filtered.vcf}
    {
        if [ ! -f SC136_KCIS8_Post/mutect2/$filtered ]
    then
        samtools index $file
        gatk Mutect2 -R Homo_sapiens_assembly38.fasta -I $file -O SC136_KCIS8_Post/mutect2/$unfiltered
        gatk FilterMutectCalls -R Homo_sapiens_assembly38.fasta -V SC136_KCIS8_Post/mutect2/$unfiltered -O SC136_KCIS8_Post/mutect2/$filtered
        rm SC136_KCIS8_Post/mutect2/$unfiltered
    else
        echo $name ' ' 'mutect2 variant calling finished'
    fi
    }
done
