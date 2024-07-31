#!/bin/bash
#SBATCH --job bamtobed    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail    
#SBATCH --ntasks=32                   # Run on a single CPU
#SBATCH --mem=48gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load singularity
module load samtools

#CONVERT TO BED FILE (ORGANIZES ALIGNMENT FILE BY GENOME POSITION AND ANNOTATES REGIONS)
for file in *final.bam
do
output=${file/_final.bam/.bed}
samtools view -b $file | bedtools bamtobed -i /dev/stdin > $output
done

