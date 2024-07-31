#!/bin/bash
#SBATCH --job-name=kcis02_freebayes    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=150:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

eval "$(conda shell.bash hook)"
module load samtools

conda activate freebayes

freebayes -f Homo_sapiens_assembly38.fasta KCIS02_Pre/bam/*.bam | vcffilter -f "QUAL > 20" > KCIS02_Pre/kcis02_freebayes.vcf