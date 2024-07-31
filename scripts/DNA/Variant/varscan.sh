#!/bin/bash
#SBATCH --job-name=kcis02_varscan    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=75gb                     # Job memory request
#SBATCH --time=750:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log

module load samtools
module load java
# USE MPILEUP TO GENERATE LIKELIHOOD SCORES AND USE VARSCAN TO GENERATE SNP CALLS (VARIANTS)
samtools mpileup -B -f Homo_sapiens_assembly38.fasta KCIS02_Pre/bam/*final.bam | java -jar VarScan.v2.3.9.jar mpileup2snp > KCIS02_Pre/kcis02_varscan.vcf
