#!/bin/bash
#!/bin/bash
#SBATCH --job-name=mutect2_filter    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=u247529@bcm.edu     # Where to send mail
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --mem=30gb                     # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
module load bcftools
for folder in KCIS*;do #loop through folder
    sample=$(basename $folder)
    sample=${sample/_Pre/}
    sample=${sample/_Post/}
    mkdir -p $sample/mutect2
    for cell in "$folder/mutect2/*.vcf";do #loop through sample vcf
        name=$(basename $cell .vcf)
        bcftools view --threads 8 -i "Qual>20" -Oz -o "$sample/mutect2/${name}.vcf.gz" $cell #filter variants QUAL > 20
        bcftools index --threads 8 "$sample/mutect2/${name}.vcf.gz" #index vcf file
        echo "Finished Filtering and Indexing Cell ID: $name"
    done
    echo "Finished Filtering Sample: $sample"
done
