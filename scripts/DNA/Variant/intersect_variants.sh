#!/bin/bash
#THIS CAN BE USED WHEN YOU HAVE RUN MULTIPLE VARIANT CALLING PIPELINES AND WANT TO FIND VARIANTS CALLED IN AT LEAST TWO METHODS
for sample in KCIS*;do
    # Define the directories
    dir1="$sample/bcftools"
    dir2="$sample/mutect2"
    dir3="$sample/freebayes"
    # Output directory for filtered VCFs
    output_dir="$sample/filtered_vcfs"
    mkdir -p "$output_dir"
    
    # Loop over files in the first directory
    for file1 in "$dir1"/*.vcf.gz; do
        base1=$(basename "$file1" .vcf)

        # Check for matching files in directory 2 and 3
        file2="$dir2/$base1.vcf"
        file3="$dir3/$base1.vcf"

        if [[ -f "$file2" && -f "$file3" ]]; then
            # Use bcftools isec to find variants present in at least two of the three files
            bcftools isec -n=+2 -p "$output_dir" -O z $file1 $file2 $file3
            bcftools index "$output_dir/*vcf.gz"
        fi
        echo "Found Intersection For Cell $base1"
    done
    echo "Finished Processing Sample $sample" 
done
echo "All processes completed."
